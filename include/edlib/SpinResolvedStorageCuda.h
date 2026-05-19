#ifndef EDLIB_SPINRESOLVEDSTORAGECUDA_H
#define EDLIB_SPINRESOLVEDSTORAGECUDA_H

// GPU spin-resolved storage.
//
// A standalone storage that mirrors the serial (non-MPI) numerics of
// edlib::SpinResolvedStorage but evaluates the matrix-vector product with a
// CUDA kernel and drives the eigensolve through cpp-arnoldi's CudaBackend.
//
// State vector layout (Kronecker / tensor product):
//
//     idx = i_up * down_size + i_down
//
//     H = D  +  (I_up (x) H_down)  +  (H_up (x) I_down)  +  H_loc
//
// where D is the diagonal, H_up / H_down are the per-spin hopping CRS
// matrices and H_loc the (sparse) off-diagonal interaction. The host fill
// is identical in spirit to SpinResolvedStorage::fill (it is reproduced
// here so this class has no dependency on that one); the CRS arrays and
// the diagonal are uploaded to the device once per sector and the four
// contributions are applied by CUDA kernels.
//
// Compilation contract (same as cpp-arnoldi's cuda.hpp):
//   - This header is only active when EDLIB_USE_CUDA is defined.
//   - The translation unit that instantiates a Hamiltonian over this
//     storage MUST be compiled by nvcc (the matvec launches `<<<>>>`).
//
// Scope: single process. There is no MPI domain decomposition here; a
// runtime guard rejects communicators with more than one rank. A concrete
// MPI port plan (decomposition, device halo exchange, distributed a_adag,
// effort tiers) is documented in docs/cuda-mpi.md.
//
// The class is only declared when EDLIB_USE_CUDA is set AND the current
// translation unit is being compiled by a CUDA compiler (__CUDACC__).
// Other TUs in a CUDA-enabled build (plain .cpp) include this header
// transitively via Hamiltonian.h and must see an empty file.

#if defined(EDLIB_USE_CUDA) && defined(__CUDACC__)

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <cuda_runtime.h>

#include <arnoldi/arnoldi.hpp>
#include <arnoldi/cuda.hpp>

#include "edlib/CRSMatrix.h"
#include "edlib/CudaKernel.h"   // cuda_detail utils + cuda::DeviceVec + CudaKernel
#include "edlib/NSymmetry.h"
#include "edlib/Parameters.h"
#include "edlib/Storage.h"
#include "edlib/SzSymmetry.h"

namespace edlib {

  // Storage-specific matvec kernels. The shared CUDA utilities (ck / grid)
  // and the generic vector-op kernels live in CudaKernel.h; cuda_detail is
  // reopened here for these four.
  namespace cuda_detail {

    // w[i] = diag[i]*v[i]              (clear)
    // w[i] += diag[i]*v[i]             (!clear, Lanczos accumulation)
    template <class Prec>
    __global__ void k_diag(long long n, const Prec* diag, const Prec* v, Prec* w, int clear) {
      long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x;
      if (i >= n) return;
      Prec d = diag[i] * v[i];
      w[i] = clear ? d : w[i] + d;
    }

    // (I_up (x) H_down): for output idx = k*down + i, accumulate row i of
    // H_down over the down-subvector of up-block k.
    template <class Prec>
    __global__ void k_spmv_down(long long up_size, long long down_size, const int* rptr, const int* cind,
                                const Prec* val, const Prec* v, Prec* w) {
      long long idx = blockIdx.x * (long long)blockDim.x + threadIdx.x;
      if (idx >= up_size * down_size) return;
      long long k = idx / down_size;
      long long i = idx % down_size;
      Prec acc = Prec(0);
      for (int j = rptr[i]; j < rptr[i + 1]; ++j)
        acc += val[j] * v[cind[j] + k * down_size];
      w[idx] += acc;
    }

    // (H_up (x) I_down): for output idx = i*down + kk, accumulate row
    // (i + up_shift) of H_up, stride down_size along the up index.
    template <class Prec>
    __global__ void k_spmv_up(long long up_size, long long down_size, long long up_shift, const int* rptr,
                              const int* cind, const Prec* val, const Prec* v, Prec* w) {
      long long idx = blockIdx.x * (long long)blockDim.x + threadIdx.x;
      if (idx >= up_size * down_size) return;
      long long i  = idx / down_size;
      long long kk = idx % down_size;
      Prec acc = Prec(0);
      for (int j = rptr[i + up_shift]; j < rptr[i + up_shift + 1]; ++j)
        acc += val[j] * v[(long long)cind[j] * down_size + kk];
      w[i * down_size + kk] += acc;
    }

    // Off-diagonal interaction: full-index CRS over rows [int_start, n).
    template <class Prec>
    __global__ void k_spmv_loc(long long int_start, long long n, const int* rptr, const int* cind,
                               const Prec* val, const Prec* v, Prec* w) {
      long long i = int_start + blockIdx.x * (long long)blockDim.x + threadIdx.x;
      if (i >= n) return;
      Prec acc = Prec(0);
      for (int j = rptr[i]; j < rptr[i + 1]; ++j)
        acc += val[j] * v[cind[j]];
      w[i] += acc;
    }

  }  // namespace cuda_detail

  /**
   * Single-process spin-resolved storage with a CUDA matvec.
   *
   * Uses the cpp-arnoldi CudaBackend solver and keeps eigenvectors on the
   * device. The Green's-function / Lanczos path is fully device-resident via
   * the CudaKernel.
   */
  template <class ModelType>
  class SpinResolvedStorageCuda : public Storage<typename ModelType::precision> {
    static_assert(std::is_base_of<SzSymmetry, typename ModelType::SYMMETRY>::value,
                  "SpinResolvedStorageCuda: model must use SzSymmetry");
  public:
    using Model  = ModelType;
    using prec   = typename ModelType::precision;
    using Matrix = CRSMatrix<prec>;
    using eigenvector_type = cuda::DeviceVec<prec>;
    using kernel_type      = CudaKernel<SpinResolvedStorageCuda>;
    using Sector           = typename Model::Sector;
    using Storage<prec>::n;
    using Storage<prec>::ntot;

#ifdef USE_MPI
    SpinResolvedStorageCuda(const Parameters& p, Model& m, MPI_Comm comm)
        : Storage<prec>(p, comm), _model(m),
          _up_symmetry(p.nsites), _down_symmetry(p.nsites),
          _Ns(p.nsites) {
      int size = 1;
      MPI_Comm_size(comm, &size);
      if (size > 1)
        throw std::runtime_error("SpinResolvedStorageCuda is single-process only "
                                 "(launch with one MPI rank)");
      cache_solver_params(p);
    }
#endif
    SpinResolvedStorageCuda(const Parameters& p, Model& m)
        : Storage<prec>(p), _model(m),
          _up_symmetry(p.nsites), _down_symmetry(p.nsites),
          _Ns(p.nsites) {
      cache_solver_params(p);
    }

    void zero_eigenapair() override {
      _ev_vals.assign(1, _h_diagonal[0]);
      cuda::DeviceVec<prec> v(1);
      prec one = prec(1);
      cuda_detail::ck(cudaMemcpy(v.data(), &one, sizeof(prec), cudaMemcpyHostToDevice),
                      "zero_eigenapair H2D");
      _ev_vecs.clear();
      _ev_vecs.push_back(std::move(v));
    }

    // Host-facing matvec required by the Storage
    void av(prec* v, prec* w, int n_local, bool clear = true) override {
      ensure_io_buffers(n_local);
      cuda_detail::ck(cudaMemcpy(_d_v.data(), v, n_local * sizeof(prec), cudaMemcpyHostToDevice),
                      "av H2D v");
      if (!clear)
        cuda_detail::ck(cudaMemcpy(_d_w.data(), w, n_local * sizeof(prec), cudaMemcpyHostToDevice),
                        "av H2D w");
      matvec_device(_d_v.data(), _d_w.data(), clear, /*stream=*/0);
      cuda_detail::ck(cudaMemcpy(w, _d_w.data(), n_local * sizeof(prec), cudaMemcpyDeviceToHost),
                      "av D2H w");
    }

    // Device matvec entry point used by the CudaKernel
    void device_matvec(const prec* dv, prec* dw, bool clear, cudaStream_t s) {
      matvec_device(dv, dw, clear, s);
    }

    int  num_eigenpairs() const { return static_cast<int>(_ev_vals.size()); }
    const prec& eigenpair_value(int i) const { return _ev_vals[i]; }
    const cuda::DeviceVec<prec>& eigenpair_vector(int i) const { return _ev_vecs[i]; }

    // Host enumeration of the c / c+ scatter for the current sector:
    // idx[ind] = target row in next_sec (or -1), sgn[ind] = fermionic sign.
    void build_adag_map(int i, const Sector& next_sec, bool a, std::size_t locsize,
                        std::vector<int>& idx, std::vector<int>& sgn) {
      idx.assign(locsize, -1);
      sgn.assign(locsize, 0);
      long long k;
      int sign;
      for (std::size_t ind = 0; ind < locsize; ++ind) {
        _model.symmetry().next_state();
        long long nst = _model.symmetry().state();
        if (_model.checkState(nst, i, _model.max_total_electrons()) == (a ? 1 : 0)) {
          if (a) _model.a   (i, nst, k, sign);
          else   _model.adag(i, nst, k, sign);
          idx[ind] = _model.symmetry().index(k, next_sec);
          sgn[ind] = sign;
        }
      }
    }

    void prepare_work_arrays(prec* /*data*/, std::size_t /*shift*/ = 0) override {}

    int finalize(int info, bool /*bcast*/ = true, bool /*empty*/ = true) override { return info; }


    /**
     *  CUDA based diagonalization
     */
    int diag() {
      if (n() == 0) return finalize(0, true, true);
      if (ntot() == 1) {
        zero_eigenapair();
        return finalize(0);
      }

      const int ncv = std::min(_ncv, ntot());
      const int nev = std::min(_nev, ncv - 1);

      arnoldi::Arnoldi<arnoldi::Kind::Sym, prec, arnoldi::SerialComm, arnoldi::CudaBackend>
          solver("I", n(), "SA", nev, ncv);
      solver.tol(prec(1e-14)).maxiter(1000).mode(1).ishift(1);

      cudaStream_t stream = solver.backend().stream();
      solver.solve([this, stream](const prec* x, prec* y) {
        this->matvec_device(x, y, /*clear=*/true, stream);
      });

      const int info = solver.info();
      if (info < 0) return finalize(info);

      const int nconv = solver.num_converged();
      // Eigenvectors stay on the device. No device->host eigenvector copy.
      auto dr = solver.eigenpairs_device(/*compute_vectors=*/_eval_only == 0, prec(0));

      _ev_vals.assign(dr.values.begin(), dr.values.begin() + nconv);
      _ev_vecs.clear();
      if (_eval_only == 0) {
        _ev_vecs.reserve(nconv);
        for (int i = 0; i < nconv; ++i) {
          cuda::DeviceVec<prec> col(n());
          cuda_detail::ck(cudaMemcpy(col.data(),
                                     dr.vectors.data() + std::size_t(i) * n(),
                                     n() * sizeof(prec), cudaMemcpyDeviceToDevice),
                          "diag eigenvector slice");
          _ev_vecs.push_back(std::move(col));
        }
      } else {
        for (int i = 0; i < nconv; ++i) _ev_vecs.emplace_back(0);
      }

      std::cout << "Here is eigenvalues" << std::endl;                                                                                                                                     
      for (auto e : _ev_vals) std::cout << e << "\n";                                                                                                                                      
      std::cout << " ========================= \n"                                                                                                                                         
                << " Size of the matrix is " << ntot() << "\n"                                                                                                                             
                << " The number of converged Ritz values is:  " << nconv << "\n"                                                                                                           
                << " The number of OP*x is: " << solver.num_op_applies() << "\n"                                                                                                           
                << " ========================= " << std::endl;

      return finalize(info);
    }

    // ---- sector lifecycle (serial, host enumeration) ----------------------

    void init() { _model.symmetry().init(); }

    void reset(int /*t*/ = 0) {
      _model.symmetry().init();
      const SzSymmetry& symmetry = static_cast<SzSymmetry&>(_model.symmetry());
      const typename SzSymmetry::Sector& sector = symmetry.sector();
      _up_symmetry  .set_sector(NSymmetry::Sector(sector.nup(),   symmetry.comb().c_n_k(_Ns, sector.nup())));
      _down_symmetry.set_sector(NSymmetry::Sector(sector.ndown(), symmetry.comb().c_n_k(_Ns, sector.ndown())));
      const std::size_t up_size   = _up_symmetry  .sector().size();
      const std::size_t down_size = _down_symmetry.sector().size();
      H_up  .init(up_size,   100);
      H_down.init(down_size, 100);
      _locsize  = up_size * down_size;
      _up_size  = up_size;
      _down_size = down_size;
      _up_shift = 0;
      _diagonal.assign(_locsize, prec(0));
      if (_model.V_states().size() > 0) H_loc.init(_locsize, 3);
      n()    = static_cast<int>(_locsize);
      ntot() = static_cast<int>(sector.size());
    }

    void fill() {
      reset();
      if (n() == 0) return;

      fill_spin(_up_symmetry,   _Ns, H_up);
      fill_spin(_down_symmetry, 0,   H_down);

      int isign;
      long long k;
      _int_start = _locsize;
      for (std::size_t i = 0; i < _locsize; ++i) {
        _model.symmetry().next_state();
        long long nst = _model.symmetry().state();
        _diagonal[i] = _model.diagonal(nst);
        if (_model.V_states().size() > 0) {
          for (std::size_t kkk = 0; kkk < _model.V_states().size(); ++kkk) {
            if (_model.valid(_model.V_states()[kkk], nst)) {
              _int_start = std::min(i, _int_start);
              _model.set(_model.V_states()[kkk], nst, k, isign);
              int j = _model.symmetry().index(k);
              H_loc.addElement(i, j, _model.V_states()[kkk].value(), isign);
            }
          }
          H_loc.endLine(i);
        }
      }
      _h_diagonal = _diagonal;
      upload_to_device();
    }

    std::size_t vector_size(typename Model::Sector sector) const { return sector.size(); }

    void constant_shift(prec shift) {
      std::transform(_h_diagonal.begin(), _h_diagonal.end(), _h_diagonal.begin(),
                     [shift](prec x) { return x + shift; });
      if (_d_diag.size() == _h_diagonal.size() && !_h_diagonal.empty())
        cuda_detail::ck(cudaMemcpy(_d_diag.data(), _h_diagonal.data(),
                                   _h_diagonal.size() * sizeof(prec), cudaMemcpyHostToDevice),
                        "constant_shift H2D");
    }

#ifdef USE_MPI
    MPI_Comm comm() override { return Storage<prec>::comm(); }
    std::size_t offset() const { return 0; }
#endif

  private:
    void cache_solver_params(const Parameters& p) {
      _nev       = p.arpack_nev;
      _ncv       = p.arpack_ncv > 0 ? p.arpack_ncv : 2 * p.arpack_nev + 3;
      _eval_only = p.eigenvalues_only ? 1 : 0;
    }

    void fill_spin(NSymmetry& spin_symmetry, int shift, Matrix& spin_matrix) {
      long long k = 0;
      int isign = 0;
      int i = 0;
      while (spin_symmetry.next_state()) {
        long long nst = spin_symmetry.state();
        for (std::size_t kkk = 0; kkk < _model.T_states().size(); ++kkk) {
          if (_model.valid(_model.T_states()[kkk], nst << shift)) {
            _model.set(_model.T_states()[kkk], nst << shift, k, isign);
            int j = spin_symmetry.index(k >> shift);
            spin_matrix.addElement(i, j, _model.T_states()[kkk].value(), isign);
          }
        }
        spin_matrix.endLine(i);
        ++i;
      }
    }

    template <class T>
    static void up(arnoldi::detail::device_vector<T>& d, const std::vector<T>& h) {
      d.assign(h.size(), T{});
      if (!h.empty())
        cuda_detail::ck(cudaMemcpy(d.data(), h.data(), h.size() * sizeof(T), cudaMemcpyHostToDevice),
                        "upload H2D");
    }

    void upload_to_device() {
      up(_d_diag,      _h_diagonal);
      up(_d_up_rptr,   H_up.row_ptr());
      up(_d_up_cind,   H_up.col_ind());
      up(_d_up_val,    H_up.values());
      up(_d_down_rptr, H_down.row_ptr());
      up(_d_down_cind, H_down.col_ind());
      up(_d_down_val,  H_down.values());
      _has_loc = _model.V_states().size() > 0 && !H_loc.row_ptr().empty();
      if (_has_loc) {
        up(_d_loc_rptr, H_loc.row_ptr());
        up(_d_loc_cind, H_loc.col_ind());
        up(_d_loc_val,  H_loc.values());
      }
    }

    void ensure_io_buffers(int n_local) {
      if ((int)_d_v.size() != n_local) {
        _d_v.assign(n_local, prec(0));
        _d_w.assign(n_local, prec(0));
      }
    }

    // Launch the four contributions on `stream` (0 = default/Lanczos path).
    // Same-stream ordering serialises the kernels, so the diagonal write
    // followed by the three accumulations needs no atomics.
    void matvec_device(const prec* dv, prec* dw, bool clear, cudaStream_t stream) {
      const int block = 256;
      const long long total = (long long)_up_size * _down_size;

      cuda_detail::k_diag<prec><<<cuda_detail::grid(total, block), block, 0, stream>>>(
          total, _d_diag.data(), dv, dw, clear ? 1 : 0);

      cuda_detail::k_spmv_down<prec><<<cuda_detail::grid(total, block), block, 0, stream>>>(
          (long long)_up_size, (long long)_down_size,
          _d_down_rptr.data(), _d_down_cind.data(), _d_down_val.data(), dv, dw);

      cuda_detail::k_spmv_up<prec><<<cuda_detail::grid(total, block), block, 0, stream>>>(
          (long long)_up_size, (long long)_down_size, (long long)_up_shift,
          _d_up_rptr.data(), _d_up_cind.data(), _d_up_val.data(), dv, dw);

      if (_has_loc) {
        const long long loc_rows = (long long)_locsize - (long long)_int_start;
        if (loc_rows > 0)
          cuda_detail::k_spmv_loc<prec><<<cuda_detail::grid(loc_rows, block), block, 0, stream>>>(
              (long long)_int_start, (long long)_locsize,
              _d_loc_rptr.data(), _d_loc_cind.data(), _d_loc_val.data(), dv, dw);
      }
      cuda_detail::ck(cudaGetLastError(), "matvec kernel launch");
      if (stream == 0) cuda_detail::ck(cudaStreamSynchronize(0), "matvec sync");
    }

    Model&  _model;
    Matrix  H_loc, H_up, H_down;

    std::vector<prec> _diagonal;     // working diagonal during fill
    std::vector<prec> _h_diagonal;   // host master copy (constant_shift, zero pair)

    NSymmetry _up_symmetry;
    NSymmetry _down_symmetry;
    int       _Ns;

    std::size_t _up_size   = 0;
    std::size_t _down_size = 0;
    std::size_t _up_shift  = 0;
    std::size_t _locsize   = 0;
    std::size_t _int_start = 0;
    bool        _has_loc   = false;

    int _nev = 1, _ncv = 0, _eval_only = 0;

    // Device-resident matrix + I/O scratch.
    arnoldi::detail::device_vector<prec> _d_diag;
    arnoldi::detail::device_vector<int>  _d_up_rptr,   _d_up_cind;
    arnoldi::detail::device_vector<prec> _d_up_val;
    arnoldi::detail::device_vector<int>  _d_down_rptr, _d_down_cind;
    arnoldi::detail::device_vector<prec> _d_down_val;
    arnoldi::detail::device_vector<int>  _d_loc_rptr,  _d_loc_cind;
    arnoldi::detail::device_vector<prec> _d_loc_val;
    arnoldi::detail::device_vector<prec> _d_v, _d_w;

    // Device-resident eigenpairs (values on host, vectors on device).
    std::vector<prec>                  _ev_vals;
    std::vector<cuda::DeviceVec<prec>> _ev_vecs;
  };

}  // namespace edlib

#endif  // EDLIB_USE_CUDA && __CUDACC__

#endif  // EDLIB_SPINRESOLVEDSTORAGECUDA_H
