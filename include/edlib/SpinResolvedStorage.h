#ifndef EDLIB_SPINRESOLVEDSTORAGE_H
#define EDLIB_SPINRESOLVEDSTORAGE_H

#include <algorithm>
#include <bitset>
#include <climits>
#include <cstddef>
#include <iomanip>
#include <type_traits>
#include <vector>

#include "edlib/CRSMatrix.h"
#include "edlib/HostKernel.h"
#include "edlib/MpiTypes.h"
#include "edlib/NSymmetry.h"
#include "edlib/Parameters.h"
#include "edlib/Storage.h"
#include "edlib/SzSymmetry.h"

namespace edlib {

  /**
   * Spin-resolved storage: keeps hopping matrices for each spin channel
   * separately, and the (smaller) off-diagonal interaction part. MPI-aware
   * via RMA on the spin-up channel.
   */
  template <class ModelType>
  class SpinResolvedStorage : public Storage<typename ModelType::precision> {
    static_assert(std::is_base_of<SzSymmetry, typename ModelType::SYMMETRY>::value,
                  "SpinResolvedStorage: model must use SzSymmetry");
  public:
    using Model = ModelType;
    using prec  = typename ModelType::precision;
    using kernel_type = HostKernel<SpinResolvedStorage>;
    using Matrix = CRSMatrix<prec>;
    using Storage<prec>::n;
    using Storage<prec>::ntot;
#ifdef USE_MPI
    using Storage<prec>::comm;
    using Storage<prec>::broadcast_evals;
#endif
    using Storage<prec>::prepare_work_arrays;
    using Storage<prec>::finalize;

#ifdef USE_MPI
    SpinResolvedStorage(const Parameters& p, Model& m, MPI_Comm comm)
        : Storage<prec>(p, comm),
          _model(m),
          _up_symmetry  (p.nsites),
          _down_symmetry(p.nsites),
          _interaction_size(m.interacting_orbitals()),
          _Ns(p.nsites), _ms(p.nspins),
          _comm(comm), _run_comm(MPI_COMM_NULL), _win(MPI_WIN_NULL) {
      MPI_Comm_size(_comm, &_nprocs);
      MPI_Comm_rank(_comm, &_myid);
    }
#else
    SpinResolvedStorage(const Parameters& p, Model& m)
        : Storage<prec>(p),
          _model(m),
          _up_symmetry  (p.nsites),
          _down_symmetry(p.nsites),
          _interaction_size(m.interacting_orbitals()),
          _Ns(p.nsites), _ms(p.nspins) {}
#endif

    void zero_eigenapair() override {
      this->eigenvalues().resize(1);
      this->eigenvalues()[0] = _diagonal[0];
      this->eigenvectors().assign(1, std::vector<prec>(1, prec(1)));
    }

    void av(prec* v, prec* w, int n_local, bool clear = true) override {
#ifdef USE_MPI
      MPI_Win_fence(MPI_MODE_NOPRECEDE, _win);
      for (std::size_t i = 0; i < _procs.size(); ++i) {
        if (_procs[i] != 0) {
          MPI_Get(&_vecval[_proc_offset[i]], _proc_size[i], mpi_type<prec>(),
                  static_cast<int>(i), _loc_min[i], _proc_size[i], mpi_type<prec>(), _win);
        }
      }
#endif
      for (int i = 0; i < n_local; ++i) {
        w[i] = _diagonal[i] * v[i] + (clear ? prec(0) : w[i]);
      }
      // spin-down hopping
      for (int k = 0; k < (int)_up_size; ++k) {
        for (int i = 0; i < (int)_down_symmetry.sector().size(); ++i) {
          for (int j = H_down.row_ptr()[i]; j < H_down.row_ptr()[i + 1]; ++j) {
            w[i + k * _down_symmetry.sector().size()] +=
                H_down.values()[j] * v[H_down.col_ind()[j] + k * _down_symmetry.sector().size()];
          }
        }
      }
#ifdef USE_MPI
      MPI_Win_fence(MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE, _win);
#endif
      // spin-up hopping
      for (int i = 0; i < (int)_up_size; ++i) {
        for (int j = H_up.row_ptr()[i + _up_shift]; j < H_up.row_ptr()[i + _up_shift + 1]; ++j) {
          for (int k = 0; k < (int)_down_symmetry.sector().size(); ++k) {
#ifdef USE_MPI
            w[i * _down_symmetry.sector().size() + k] +=
                H_up.values()[j] * _vecval[H_up.col_ind()[j] * _down_symmetry.sector().size() + k];
#else
            w[i * _down_symmetry.sector().size() + k] +=
                H_up.values()[j] * v[H_up.col_ind()[j] * _down_symmetry.sector().size() + k];
#endif
          }
        }
      }
      // off-diagonal interaction
      if (!H_loc.row_ptr().empty()) {
        for (std::size_t i = _int_start; i < (std::size_t)n_local; ++i) {
          for (int j = H_loc.row_ptr()[i]; j < H_loc.row_ptr()[i + 1]; ++j) {
#ifdef USE_MPI
            w[i] += H_loc.values()[j] * _vecval[H_loc.col_ind()[j]];
#else
            w[i] += H_loc.values()[j] * v[H_loc.col_ind()[j]];
#endif
          }
        }
      }
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
#ifdef USE_MPI
      find_neighbours();
#endif
    }

    void init() {
      _model.symmetry().init();
#ifdef USE_MPI
      _model.symmetry().set_offset(_offset);
#endif
    }

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
#ifdef USE_MPI
      if (_comm != _run_comm && _run_comm != MPI_COMM_NULL) MPI_Comm_free(&_run_comm);
      int color = _myid < (int)up_size ? 1 : 0;
      MPI_Comm_split(_comm, color, (color == 1 ? _myid : 0), &_run_comm);
      if (color == 1) {
        int myid; MPI_Comm_rank(_run_comm, &myid);
        int size; MPI_Comm_size(_run_comm, &size);
        std::size_t locsize = up_size / size;
        if ((up_size % size) > (std::size_t)myid) {
          locsize += 1;
          _offset = myid * locsize * _down_symmetry.sector().size();
        } else {
          _offset = (myid * locsize + (up_size % size)) * _down_symmetry.sector().size();
        }
        _up_size  = locsize;
        _up_shift = _offset / down_size;
        _locsize  = locsize * down_size;
        _model.symmetry().set_offset(_offset);
        _procs       .assign(size, 0);
        _proc_offset .assign(size, 0);
        _proc_size   .assign(size, 0);
        _loc_min     .assign(size, 0);
      } else {
        _up_size = 0;
        _locsize = 0;
      }
#else
      _locsize  = up_size * down_size;
      _up_size  = up_size;
      _up_shift = 0;
#endif
      _diagonal.assign(_locsize, prec(0));
      if (_model.V_states().size() > 0) H_loc.init(_locsize, 3);
      n()    = static_cast<int>(_locsize);
      ntot() = static_cast<int>(sector.size());
    }

    std::size_t vector_size(typename Model::Sector sector) const {
      const std::size_t sector_size = sector.size();
#ifdef USE_MPI
      int myid, size;
      MPI_Comm_rank(_comm, &myid);
      MPI_Comm_size(_comm, &size);
      std::size_t up_size   = _model.symmetry().comb().c_n_k(_Ns, sector.nup());
      std::size_t down_size = sector_size / up_size;
      size = (int)(up_size > (std::size_t)size ? size : up_size);
      if (myid >= size) return 0;
      std::size_t locsize = up_size / size;
      if ((up_size % size) > (std::size_t)myid) locsize += 1;
      return locsize * down_size;
#else
      return sector_size;
#endif
    }

    void a_adag(int i, const std::vector<prec>& invec, std::vector<prec>& outvec,
                const typename Model::Sector& next_sec, bool a) {
      std::size_t locsize     = invec.size();
      std::size_t locsize_max = locsize;
      std::size_t next_size   = next_sec.size();
      std::size_t up_size     = _model.symmetry().comb().c_n_k(_Ns, next_sec.nup());
      std::size_t down_size   = next_size / up_size;
      long long k;
      int sign;
#ifdef USE_MPI
      MPI_Allreduce(MPI_IN_PLACE, &locsize_max, 1, mpi_type<std::size_t>(), MPI_MAX, _comm);
      int ci, cid;
      int myid; MPI_Comm_rank(_comm, &myid);
      int size; MPI_Comm_size(_comm, &size);
      int t = 0;
      bool fence = false;
      if (_up_symmetry.sector().size() % size != 0) locsize_max += _down_symmetry.sector().size();
      std::vector<prec> buff(1000, prec(0));
      int rank;
      MPI_Comm_rank(_run_comm, &rank);
      MPI_Comm_size(_run_comm, &size);
      size = outvec.size() > 0 ? 1 : 0;
      MPI_Allreduce(MPI_IN_PLACE, &size, 1, mpi_type<int>(), MPI_SUM, _comm);
      MPI_Win eigwin;
      MPI_Win_create(outvec.data(), sizeof(prec) * outvec.size(), sizeof(prec),
                     MPI_INFO_NULL, _comm, &eigwin);
      MPI_Win_fence(MPI_MODE_NOPRECEDE, eigwin);
#endif
      for (std::size_t ind = 0; ind < locsize_max; ++ind) {
#ifdef USE_MPI
        if (fence) MPI_Win_fence(MPI_MODE_NOPRECEDE, eigwin);
        fence = false;
#endif
        if (ind < locsize) {
          _model.symmetry().next_state();
          long long nst = _model.symmetry().state();
          if (_model.checkState(nst, i, _model.max_total_electrons()) == (a ? 1 : 0)) {
            if (a) _model.a   (i, nst, k, sign);
            else   _model.adag(i, nst, k, sign);
            int i1 = _model.symmetry().index(k, next_sec);
#ifdef USE_MPI
            calcIndex(ci, cid, i1, up_size, down_size, size);
            if (myid == cid) {
              outvec[ci] = sign * invec[ind];
            } else {
              buff[t] = sign * invec[ind];
              MPI_Put(&buff[t], 1, mpi_type<prec>(), cid, ci, 1, mpi_type<prec>(), eigwin);
            }
#else
            outvec[i1] = sign * invec[ind];
#endif
          }
        }
#ifdef USE_MPI
        if ((++t) == (int)buff.size()) { fence = true; t = 0; }
        if (fence) MPI_Win_fence(MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE, eigwin);
#endif
      }
#ifdef USE_MPI
      if (!fence) MPI_Win_fence(MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE, eigwin);
      MPI_Win_free(&eigwin);
#endif
    }

    prec vv(const std::vector<prec>& v, const std::vector<prec>& w) const {
#ifdef USE_MPI
      return vv(v, w, const_cast<SpinResolvedStorage*>(this)->comm());
#else
      prec alf = prec(0);
      for (std::size_t k = 0; k < v.size(); ++k) alf += w[k] * v[k];
      return alf;
#endif
    }

#ifdef USE_MPI
    prec vv(const std::vector<prec>& v, const std::vector<prec>& w, MPI_Comm com) const {
      prec temp = prec(0), alf = prec(0);
      for (std::size_t k = 0; k < v.size(); ++k) temp += w[k] * v[k];
      MPI_Allreduce(&temp, &alf, 1, mpi_type<prec>(), MPI_SUM, com);
      return alf;
    }

    void prepare_work_arrays(prec* data, std::size_t shift = 0) override {
      MPI_Win_create(&data[shift], (MPI_Aint)(n() * sizeof(prec)),
                     (int)sizeof(prec), MPI_INFO_NULL, _run_comm, &_win);
    }

    MPI_Comm comm() override { return _run_comm; }

    int finalize(int info, bool bcast = true, bool empty = true) override {
      MPI_Bcast(&info, 1, MPI_INT, 0, Storage<prec>::comm());
      if (info >= 0 && bcast) broadcast_evals(empty);
      if (ntot() > 1 && n() > 0) MPI_Win_free(&_win);
      if (_run_comm != MPI_COMM_WORLD && _run_comm != MPI_COMM_NULL) MPI_Comm_free(&_run_comm);
      _run_comm = Storage<prec>::comm();
      return info;
    }

    std::size_t offset() const { return _offset; }
#endif

    void constant_shift(prec shift) {
      std::transform(_diagonal.begin(), _diagonal.end(), _diagonal.begin(),
                     [shift](prec x) { return x + shift; });
    }

  private:
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

#ifdef USE_MPI
    void find_neighbours() {
      int ci, cid;
      int nprocs; MPI_Comm_size(_run_comm, &nprocs);
      std::vector<int> loc_offset(nprocs, 0);
      std::vector<int> l_loc_max(_loc_min.size(), INT_MIN);
      std::vector<int> l_loc_min(_loc_min.size(), INT_MAX);
      for (int i = 0; i < (int)_up_size; ++i) {
        for (int j = H_up.row_ptr()[i + _up_shift]; j < H_up.row_ptr()[i + _up_shift + 1]; ++j) {
          calcIndex(ci, cid,
                    H_up.col_ind()[j] * _down_symmetry.sector().size(),
                    _up_symmetry.sector().size(), _down_symmetry.sector().size(), nprocs);
          l_loc_max[cid] = std::max(ci, l_loc_max[cid]);
          l_loc_min[cid] = std::min(ci, l_loc_min[cid]);
          if (_procs[cid] == 0) _procs[cid] = 1;
        }
      }
      if (!H_loc.row_ptr().empty()) {
        for (std::size_t i = _int_start; i < _locsize; ++i) {
          for (int j = H_loc.row_ptr()[i]; j < H_loc.row_ptr()[i + 1]; ++j) {
            calcIndex(ci, cid,
                      _down_symmetry.sector().size() *
                          (H_loc.col_ind()[j] / _down_symmetry.sector().size()),
                      _up_symmetry.sector().size(), _down_symmetry.sector().size(), nprocs);
            l_loc_max[cid] = std::max(ci, l_loc_max[cid]);
            l_loc_min[cid] = std::min(ci, l_loc_min[cid]);
            if (_procs[cid] == 0) _procs[cid] = 1;
          }
        }
      }
      int oset = 0;
      for (int i = 0; i < nprocs; ++i) {
        if (_procs[i]) {
          _procs[i] = 1;
          _proc_offset[i] = oset * _down_symmetry.sector().size() + l_loc_min[i];
          _loc_min[i]     = l_loc_min[i];
          int ls = _up_symmetry.sector().size() / nprocs;
          if ((_up_symmetry.sector().size() % nprocs) > (std::size_t)i) {
            ls++;
            loc_offset[i] = (i * ls) - oset;
          } else {
            loc_offset[i] = i * ls + (_up_symmetry.sector().size() % nprocs) - oset;
          }
          _proc_size[i] = l_loc_max[i] - l_loc_min[i] + _down_symmetry.sector().size();
          oset += ls;
        }
      }
      _vecval.assign(oset * _down_symmetry.sector().size(), prec(0));
      for (int i = 0; i < (int)_up_size; ++i) {
        for (int j = H_up.row_ptr()[i + _up_shift]; j < H_up.row_ptr()[i + _up_shift + 1]; ++j) {
          calcIndex(ci, cid,
                    H_up.col_ind()[j] * _down_symmetry.sector().size(),
                    _up_symmetry.sector().size(), _down_symmetry.sector().size(), nprocs);
          H_up.col_ind()[j] -= loc_offset[cid];
        }
      }
      if (!H_loc.row_ptr().empty()) {
        for (std::size_t i = _int_start; i < _locsize; ++i) {
          for (int j = H_loc.row_ptr()[i]; j < H_loc.row_ptr()[i + 1]; ++j) {
            calcIndex(ci, cid, H_loc.col_ind()[j],
                      _up_symmetry.sector().size(), _down_symmetry.sector().size(), nprocs);
            H_loc.col_ind()[j] -= loc_offset[cid] * _down_symmetry.sector().size();
          }
        }
      }
    }

    void calcIndex(int& ci, int& cid, int i, std::size_t u_s, std::size_t d_s, int nprocs) const {
      int i_rest = i % d_s;
      int i_up   = i / d_s;
      int tmp1 = u_s / nprocs + 1;
      int tmp2 = u_s % nprocs;
      int tmp3 = u_s / nprocs;
      int tmp4 = i_up - (tmp1 * tmp2);
      if (i_up > (tmp1 * tmp2)) {
        ci  = (tmp4 % tmp3) * d_s + i_rest;
        cid = (i_up - tmp2) / tmp3;
      } else {
        ci  = (i_up % (tmp3 + 1)) * d_s + i_rest;
        cid = i_up / (tmp3 + 1);
      }
    }
#endif

    Model& _model;
    Matrix H_loc;
    Matrix H_up;
    Matrix H_down;

    std::vector<prec> _diagonal;
    std::vector<prec> _vecval;

    NSymmetry _up_symmetry;
    NSymmetry _down_symmetry;

    int _interaction_size;
    int _Ns;
    int _ms;

    std::size_t _up_size  = 0;
    std::size_t _up_shift = 0;
    std::size_t _locsize  = 0;
    std::size_t _int_start = 0;

#ifdef USE_MPI
    MPI_Comm    _comm;
    MPI_Comm    _run_comm;
    std::size_t _offset = 0;
    int         _myid = 0;
    int         _nprocs = 0;
    std::vector<int> _proc_offset;
    std::vector<int> _procs;
    std::vector<int> _loc_min;
    std::vector<int> _proc_size;
    MPI_Win     _win;
#endif
  };

}

#endif
