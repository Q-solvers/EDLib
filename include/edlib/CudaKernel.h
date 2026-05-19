#ifndef EDLIB_CUDAKERNEL_H
#define EDLIB_CUDAKERNEL_H

// Device-resident Lanczos kernel + its CUDA support layer.
//
// This is the GPU counterpart of HostKernel.h. A storage exposes
// `using kernel_type = CudaKernel<ThatStorage>;` and Lanczos pulls it from
// there (no template parameter is threaded through Lanczos / GreensFunction
// / ChiLoc).

#if defined(EDLIB_USE_CUDA) && defined(__CUDACC__)

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <cuda_runtime.h>

#include <arnoldi/arnoldi.hpp>
#include <arnoldi/cuda.hpp>

namespace edlib {

  namespace cuda_detail {

    inline void ck(cudaError_t e, const char* what) {
      if (e != cudaSuccess)
        throw std::runtime_error(std::string("edlib/cuda: ") + what + ": " + cudaGetErrorString(e));
    }

    inline int grid(long long total, int block) {
      return static_cast<int>((total + block - 1) / block);
    }

    // ---- generic device-resident Lanczos vector ops -----------------------

    // Lanczos three-term recurrence swap: d=v; v=w/bet; w=-bet*d.
    template <class Prec>
    __global__ void k_recurrence(long long n, Prec* v, Prec* w, Prec bet) {
      long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x;
      if (i >= n) return;
      Prec d = v[i];
      v[i]   = w[i] / bet;
      w[i]   = -bet * d;
    }

    template <class Prec>
    __global__ void k_axpy(long long n, Prec a, const Prec* x, Prec* y) {
      long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x;
      if (i < n) y[i] += a * x[i];
    }

    template <class Prec>
    __global__ void k_scale(long long n, Prec* v, Prec s) {
      long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x;
      if (i < n) v[i] *= s;
    }

    // c / c+ scatter: out[idx[ind]] = sgn[ind] * in[ind]  (idx<0 => skip).
    // idx is a permutation, so plain assignment is race-free.
    template <class Prec>
    __global__ void k_scatter(long long nind, const int* idx, const int* sgn,
                              const Prec* in, Prec* out) {
      long long ind = blockIdx.x * (long long)blockDim.x + threadIdx.x;
      if (ind >= nind) return;
      int t = idx[ind];
      if (t >= 0) out[t] = static_cast<Prec>(sgn[ind]) * in[ind];
    }

  }  // namespace cuda_detail

  namespace cuda {

    // Reference-counted RAII device buffer. Copyable (shared ownership) so
    // it can live inside a value-semantic EigenPair / std::set; the actual
    // cudaMalloc is shared, never deep-copied. Used both as the storage's
    // eigenvector_type and as the CudaKernel working-vector type, so the
    // eigenvector, the c/c+ output, and the Krylov vectors all share one
    // device-resident representation.
    template <class T>
    class DeviceVec {
    public:
      DeviceVec() = default;
      explicit DeviceVec(std::size_t n) : _n(n) {
        if (n) {
          T* d = nullptr;
          cuda_detail::ck(cudaMalloc(&d, n * sizeof(T)), "DeviceVec cudaMalloc");
          _p = std::shared_ptr<T>(d, [](T* q) { cudaFree(q); });
        }
      }
      T*          data()  const { return _p.get(); }
      std::size_t size()  const { return _n; }
    private:
      std::shared_ptr<T> _p;
      std::size_t        _n = 0;
    };

  }  // namespace cuda

  /**
   * Fully device-resident Lanczos kernel.
   *
   * All n-size vectors live in cuda::DeviceVec (GPU). There is NO per-iteration host<->device staging:
   * av is the storage device matvec on a private stream, vv is cuBLAS dot,
   * the three-term recurrence / axpy / scale are device kernels.
   * The stream is synchronised only when a host scalar is actually needed
   * (inside dot), i.e. ~twice per Lanczos iteration instead of every av.
   */
  template <class Storage>
  class CudaKernel {
  public:
    using prec   = typename Storage::prec;
    using Model  = typename Storage::Model;
    using Sector = typename Model::Sector;
    using Vector = cuda::DeviceVec<prec>;

    explicit CudaKernel(Storage& s) : _s(s) {
      cuda_detail::ck(cudaStreamCreate(&_stream), "CudaKernel stream");
      arnoldi::detail::cublas_check(cublasCreate(&_h), "CudaKernel cublasCreate");
      arnoldi::detail::cublas_check(cublasSetStream(_h, _stream), "CudaKernel setStream");
      arnoldi::detail::cublas_check(cublasSetPointerMode(_h, CUBLAS_POINTER_MODE_HOST),
                                    "CudaKernel pointerMode");
    }
    CudaKernel(const CudaKernel&)            = delete;
    CudaKernel& operator=(const CudaKernel&) = delete;
    CudaKernel(CudaKernel&& o) noexcept
        : _s(o._s), _h(o._h), _stream(o._stream),
          _d_idx(std::move(o._d_idx)), _d_sgn(std::move(o._d_sgn)) {
      o._h = nullptr; o._stream = nullptr;
    }
    ~CudaKernel() {
      if (_h) cublasDestroy(_h);
      if (_stream) cudaStreamDestroy(_stream);
    }

    Vector      make_vector(std::size_t n) const {
      Vector v(n);
      if (n) cuda_detail::ck(cudaMemsetAsync(v.data(), 0, n * sizeof(prec), _stream),
                             "make_vector memset");
      return v;
    }
    std::size_t size(const Vector& v) const { return v.size(); }
    void        prepare(Vector&) {}

    // w <- w + H v   (no host staging, no per-call sync)
    void av(Vector& v, Vector& w) {
      _s.device_matvec(v.data(), w.data(), /*clear=*/false, _stream);
    }

    prec dot(const Vector& a, const Vector& b) {
      prec r = prec(0);
      int  n = static_cast<int>(a.size());
      if constexpr (std::is_same_v<prec, double>)
        arnoldi::detail::cublas_check(cublasDdot(_h, n, a.data(), 1, b.data(), 1, &r), "cublasDdot");
      else
        arnoldi::detail::cublas_check(cublasSdot(_h, n, a.data(), 1, b.data(), 1, &r), "cublasSdot");
      cuda_detail::ck(cudaStreamSynchronize(_stream), "dot sync");
      return r;
    }

    void recurrence(Vector& v, Vector& w, prec bet) {
      const std::size_t n = v.size();
      if (!n) return;
      cuda_detail::k_recurrence<prec><<<g(n), B, 0, _stream>>>(
          (long long)n, v.data(), w.data(), bet);
      cuda_detail::ck(cudaGetLastError(), "k_recurrence");
    }
    void axpy(prec a, const Vector& x, Vector& y) {
      const std::size_t n = y.size();
      if (!n) return;
      cuda_detail::k_axpy<prec><<<g(n), B, 0, _stream>>>(
          (long long)n, a, x.data(), y.data());
      cuda_detail::ck(cudaGetLastError(), "k_axpy");
    }
    void scale(Vector& v, prec s) {
      const std::size_t n = v.size();
      if (!n) return;
      cuda_detail::k_scale<prec><<<g(n), B, 0, _stream>>>((long long)n, v.data(), s);
      cuda_detail::ck(cudaGetLastError(), "k_scale");
    }
    void add(Vector& acc, const Vector& t) {
      const std::size_t n = acc.size();
      if (!n) return;
      cuda_detail::k_axpy<prec><<<g(n), B, 0, _stream>>>(
          (long long)n, prec(1), t.data(), acc.data());
      cuda_detail::ck(cudaGetLastError(), "k_add");
    }

    // c / c+ application: host builds the (idx, sign) map for the current
    // sector; the device scatters the (device-resident) eigenvector.
    void a_adag(int op, const Vector& invec, Vector& out,
                const Sector& next, bool annihilate) {
      const std::size_t loc = invec.size();
      _s.build_adag_map(op, next, annihilate, loc, _h_idx, _h_sgn);
      grow(_d_idx, loc);
      grow(_d_sgn, loc);
      cuda_detail::ck(cudaMemcpyAsync(_d_idx.data(), _h_idx.data(), loc * sizeof(int),
                                      cudaMemcpyHostToDevice, _stream), "a_adag idx H2D");
      cuda_detail::ck(cudaMemcpyAsync(_d_sgn.data(), _h_sgn.data(), loc * sizeof(int),
                                      cudaMemcpyHostToDevice, _stream), "a_adag sgn H2D");
      const int block = 256;
      cuda_detail::k_scatter<prec><<<cuda_detail::grid((long long)loc, block), block, 0, _stream>>>(
          (long long)loc, _d_idx.data(), _d_sgn.data(), invec.data(), out.data());
      cuda_detail::ck(cudaGetLastError(), "a_adag scatter launch");
    }

  private:
    static constexpr int B = 256;
    static int g(std::size_t n) { return cuda_detail::grid((long long)n, B); }
    static void grow(arnoldi::detail::device_vector<int>& d, std::size_t n) {
      if (d.size() < n) d.assign(n, 0);
    }

    Storage&       _s;
    cublasHandle_t _h      = nullptr;
    cudaStream_t   _stream = nullptr;
    std::vector<int>                    _h_idx, _h_sgn;
    arnoldi::detail::device_vector<int> _d_idx, _d_sgn;
  };

}  // namespace edlib

#endif  // EDLIB_USE_CUDA && __CUDACC__

#endif  // EDLIB_CUDAKERNEL_H
