#ifndef EDLIB_HOSTKERNEL_H
#define EDLIB_HOSTKERNEL_H

#include <cstddef>
#include <vector>

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace edlib {

  /**
   * Host (CPU) Lanczos kernel.
   *
   * A kernel owns the Lanczos working vectors and implements every
   * vector-level operation the Lanczos / Green's-function needs.
   *
   * Storages expose their kernel via `using kernel_type = ...;` and Lanczos
   * pulls it from there -- no template parameter is threaded through
   * Lanczos / GreensFunction / ChiLoc.
   *
   */
  template <class Storage>
  class HostKernel {
  public:
    using prec   = typename Storage::prec;
    using Vector = std::vector<prec>;
    using Model  = typename Storage::Model;
    using Sector = typename Model::Sector;

    explicit HostKernel(Storage& s) : _s(s) {}

    Vector      make_vector(std::size_t n) const { return Vector(n, prec(0)); }
    std::size_t size(const Vector& v)      const { return v.size(); }

    void prepare(Vector& v) { _s.prepare_work_arrays(v.data()); }

    // w <- w + H v   (storage av with clear == false)
    void av(Vector& v, Vector& w) {
      _s.av(v.data(), w.data(), static_cast<int>(v.size()), /*clear=*/false);
    }

    prec dot(const Vector& a, const Vector& b) const { return _s.vv(a, b); }

#ifdef USE_MPI
    prec dot(const Vector& a, const Vector& b, MPI_Comm comm) const {
      return _s.vv(a, b, comm);
    }
#endif

    // Lanczos three-term recurrence vector swap:
    //   dummy = v;  v = w / bet;  w = -bet * dummy
    void recurrence(Vector& v, Vector& w, prec bet) const {
      for (std::size_t j = 0; j < v.size(); ++j) {
        prec d = v[j];
        v[j]   = w[j] / bet;
        w[j]   = -bet * d;
      }
    }

    void axpy(prec a, const Vector& x, Vector& y) const {
      for (std::size_t j = 0; j < y.size(); ++j) y[j] += a * x[j];
    }

    void scale(Vector& v, prec s) const {
      for (auto& x : v) x *= s;
    }

    void add(Vector& acc, const Vector& t) const {
      for (std::size_t j = 0; j < acc.size(); ++j) acc[j] += t[j];
    }

    // c / c+ application producing (part of) the Lanczos start vector.
    void a_adag(int op, const Vector& invec, Vector& out,
                const Sector& next, bool annihilate) {
      _s.a_adag(op, invec, out, next, annihilate);
    }

  private:
    Storage& _s;
  };

}  // namespace edlib

#endif  // EDLIB_HOSTKERNEL_H
