#ifndef EDLIB_GF_H
#define EDLIB_GF_H

#include <array>
#include <complex>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace edlib {

  namespace detail {
    template <int N>
    std::size_t linear_index(const std::array<int, N>& shape,
                             const std::array<int, N>& idx) {
      std::size_t off = 0;
      for (int d = 0; d < N; ++d) {
        off = off * static_cast<std::size_t>(shape[d])
              + static_cast<std::size_t>(idx[d]);
      }
      return off;
    }

    template <int N>
    std::size_t total_size(const std::array<int, N>& shape) {
      std::size_t s = 1;
      for (int d = 0; d < N; ++d) s *= static_cast<std::size_t>(shape[d]);
      return s;
    }
  }

  /**
   * Multi-index Green's function container.
   *
   * Row-major contiguous storage in `data`. The leading axis is the frequency
   * axis (length matches the associated mesh). Remaining axes are integer
   * index meshes (orbitals, spins, real-space, etc.).
   */
  template <class T, int N>
  class Gf {
    static_assert(N >= 1, "Gf rank must be >= 1");

  public:
    using value_type = T;
    static constexpr int rank = N;

    Gf() { _shape.fill(0); }

    explicit Gf(const std::array<int, N>& shape) : _shape(shape) {
      _data.assign(detail::total_size<N>(_shape), T{});
    }

    const std::array<int, N>& shape() const { return _shape; }
    int                       shape(int d) const { return _shape[d]; }
    const std::vector<T>&     data() const { return _data; }
    std::vector<T>&           data()       { return _data; }

    template <class... Idx>
    T& operator()(Idx... ix) {
      static_assert(sizeof...(Idx) == N, "operator(): wrong number of indices");
      return _data[detail::linear_index<N>(_shape, {static_cast<int>(ix)...})];
    }

    template <class... Idx>
    const T& operator()(Idx... ix) const {
      static_assert(sizeof...(Idx) == N, "operator(): wrong number of indices");
      return _data[detail::linear_index<N>(_shape, {static_cast<int>(ix)...})];
    }

    Gf& operator+=(const Gf& o) {
      check_shape(o);
      for (std::size_t i = 0; i < _data.size(); ++i) _data[i] += o._data[i];
      return *this;
    }
    Gf& operator-=(const Gf& o) {
      check_shape(o);
      for (std::size_t i = 0; i < _data.size(); ++i) _data[i] -= o._data[i];
      return *this;
    }
    Gf& operator*=(const T& s) {
      for (auto& v : _data) v *= s;
      return *this;
    }
    Gf& operator/=(const T& s) {
      for (auto& v : _data) v /= s;
      return *this;
    }

    friend Gf operator+(Gf a, const Gf& b) { a += b; return a; }
    friend Gf operator-(Gf a, const Gf& b) { a -= b; return a; }

  private:
    void check_shape(const Gf& o) const {
      if (_shape != o._shape) {
        throw std::invalid_argument("Gf: shape mismatch in element-wise op");
      }
    }

    std::array<int, N> _shape;
    std::vector<T>     _data;
  };

  using GF2 = Gf<std::complex<double>, 2>;
  using GF3 = Gf<std::complex<double>, 3>;
  using GF4 = Gf<std::complex<double>, 4>;

}

#endif
