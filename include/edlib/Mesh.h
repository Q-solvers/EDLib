#ifndef EDLIB_MESH_H
#define EDLIB_MESH_H

#include <cmath>
#include <stdexcept>
#include <vector>

namespace edlib {

  enum class Statistics { Fermionic, Bosonic };

  class MatsubaraMesh {
  public:
    MatsubaraMesh() = default;
    MatsubaraMesh(double beta, int n, Statistics stat)
        : _beta(beta), _n(n), _stat(stat) {}

    double beta()             const { return _beta; }
    int    extent()           const { return _n; }
    Statistics statistics()   const { return _stat; }

    std::vector<double> points() const {
      if (_beta <= 0.0) {
        throw std::invalid_argument("MatsubaraMesh: beta must be positive");
      }
      std::vector<double> w(_n);
      const double pi = 3.14159265358979323846;
      const int shift = (_stat == Statistics::Fermionic) ? 1 : 0;
      for (int i = 0; i < _n; ++i) {
        w[i] = pi * (2 * i + shift) / _beta;
      }
      return w;
    }

  private:
    double _beta = 0.0;
    int    _n = 0;
    Statistics _stat = Statistics::Fermionic;
  };

  class RealFreqMesh {
  public:
    RealFreqMesh() = default;
    RealFreqMesh(double emin, double emax, int n)
        : _emin(emin), _emax(emax), _n(n) {}

    double emin()   const { return _emin; }
    double emax()   const { return _emax; }
    int    extent() const { return _n; }

    std::vector<double> points() const {
      std::vector<double> w(_n);
      if (_n == 1) {
        w[0] = _emin;
        return w;
      }
      const double step = (_emax - _emin) / static_cast<double>(_n - 1);
      for (int i = 0; i < _n; ++i) {
        w[i] = _emin + step * i;
      }
      return w;
    }

  private:
    double _emin = 0.0;
    double _emax = 0.0;
    int    _n = 0;
  };

}

#endif
