#ifndef EDLIB_CRSMATRIX_H
#define EDLIB_CRSMATRIX_H

#include <cmath>
#include <cstddef>
#include <vector>

namespace edlib {

  /**
   * Simple compressed row sparse matrix used to store hopping matrices and
   * off-diagonal interaction terms.
   */
  template <class Prec>
  class CRSMatrix {
  public:
    CRSMatrix() = default;

    void init(std::size_t N, std::size_t nnzl = 100) {
      _nnz = N * nnzl;
      _values.assign(_nnz, Prec(0));
      _col_ind.assign(_nnz, 0);
      _row_ptr.assign(N + 1, 0);
      _vind = 0;
    }

    inline void addElement(int i, int j, Prec t, int sign) {
      if (std::abs(t) == Prec(0)) return;

      bool        hasstate = false;
      std::size_t foundstate = 0;
      for (int k = _row_ptr[i]; k < _vind; ++k) {
        if (_col_ind[k] == j) { hasstate = true; foundstate = k; }
      }
      if (hasstate) {
        _values[foundstate] += static_cast<Prec>(sign) * t;
      } else {
        _col_ind[_vind] = j;
        _values[_vind]  = static_cast<Prec>(sign) * t;
        ++_vind;
        if (static_cast<std::size_t>(_vind) == _nnz) {
          _nnz *= 2;
          _values.resize(_nnz);
          _col_ind.resize(_nnz);
        }
      }
    }

    inline void compress(int i) {
      for (int k = _row_ptr[i]; k < _vind; ++k) {
        if (std::abs(_values[k]) < 1e-15) {
          for (int m = k; m < _vind - 1; ++m) {
            _values[m]  = _values[m + 1];
            _col_ind[m] = _col_ind[m + 1];
          }
          --_vind;
        }
      }
    }

    inline void endLine(int i) {
      compress(i);
      _row_ptr[i + 1] = _vind;
    }

    std::vector<int>&  row_ptr() { return _row_ptr; }
    std::vector<int>&  col_ind() { return _col_ind; }
    std::vector<Prec>& values()  { return _values; }

  private:
    std::vector<Prec> _values;
    std::vector<int>  _row_ptr;
    std::vector<int>  _col_ind;
    int               _vind = 0;
    std::size_t       _nnz  = 0;
  };

}

#endif
