//
// Created by iskakoff on 20/04/17.
//

#ifndef HUBBARD_CSRMATRIX_H
#define HUBBARD_CSRMATRIX_H

namespace EDLib {
  namespace Storage {
    /**
     * @brief CSRMatrix class
     *
     * @author iskakoff
     */
    /**
       * Simple CRS matrix class. This class is used to store hopping matrices and off-diagonal interactions
       * @tparam p
       */
    template<typename prec>
    class CRSMatrix {
    public:
      CRSMatrix() {
      }

      /**
       * init matrix arrays
       * @param N -- leading dimension
       * @param nnzl -- average number of non-zero elements per line
       */
      void init(size_t N, size_t nnzl = 100) {
        _nnz = N * nnzl;
        _values.assign(_nnz, prec(0));
        _col_ind.assign(_nnz, 0);
        _row_ptr.assign(N + 1, 0);
        _vind = 0;
      }

      /**
       * Add off-diagonal matrix element at the position (i,j)
       *
       * @param i - row number
       * @param j - column number
       * @param t - value
       * @param sign - fermionic sign
       */
      void inline addElement(int i, int j, prec t, int sign) {
        /// flag that we already have data for (i,j)
        bool hasstate = false;
        /// index of (i,j) element in spare storage
        size_t foundstate = 0;
        if (std::abs(t) == 0) {
          return;
        }
        /// check that there is no any data on the k state
        /// In case of off-diagonal interaction there can be multiple possible transition from i-state to j-state
        for (int iii = _row_ptr[i]; iii < _vind; ++iii) {
          if (_col_ind[iii] == j) {
            hasstate = true;
            foundstate = iii;
          }
        }
        if(hasstate) {
          /// update existing value
          _values[foundstate] += sign * t;
        }else {
          /// create new element in CRS arrays
          _col_ind[_vind] = j;
          _values[_vind] = sign * t;
          ++_vind;
          /// check that we have exceed the upper bound
          if(_vind == _nnz) {
            /// resize storage
            _nnz *= 2;
            _values.resize(_nnz);
            _col_ind.resize(_nnz);
          }
        }
      }

      /// some of the interaction terms can compensate each other
      /// in this case we need to remove zero elements from storage to reduce required memory and communications
      void inline compress(int i) {
        int shift = 0;
        for (int iii = _row_ptr[i]; iii < _vind; ++iii){
          if(std::abs(_values[iii])<1e-15) {
            for(int kkk = iii; kkk<_vind-1; ++kkk) {
              _values[kkk] = _values[kkk+1];
              _col_ind[kkk] = _col_ind[kkk+1];
            }
            --_vind;
          }
        }
      }

      void inline endLine(int i) {
        compress(i);
        _row_ptr[i + 1] = _vind;
      }

      std::vector < int > &row_ptr() {
        return _row_ptr;
      };

      std::vector < int > &col_ind() {
        return _col_ind;
      }

      std::vector < prec > &values() {
        return _values;
      }

    private:
      /// matrix values
      std::vector < prec > _values;
      /// pointer to a row
      std::vector < int > _row_ptr;
      /// column indices
      std::vector < int > _col_ind;
      /// internal index of non-zero values
      int _vind;
      /// number of non-zero elements allocated in memory
      size_t _nnz;
    };
  }
}

#endif //HUBBARD_CSRMATRIX_H
