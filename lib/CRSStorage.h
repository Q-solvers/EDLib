//
// Created by iskakoff on 20/07/16.
//

#ifndef HUBBARD_CRSSTORAGE_H
#define HUBBARD_CRSSTORAGE_H


#include <vector>
#include <alps/params.hpp>

template<typename prec>
class CRSStorage {
public:
  CRSStorage(size_t max_size, size_t max_dim, alps::params & p): _vind(0), _max_size(max_size), _max_dim(max_dim) {
    // init what you need from parameters
  };

  void reset(){
    _vind = 0;
    row_ptr.assign(_max_dim+1, 0);
    col_ind.assign(_max_size, 0);
    values.assign(_max_size, prec(0.0));
  }

  void inline addDiagonal(const int &i, prec v) {
    row_ptr[i] = _vind;
    col_ind[_vind] = i + 1;
    values[_vind] = v;
    _vind++;
  }

  /**
   * Add off-diagonal H(i,j) element
   */
  void inline addElement(const int &i, int j, prec t) {
    int findedstate = 0;
    bool hasstate = false;
    // check that there is no any data on the k state
    for (int iii = row_ptr[i]; iii <= _vind; iii++) {
      if (col_ind[iii] == (j + 1)) {
        hasstate = true;
        findedstate = iii;
        break;
      }
    }
    // if there is data add value
    if (hasstate) {
      values[findedstate] += t;
    } else {
      // create new element in CRS arrays
      col_ind[_vind] = j + 1;
      values[_vind] = t;
      _vind++;
    }
  }

  void av(){

  }

private:
  std::vector<prec> values;
  std::vector<int> row_ptr;
  std::vector<int> col_ind;
  size_t _max_size;
  size_t _max_dim;

  int _vind;
};


#endif //HUBBARD_CRSSTORAGE_H
