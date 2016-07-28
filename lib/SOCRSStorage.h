//
// Created by iskakoff on 28/07/16.
//

#ifndef HUBBARD_SOCRSSTORAGE_H
#define HUBBARD_SOCRSSTORAGE_H

#include <vector>
#include "Storage.h"

template <typename prec>
class SOCRSStorage : public Storage{
  SOCRSStorage(size_t max_size, size_t max_dim, alps::params & p):  Storage<prec>(max_dim, p),
                                                                    _vind(0), _max_size(max_size),
                                                                    _max_dim(max_dim) {/** init what you need from parameters*/};
public:
  virtual void av(prec *v, prec *w, int n) override {

  }

  void reset(){
    _vind = 0;
    col_ind.assign(_max_size, 0);
    dvalues.assign(_max_dim, prec(0.0));
    n()= 0;
  }

  void inline addDiagonal(const int &i, prec v) {
    dvalues[i] = v;
    n()++;
    _vind_start = _vind;
  }

  /**
   * Add off-diagonal H(i,j) element
   */
  void inline addElement(const int &i, int j, prec t) {
    char sign = t/std::abs(t);
    int findedstate = 0;
    bool hasstate = false;
    // check that there is no any data on the k state
    for (int iii = _vind_start; iii <= _vind; iii++) {
      if (col_ind[iii] == (j + 1)) {
        hasstate = true;
        findedstate = iii;
        throw std::logic_error("");
      }
    }
    // create new element in CRS arrays
    col_ind[_vind] = j + 1;
    dvalues[_vind] = -sign;
    _vind++;
  }

  void endMatrix() {
//    row_ptr[n()] = _vind;
  }

  virtual void zero_eigenapair() override {
    Storage<prec>::eigenvalues().resize(1);
    Storage<prec>::eigenvalues()[0] = dvalues[0];
    Storage<prec>::eigenvectors().assign(1, std::vector<prec>(1, prec(1.0)));
  }
private:
  //
  std::vector<prec> dvalues;
  std::vector<size_t> col_ind;
  std::vector<char> signs;

  size_t _max_size;
  size_t _max_dim;

  size_t _vind;
  size_t _vind_start;
};


#endif //HUBBARD_SOCRSSTORAGE_H
