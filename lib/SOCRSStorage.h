//
// Created by iskakoff on 28/07/16.
//

#ifndef HUBBARD_SOCRSSTORAGE_H
#define HUBBARD_SOCRSSTORAGE_H

#include <vector>
#include "Storage.h"

template <typename prec, class Symmetry, class Model>
class SOCRSStorage : public Storage<prec>{
public:
  using Storage<prec>::n;
  SOCRSStorage(alps::params & p):  Storage<prec>(p), _vind(0), _max_size(p["storage.MAX_SIZE"]),
                                   _max_dim(p["storage.MAX_DIM"]), symmetry(p), model(p) {
    /** init what you need from parameters*/
    col_ind.assign(_max_size, 0);
    signs.assign(_max_size, 0);
    dvalues.assign(_max_dim, prec(0.0));
  };

  virtual void av(prec *v, prec *w, int n) override {
    long long k;
    int isign;
    symmetry.init();
    _vind = 0;
    _vind_byte = 0;
    _vind_bit =0;
    for(int i = 0; i<n; ++i){
      symmetry.next_state();
      long long nst = symmetry.state();
      w[i] = dvalues[i] * v[i];
      for(auto & state: model.states()) {
        int test = model.valid(state, nst);
        w[i] += test * state.value() * (1- 2* ((signs[_vind_byte]>>_vind_bit)&1)) * v[col_ind[_vind] - 1];
        _vind_bit+=test;
        _vind_byte+=_vind_bit/sizeof(char);
        _vind_bit%= sizeof(char);
        _vind+= test;
      }
    }
  }

  void reset(size_t sector_size){
    if(sector_size>_max_dim) {
      std::stringstream s;
      s<<"Current sector request more memory than allocated. Increase MAX_DIM parameter. Requested "<<sector_size<<", allocated "<<_max_dim<<".";
      throw std::runtime_error(s.str().c_str());
    }
    symmetry.next_sector();
    _vind = 0;
    _vind_byte = 0;
    _vind_bit =0;
    n()= 0;
  }

  void inline addDiagonal(const int &i, prec v) {
    dvalues[i] = v;
    ++n();
    _vind_start = _vind;
  }

  /**
   * Add off-diagonal H(i,j) element
   */
  void inline addElement(const int &i, int j, prec t) {
    char sign = (prec(0) < t) - (t < prec(0));
    int findedstate = 0;
    bool hasstate = false;
    // check that there is no any data on the k state
    for (int iii = _vind_start; iii < _vind; iii++) {
      if (col_ind[iii] == (j + 1)) {
        throw std::logic_error("Collision. Check a, adag, numState, ninv_value!");
      }
    }
    // create new element in CRS arrays
    col_ind[_vind] = j + 1;
    signs[_vind_byte] &= ~(1ll << _vind_bit);
    signs[_vind_byte] |= sign<0 ? 1ll<<_vind_bit:0;
    ++_vind_bit;
    ++_vind;
    _vind_byte += _vind_bit/sizeof(char);
    _vind_bit%=sizeof(char);
    if(_vind>_max_size || _vind_byte>_max_size) {
      std::stringstream s;
      s<<"Current sector request more memory than allocated. Increase MAX_SIZE parameter. Requested "<<_vind<<", allocated "<<_max_size<<".";
      throw std::runtime_error(s.str().c_str());
    }
  }

  void endMatrix() {
    // Nothing to do
  }

  virtual void zero_eigenapair() override {
    Storage<prec>::eigenvalues().resize(1);
    Storage<prec>::eigenvalues()[0] = dvalues[0];
    Storage<prec>::eigenvectors().assign(1, std::vector<prec>(1, prec(1.0)));
  }
private:
  // System symmetry to walk through all states
  Symmetry symmetry;
  // Internal storage structure
  std::vector<prec> dvalues;
  std::vector<size_t> col_ind;
  std::vector<char> signs;

  // the maximum sizes of all the objects
  size_t _max_size;
  size_t _max_dim;

  // internal indicies
  size_t _vind;
  size_t _vind_start;
  size_t _vind_bit;
  size_t _vind_byte;

  // Hubbard model parameters
  Model model;

};


#endif //HUBBARD_SOCRSSTORAGE_H
