//
// Created by iskakoff on 28/07/16.
//

#ifndef HUBBARD_SOCRSSTORAGE_H
#define HUBBARD_SOCRSSTORAGE_H

#include <vector>
#include <Symmetry.h>
#include "Storage.h"

template <typename prec, class Symmetry=Symmetry>
class SOCRSStorage : public Storage<prec>{
public:
  using Storage<prec>::n;
  SOCRSStorage(size_t max_size, size_t max_dim, alps::params & p):  Storage<prec>(max_dim, p),
                                                                    _vind(0), _max_size(max_size),
                                                                    _max_dim(max_dim), symmetry(p),
                                                                    _Ns(p["NSITES"]), _Ip(2*_Ns), _ms(p["NSPINS"]) {
    /** init what you need from parameters*/
    std::string input = p["INPUT_FILE"];
    alps::hdf5::archive input_data(input.c_str(), "r");
    input_data>>alps::make_pvp("hopping/values", t);
    input_data.close();
  };

  virtual void av(prec *v, prec *w, int n) override {
    symmetry.init();
    _vind = 0;
    _vind_byte = 0;
    _vind_bit =0;
    for(int i = 0; i<n; ++i){
      symmetry.next_state();
      long long nst = symmetry.state();
      w[i] = dvalues[i] * v[i];
      for(int ii = 0; ii< _Ns; ++ii) {
        for(int jj = 0; jj< _Ns; ++jj) {
          if(ii!=jj && std::abs(t[ii][jj])) {
            for(int ispin = 0; ispin<_ms ; ++ispin) {
              int test = checkState(nst, ii + ispin * _Ns) * (1 - checkState(nst, jj + ispin * _Ns));
              if(test) {
                w[i] += test * t[ii][jj] * (1- 2* ((signs[_vind_byte]>>_vind_bit)&1)) * v[col_ind[_vind] - 1];
              }
              _vind_bit+=test;
              _vind_byte+=_vind_bit/sizeof(char);
              _vind_bit%= sizeof(char);
              _vind+= test;
            }
          }
        }
      }
    }
  }

  void reset(){
    symmetry.next_sector();
    _vind = 0;
    _vind_byte = 0;
    _vind_bit =0;
    col_ind.assign(_max_size, 0);
    signs.assign(_max_size, 0);
    dvalues.assign(_max_dim, prec(0.0));
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
    for (int iii = _vind_start; iii <= _vind; iii++) {
      if (col_ind[iii] == (j + 1)) {
        throw std::logic_error("Collision. Check a, adag, numState, ninv_value!");
      }
    }
    // create new element in CRS arrays
//    vind_bit = vind_bit + 1;
//    vind_byte = vind_byte + vind_bit/8;
//    vind_bit = modulo(vind_bit, 8);
    col_ind[_vind] = j + 1;
    signs[_vind_byte] |= sign<0 ? 1ll<<_vind_bit:0;
    ++_vind_bit;
    ++_vind;
    _vind_byte += _vind_bit/sizeof(char);
    _vind_bit%=sizeof(char);
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
  Symmetry symmetry;
  //
  std::vector<prec> dvalues;
  std::vector<size_t> col_ind;
  std::vector<char> signs;

  size_t _max_size;
  size_t _max_dim;

  size_t _vind;
  size_t _vind_start;
  size_t _vind_bit;
  size_t _vind_byte;

  // Hubbard model parameters
  int _Ns;
  int _Ip;
  int _ms;

  std::vector<std::vector<prec> > t;

  // TODO: think about dependency of Hamiltonian
  int inline checkState(const long long& nst, const int& im) {
    return (int) ((nst & (1ll << (_Ip - 1 - im))) >> (_Ip - 1 - im));
  }
};


#endif //HUBBARD_SOCRSSTORAGE_H
