//
// Created by iskakoff on 28/07/16.
//

#ifndef HUBBARD_SOCRSSTORAGE_H
#define HUBBARD_SOCRSSTORAGE_H

#include <vector>
#include <iomanip>
#include "Storage.h"

template <typename prec, class Model>
class SOCRSStorage : public Storage<prec>{
public:
  using Storage<prec>::n;
  template<class ModelType>
  SOCRSStorage(EDParams& p, ModelType &m):  Storage<prec>(p), _vind(0), _max_size(p["storage.MAX_SIZE"]),
                                   _max_dim(p["storage.MAX_DIM"]), model(m) {
    /** init what you need from parameters*/
    col_ind.assign(_max_size, 0);
    signs.assign(_max_size, 0);
    dvalues.assign(_max_dim, prec(0.0));
  };

  virtual void av(prec *v, prec *w, int n, bool clear=true) override {
    model.symmetry().init();
    _vind = 0;
    _vind_byte = 0;
    _vind_bit =0;
    // Iteration over rows.
    for(int i = 0; i<n; ++i){
      model.symmetry().next_state();
      long long nst = model.symmetry().state();
      // Diagonal contribution.
      w[i] = dvalues[i] * v[i] + (clear? 0.0: w[i]);
      // Offdiagonal contribution.
      // Iteration over columns(unordered).
      for(auto & state: model.states()) {
        int test = model.valid(state, nst);
        // If transition between states corresponding to row and column is possible, calculate the offdiagonal element.
        w[i] += test * state.value() * (1 - 2* ((signs[_vind_byte]>>_vind_bit)&1)) * v[col_ind[_vind]];
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
    _vind = 0;
    _vind_byte = 0;
    _vind_bit =0;
    n()= 0;
  }

  /**
   * Add diagonal H(i,i) element with value v.
   */
  void inline addDiagonal(const int &i, prec v) {
    dvalues[i] = v;
    ++n();
    _vind_start = _vind;
  }

  /**
   * Add off-diagonal H(i,j) element with value t (discarded here, restored in av) and Fermi sign.
   */
  void inline addElement(const int &i, int j, prec t, int sign) {
    if (i == j) {
      throw std::logic_error("Attempt to use addElement() to add diagonal element. Use addDiagonal() instead!");
    }
    // It is an error to add the element (i, j) twice.
    for (size_t iii = _vind_start; iii < _vind; iii++) {
      if (col_ind[iii] == j) {
        throw std::logic_error("Collision. Check a, adag, numState, ninv_value!");
      }
    }
    // Store sign in CRS-like array, one bit per sign.
    col_ind[_vind] = j;
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

  void print() {
    // See: av().
    // Each row of the matrix is first restored from the arrays.
    std::vector<prec> line(n(), prec(0.0));
    model.symmetry().init();
    _vind = 0;
    _vind_byte = 0;
    _vind_bit =0;
    std::cout<< std::setprecision(2)<<std::fixed;
    std::cout<<"[";
    for(int i = 0; i<n(); ++i){
      model.symmetry().next_state();
      long long nst = model.symmetry().state();
      std::fill(line.begin(), line.end(), prec(0.0));
      line[i] = dvalues[i];
      for(auto & state: model.states()) {
        int test = model.valid(state, nst);
        line[col_ind[_vind]] += test * state.value() * (1 - 2* ((signs[_vind_byte]>>_vind_bit)&1));
        _vind_bit+=test;
        _vind_byte+=_vind_bit/sizeof(char);
        _vind_bit%= sizeof(char);
        _vind+= test;
      }
      std::cout<<"[";
      for(int j = 0; j<n(); ++j){
        std::cout<<std::setw(6)<<line[j]<<(j==n()-1? "" : ", ");
      }
      std::cout<<"]"<<(i==n()-1? "" : ", \n");
    }
    std::cout<<"]"<<std::endl;
  }

  virtual void zero_eigenapair() override {
    Storage<prec>::eigenvalues().resize(1);
    Storage<prec>::eigenvalues()[0] = dvalues[0];
    Storage<prec>::eigenvectors().assign(1, std::vector<prec>(1, prec(1.0)));
  }
private:
  // Internal storage structure
  std::vector<prec> dvalues;
  std::vector<size_t> col_ind;
  std::vector<char> signs;

  // the maximum sizes of all the objects
  size_t _max_size;
  size_t _max_dim;

  // internal indicies
  size_t _vind;
  // start of current row, used for checks
  size_t _vind_start;
  // bit and byte of the bitmap corresponding to CRS index
  size_t _vind_bit;
  size_t _vind_byte;

  // Hubbard model parameters
  Model &model;

};


#endif //HUBBARD_SOCRSSTORAGE_H
