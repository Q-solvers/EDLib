//
// Created by iskakoff on 20/07/16.
//

#ifndef HUBBARD_CRSSTORAGE_H
#define HUBBARD_CRSSTORAGE_H


#include <vector>
#include <alps/params.hpp>
#include <iomanip>
#include "fortranbinding.h"
#include "Storage.h"

template<typename prec, class Symmetry>
class CRSStorage: public Storage<prec> {
  using Storage<prec>::n;
public:
  CRSStorage(alps::params & p, Symmetry& s):  Storage<prec>(p), _vind(0), _max_size(p["storage.MAX_SIZE"]), _max_dim(p["storage.MAX_DIM"]) {
    // init what you need from parameters
  };

  void reset(size_t sector_size){
    if(sector_size>_max_dim) {
      std::stringstream s;
      s<<"Current sector request more memory than allocated. Increase MAX_DIM parameter. Requested "<<sector_size<<", allocated "<<_max_dim<<".";
      throw std::runtime_error(s.str().c_str());
    }
    _vind = 0;
    row_ptr.assign(_max_dim+1, 0);
    col_ind.assign(_max_size, 0);
    values.assign(_max_size, prec(0.0));
    n()= 0;
  }

  void inline addDiagonal(const int &i, prec v) {
    row_ptr[i] = _vind;
    col_ind[_vind] = i + 1;
    values[_vind] = v;
    ++_vind;
    ++n();
  }

  /**
   * Add off-diagonal H(i,j) element
   */
  void inline addElement(const int &i, int j, prec t) {
    int findedstate = 0;
    bool hasstate = false;
    // check that there is no any data on the k state
    for (int iii = row_ptr[i]; iii <= _vind; ++iii) {
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
      ++_vind;
      if(_vind>_max_size) {
        std::stringstream s;
        s<<"Current sector request more memory than allocated. Increase MAX_SIZE parameter. Requested "<<_vind<<", allocated "<<_max_size<<".";
        throw std::runtime_error(s.str().c_str());
      }
    }
  }

  /**
   * Simple Compressed-Row-Storage Matrix-Vector product
   */
  virtual void av(prec* v, prec* w, int n, bool clear=true) override {
    for (int i = 0; i < n; ++i) {
      w[i] = clear? 0.0: w[i];
      for(int j = row_ptr[i]; j<row_ptr[i+1];++j){
        w[i] = w[i] + values[j] * v[col_ind[j]-1];
      }
    }
  }

  void endMatrix() {
    row_ptr[n()] = _vind;
  }

  void print() {
    std::cout<< std::setprecision(2)<<std::fixed;
    std::cout<<"[";
    for(int i = 0; i<n(); ++i) {
      std::cout<<"[";
      for(int j = 0; j<n(); ++j) {
        bool f = true;
        for(int k = row_ptr[i]; k<row_ptr[i+1]; ++k) {
          if((col_ind[k]-1) == j) {
            std::cout<<std::setw(6)<<values[k]<<(j==n()-1? "" : ", ");
            f = false;
          } /*else {
            std::cout<<"0.0 ";
          }*/
        }
        if(f) {
          std::cout<<std::setw(6)<<0.0<<(j==n()-1? "" : ", ");
        }
      }
      std::cout<<"]"<<(i==n()-1? "" : ", \n");
    }
    std::cout<<"]"<<std::endl;
  }

  virtual void zero_eigenapair() override {
    Storage<prec>::eigenvalues().resize(1);
    Storage<prec>::eigenvalues()[0] = values[0];
    Storage<prec>::eigenvectors().assign(1, std::vector<prec>(1, prec(1.0)));
  }


private:
  std::vector<prec> values;
  std::vector<int> row_ptr;
  std::vector<int> col_ind;
  size_t _max_size;
  size_t _max_dim;


  size_t _vind;
};


#endif //HUBBARD_CRSSTORAGE_H
