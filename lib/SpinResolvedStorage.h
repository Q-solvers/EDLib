//
// Created by iskakoff on 23/08/16.
//

#ifndef HUBBARD_SPINRESOLVEDSTORAGE_H
#define HUBBARD_SPINRESOLVEDSTORAGE_H

#include <type_traits>
#include <bitset>

#include "Storage.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <SzSymmetry.h>
#include <NSymmetry.h>


template<typename prec, class Model>
class SpinResolvedStorage : public Storage<prec> {
//  static_assert(std::is_base_of<SzSymmetry,typename std::result_of<decltype(&Model::symmetry)(Model)>::type>::value, "Wrong symmetry");
public:
  using Storage<prec>::n;

  template<typename p>
  class CRSMatrix {
  public:
    CRSMatrix() {
    }

    void init(int N) {
      _values.assign(N*N, p(0));
      _col_ind.assign(N*N, 0);
      _row_ptr.assign(N+1, 0);
      _vind = 0;
    }

    void inline addElement(int i, int j, prec t, int sign) {
      if (i == j) {
        throw std::logic_error("Attempt to use addElement() to add diagonal element. Use addDiagonal() instead!");
      }
      // check that there is no any data on the k state
      for (int iii = _row_ptr[i]; iii < _vind; ++iii) {
        if (_col_ind[iii] == j) {
          throw std::logic_error("Collision. Check a, adag, numState, ninv_value!");
        }
      }
      // create new element in CRS arrays
      _col_ind[_vind] = j;
      _values[_vind] = sign * t;
      ++_vind;
    }

    void inline endLine(int i) {
      _row_ptr[i+1] = _vind;
    }
    std::vector<int> & row_ptr() {
      return _row_ptr;
    };
    std::vector<int>& col_ind() {
      return _col_ind;
    }
    std::vector<prec> & values() {
      return _values;
    }

  private:
    std::vector<prec> _values;
    std::vector<int> _row_ptr;
    std::vector<int> _col_ind;
    int _vind;
  };

  typedef CRSMatrix<prec> Matrix;

  SpinResolvedStorage(EDParams &p, Model &m) : Storage<prec>(p), _model(m), _interaction_size(int(p["spinstorage.ORBITAL_NUMBER"]) * int(p["NSPINS"])),
                                                           _Ns(p["NSITES"]), _ms(p["NSPINS"]), _up_symmetry(int(p["NSITES"])), _down_symmetry(int(p["NSITES"])) {
    H_loc.resize(_interaction_size);
  }

  virtual void zero_eigenapair() override {

  }

  virtual void av(prec *v, prec *w, int n, bool clear = true) override {
    // Iteration over rows.
    for(int i = 0; i<n; ++i){
      // Diagonal contribution.
      w[i] = _diagonal[i] * v[i] + (clear? 0.0: w[i]);
    }
    std::vector<prec> ww_u(_up_symmetry.sector().size(), 0.0);
    std::vector<prec> ww_d(_down_symmetry.sector().size(), 0.0);
    // Offdiagonal contribution.
    // Iteration over columns(unordered).
    for(int k = 0; k<_up_symmetry.sector().size(); ++k) {
      for(int i = 0; i<_down_symmetry.sector().size(); ++i) {
        for (int j = H_down.row_ptr()[k]; j < H_down.row_ptr()[k + 1]; ++j) {
          w[i + k*_up_symmetry.sector().size()] += H_down.values()[j] * v[H_down.col_ind()[j] + k*_up_symmetry.sector().size()];
        }
      }
    }
    for(int k = 0; k<_down_symmetry.sector().size(); ++k) {
      for(int i = 0; i<_up_symmetry.sector().size(); ++i) {
        for (int j = H_up.row_ptr()[k]; j < H_up.row_ptr()[k + 1]; ++j) {
          w[i + k*_down_symmetry.sector().size()] += H_up.values()[j] * v[H_up.col_ind()[j] + k*_down_symmetry.sector().size()];
        }
      }
    }
  }

  void reset() {
    const SzSymmetry &symmetry = static_cast<SzSymmetry>(_model.symmetry());
    const SzSymmetry::Sector &sector = symmetry.sector();
    _up_symmetry.set_sector(NSymmetry::Sector(sector.nup(), symmetry.comb().c_n_k(_Ns, sector.nup())));
    _down_symmetry.set_sector(NSymmetry::Sector(sector.ndown(), symmetry.comb().c_n_k(_Ns, sector.ndown())));
    H_up.init(_up_symmetry.sector().size());
    H_down.init(_down_symmetry.sector().size());
    _diagonal.assign(sector.size(), prec(0.0));
  }

  void fill() {
    _model.symmetry().init();
    reset();
    // fill off-diagonal matrix for each spin
    fill_spin(_up_symmetry, _Ns, H_up);
    int sign = (_up_symmetry.sector().n() %2) == 0 ? 1:-1;
    fill_spin(_down_symmetry, 0, H_down, sign);
    // fill diagonal;
    int i =0;
    while (_model.symmetry().next_state()) {
      long long nst = _model.symmetry().state();
      _diagonal[i] = _model.diagonal(nst);
      ++i;
    }
    n() = _model.symmetry().sector().size();
  }

  void fill_spin(NSymmetry &spin_symmetry, int shift, Matrix &spin_matrix, int sign = 1) {
    long long k = 0;
    int isign = 0;
    int i =0;
    while (spin_symmetry.next_state()) {
      long long nst = spin_symmetry.state();
      for(auto & state: _model.states()) {
        if(_model.valid(state, nst << shift)) {
          _model.set(state, nst << shift, k, isign);
          int j = spin_symmetry.index(k >> shift);
//          std::cout<<"st:"<<std::bitset<4>(nst).to_string()<<" -> "<<std::bitset<4>(k>>shift).to_string()<<" "<<i<<" "<<j<<std::endl;
          spin_matrix.addElement(i, j, state.value(), isign * sign);
        }
      }
      spin_matrix.endLine(i);
      ++i;
    }
  }

private:
  std::vector<Eigen::Matrix<prec, -1, -1> > H_loc;
  Model & _model;
  Matrix H_up;
  Matrix H_down;

  std::vector<prec> _diagonal;

  NSymmetry _up_symmetry;
  NSymmetry _down_symmetry;


  int _interaction_size;
  int _Ns;
  int _ms;
};


#endif //HUBBARD_SPINRESOLVEDSTORAGE_H
