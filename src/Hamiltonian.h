//
// Created by iskakoff on 19/07/16.
//

#ifndef HUBBARD_HAMILTONIAN_H
#define HUBBARD_HAMILTONIAN_H

#include <vector>
#include <alps/params.hpp>
#include <type_traits>

#include <fstream>
#include "Symmetry.h"
#include "EigenPair.h"

template<typename prec, class Symmetry, class Storage, class Model>
class Hamiltonian {
  static_assert(std::is_base_of<Symmetry, Symmetry>::value, "Symmetry should extend base Symmetry class");
public:
  /*
   * Allocate space for Hamiltonian matrix
   * \param [in] max_size - the maximum size of values array
   * \param [in] max_dim - the maximum dimension of Hamiltonian matrix
   * \param [in] p - alps::parameters
   */
  Hamiltonian(alps::params& p) :
    _symmetry(p),
    _storage(p, _symmetry),
    model(p) {};

  /**
   * fill current sector
   */
  void fill() {
    _symmetry.init();
    _storage.reset(_symmetry.sector().size());
    int i =0;
    long long k = 0;
    int isign = 0;
    while (_symmetry.next_state()) {
      long long nst = _symmetry.state();
      // Compute diagonal element for current i state

      _storage.addDiagonal(i, model.diagonal(nst));
      // non-diagonal terms calculation
      for(auto & state: model.states()) {
        if(model.valid(state, nst)) {
          model.set(state, nst, k, isign);
          hopping(i, nst, k, state.value(), isign);
        }
      }
      i++;
    }
    // additional steps after all data
    _storage.endMatrix();
  }
  /**
   * perform Hamiltonian diagonalization
   * result will be stored in evals and evecs
   */
  void diag() {
    while(_symmetry.next_sector()) {
      fill();
      /**
       * perform ARPACK call
       */
      int info = _storage.diag();
      if(info != 0) {

      } else {
        const std::vector<prec>& evals = _storage.eigenvalues();
        const std::vector<std::vector<prec> >& evecs = _storage.eigenvectors();
        for(int i = 0; i<evals.size(); ++i) {
          _eigenpairs.push_back(EigenPair<prec, typename Symmetry::Sector>(evals[i], evecs[i], _symmetry.sector()));
        }
      }
    }
    std::sort(_eigenpairs.begin(), _eigenpairs.end());
    std::cout<<"Here is the list of eigenvalues:"<<std::endl;
    for(auto& eigenpair : _eigenpairs) {
      std::cout<<eigenpair.eigenvalue()<<" ";
      eigenpair.sector().print();
      std::cout<<std::endl;
    }
  }

  const Symmetry & symmetry() const {
    return _symmetry;
  }

  Symmetry & symmetry() {
    return _symmetry;
  }

  Storage& storage() {
    return _storage;
  }

  const std::vector<EigenPair<prec, typename Symmetry::Sector> > & eigenpairs() const {
    return _eigenpairs;
  };
private:
  // CSR format Hamiltonian matrix storage
  Storage _storage;
  Symmetry _symmetry;

  // Eigen-pairs
  std::vector<EigenPair<prec, typename Symmetry::Sector> > _eigenpairs;

  /**
   * Model to diagonalize
   */
  Model model;


  /**
   * \param i - current Hamiltonian matrix line
   * \param nst - current state
   * \param k - state to hop between
   * \param v - hopping value
   * \param sector - current conservation law sector
   */
  void inline hopping(const int& i, const long long& nst, const long long& k, const prec &v, const int &sign) {
    int k_index = _symmetry.index(k);
    _storage.addElement(i, k_index, v, sign);
  }

};

#endif //HUBBARD_HAMILTONIAN_H
