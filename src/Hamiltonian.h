//
// Created by iskakoff on 19/07/16.
//

#ifndef HUBBARD_HAMILTONIAN_H
#define HUBBARD_HAMILTONIAN_H

#include <vector>
#include <alps/params.hpp>
#include <type_traits>

#include <CRSStorage.h>
#include <fstream>
#include "Symmetry.h"
#include "EigenPair.h"
#include "HubbardModel.h"

template<typename prec, class Symmetry=Symmetry, class Storage=CRSStorage<double>, class Model=HubbardModel<double> >
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
    storage(p),
    symmetry(p),
    model(p) {
    std::string input = p["INPUT_FILE"];
    alps::hdf5::archive input_data(input.c_str(), "r");
    input_data.close();
  };

  /**
   * fill current sector
   */
  void fill() {
    symmetry.init();
    storage.reset(symmetry.sector().size());
    int i =0;
    long long k;
    int isign;
    while (symmetry.next_state()) {
      long long nst = symmetry.state();
      // Compute diagonal element for current i state

      storage.addDiagonal(i, model.diagonal(nst));
      // non-diagonal terms calculation
      for(auto & state: model.states()) {
        if(model.valid(state, nst)) {
          model.set(state, nst, k, isign);
          hopping(i, nst, k, isign * state.value());
        }
      }
      i++;
    }
    // additional steps after all data
    storage.endMatrix();
  }
  /**
   * perform Hamiltonian diagonalization
   * result will be stored in evals and evecs
   */
  void diag() {
    while(symmetry.next_sector()) {
      fill();
      /**
       * perform ARPACK call
       */
      int info = storage.diag();
      if(info != 0) {

      } else {
        const std::vector<prec>& evals = storage.eigenvalues();
        const std::vector<std::vector<prec> >& evecs = storage.eigenvectors();
        for(int i = 0; i<evals.size(); ++i) {
//          std::vector<prec> evec(evecs[i]);
          eigenpairs.push_back(EigenPair<prec, typename Symmetry::Sector>(evals[i], evecs[i], symmetry.sector()));
        }
      }
    }
    std::sort(eigenpairs.begin(), eigenpairs.end());
    std::cout<<"Here is the list of eigenvalues:"<<std::endl;
    for(auto& eigenpair : eigenpairs) {
      std::cout<<eigenpair.eigenvalue()<<" ";
      eigenpair.sector().print();
      std::cout<<std::endl;
    }
  }

private:
  // CSR format Hamiltonian matrix storage
  Storage storage;
  Symmetry symmetry;

  // Eigen-pairs
  std::vector<EigenPair<prec, typename Symmetry::Sector> > eigenpairs;

  /**
   * Model to diagonalize
   */
  Model model;


  /**
   * \param i - current index
   * \param nst - curent state
   * \param k - state to hop between
   * \param v - hopping value
   * \param sector - current conservation law sector
   */
  void inline hopping(const int& i, const long long& nst, const long long& k, const prec &v) {
    int k_index = symmetry.index(k) - 1;
    storage.addElement(i, k_index, v);
  }

};

#endif //HUBBARD_HAMILTONIAN_H
