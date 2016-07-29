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
    model(p),
    Eps(p["NSITES"], std::vector<double>(p["NSPINS"], 0.0)),
    t(p["NSITES"], std::vector<double>(p["NSITES"], 0.0)),
    U(p["NSITES"], 0.0),
    Ns(p["NSITES"]),
    ms(p["NSPINS"]),
    Ip(Ns * ms) {
    std::string input = p["INPUT_FILE"];
    alps::hdf5::archive input_data(input.c_str(), "r");
    input_data>>alps::make_pvp("BETA", _beta);
    input_data>>alps::make_pvp("hopping/values", t);
    input_data>>alps::make_pvp("interaction/values", U);
    input_data.close();
    for(int i = 0; i< Ns; ++i ) {
      // HARDCODED Half-filling
      Eps[i][0] = Eps[i][1] = -U[i]/2.0;
    }
    // TODO: move to input file and make site-dependent
    xmu = 0;
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
   * xmu - chemical potential
   * _beta - inverse temperature
   * Eps - level shift for each spin
   * t - hoppings
   * U - onsite Coulomb interaction
   * Ns - number of orbitals
   * ms - number of spins
   * Ip - maximum total number of electrons ms*Ns
   */

  Model model;

  double xmu;
  double _beta;
  std::vector<std::vector<double> > Eps;
  std::vector<std::vector<double> > t;
  std::vector<double> U;
  int Ns;
  int ms;
  int Ip;


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
