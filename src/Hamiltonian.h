//
// Created by iskakoff on 19/07/16.
//

#ifndef HUBBARD_HAMILTONIAN_H
#define HUBBARD_HAMILTONIAN_H

#include <vector>
#include <type_traits>

#include <fstream>
#include "Symmetry.h"
#include "EigenPair.h"
#include "EDParams.h"
#include "HubbardModel.h"
#include "CRSStorage.h"
#include "SOCRSStorage.h"

template<typename prec, class Storage, class Model>
class Hamiltonian {
public:
  typedef Model ModelType;
  /*
   * Allocate space for Hamiltonian matrix
   * \param [in] max_size - the maximum size of values array
   * \param [in] max_dim - the maximum dimension of Hamiltonian matrix
   * \param [in] p - alps::parameters
   */
  Hamiltonian(EDParams& p) :
    _model(p),
    _storage(p, _model){};

  /**
   * fill current sector
   */
  void fill() {
    _storage.fill();
  }
  /**
   * perform Hamiltonian diagonalization
   * result will be stored in evals and evecs
   */
  void diag() {
    while(_model.symmetry().next_sector()) {
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
          _eigenpairs.push_back(EigenPair<prec, typename Model::Sector>(evals[i], evecs[i], _model.symmetry().sector()));
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

  Storage& storage() {
    return _storage;
  }

  const std::vector<EigenPair<prec, typename Model::Sector> > & eigenpairs() const {
    return _eigenpairs;
  };

  Model &model(){
    return _model;
  }
private:
  // CSR format Hamiltonian matrix storage
  Storage _storage;

  // Eigen-pairs
  std::vector<EigenPair<prec, typename Model::Sector> > _eigenpairs;

  /**
   * Model to diagonalize
   */
  Model _model;

};

typedef Hamiltonian<double, CRSStorage<double, HubbardModel<double> > , HubbardModel<double> > CSRHubbardHamiltonian;
typedef Hamiltonian<double, SOCRSStorage<double, HubbardModel<double> > , HubbardModel<double> > SOCSRHubbardHamiltonian;

typedef Hamiltonian<float, CRSStorage<float, HubbardModel<float> > , HubbardModel<float> > CSRHubbardHamiltonian_float;
typedef Hamiltonian<float, SOCRSStorage<float, HubbardModel<float> > , HubbardModel<float> > SOCSRHubbardHamiltonian_float;

#endif //HUBBARD_HAMILTONIAN_H
