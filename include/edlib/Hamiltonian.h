//
// Created by iskakoff on 19/07/16.
//

#ifndef HUBBARD_HAMILTONIAN_H
#define HUBBARD_HAMILTONIAN_H

#include <set>
#include <type_traits>

#include <iomanip>
#include <fstream>
#include "SpinResolvedStorage.h"
#include "Symmetry.h"
#include "EigenPair.h"
#include "HubbardModel.h"
#include "CRSStorage.h"
#include "SOCRSStorage.h"
#include "SingleImpurityAndersonModel.h"

namespace EDLib {
  template<class Storage, class Model>
  class Hamiltonian {
  public:
    typedef Model ModelType;
    typedef Storage StorageType;
    typedef typename Model::precision prec;

    /*
     * Initialize Hamiltonian for specific model and allocate storage
     * \param [in] p - alps::parameters
     */
#ifdef USE_MPI
    Hamiltonian(alps::params &p, MPI_Comm comm) :
      _comm(comm),
      _model(p),
      _storage(p, _model, comm) {};
#endif
    Hamiltonian(alps::params &p) :
      _model(p),
      _storage(p, _model) {};
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
      int k =0;
      while (_model.symmetry().next_sector()) {
        fill();
        /**
         * perform ARPACK call
         */
        int info = _storage.diag();
        if (info != 0) {

        } else {
          const std::vector < prec > &evals = _storage.eigenvalues();
          const std::vector < std::vector < prec > > &evecs = _storage.eigenvectors();
          for (int i = 0; i < evals.size(); ++i, ++k) {
            _eigenpairs.insert(EigenPair < prec, typename Model::Sector >(evals[i], evecs[i], k, _model.symmetry().sector()));
          }
        }
      }
#ifdef USE_MPI
      int rank;
      MPI_Comm_rank(_comm, &rank);
      if (rank == 0){
#endif
        std::cout << "Here is the list of eigenvalues:" << std::endl;
        for (typename std::set<EigenPair<prec, typename Model::Sector> >::iterator kkk = _eigenpairs.begin(); kkk != _eigenpairs.end(); kkk++) {
	  std::cout<<std::setprecision(14);
          std::cout << kkk->eigenvalue() << " ";
          kkk->sector().print();
          std::cout << std::endl;
        }
#ifdef USE_MPI
      }
#endif
    }

    Storage &storage() {
      return _storage;
    }

    const std::set < EigenPair < prec, typename Model::Sector > > &eigenpairs() const {
      return _eigenpairs;
    };

    Model &model() {
      return _model;
    }

#ifdef USE_MPI
    MPI_Comm& comm() {
      return _comm;
    }
#endif

  private:
    // CSR format Hamiltonian matrix storage
    Storage _storage;

    // Eigen-pairs
    std::set < EigenPair < prec, typename Model::Sector > > _eigenpairs;

    /**
     * Model to diagonalize
     */
    Model _model;

#ifdef USE_MPI
    MPI_Comm _comm;
#endif

  };

  typedef Hamiltonian < Storage::CRSStorage < Model::HubbardModel < double > >, Model::HubbardModel < double > > CSRHubbardHamiltonian;
  typedef Hamiltonian < Storage::SpinResolvedStorage < Model::HubbardModel < double > >, Model::HubbardModel < double > > SRSHubbardHamiltonian;
  typedef Hamiltonian < Storage::SOCRSStorage < Model::HubbardModel < double > >, Model::HubbardModel < double > > SOCSRHubbardHamiltonian;

  typedef Hamiltonian < Storage::CRSStorage < Model::HubbardModel < float > >, Model::HubbardModel < float > > CSRHubbardHamiltonian_float;
  typedef Hamiltonian < Storage::SpinResolvedStorage < Model::HubbardModel < float > >, Model::HubbardModel < float > > SRSHubbardHamiltonian_float;
  typedef Hamiltonian < Storage::SOCRSStorage < Model::HubbardModel < float > >, Model::HubbardModel < float > > SOCSRHubbardHamiltonian_float;

  typedef Hamiltonian < Storage::CRSStorage < Model::SingleImpurityAndersonModel < double > >, Model::SingleImpurityAndersonModel < double > > CSRSIAMHamiltonian;
  typedef Hamiltonian < Storage::CRSStorage < Model::SingleImpurityAndersonModel < float > >, Model::SingleImpurityAndersonModel < float > > CSRSIAMHamiltonian_float;

  typedef Hamiltonian < Storage::SpinResolvedStorage < Model::SingleImpurityAndersonModel < double > >, Model::SingleImpurityAndersonModel < double > > SRSSIAMHamiltonian;
  typedef Hamiltonian < Storage::SpinResolvedStorage < Model::SingleImpurityAndersonModel < float > >, Model::SingleImpurityAndersonModel < float > > SRSSIAMHamiltonian_float;
}
#endif //HUBBARD_HAMILTONIAN_H
