#ifndef EDLIB_HAMILTONIAN_H
#define EDLIB_HAMILTONIAN_H

#include <iomanip>
#include <iostream>
#include <set>

#include "edlib/CRSStorage.h"
#include "edlib/EigenPair.h"
#include "edlib/HubbardModel.h"
#include "edlib/Parameters.h"
#include "edlib/SOCRSStorage.h"
#include "edlib/SingleImpurityAndersonModel.h"
#include "edlib/SpinResolvedStorage.h"

namespace edlib {

  /**
   * Hamiltonian == Model + Storage. Diagonalisation iterates over the model's
   * symmetry sectors, lets the storage build the per-sector matrix, and
   * collects EigenPairs across sectors.
   */
  template <class Storage>
  class Hamiltonian {
  public:
    using Model      = typename Storage::Model;
    using ModelType  = typename Storage::Model;
    using StorageType = Storage;
    using prec       = typename Model::precision;

#ifdef USE_MPI
    Hamiltonian(const Parameters& p, const typename Model::ModelData& bath, MPI_Comm comm)
        : _comm(comm),
          _model(p, bath),
          _storage(p, _model, comm) {}
#endif
    Hamiltonian(const Parameters& p, const typename Model::ModelData& bath)
        : _model(p, bath),
          _storage(p, _model) {}

    void fill() { _storage.fill(); }

    void diag() {
#ifdef USE_MPI
      int rank; MPI_Comm_rank(_comm, &rank);
#endif
      int k = 0;
      while (_model.symmetry().next_sector()) {
#ifdef USE_MPI
        if (rank == 0)
#endif
          std::cout << "Diagonalize sector " << _model.symmetry().sector() << std::endl;
        fill();
        int info = _storage.diag();
        if (info != 0) {
#ifdef USE_MPI
          if (rank == 0)
#endif
            std::cerr << "Eigenvalue have not been computed." << std::endl;
        } else {
          const auto& evals = _storage.eigenvalues();
          const auto& evecs = _storage.eigenvectors();
          for (std::size_t i = 0; i < evals.size(); ++i, ++k) {
            _eigenpairs.insert(EigenPair<prec, typename Model::Sector>(
                evals[i], evecs[i], k, _model.symmetry().sector()));
          }
        }
      }
#ifdef USE_MPI
      if (rank == 0) {
#endif
        std::cout << "Here is the list of eigenvalues:" << std::endl;
        std::streamsize old_p = std::cout.precision();
        std::cout << std::setprecision(14);
        for (auto it = _eigenpairs.begin(); it != _eigenpairs.end(); ++it) {
          std::cout << it->eigenvalue() << " ";
          it->sector().print();
          std::cout << std::endl;
        }
        std::cout << std::setprecision(old_p);
#ifdef USE_MPI
      }
#endif
    }

    Storage&       storage()       { return _storage; }
    Model&         model()         { return _model;   }

    const std::set<EigenPair<prec, typename Model::Sector>>& eigenpairs() const {
      return _eigenpairs;
    }

    void constant_shift(prec shift) { _storage.constant_shift(shift); }

#ifdef USE_MPI
    const MPI_Comm& comm() const { return _comm; }
#endif

  private:
    // _model must be declared (and constructed) before _storage since the
    // storage ctor takes a Model& reference.
#ifdef USE_MPI
    MPI_Comm _comm;
#endif
    Model   _model;
    Storage _storage;
    std::set<EigenPair<prec, typename Model::Sector>> _eigenpairs;
  };

  using CSRHubbardHamiltonian          = Hamiltonian<CRSStorage<HubbardModel<double>>>;
  using SRSHubbardHamiltonian          = Hamiltonian<SpinResolvedStorage<HubbardModel<double>>>;
  using SOCSRHubbardHamiltonian        = Hamiltonian<SOCRSStorage<HubbardModel<double>>>;

  using CSRHubbardHamiltonian_float    = Hamiltonian<CRSStorage<HubbardModel<float>>>;
  using SRSHubbardHamiltonian_float    = Hamiltonian<SpinResolvedStorage<HubbardModel<float>>>;
  using SOCSRHubbardHamiltonian_float  = Hamiltonian<SOCRSStorage<HubbardModel<float>>>;

  using CSRSIAMHamiltonian             = Hamiltonian<CRSStorage<SingleImpurityAndersonModel<double>>>;
  using CSRSIAMHamiltonian_float       = Hamiltonian<CRSStorage<SingleImpurityAndersonModel<float>>>;
  using SRSSIAMHamiltonian             = Hamiltonian<SpinResolvedStorage<SingleImpurityAndersonModel<double>>>;
  using SRSSIAMHamiltonian_float       = Hamiltonian<SpinResolvedStorage<SingleImpurityAndersonModel<float>>>;

}

#endif
