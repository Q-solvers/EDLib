//
// Created by iskakoff on 18/08/16.
//

#ifndef HUBBARD_EDPARAMS_H
#define HUBBARD_EDPARAMS_H


#include <alps/params.hpp>
namespace EDLib {
  class EDParams : public alps::params {
  public:
    EDParams() : params() {
      define_parameters();
    }

    EDParams(int argc, const char **argv) : alps::params(argc, argv) {
      define_parameters();
    }

  private:
    /**
   * Define parameters used in programm
   */
    void define_parameters() {
      define < int >("NSITES", 4, "Number of sites");
      define < int >("siam.NORBITALS", 1, "Number of orbitals in single impurity Anderson Model.");
      define < int >("NSITES_BOSE", 2, "Number of bosonic sites. Should be 2^K - 1");
      define < int >("NSPINS", 2, "Number of spins");
      define < int >("arpack.NEV", 2, "Number of eigenvalues to find");
      define < int >("arpack.NCV", "Number of convergent values");
      define < bool >("arpack.SECTOR", "Read symmetry sectors from file");
      define < size_t >("storage.MAX_SIZE", 70000, "Number of eigenvalues to find");
      define < size_t >("storage.MAX_DIM", 5000, "Number of eigenvalues to find");
      define < int >("spinstorage.ORBITAL_NUMBER", 1, "Number of orbitals with interaction");
      define < std::string >("INPUT_FILE", "input.h5", "File with initial data");
      define < int >("lanc.NOMEGA", 32, "Number of fermionic frequencies");
      define < int >("lanc.NLANC", 100, "Number of Lanczos iterations");
      define < double >("lanc.BETA", 10.0, "Inverse temperature");
      define < double >("lanc.BOLTZMANN_CUTOFF", 1e-12, "Cutoff for Boltsman factor");
    }

  };

}
#endif //HUBBARD_EDPARAMS_H
