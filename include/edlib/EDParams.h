//
// Created by iskakoff on 18/08/16.
//

#ifndef HUBBARD_EDPARAMS_H
#define HUBBARD_EDPARAMS_H


#include <alps/params.hpp>
namespace EDLib {
  void define_parameters(alps::params &params) {
    // General parameters
    params.define < int >("NSITES", 4, "Number of sites");
    params.define < int >("NSPINS", 2, "Number of spins");
    params.define < std::string >("INPUT_FILE", "input.h5", "File with initial data");
    params.define < std::string >("OUTPUT_FILE", "sim.h5", "File with results");
    // Symmetry parameters
    params.define < bool >("arpack.SECTOR", "Read symmetry sectors from file");
    // Storage parameters
    params.define < size_t >("storage.MAX_DIM", 5000, "Maximum dimension of the Hamiltonian matrix.");
    params.define < size_t >("storage.MAX_SIZE", 70000, "Maximum size of the matrix arrays. Must be between MAX_DIM and MAX_DIM^2.");
    params.define < int >("storage.EIGENVALUES_ONLY", 0, "Compute only eigenvalues.");
    params.define < int >("spinstorage.ORBITAL_NUMBER", 1, "Number of orbitals with interaction");
    // ARPACK parameters
    params.define < int >("arpack.NEV", 2, "Number of eigenvalues to find");
    params.define < int >("arpack.NCV", "Number of convergent values");
    // Lanczos parameters
    params.define < int >("lanc.NOMEGA", 32, "Number of fermionic frequencies");
    params.define < double >("lanc.EMIN", -3, "Lowest real frequency value");
    params.define < double >("lanc.EMAX", 3, "Largest real frequency value");
    params.define < int >("lanc.NLANC", 100, "Number of Lanczos iterations");
    params.define < double >("lanc.BETA", 10.0, "Inverse temperature");
    params.define < double >("lanc.BOLTZMANN_CUTOFF", 1e-12, "Cutoff for Boltzmann factor");

    // Anderson model
    params.define < int >("siam.NORBITALS", 1, "Number of orbitals in single impurity Anderson Model.");
  }
}
#endif //HUBBARD_EDPARAMS_H
