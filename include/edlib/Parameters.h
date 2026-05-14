#ifndef EDLIB_PARAMETERS_H
#define EDLIB_PARAMETERS_H

#include <cstddef>
#include <string>

namespace edlib {

  struct Parameters {
    int    nsites                      = 4;
    int    nspins                      = 2;

    bool   arpack_sector               = false;
    int    arpack_nev                  = 2;
    int    arpack_ncv                  = 0;

    std::size_t storage_max_dim        = 5000;
    std::size_t storage_max_size       = 70000;
    bool        eigenvalues_only       = false;
    int         spinstorage_orbital_number = 1;

    int    lanc_nomega                 = 32;
    double lanc_emin                   = -3.0;
    double lanc_emax                   =  3.0;
    int    lanc_nlanc                  = 100;
    double lanc_beta                   = 10.0;
    double lanc_boltzmann_cutoff       = 1e-12;

    int    siam_norbitals              = 1;
  };

}

#endif
