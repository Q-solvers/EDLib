//
// Created by iskakoff on 22/08/16.
//

#include <gtest/gtest.h>
#include "edlib/Hamiltonian.h"
#include "edlib/HubbardModel.h"
#include "edlib/Storage.h"

TEST(HubbardModelTest, Filling) {
  const char *string = "--NSITES=16 --INPUT_FILE=../input/input.h5";
  alps::params p(1, &string);
  EDLib::SRSHubbardHamiltonian_float ham(p
#ifdef USE_MPI
  , MPI_COMM_WORLD
#endif
);
  while(ham.model().symmetry().next_sector()) {
    std::cout<<ham.model().symmetry().sector().size()<<std::endl;
  }
}
