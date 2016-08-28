//
// Created by iskakoff on 22/08/16.
//

#include <gtest/gtest.h>
#include "Hamiltonian.h"
#include "HubbardModel.h"
#include "Storage.h"

TEST(HubbardModelTest, Filling) {
  const char *string = "--NSITES=16 --INPUT_FILE=../input/input.h5";
  EDLib::EDParams p(1, &string);
  EDLib::SOCSRHubbardHamiltonian_float ham(p);
  while(ham.model().symmetry().next_sector()) {
    std::cout<<ham.model().symmetry().sector().size()<<std::endl;
  }
}