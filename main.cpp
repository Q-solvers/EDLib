#include <iostream>
#include <alps/params.hpp>


#include "Hamiltonian.h"
#include <SzCombination.h>
#include <CRSStorage.h>

int main(int argc, const char ** argv) {
  alps::params params(argc, argv);
  Hamiltonian<double, SzCombination, CRSStorage<double> > ham(1000, 10000, params);
  return 0;
}
