#include <iostream>
#include <alps/params.hpp>


#include "Hamiltonian.h"
#include <SzCombination.h>

int main(int argc, const char ** argv) {
  alps::params params(argc, argv);
  Hamiltonian<double, SzCombination> ham(1000, 10000, params);
  return 0;
}
