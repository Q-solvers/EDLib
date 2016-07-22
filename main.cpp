#include <iostream>
#include <alps/params.hpp>


#include <Hamiltonian.h>
#include <SzSymmetry.h>


/**
 * Define parameters used in programm
 * \param p - alps paramters
 */
void define_parameters(alps::params &p) {
  p.define<int>("NSITES", 4, "Number of sites");
  p.define<int>("NSPINS", 2, "Number of spins");
}


int main(int argc, const char ** argv) {
  alps::params params(argc, argv);
  define_parameters(params);
  Hamiltonian<double, SzSymmetry, CRSStorage<double> > ham(1000, 10000, params);
  ham.diag();
  return 0;
}
