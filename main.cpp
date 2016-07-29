#include <iostream>
#include <alps/params.hpp>


#include <Hamiltonian.h>
#include <SzSymmetry.h>
#include <SOCRSStorage.h>


/**
 * Define parameters used in programm
 * \param p - alps paramters
 */
void define_parameters(alps::params &p) {
  p.define<int>("NSITES", 4, "Number of sites");
  p.define<int>("NSPINS", 2, "Number of spins");
  p.define<int>("arpack.NEV", 2, "Number of eigenvalues to find");
  p.define<int>("arpack.NCV", "Number of convergent values");
  p.define<bool>("arpack.SECTOR", "Read symmetry sectors from file");
  p.define<size_t>("storage.MAX_SIZE", 70000, "Number of eigenvalues to find");
  p.define<size_t>("storage.MAX_DIM", 5000, "Number of eigenvalues to find");
  p.define<std::string>("INPUT_FILE", "input.h5", "File with initial data");
}


int main(int argc, const char ** argv) {
  alps::params params(argc, argv);
  define_parameters(params);
  if(params.help_requested(std::cout)) {
    exit(0);
  }
  Hamiltonian<double, SzSymmetry, SOCRSStorage<double, SzSymmetry, HubbardModel<double> > , HubbardModel<double> > ham(params);
//  Hamiltonian<double, SzSymmetry, CRSStorage<double> > ham(100000, 100000, params);
  ham.diag();
  return 0;
}
