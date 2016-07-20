#include <iostream>
#include <alps/params.hpp>


#include "Hamiltonian.h"
#include <SzCombination.h>
#include <CRSStorage.h>
#include <fortranbinding.h>


/**
 * Define parameters used in programm
 * \param p - alps paramters
 */
void define_parameters(alps::params &p) {
  p.define<int>("NS", 3, "Number of sites");
  p.define<int>("SPIN", 2, "Number of spins");
}


int main(int argc, const char ** argv) {
  alps::params params(argc, argv);
  define_parameters(params);
  Hamiltonian<double, SzCombination, CRSStorage<double> > ham(1000, 10000, params);
  int nloc = 1;
  std::vector<double> vout(10, 0.0);
  std::vector<double> eout(1, 0.0);
  int ncv = 1;
  double Hstate0 = 0.0;
  int nev = 1;
  int ierr = 0;
  int info = 0;
  std::cout<<"Before Arpack test"<<std::endl;
  darnoldi(nloc,vout,eout, ncv,Hstate0, nev,ierr,info);
  std::cout<<"After Arpack test"<<std::endl;
  return 0;
}
