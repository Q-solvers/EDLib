#include <iostream>

#include <edlib/EDParams.h>
#include "edlib/Hamiltonian.h"
#include "edlib/SzSymmetry.h"
#include "edlib/SOCRSStorage.h"
#include "edlib/CRSStorage.h"
#include "edlib/HubbardModel.h"
#include "edlib/GreensFunction.h"
#include "edlib/SpinResolvedStorage.h"


int main(int argc, const char ** argv) {
#ifdef USE_MPI
  MPI_Init(&argc, (char ***) &argv);
  alps::mpi::communicator comm;
#endif
  alps::params params(argc, argv);
  if(params.help_requested(std::cout)) {
    exit(0);
  }
  EDLib::define_parameters(params);
//  try {
//  CSRHubbardHamiltonian_float ham(params);
#ifdef USE_MPI
    EDLib::SRSHubbardHamiltonian ham(params, comm);
#else
    EDLib::SRSHubbardHamiltonian ham(params);
#endif
//  SOCSRHubbardHamiltonian_float ham(params);
    ham.diag();
//  GreensFunction<float, CSRHubbardHamiltonian_float > greensFunction(params, ham);
    EDLib::gf::GreensFunction < double, EDLib::SRSHubbardHamiltonian > greensFunction(params, ham);
    greensFunction.compute();
//    EDLib::CSRSIAMHamiltonian ham2(params);
  /*} catch (std::exception & e) {
#ifdef USE_MPI
    if(comm.rank() == 0) std::cerr<<e.what();
#else
    std::cerr<<e.what();
#endif
  }*/
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
