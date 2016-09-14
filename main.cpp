#include <iostream>


#include <EDParams.h>
#include <Hamiltonian.h>
#include <SzSymmetry.h>
#include <SOCRSStorage.h>
#include <CRSStorage.h>
#include <HubbardModel.h>
#include <GreensFunction.h>
#include <SpinResolvedStorage.h>
#include "HolsteinAndersonModel.h"


int main(int argc, const char ** argv) {
#ifdef ALPS_HAVE_MPI
  MPI_Init(&argc, (char ***) &argv);
  alps::mpi::communicator comm;
#endif
  EDLib::EDParams params(argc, argv);
  if(params.help_requested(std::cout)) {
    exit(0);
  }
  try {
//  CSRHubbardHamiltonian_float ham(params);
#ifdef ALPS_HAVE_MPI
    EDLib::SRSHubbardHamiltonian ham(params, comm);
#else
    EDLib::SRSHubbardHamiltonian ham(params);
#endif
//  SOCSRHubbardHamiltonian_float ham(params);
    ham.diag();
//  GreensFunction<float, CSRHubbardHamiltonian_float > greensFunction(params, ham);
//    EDLib::gf::GreensFunction < double, EDLib::SRSHubbardHamiltonian > greensFunction(params, ham);
//    greensFunction.compute();
//    EDLib::CSRSIAMHamiltonian ham2(params);
  } catch (std::exception & e) {
#ifdef ALPS_HAVE_MPI
    if(comm.rank() == 0) std::cerr<<e.what();
#else
    std::cerr<<e.what();
#endif
  }
#ifdef ALPS_HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
