//
// Created by iskakoff on 19/07/16.
//

#include "edlib/Hamiltonian.h"
#include "edlib/GreensFunction.h"

int main(int argc, const char ** argv) {
#ifdef USE_MPI
  MPI_Init(&argc, (char ***) &argv);
  alps::mpi::communicator comm;
#endif
  EDLib::EDParams params(argc, argv);
  if(params.help_requested(std::cout)) {
    exit(0);
  }
  try {
#ifdef USE_MPI
    EDLib::SRSSIAMHamiltonian ham(params, comm);
#else
    EDLib::SRSSIAMHamiltonian ham(params);
#endif
    ham.diag();
    EDLib::gf::GreensFunction < double, EDLib::SRSSIAMHamiltonian > greensFunction(params, ham);
    greensFunction.compute();
  } catch (std::exception & e) {
#ifdef USE_MPI
    if(comm.rank() == 0) std::cerr<<e.what();
#else
    std::cerr<<e.what();
#endif
  }
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}