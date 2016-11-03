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
  typedef EDLib::SRSHubbardHamiltonian HamType;
#ifdef USE_MPI
  MPI_Init(&argc, (char ***) &argv);
#endif
  alps::params params(argc, argv);
  if(params.help_requested(std::cout)) {
    exit(0);
  }
  EDLib::define_parameters(params);
  try {
//  CSRHubbardHamiltonian_float ham(params);
#ifdef USE_MPI
//    EDLib::SRSHubbardHamiltonian ham(params, comm);
    HamType ham(params, MPI_COMM_WORLD);
#else
    EDLib::SRSHubbardHamiltonian ham(params);
#endif
    ham.diag();
//  GreensFunction<float, CSRHubbardHamiltonian_float > greensFunction(params, ham);
    EDLib::gf::GreensFunction < double, HamType > greensFunction(params, ham);
    greensFunction.compute();
//    EDLib::CSRSIAMHamiltonian ham2(params);
  } catch (std::exception & e) {
#ifdef USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0) std::cerr<<e.what()<<std::endl;
#else
    std::cerr<<e.what();
#endif
  }
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
