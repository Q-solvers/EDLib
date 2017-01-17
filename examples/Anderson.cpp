//
// Created by iskakoff on 19/07/16.
//

#include <edlib/EDParams.h>
#include <edlib/HDF5Utils.h>
#include "edlib/Hamiltonian.h"
#include "edlib/GreensFunction.h"
#include "edlib/StateDescription.h"

int main(int argc, const char ** argv) {
#ifdef USE_MPI
  MPI_Init(&argc, (char ***) &argv);
  alps::mpi::communicator comm;
#endif
  alps::params params(argc, argv);
  EDLib::define_parameters(params);
  if(params.help_requested(std::cout)) {
    exit(0);
  }
  alps::hdf5::archive ar(params["OUTPUT_FILE"], alps::hdf5::archive::WRITE);
  try {
#ifdef USE_MPI
    EDLib::SRSSIAMHamiltonian ham(params, comm);
#else
    EDLib::SRSSIAMHamiltonian ham(params);
#endif
    ham.diag();
    EDLib::hdf5::save_eigen_pairs(ham, ar, "results");
    EDLib::StateDescription<HamType> sd(ham);
    for(auto & pair:ham.eigenpairs())
      sd.print(pair, 10, 1e-5);
    //EDLib::gf::GreensFunction < EDLib::SRSSIAMHamiltonian > greensFunction(params, ham);
    //greensFunction.compute();
    //greensFunction.save(ar, "results");
  } catch (std::exception & e) {
#ifdef USE_MPI
    if(comm.rank() == 0) std::cerr<<e.what();
#else
    std::cerr<<e.what();
#endif
  }
  ar.close();
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
