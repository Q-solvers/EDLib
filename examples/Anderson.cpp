//
// Created by iskakoff on 19/07/16.
//
#include <iostream>
#include <iomanip>

#include <edlib/EDParams.h>
#include <edlib/HDF5Utils.h>
#include <edlib/ChiLoc.h>
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
#ifdef USE_MPI
  if(comm.rank())
    ar.close();
#endif
  try {
#ifdef USE_MPI
    EDLib::SRSSIAMHamiltonian ham(params, comm);
#else
    EDLib::SRSSIAMHamiltonian ham(params);
#endif
    ham.diag();
    EDLib::hdf5::save_eigen_pairs(ham, ar, "results");
    EDLib::gf::GreensFunction < EDLib::SRSSIAMHamiltonian, alps::gf::matsubara_positive_mesh, alps::gf::statistics::statistics_type> greensFunction(params, ham,alps::gf::statistics::statistics_type::FERMIONIC);
    greensFunction.compute();
    greensFunction.save(ar, "results");
    EDLib::gf::ChiLoc<EDLib::SRSSIAMHamiltonian, alps::gf::matsubara_positive_mesh, alps::gf::statistics::statistics_type> susc(params, ham, alps::gf::statistics::statistics_type::BOSONIC);
    susc.compute();
    susc.save(ar, "results");
    susc.compute<EDLib::gf::NOperator<double> >();
    susc.save(ar, "results");
  } catch (std::exception & e) {
#ifdef USE_MPI
    if(comm.rank() == 0) std::cerr<<e.what();
#else
    std::cerr<<e.what();
#endif
  }
#ifdef USE_MPI
  if(!comm.rank())
#endif
  ar.close();
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
