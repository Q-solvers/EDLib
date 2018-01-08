//
// Created by iskakoff on 19/07/16.
//
#include <iostream>

#include <edlib/EDParams.h>
#include "edlib/Hamiltonian.h"
#include "edlib/SzSymmetry.h"
#include "edlib/SOCRSStorage.h"
#include "edlib/CRSStorage.h"
#include "edlib/HubbardModel.h"
#include "edlib/GreensFunction.h"
#include "edlib/ChiLoc.h"
#include "edlib/HDF5Utils.h"
#include "edlib/SpinResolvedStorage.h"
#include "edlib/StaticObservables.h"
#include "edlib/MeshFactory.h"
#include "edlib/ExecutionStatistic.h"



int main(int argc, const char ** argv) {
  // Initialize MPI environment if was enabled
#ifdef USE_MPI
  MPI_Init(&argc, (char ***) &argv);
  alps::mpi::communicator comm;
#endif
  // Define and read parameters
  alps::params params(argc, argv);
  EDLib::define_parameters(params);
  if(params.help_requested(std::cout)) {
    exit(0);
  }
  // open output file
  alps::hdf5::archive ar;
#ifdef USE_MPI
  if(!comm.rank())
#endif
    ar.open(params["OUTPUT_FILE"].as<std::string>(), "w");;

  try {
    // Construct Hamiltonian object
#ifdef USE_MPI
    EDLib::SRSHubbardHamiltonian ham(params, comm);
#else
    EDLib::SRSHubbardHamiltonian ham(params);
#endif
    // Perform Hamiltonian matrix diagonalization
    ham.diag();
    // Save eigenvalues to HDF5 file
    EDLib::hdf5::save_eigen_pairs(ham, ar, "results");
    // Construct single-particle local Green's function object
    EDLib::gf::GreensFunction < EDLib::SRSHubbardHamiltonian, alps::gf::matsubara_positive_mesh, alps::gf::statistics::statistics_type> greensFunction(params, ham,alps::gf::statistics::statistics_type::FERMIONIC);
    // Compute and store local Green's function
    greensFunction.compute();
    greensFunction.save(ar, "results");
    // Construct two particle Green's function object
    EDLib::gf::ChiLoc<EDLib::SRSHubbardHamiltonian, alps::gf::matsubara_positive_mesh, alps::gf::statistics::statistics_type> susc(params, ham, alps::gf::statistics::statistics_type::BOSONIC);
    // Compute and store spin susceptibility
    susc.compute();
    susc.save(ar, "results");
    // Compute and store charge susceptibility
    susc.compute<EDLib::gf::NOperator<double> >();
    susc.save(ar, "results");
    EDLib::common::statistics.print();
  } catch (std::exception & e) {
#ifdef USE_MPI
    if(comm.rank() == 0) {
      std::cerr<<e.what();
      ar.close();
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
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
