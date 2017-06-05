//
// Created by iskakoff on 01/06/17.
//

#include <iostream>

#include <edlib/EDParams.h>
#include "edlib/Hamiltonian.h"
#include "edlib/SzSymmetry.h"
#include "edlib/SOCRSStorage.h"
#include "edlib/CRSStorage.h"
#include "edlib/GreensFunction.h"
#include "edlib/ChiLoc.h"
#include "edlib/HDF5Utils.h"
#include "edlib/SpinResolvedStorage.h"
#include "edlib/StaticObervables.h"
#include "edlib/MeshFactory.h"
#include "ext/HolsteinAndersonModel.h"
#include "ext/HolsteinAndersonParameter.h"

int main(int argc, const char ** argv) {
  typedef EDLib::Hamiltonian<EDLib::Storage::CRSStorage<EDLib::Ext::Model::HolsteinAndersonModel<double> >, EDLib::Ext::Model::HolsteinAndersonModel<double> > HamType;
#ifdef USE_MPI
  MPI_Init(&argc, (char ***) &argv);
  alps::mpi::communicator comm;
#endif
  alps::params params(argc, argv);
  EDLib::define_parameters(params);
  EDLib::Ext::define_parameters(params);
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
    HamType ham(params, comm);
#else
    HamType ham(params);
#endif
    ham.diag();
    EDLib::common::statistics.updateEvent("diag");
    EDLib::common::statistics.registerEvent("GF");
    EDLib::hdf5::save_eigen_pairs(ham, ar, "results");
    EDLib::gf::GreensFunction < HamType, alps::gf::real_frequency_mesh> greensFunction(params, ham);
    greensFunction.compute();
    greensFunction.save(ar, "results");
    EDLib::gf::ChiLoc<HamType, alps::gf::real_frequency_mesh> susc(params, ham);
    susc.compute();
    susc.save(ar, "results");
    susc.compute<EDLib::gf::NOperator<double> >();
    susc.save(ar, "results");
    EDLib::common::statistics.updateEvent("GF");
    EDLib::StaticObervables<HamType> sd(params);
    std::map < std::string, std::vector < double>> observables = sd.calculate_static_observables(ham);
    EDLib::hdf5::save_static_observables(observables, ar, "results");
    
    EDLib::common::statistics.updateEvent("total");
    EDLib::common::statistics.print();
  } catch (std::exception & e) {
#ifdef USE_MPI
    if(comm.rank() == 0) std::cerr<<e.what();
    MPI_Abort(comm, -1);
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
