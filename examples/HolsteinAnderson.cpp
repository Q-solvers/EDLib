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
    EDLib::gf::GreensFunction < HamType, alps::gf::matsubara_positive_mesh, alps::gf::statistics::statistics_type> greensFunction(params, ham, alps::gf::statistics::FERMIONIC);
    greensFunction.compute();
    greensFunction.save(ar, "results");
    EDLib::StaticObervables<HamType> sd(params);
    // compute static observables
    std::map < std::string, std::vector < double>> observables = sd.calculate_static_observables(ham);
    EDLib::gf::ChiLoc<HamType, alps::gf::matsubara_positive_mesh, alps::gf::statistics::statistics_type > susc(params, ham, alps::gf::statistics::BOSONIC);
    // compute average magnetic moment
    double avg = 0;
    for(auto x : observables[sd._M_]) {
      avg += x / (2.0*observables[sd._M_].size());
    }
    susc.compute<EDLib::gf::SzOperator<double>>(&avg);
    susc.save(ar, "results");
    // compute average occupancy moment
    avg = 0;
    for(auto x : observables[sd._N_]) {
      avg += x / double(observables[sd._N_].size());
    }
    susc.compute<EDLib::gf::NOperator<double>>(&avg);
    susc.save(ar, "results");
    EDLib::common::statistics.updateEvent("GF");
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
