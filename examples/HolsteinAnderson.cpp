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
#include "edlib/StaticObservables.h"
#include "edlib/MeshFactory.h"
#include "ext/HolsteinAndersonModel.h"
#include "ext/HolsteinAndersonParameter.h"

int main(int argc, const char ** argv) {
  // Define Hamiltonian object type
  typedef EDLib::Hamiltonian<EDLib::Storage::CRSStorage<EDLib::Ext::Model::HolsteinAndersonModel<double> > > HamType;
  typedef typename HamType::ModelType::Sector sector;
  // Initialize MPI environment if enabled
#ifdef USE_MPI
  MPI_Init(&argc, (char ***) &argv);
  alps::mpi::communicator comm;
#endif
  // Define and read parameters
  alps::params params(argc, argv);
  EDLib::define_parameters(params);
  EDLib::Ext::define_parameters(params);
  params.define<int>("NORBITALS", 1, "");
  params.define<bool>("COMPUTE_REAL", false, "");
  if(params.help_requested(std::cout)) {
    exit(0);
  }
  // Open output HDF5 archive
  alps::hdf5::archive ar;
#ifdef USE_MPI
  if(!comm.rank())
#endif
    ar.open(params["OUTPUT_FILE"].as<std::string>(), "w");
  try {
    // Construct Hamiltonian object
#ifdef USE_MPI
    HamType ham(params, comm);
#else
    HamType ham(params);
#endif
    // Perform Hamiltonian matrix diagonalization
    ham.diag();
    // Store eigenvalues in HDF5 archive
    EDLib::hdf5::save_eigen_pairs(ham, ar, "results");
    // Construct single-particle Green's function object
    EDLib::gf::GreensFunction < HamType, alps::gf::matsubara_positive_mesh, alps::gf::statistics::statistics_type> greensFunction(params, ham, alps::gf::statistics::FERMIONIC);
    // Compute and store single particle Green's function
    greensFunction.compute();
    greensFunction.save(ar, "results");
    if(params["COMPUTE_REAL"].as<bool>()) {
      // Construct Green's function object
      EDLib::gf::GreensFunction < HamType, alps::gf::real_frequency_mesh> greens(params, ham);
      // Compute and save Green's function
      greens.compute();
      greens.save(ar, "results");
    }
    // Construct static observables object
    EDLib::StaticObservables<HamType> sd(params);
    // compute static observables
    std::map < std::string, std::vector < double>> observables = sd.calculate_static_observables(ham);
    EDLib::gf::ChiLoc<HamType, alps::gf::matsubara_positive_mesh, alps::gf::statistics::statistics_type > susc(params, ham, alps::gf::statistics::BOSONIC);
    // compute average magnetic moment
    double avg = 0;
    for(auto x : observables[sd._M_]) {
      avg += x / (2.0*observables[sd._M_].size());
    }
    // compute spin susceptibility
    susc.compute<EDLib::gf::SzOperator<double>>(&avg);
    susc.save(ar, "results");
    // compute average occupancy moment
    avg = 0;
    for(auto x : observables[sd._N_]) {
      avg += x / double(observables[sd._N_].size());
    }
    // Compute sharge susceptibility
    susc.compute<EDLib::gf::NOperator<double>>(&avg);
    susc.save(ar, "results");
    EDLib::common::statistics.updateEvent("GF");
    EDLib::hdf5::save_static_observables(observables, ar, "results");

    EDLib::common::statistics.updateEvent("total");
    EDLib::common::statistics.print();

    double sum = 0.0;
    double beta = params["lanc.BETA"].as<double>();
    int bos_dim = ham.model().bos_dim();
    std::vector<double> bosonic_occ(ham.model().bos_dim(), 0.0);
    std::vector<std::vector<double> > bosonic_occd(ham.model().bos_dim(), std::vector<double>(100, 0.0) );
    const EDLib::EigenPair<double, sector> &groundstate =  *ham.eigenpairs().begin();
    for(auto ipair = ham.eigenpairs().begin(); ipair != ham.eigenpairs().end(); ++ipair){
        const EDLib::EigenPair<double, sector>& pair = *ipair;
        // Calculate Boltzmann factor, skip the states with trivial contribution.
        double boltzmann_f = std::exp(
         -(pair.eigenvalue() - groundstate.eigenvalue()) * beta
        );
        if(boltzmann_f < 1e-9){
//          std::cout<<"Skipped by Boltzmann factor."<<std::endl;
          continue;
        }
        ham.model().symmetry().set_sector(pair.sector());
        ham.storage().reset();
        for(int i = 0; i < pair.eigenvector().size(); ++i){
          double weight = pair.eigenvector()[i] * pair.eigenvector()[i];
          // Calculate static variables for each orbital.
          ham.model().symmetry().next_state();
          long long nst = ham.model().symmetry().state();
          for(int orb = 0; orb < bos_dim; ++orb){
            int nbos = ham.model().number_of_bosons(nst, orb);
            //std::cout<<nst<<" "<<weight<<" "<<nbos<<std::endl;
            bosonic_occ[orb] += nbos * weight * boltzmann_f;
            bosonic_occd[orb][nbos] += 1.0 * weight * boltzmann_f;
          }
        }
        sum += boltzmann_f;
    }
    for(int orb = 0; orb < bos_dim; ++orb){
      bosonic_occ[orb]/= sum;
      for(int nbos= 0; nbos<100;++nbos)
      bosonic_occd[orb][nbos] /= sum;
    }
    std::cout<<"statsum: "<<sum<<std::endl;
    std::cout<<"=========================\nBosonic occupancy:"<<std::endl;
    for(int orb = 0; orb < bos_dim; ++orb){
      std::cout<<bosonic_occ[orb]<<std::endl;
      for(int i = 0; i< 100; ++i) {
        std::cout<<bosonic_occd[orb][i]<<"\n";
      }
      std::cout<<std::endl;
    }
    std::cout<<"========================="<<std::endl;
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
