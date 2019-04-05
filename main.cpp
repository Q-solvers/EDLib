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

// freq, I, J, spin
template<typename freq_grid>
using gf_type = alps::gf::four_index_gf<std::complex<double>, freq_grid, alps::gf::real_space_index_mesh, alps::gf::real_space_index_mesh, alps::gf::index_mesh >;

using r_mesh_t = alps::gf::real_space_index_mesh;
using s_mesh_t = alps::gf::index_mesh;
using mat_mesh_t = alps::gf::matsubara_positive_mesh;

int main(int argc, const char ** argv) {
#ifdef USE_MPI
  typedef EDLib::SRSHubbardHamiltonian HamType;
#else
  typedef EDLib::SRSHubbardHamiltonian HamType;
#endif
#ifdef USE_MPI
  MPI_Init(&argc, (char ***) &argv);
#endif
  alps::params params(argc, argv);
  if(params.help_requested(std::cout)) {
    exit(0);
  }
  EDLib::define_parameters(params);
  params.define<std::string>("CLUSTER_DATA", "Cluster_results.h5", "Greens Function and Self-energy of cluster");
  alps::hdf5::archive ar;
#ifdef USE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(!rank)
#endif
  ar.open(params["OUTPUT_FILE"].as<std::string>().c_str(), "w");
  try {
#ifdef USE_MPI
    HamType ham(params, MPI_COMM_WORLD);
#else
    HamType ham(params);
#endif
    ham.diag();
    EDLib::StaticObservables<HamType> so(params);
    so.print_static_observables(ham);
    for (const auto& pair :ham.eigenpairs()) {
      so.print_major_electronic_configuration(ham, pair, 256, 1e-5);
      so.print_class_contrib(ham, pair, 256, 1e-5, true);
    }
    EDLib::hdf5::save_eigen_pairs(ham, ar, "results");
    EDLib::gf::GreensFunction < HamType, alps::gf::matsubara_positive_mesh, alps::gf::statistics::statistics_type> greensFunction(params, ham,alps::gf::statistics::statistics_type::FERMIONIC);
    //EDLib::gf::GreensFunction < HamType, alps::gf::real_frequency_mesh> greensFunction(params, ham);
    greensFunction.compute();
    greensFunction.save(ar, "results");
    EDLib::gf::ChiLoc<HamType, alps::gf::matsubara_positive_mesh, alps::gf::statistics::statistics_type> susc(params, ham, alps::gf::statistics::statistics_type::BOSONIC);
    //EDLib::gf::ChiLoc< HamType, alps::gf::real_frequency_mesh> susc(params, ham);
    susc.compute();
    susc.save(ar, "results");
    susc.compute<EDLib::gf::NOperator<double> >();
    susc.save(ar, "results");
    
    auto & G_ij = greensFunction.G_ij();
    auto sigma= greensFunction.compute_selfenergy(ar, "results");
    alps::gf::real_space_index_mesh r_mesh(ham.model().interacting_orbitals(), 2);
    gf_type<alps::gf::matsubara_positive_mesh> cluster_gf(G_ij.mesh1(), r_mesh, r_mesh, G_ij.mesh3());
    gf_type<alps::gf::matsubara_positive_mesh> cluster_sigma(G_ij.mesh1(), r_mesh, r_mesh, G_ij.mesh3());
    for(int iw = 0 ; iw<G_ij.mesh1().extent(); ++iw) {
      for(int is = 0; is<G_ij.mesh3().extent(); ++is) {
        for(int i = 0; i < G_ij.mesh2().extent(); ++i) {
          int I = i / r_mesh.extent();
          int J = i % r_mesh.extent();
          cluster_gf(mat_mesh_t::index_type(iw), r_mesh_t::index_type(I), r_mesh_t::index_type(J), s_mesh_t::index_type(is)) = 
              G_ij(mat_mesh_t::index_type(iw), s_mesh_t::index_type(i), s_mesh_t::index_type(is));
          cluster_sigma(mat_mesh_t::index_type(iw), r_mesh_t::index_type(I), r_mesh_t::index_type(J), s_mesh_t::index_type(is)) = 
              sigma(mat_mesh_t::index_type(iw), s_mesh_t::index_type(i), s_mesh_t::index_type(is));
        }
      }
    }
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(!rank) 
#endif
    ar.close();
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(!rank) 
#endif
    ar.open(params["CLUSTER_DATA"].as<std::string>().c_str(), "w");
    ar["G_ij"]<<cluster_gf;
    ar["Sigma_ij"]<<cluster_sigma;
//    EDLib::CSRSIAMHamiltonian ham2(params);
  } catch (std::exception & e) {
#ifdef USE_MPI
    if(!rank) std::cerr<<e.what()<<std::endl;
#else
    std::cerr<<e.what();
#endif
  }
#ifdef USE_MPI
  if(!rank)
#endif
  ar.close();
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
