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
#include "edlib/EigenPair.h"

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
    EDLib::gf::GreensFunction < HamType, EDLib::MatsubaraMeshFactory, alps::gf::statistics::statistics_type> greensFunction(params, ham,alps::gf::statistics::statistics_type::FERMIONIC);
    size_t count = 0;
    std::vector<decltype(ham.eigenpairs().begin())> pairs(ham.eigenpairs().size());
    for(auto ipair = ham.eigenpairs().begin(); ipair != ham.eigenpairs().end(); ++ipair){
      pairs[count] = ipair;
      ++count;
    }
    if((pairs[0]->sector().nup() != 7) || (pairs[0]->sector().ndown() != 7)){
     std::cout << "Wrong sector 0 " << pairs[0]->sector().nup() << " " << pairs[0]->sector().ndown() << std::endl;
    }
    if((pairs[1]->sector().nup() != 7) || (pairs[1]->sector().ndown() != 7)){
     std::cout << "Wrong sector 1 " << pairs[1]->sector().nup() << " " << pairs[1]->sector().ndown() << std::endl;
    }
    if((pairs[2]->sector().nup() != 7) || (pairs[2]->sector().ndown() != 7)){
     std::cout << "Wrong sector 1 " << pairs[2]->sector().nup() << " " << pairs[2]->sector().ndown() << std::endl;
    }
    if((pairs[32]->sector().nup() != 8) || (pairs[32]->sector().ndown() != 8)){
     std::cout << "Wrong sector 2 " << pairs[32]->sector().nup() << " " << pairs[32]->sector().ndown() << std::endl;
    }
    std::vector<double> outvec(1, double(0.0));
    std::vector<double> outvec2(1, double(0.0));
    double expectation_value = 0.0;
    for(size_t ipair = 0; ipair < 3; ++ipair){
      for(size_t ispin = 0; ispin < ham.model().spins(); ++ispin){
        std::ostringstream cc_name;
        cc_name << "cc_" << ispin << "_pair" << ipair << ".txt";
        std::ofstream cc_out(cc_name.str().c_str());
        for(size_t orb1 = 0; orb1 < ham.model().interacting_orbitals(); ++orb1){
          for(size_t orb2 = 0; orb2 < ham.model().interacting_orbitals(); ++orb2){
            double product = 0.0;
            ham.model().symmetry().set_sector(pairs[32]->sector());
            if(greensFunction.annihilate_particles(std::array<size_t, 1>{{size_t(orb1)}}, ispin, pairs[32]->eigenvector(), outvec2, expectation_value)) {
              if(greensFunction.annihilate_particles(std::array<size_t, 1>{{size_t(orb2)}}, (1 - ispin), outvec2, outvec, expectation_value)) {
                product = ham.storage().vv(pairs[ipair]->eigenvector(), outvec
#ifdef USE_MPI
                  , ham.comm()
#endif
                );
              }
            }
            cc_out << product << "\t";
          }
          cc_out << std::endl;
        }
        cc_out.close();
      }
    }
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
