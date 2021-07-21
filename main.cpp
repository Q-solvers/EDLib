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
    std::vector<double> outvec(1, double(0.0));
    std::vector<double> outvec2(1, double(0.0));
    double expectation_value = 0.0;
    std::vector<std::vector<std::vector<double>>> cc(ham.model().spins(), std::vector<std::vector<double>>(ham.model().interacting_orbitals(), std::vector<double>(ham.model().interacting_orbitals(), double(0.0))));
    double Z = 0.0;
    double cutoff = params["lanc.BOLTZMANN_CUTOFF"];
    double beta = params["lanc.BETA"].as<double>();

    for(auto ipair = ham.eigenpairs().begin(); ipair != ham.eigenpairs().end(); ++ipair){
      double boltzmann_f = std::exp(-(ipair->eigenvalue() - ham.eigenpairs().begin()->eigenvalue()) * beta);
      if (std::abs(cutoff - boltzmann_f) > std::numeric_limits<double>::epsilon() && boltzmann_f < cutoff ) {
        continue;
      }
#ifdef USE_MPI
      if(!rank)
#endif
      std::cout << "Compute c+c contribution for eigenvalue E=" << ipair->eigenvalue() << " with Boltzmann factor = " << boltzmann_f << "; for sector" << ipair->sector() << std::endl;
      Z += boltzmann_f;
      for(size_t ispin = 0; ispin < ham.model().spins(); ++ispin){
        for(size_t orb1 = 0; orb1 < ham.model().interacting_orbitals(); ++orb1){
          for(size_t orb2 = 0; orb2 < ham.model().interacting_orbitals(); ++orb2){
            ham.model().symmetry().set_sector(ipair->sector());
            if(greensFunction.annihilate_particles(std::array<size_t, 1>{{size_t(orb1)}}, ispin, ipair->eigenvector(), outvec, expectation_value)) {
              ham.model().symmetry().set_sector(ipair->sector());
              if(greensFunction.annihilate_particles(std::array<size_t, 1>{{size_t(orb2)}}, ispin, ipair->eigenvector(), outvec2, expectation_value)) {
                cc[ispin][orb1][orb2] += boltzmann_f * ham.storage().vv(outvec, outvec2
#ifdef USE_MPI
                  , ham.comm()
#endif
                );
              }
            }
          }
        }
      }
    }
    for(size_t ispin = 0; ispin < ham.model().spins(); ++ispin){
      std::ostringstream cc_name;
      cc_name << "cdc_" << ispin << ".txt";
      std::ofstream cc_out(cc_name.str().c_str());
      for(size_t orb1 = 0; orb1 < ham.model().interacting_orbitals(); ++orb1){
        for(size_t orb2 = 0; orb2 < ham.model().interacting_orbitals(); ++orb2){
          cc_out << cc[ispin][orb1][orb2] / Z << "\t";
        }
        cc_out << std::endl;
      }
      cc_out.close();
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
