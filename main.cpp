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
    for(auto ipair = ham.eigenpairs().begin(); ipair != ham.eigenpairs().end(); ++ipair){
      pairs[count] = ipair;
      ++count;
    }
    int nev = params["arpack.NEV"].as<int>();
    if(
      ((pairs[0]->sector().nup() + 1) != pairs[nev]->sector().nup()) ||
      ((pairs[0]->sector().ndown() + 1) != pairs[nev]->sector().ndown())
    ){
      std::stringstream s;
      s << "Wrong sectors " << pairs[0]->sector().nup()   << " " << pairs[0]->sector().ndown()
        <<          " and " << pairs[nev]->sector().nup() << " " << pairs[nev]->sector().ndown()
        << ". Both number of electrons up and down must differ by 1." << std::endl;
      throw std::invalid_argument(s.str().c_str());
    }
    std::vector<double> outvec(1, double(0.0));
    std::vector<double> outvec2(1, double(0.0));
    double expectation_value = 0.0;
    std::vector<std::vector<std::vector<double>>> cc(ham.model().spins(), std::vector<std::vector<double>>(ham.model().interacting_orbitals(), std::vector<double>(ham.model().interacting_orbitals(), double(0.0))));
    double Z = 0.0;
    double cutoff = params["lanc.BOLTZMANN_CUTOFF"];
    double beta = params["lanc.BETA"].as<double>();
    for(size_t ipair1 = 0; ipair1 < nev; ++ipair1){
      double boltzmann_f_1 = std::sqrt(std::exp(-(pairs[ipair1]->eigenvalue() - pairs[0]->eigenvalue()) * beta));
      for(size_t ipair2 = nev; ipair2 < 2 * nev; ++ipair2){
        double boltzmann_f_2 = std::sqrt(std::exp(-(pairs[ipair2]->eigenvalue() - pairs[nev]->eigenvalue()) * beta));
        double boltzmann_f = boltzmann_f_1 * boltzmann_f_2;
        if (std::abs(cutoff - boltzmann_f) > std::numeric_limits<double>::epsilon() && boltzmann_f < cutoff ) {
          continue;
        }
#ifdef USE_MPI
        if(!rank)
#endif
        std::cout << "Compute cc contribution for eigenvalues E=" << pairs[ipair1]->eigenvalue() << ", " << pairs[ipair2]->eigenvalue()
                  << " with Boltzmann factor = " << boltzmann_f
                  << " for sectors" << pairs[ipair1]->sector() << " and" << pairs[ipair2]->sector() << std::endl;
        Z += boltzmann_f;
        for(size_t ispin = 0; ispin < ham.model().spins(); ++ispin){
          std::ostringstream pair_name;
          pair_name << "cc_" << ispin << "_pairs_" << ipair1 << "_" << ipair2 - nev << ".txt";
          std::ofstream pair_out(pair_name.str().c_str());
          for(size_t orb1 = 0; orb1 < ham.model().interacting_orbitals(); ++orb1){
            for(size_t orb2 = 0; orb2 < ham.model().interacting_orbitals(); ++orb2){
              double product = 0.0;
              ham.model().symmetry().set_sector(pairs[ipair2]->sector());
              if(greensFunction.annihilate_particles(std::array<size_t, 1>{{size_t(orb1)}}, ispin, pairs[ipair2]->eigenvector(), outvec2, expectation_value)) {
                if(greensFunction.annihilate_particles(std::array<size_t, 1>{{size_t(orb2)}}, (1 - ispin), outvec2, outvec, expectation_value)) {
                  product = ham.storage().vv(pairs[ipair1]->eigenvector(), outvec
#ifdef USE_MPI
                    , ham.comm()
#endif
                  );
                }
              }
              cc[ispin][orb1][orb2] += boltzmann_f * product;
              pair_out << product << "\t";
            }
            pair_out << std::endl;
          }
          pair_out.close();
        }
      }
    }
    for(size_t ispin = 0; ispin < ham.model().spins(); ++ispin){
      std::ostringstream cc_name;
      cc_name << "cc_" << ispin << ".txt";
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
