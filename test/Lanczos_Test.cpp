//
// Created by iskakoff on 22/08/16.
//

#include <gtest/gtest.h>
#include "edlib/Hamiltonian.h"
#include "edlib/HubbardModel.h"
#include "edlib/Storage.h"
#include "edlib/EDParams.h"
#include "edlib/StaticObservables.h"
#include "edlib/GreensFunction.h"
#include "edlib/ChiLoc.h"


#ifdef USE_MPI

class HubbardModelTestEnv : public ::testing::Environment {
  protected:
  virtual void SetUp() {
    char** argv;
    int argc = 0;
    int mpiError = MPI_Init(&argc, &argv);
  }

  virtual void TearDown() {
    MPI_Finalize();
  }

  ~HubbardModelTestEnv(){};

};

::testing::Environment* const foo_env = AddGlobalTestEnvironment(new HubbardModelTestEnv);

#endif


TEST(HubbardModelTest, ReferenceTest) {
  alps::params p;
  EDLib::define_parameters(p);
  p["NSITES"]=4;
  p["NSPINS"]=2;
  p["INPUT_FILE"]="test/input/GF_Chi/input.h5";
  p["arpack.SECTOR"]=false;
  p["storage.MAX_SIZE"]=864;
  p["storage.MAX_DIM"]=36;
  p["storage.EIGENVALUES_ONLY"]=false;
  p["storage.ORBITAL_NUMBER"]=1;
  p["arpack.NEV"]=100;
  p["lanc.BETA"]=20;
  p["lanc.EMIN"]=-4.0;
  p["lanc.EMAX"]=4.0;
  p["lanc.NOMEGA"]=1000;

#ifdef USE_MPI
  typedef EDLib::SRSHubbardHamiltonian HamType;
#else
  typedef EDLib::SOCSRHubbardHamiltonian HamType;
#endif
  HamType ham(p
#ifdef USE_MPI
  , MPI_COMM_WORLD
#endif
  );

  ham.diag();

  // Compute our GFs for the reference model.
  EDLib::gf::GreensFunction < HamType, alps::gf::matsubara_positive_mesh, alps::gf::statistics::statistics_type> greensFunction(p, ham,alps::gf::statistics::statistics_type::FERMIONIC);
  greensFunction.compute();
  auto G = greensFunction.G();
  EDLib::StaticObservables<HamType> so(p);
  std::map<std::string, std::vector<double>> observables = so.calculate_static_observables(ham);
  EDLib::gf::ChiLoc<HamType, alps::gf::matsubara_positive_mesh, alps::gf::statistics::statistics_type> susc(p, ham, alps::gf::statistics::statistics_type::BOSONIC);
  // compute average magnetic moment
  double avg = 0.0;
  for(auto x : observables[so._M_]) {
    avg += x / (2.0*observables[so._M_].size());
  }
  // compute spin susceptibility
  susc.compute<EDLib::gf::SzOperator<double>>(&avg);
  auto ChiSz = susc.G();
  // compute average occupancy moment
  avg = 0.0;
  for(auto x : observables[so._N_]) {
    avg += x / double(observables[so._N_].size());
  }
  // Compute sharge susceptibility
  susc.compute<EDLib::gf::NOperator<double>>(&avg);
  auto ChiN = susc.G();

  // Read in the reference GFs.
  // FIXME Desperate kludges. Must take the indextypes from GFs instead.
  auto G_file = G;
  std::ifstream infile("test/input/GF_Chi/gom1.dat");
  for(size_t ii = 0; ii < 200; ++ii){
   double omega, real, imag;
   infile >> omega >> real >> imag;
   for(size_t is = 0; is < ham.model().spins(); ++is){
    G_file(alps::gf::generic_index<alps::gf::matsubara_mesh<alps::gf::mesh::POSITIVE_ONLY>>(ii), alps::gf::generic_index<alps::gf::index_mesh>(0), alps::gf::generic_index<alps::gf::index_mesh>(is)) = std::complex<double>(real, imag);
   }
  }
  infile.close();
  auto ChiSz_file = ChiSz;
  infile.open("test/input/GF_Chi/xiats.dat");
  for(size_t ii = 0; ii < 200; ++ii){
   double omega, real, imag;
   infile >> omega >> real >> imag;
   // S_z = 0.5 M, <S_z S_z> = 0.25 <M M>
   ChiSz_file(alps::gf::generic_index<alps::gf::matsubara_mesh<alps::gf::mesh::POSITIVE_ONLY>>(ii), alps::gf::generic_index<alps::gf::index_mesh>(0)) = std::complex<double>(-0.25 * real, 0.0);
  }
  infile.close();
  auto ChiN_file = ChiN;
  infile.open("test/input/GF_Chi/xiatd.dat");
  for(size_t ii = 0; ii < 200; ++ii){
   double omega, real, imag;
   infile >> omega >> real >> imag;
   ChiN_file(alps::gf::generic_index<alps::gf::matsubara_mesh<alps::gf::mesh::POSITIVE_ONLY>>(ii), alps::gf::generic_index<alps::gf::index_mesh>(0)) = std::complex<double>(-real, 0.0);
  }
  infile.close();

  // Subtract the reference GF from our result, the norm() is then the largest diff.
  G -= G_file;
  ASSERT_NEAR(G.norm(), 0.0, 1e-10);
  ChiSz -= ChiSz_file;
  ASSERT_NEAR(ChiSz.norm(), 0.0, 1e-9);
  ChiN -= ChiN_file;
  ASSERT_NEAR(ChiN.norm(), 0.0, 1e-9);

}
