//
// Created by iskakoff on 22/08/16.
//

#include <gtest/gtest.h>
#include "edlib/Hamiltonian.h"
#include "edlib/HubbardModel.h"
#include "edlib/Storage.h"
#include "edlib/EDParams.h"


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
  p["INPUT_FILE"]="test/input/4ring/input.h5";
  p["arpack.SECTOR"]=false;
  p["storage.MAX_SIZE"]=576;
  p["storage.MAX_DIM"]=36;
  p["storage.EIGENVALUES_ONLY"]=true;
  p["storage.ORBITAL_NUMBER"]=1;
  p["arpack.NEV"]=1;

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

  // [arXiv:cond-mat/0101476 [cond-mat.str-el]]
  double ref[25][3]={
    {-11.8443, 2, 2},
    {-11.5336, 3, 1},
    {-11.4936, 1, 3},
    {-10.1129, 2, 1},
    {-10.1129, 3, 2},
    {-10.0929, 2, 3},
    {-10.0929, 1, 2},
    {-10.04, 4, 0},
    {-9.96, 0, 4},
    {-9.53, 3, 0},
    {-9.53, 4, 1},
    {-9.47, 0, 3},
    {-9.47, 1, 4},
    {-8.34789, 3, 3},
    {-8.34789, 1, 1},
    {-7.02, 2, 0},
    {-7.02, 4, 2},
    {-6.98, 2, 4},
    {-6.98, 0, 2},
    {-4.51, 1, 0},
    {-4.51, 4, 3},
    {-4.49, 0, 1},
    {-4.49, 3, 4},
    {0, 0, 0},
    {0, 4, 4},
  };

  size_t i = 0;

  for(auto pair = ham.eigenpairs().begin(); pair != ham.eigenpairs().end(); ++pair){
    ASSERT_NEAR(pair->eigenvalue(), ref[i][0], 1e-3);
    ASSERT_EQ(pair->sector().nup(), ref[i][1]);
    ASSERT_EQ(pair->sector().ndown(), ref[i][2]);
    i++;
  }

  while(ham.model().symmetry().next_sector()) {
    std::cout<<ham.model().symmetry().sector().size()<<std::endl;
  }
}
