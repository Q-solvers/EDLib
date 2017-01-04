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
  ASSERT_NEAR(ham.eigenpairs().begin()->eigenvalue(), -11.8443, 1e-4);
  ASSERT_EQ(ham.eigenpairs().begin()->sector().nup(), 2);
  ASSERT_EQ(ham.eigenpairs().begin()->sector().ndown(), 2);

  while(ham.model().symmetry().next_sector()) {
    std::cout<<ham.model().symmetry().sector().size()<<std::endl;
  }
}
