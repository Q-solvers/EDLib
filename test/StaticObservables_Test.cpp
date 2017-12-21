#include <gtest/gtest.h>
#include "edlib/Hamiltonian.h"
#include "edlib/HubbardModel.h"
#include "edlib/Storage.h"
#include "edlib/EDParams.h"
#include "edlib/StaticObservables.h"


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
  p["INPUT_FILE"]="test/input/4ring_nomagfield/input.h5";
  p["arpack.SECTOR"]=false;
  p["storage.MAX_SIZE"]=576;
  p["storage.MAX_DIM"]=36;
  p["storage.EIGENVALUES_ONLY"]=false;
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

  EDLib::StaticObservables<HamType> so(p);
  std::map<std::string, std::vector<double>> result = so.calculate_static_observables(ham);

  ASSERT_EQ(result["N"].size(), 4);
  ASSERT_EQ(result["N_up"].size(), 4);
  ASSERT_EQ(result["N_dn"].size(), 4);
  ASSERT_EQ(result["M"].size(), 4);
  ASSERT_EQ(result["D_occ"].size(), 4);
  ASSERT_EQ(result["N_eff"].size(), 1);
  ASSERT_EQ(result["M_i M_j"].size(), 16);

  for(int orb = 0; orb < ham.model().interacting_orbitals(); ++orb){
   ASSERT_NEAR(result["N"][orb], 1.0, 1e-8);
   ASSERT_NEAR(result["N_up"][orb], 0.5, 1e-8);
   ASSERT_NEAR(result["N_dn"][orb], 0.5, 1e-8);
   ASSERT_NEAR(result["M"][orb], 0.0, 1e-8);
  }

  for(int orb = 1; orb < ham.model().interacting_orbitals(); ++orb){
   ASSERT_NEAR(result["D_occ"][orb], result["D_occ"][0], 1e-8);
  }
  ASSERT_GT(result["N_eff"][0], 0.0);

}
