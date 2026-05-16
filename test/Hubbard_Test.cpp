// Core-only Hubbard test. Mirrors the legacy HubbardModelTest 4-ring case
// using hand-built bath data (no HDF5).

#include <edlib/Hamiltonian.h>
#include <edlib/StaticObservables.h>

#include <gtest/gtest.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace {

edlib::Parameters make_params() {
  edlib::Parameters p;
  p.nsites             = 4;
  p.nspins             = 2;
  p.arpack_nev         = 1;
  p.storage_max_size   = 576;
  p.storage_max_dim    = 36;
  p.lanc_beta          = 10.0;
  p.lanc_boltzmann_cutoff = 1e-12;
  return p;
}

edlib::HubbardModel<double>::ModelData make_4ring_bath() {
  edlib::HubbardModel<double>::ModelData b;
  b.hopping = {
    { 0.0, -1.0,  0.0, -1.0},
    {-1.0,  0.0, -1.0,  0.0},
    { 0.0, -1.0,  0.0, -1.0},
    {-1.0,  0.0, -1.0,  0.0}
  };
  b.U              = {5.0, 5.0, 5.0, 5.0};
  b.mu             = {2.5, 2.5, 2.5, 2.5};
  b.magnetic_field = {0.01, 0.01, 0.01, 0.01};
  return b;
}

}

TEST(HubbardCore, ReferenceTest) {
  auto p    = make_params();
  auto bath = make_4ring_bath();

#ifdef USE_MPI
  using HamType = edlib::SRSHubbardHamiltonian;
  HamType ham(p, bath, MPI_COMM_WORLD);
#else
  using HamType = edlib::SOCSRHubbardHamiltonian;
  HamType ham(p, bath);
#endif

  ham.diag();

  // Ground state reference from arXiv:cond-mat/0101476
  ASSERT_NEAR(ham.eigenpairs().begin()->eigenvalue(), -11.8443, 1e-4);
  ASSERT_EQ(ham.eigenpairs().begin()->sector().nup(),   2);
  ASSERT_EQ(ham.eigenpairs().begin()->sector().ndown(), 2);

  edlib::StaticObservables<HamType> so(p);
  auto obs = so.calculate_static_observables(ham);
#ifdef USE_MPI
  for (auto& kv : obs)
    MPI_Bcast(kv.second.data(), kv.second.size(), MPI_DOUBLE, 0, ham.comm());
#endif

  for (int orb = 0; orb < ham.model().interacting_orbitals(); ++orb) {
    ASSERT_NEAR(obs[edlib::StaticObservables<HamType>::_N_][orb], 1.0, 1e-8);
    ASSERT_GT(obs[edlib::StaticObservables<HamType>::_N_UP_][orb],
              obs[edlib::StaticObservables<HamType>::_N_DN_][orb]);
    ASSERT_GT(obs[edlib::StaticObservables<HamType>::_M_][orb], 0.0);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif
  int res = RUN_ALL_TESTS();
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return res;
}
