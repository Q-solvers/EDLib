// Core-only MPI test. Exercises the alpscore-free edlib:: API on the
// MPI-distributed SpinResolvedStorage Hamiltonian and asserts the parallel
// result is consistent across all ranks and matches the serial reference.
//
// Only meaningful when built with USE_MPI and launched under mpirun with
// >= 1 rank; with 2+ ranks it catches storage-distribution / reduction bugs
// that a single-process run cannot.

#include <edlib/Hamiltonian.h>

#include <gtest/gtest.h>

#ifdef USE_MPI
#include <mpi.h>

namespace {

edlib::Parameters make_params() {
  edlib::Parameters p;
  p.nsites                = 4;
  p.nspins                = 2;
  p.arpack_nev            = 1;
  p.storage_max_size      = 576;
  p.storage_max_dim       = 36;
  p.lanc_beta             = 10.0;
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

}  // namespace

// Distributed diagonalization must give every rank the same ground state,
// equal to the serial reference (arXiv:cond-mat/0101476, 4-ring, U=5).
TEST(MpiCore, DistributedGroundStateIsRankConsistent) {
  auto p    = make_params();
  auto bath = make_4ring_bath();

  edlib::SRSHubbardHamiltonian ham(p, bath, MPI_COMM_WORLD);
  ham.diag();

  ASSERT_FALSE(ham.eigenpairs().empty());
  const auto& gp = *ham.eigenpairs().begin();

  int size = 0;
  MPI_Comm_size(ham.comm(), &size);

  // 1. Reference value, checked locally on every rank.
  EXPECT_NEAR(gp.eigenvalue(), -11.8443, 1e-4);
  EXPECT_EQ(gp.sector().nup(),   2);
  EXPECT_EQ(gp.sector().ndown(), 2);

  // 2. Cross-rank consistency: the spread of the ground-state energy across
  //    all ranks must be zero (bit-for-bit; the reduced scalar is replicated).
  const double local_gs = gp.eigenvalue();
  double gs_min = 0.0, gs_max = 0.0;
  MPI_Allreduce(&local_gs, &gs_min, 1, MPI_DOUBLE, MPI_MIN, ham.comm());
  MPI_Allreduce(&local_gs, &gs_max, 1, MPI_DOUBLE, MPI_MAX, ham.comm());
  EXPECT_DOUBLE_EQ(gs_min, gs_max)
      << "ground-state energy differs across " << size << " ranks";

  // 3. Sector identification must also agree on every rank.
  const int local_nup = gp.sector().nup();
  int nup_min = 0, nup_max = 0;
  MPI_Allreduce(&local_nup, &nup_min, 1, MPI_INT, MPI_MIN, ham.comm());
  MPI_Allreduce(&local_nup, &nup_max, 1, MPI_INT, MPI_MAX, ham.comm());
  EXPECT_EQ(nup_min, nup_max);
}

#else  // !USE_MPI

TEST(MpiCore, SkippedWithoutMpi) {
  GTEST_SKIP() << "built without USE_MPI";
}

#endif

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
