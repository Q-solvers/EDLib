// Core-only Lanczos / GF / Chi test. Mirrors legacy LanczosTest's reference
// comparison against test/input/GF_Chi/{gom1,xiats,xiatd}.dat using only the
// new edlib:: API.

#include <edlib/ChiLoc.h>
#include <edlib/GreensFunction.h>
#include <edlib/Hamiltonian.h>
#include <edlib/Mesh.h>
#include <edlib/StaticObservables.h>

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <fstream>
#include <vector>

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace {

double frob_diff(const std::vector<std::complex<double>>& a,
                 const std::vector<std::complex<double>>& b) {
  double s = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    auto d = a[i] - b[i];
    s += d.real() * d.real() + d.imag() * d.imag();
  }
  return std::sqrt(s);
}

}

TEST(LanczosCore, ReferenceGFAndChi) {
  edlib::Parameters p;
  p.nsites             = 4;
  p.nspins             = 2;
  p.arpack_nev         = 100;
  p.storage_max_size   = 864;
  p.storage_max_dim    = 36;
  p.lanc_beta          = 20.0;
  p.lanc_nomega        = 1000;
  p.lanc_nlanc         = 100;
  p.lanc_emin          = -4.0;
  p.lanc_emax          =  4.0;
  p.lanc_boltzmann_cutoff = 1e-12;

  edlib::HubbardModel<double>::ModelData bath;
  bath.hopping = {
    { 0.0,  1.0,  1.0, -0.3},
    { 1.0,  0.0, -0.3,  1.0},
    { 1.0, -0.3,  0.0,  1.0},
    {-0.3,  1.0,  1.0,  0.0}
  };
  bath.U  = {3.0, 3.0, 3.0, 3.0};
  bath.mu = {1.0, 1.0, 1.0, 1.0};

#ifdef USE_MPI
  using HamType = edlib::SRSHubbardHamiltonian;
  HamType ham(p, bath, MPI_COMM_WORLD);
#else
  using HamType = edlib::CSRHubbardHamiltonian;
  HamType ham(p, bath);
#endif
  ham.diag();

  // Single-particle GF on orbital 0
  edlib::MatsubaraMesh fmesh(p.lanc_beta, p.lanc_nomega, edlib::Statistics::Fermionic);
  edlib::GreensFunction<HamType, edlib::MatsubaraMesh> gf(p, ham, fmesh, {{0, 0}});
  gf.compute();

  edlib::StaticObservables<HamType> so(p);
  auto obs = so.calculate_static_observables(ham);
#ifdef USE_MPI
  for (auto& kv : obs)
    MPI_Bcast(kv.second.data(), kv.second.size(), MPI_DOUBLE, 0, ham.comm());
#endif

  double avg_M = 0;
  for (double m : obs[edlib::StaticObservables<HamType>::_M_]) {
    avg_M += m / (2.0 * obs[edlib::StaticObservables<HamType>::_M_].size());
  }
  double avg_N = 0;
  for (double n : obs[edlib::StaticObservables<HamType>::_N_]) {
    avg_N += n / obs[edlib::StaticObservables<HamType>::_N_].size();
  }

  edlib::MatsubaraMesh bmesh(p.lanc_beta, p.lanc_nomega, edlib::Statistics::Bosonic);
  edlib::ChiLoc<HamType, edlib::MatsubaraMesh> susc(p, ham, bmesh, {{0, 0}});
  susc.compute<edlib::SzOperator<double>>(&avg_M);
  auto chiSz = susc.G();
  susc.compute<edlib::NOperator<double>>(&avg_N);
  auto chiN = susc.G();

#ifdef USE_MPI
  int rank; MPI_Comm_rank(ham.comm(), &rank);
  if (rank != 0) return;
#endif

  const std::string root = "test/input/GF_Chi/";

  // Compare 200 frequencies × 2 spins for G
  std::vector<std::complex<double>> G_ours(200 * 2), G_ref(200 * 2);
  {
    std::ifstream f(root + "gom1.dat");
    ASSERT_TRUE(f.is_open()) << "Cannot open " << root << "gom1.dat";
    for (int ii = 0; ii < 200; ++ii) {
      double w, r, i; f >> w >> r >> i;
      for (int is = 0; is < 2; ++is) {
        G_ref [ii * 2 + is] = {r, i};
        G_ours[ii * 2 + is] = gf.G()(ii, 0, is);
      }
    }
  }
  EXPECT_LT(frob_diff(G_ours, G_ref), 5e-10);

  std::vector<std::complex<double>> Cs_ours(200), Cs_ref(200);
  {
    std::ifstream f(root + "xiats.dat");
    ASSERT_TRUE(f.is_open());
    for (int ii = 0; ii < 200; ++ii) {
      double w, r, i; f >> w >> r >> i;
      Cs_ref [ii] = {-0.25 * r, 0.0};   // legacy convention <Sz Sz> = 0.25 <M M>
      Cs_ours[ii] = chiSz(ii, 0);
    }
  }
  EXPECT_LT(frob_diff(Cs_ours, Cs_ref), 1e-9);

  std::vector<std::complex<double>> Cn_ours(200), Cn_ref(200);
  {
    std::ifstream f(root + "xiatd.dat");
    ASSERT_TRUE(f.is_open());
    for (int ii = 0; ii < 200; ++ii) {
      double w, r, i; f >> w >> r >> i;
      Cn_ref [ii] = {-r, 0.0};
      Cn_ours[ii] = chiN(ii, 0);
    }
  }
  EXPECT_LT(frob_diff(Cn_ours, Cn_ref), 1e-9);
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
