// CPU vs GPU timing benchmark for the spin-resolved storage.
//
// Builds a single-impurity Anderson model (1 correlated impurity orbital +
// Nbath bath levels, total nsites = 1 + Nbath) and diagonalises the
// half-filling sector with:
//
//   * SpinResolvedStorage      (CPU,  cpp-arnoldi CpuBackend)
//   * SpinResolvedStorageCuda  (GPU,  cpp-arnoldi CudaBackend, CUDA matvec)
//
// A window of particle-number sectors around half filling is diagonalised
// (nup, ndown in {nhalf-1, nhalf, nhalf+1}) and the `nev` lowest states of
// each are kept, so the finite-temperature Green's function is a genuine
// Boltzmann sum over many eigenpairs across sectors -- not just the ground
// state. Both backends run the identical model; eigenvalues must agree.
// For each backend the wall time of ham.diag() (eigensolve over the sector
// window) AND of GreensFunction::compute() (one Lanczos continued fraction
// per contributing eigenpair x orbital x spin x {a, a+}) is reported.
//
// Build: configure EDLib with -DEDLIB_USE_CUDA=ON  (this TU is .cu/nvcc).
// Run:   ./anderson-cuda-bench [nsites] [nev]   (defaults: nsites=12 nev=6)
//
// Hilbert space of the half-filling sector is C(nsites, nsites/2)^2:
//   nsites=10 ->     63 504      nsites=14 ->  11 778 624
//   nsites=12 ->    853 776      nsites=16 -> 165 636 900

#include <edlib/GreensFunction.h>
#include <edlib/Hamiltonian.h>
#include <edlib/Mesh.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>

using prec = double;

static edlib::SingleImpurityAndersonModel<prec>::ModelData
make_anderson(int nsites, int nhalf, prec U) {
  const int ml  = 1;             // one correlated impurity orbital
  const int Nk  = nsites - ml;   // bath levels
  const int ns  = 2;             // spins

  edlib::SingleImpurityAndersonModel<prec>::ModelData d;

  // Hybridisation Vk[ml][Nk][ns] and bath energies Epsk[Nk][ns]:
  // a flat, symmetric bath around the Fermi level.
  d.Vk.assign(ml, std::vector<std::vector<prec>>(Nk, std::vector<prec>(ns, 0.0)));
  d.Epsk.assign(Nk, std::vector<prec>(ns, 0.0));
  for (int ik = 0; ik < Nk; ++ik) {
    prec eps = -2.0 + 4.0 * (ik + 0.5) / Nk;   // bath energies in [-2, 2]
    for (int is = 0; is < ns; ++is) {
      d.Epsk[ik][is]     = eps;
      d.Vk[0][ik][is]    = 0.3;
    }
  }

  // Non-interacting impurity level H0[ml][ml][ns] = 0 (particle-hole
  // symmetric together with mu = U/2).
  d.H0.assign(ml, std::vector<std::vector<prec>>(ml, std::vector<prec>(ns, 0.0)));
  d.mu = U / 2.0;

  // Interaction tensor U[ns,ns,ml,ml,ml,ml]; only the impurity density-
  // density term is non-zero (single-orbital Hubbard U on the impurity).
  d.U = edlib::Gf<prec, 6>(std::array<int, 6>{ns, ns, ml, ml, ml, ml});
  for (int is = 0; is < ns; ++is) d.U(is, is, 0, 0, 0, 0) = U;

  // Diagonalise a window of sectors around half filling so the thermal GF
  // sum spans several particle-number sectors (not just the ground state).
  for (int du = -1; du <= 1; ++du) {
    for (int dd = -1; dd <= 1; ++dd) {
      int nu = nhalf + du, nd = nhalf + dd;
      if (nu >= 0 && nu <= nsites && nd >= 0 && nd <= nsites)
        d.sectors.push_back({nu, nd});
    }
  }
  return d;
}

// Replicate GreensFunction's Boltzmann selection to report how many of the
// computed eigenpairs actually contribute to the thermal sum.
template <class Ham>
void report_ensemble(const Ham& ham, prec beta, prec cutoff) {
  const auto& eps = ham.eigenpairs();
  if (eps.empty()) return;
  double e0 = eps.begin()->eigenvalue();
  int total = 0, active = 0;
  for (const auto& pr : eps) {
    ++total;
    if (std::exp(-(pr.eigenvalue() - e0) * beta) >= cutoff) ++active;
  }
  std::printf("    eigenpairs computed = %d, contributing to GF = %d\n",
              total, active);
}

struct Timing {
  double diag_secs;
  double gf_secs;
};

template <class Ham>
Timing run_case(const char* tag, const edlib::Parameters& p,
                const typename Ham::Model::ModelData& d) {
  Ham ham(p, d);

  auto t0 = std::chrono::steady_clock::now();
  ham.diag();
  auto t1 = std::chrono::steady_clock::now();
  double diag_secs = std::chrono::duration<double>(t1 - t0).count();
  double e0 = ham.eigenpairs().begin()->eigenvalue();

  // Single-particle Matsubara Green's function. compute() runs the Lanczos
  // continued fraction, whose matvec goes through the storage (fully
  // device-resident for the CUDA case).
  edlib::MatsubaraMesh fmesh(p.lanc_beta, p.lanc_nomega, edlib::Statistics::Fermionic);
  edlib::GreensFunction<Ham, edlib::MatsubaraMesh> gf(p, ham, fmesh);
  auto t2 = std::chrono::steady_clock::now();
  gf.compute();
  auto t3 = std::chrono::steady_clock::now();
  double gf_secs = std::chrono::duration<double>(t3 - t2).count();

  std::cout << std::setprecision(10);
  std::cout << "[" << tag << "] ground state = " << e0
            << "   G(iw_0,0,0) = " << gf.G()(0, 0, 0) << "\n";
  report_ensemble(ham, static_cast<prec>(p.lanc_beta),
                  static_cast<prec>(p.lanc_boltzmann_cutoff));
  std::cout << "    diag() wall = " << std::setprecision(4) << diag_secs
            << " s    GF compute() wall = " << gf_secs << " s\n";
  return {diag_secs, gf_secs};
}

int main(int argc, char** argv) {
  int nsites = (argc > 1) ? std::atoi(argv[1]) : 12;
  int nev    = (argc > 2) ? std::atoi(argv[2]) : 6;
  int nhalf  = nsites / 2;
  const prec U = 4.0;

  edlib::Parameters p;
  p.nsites          = nsites;
  p.nspins          = 2;
  p.siam_norbitals  = 1;
  p.arpack_nev      = nev;     // several lowest states per sector
  p.arpack_sector   = true;
  p.lanc_beta       = 2.0;     // finite T: many eigenpairs clear the cutoff

  auto d = make_anderson(nsites, nhalf, U);

  std::printf("Single-impurity Anderson: nsites=%d (1 imp + %d bath), "
              "sector window around (%d,%d), nev=%d per sector, beta=%.1f\n",
              nsites, nsites - 1, nhalf, nhalf, nev, p.lanc_beta);

  Timing cpu = run_case<edlib::SRSSIAMHamiltonian>("CPU SpinResolvedStorage", p, d);
  Timing gpu = run_case<edlib::SRSCudaSIAMHamiltonian>("GPU SpinResolvedStorageCuda", p, d);

  std::printf("speedup  diag()       = %.2fx\n", cpu.diag_secs / gpu.diag_secs);
  std::printf("speedup  GF compute() = %.2fx\n", cpu.gf_secs   / gpu.gf_secs);
  return 0;
}
