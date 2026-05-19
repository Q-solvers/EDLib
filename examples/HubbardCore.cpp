// Header-only edlib:: example: 4-site Hubbard ring at half-filling.
//
// Uses ONLY the new edlib:: API (no ALPSCore). Mirrors the parameters of the
// legacy test/input/4ring/ HDF5 input, but constructs the bath in code so we
// don't need any file I/O.

#include <edlib/ChiLoc.h>
#include <edlib/GreensFunction.h>
#include <edlib/Hamiltonian.h>
#include <edlib/Mesh.h>
#include <edlib/StaticObservables.h>

#include <iomanip>
#include <iostream>

#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char** argv) {
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  int rank = 0;
#endif

  edlib::Parameters p;
  p.nsites             = 4;
  p.nspins             = 2;
  p.arpack_nev         = 1;
  p.storage_max_size   = 576;
  p.storage_max_dim    = 36;
  p.lanc_beta          = 10.0;
  p.lanc_nomega        = 32;
  p.lanc_nlanc         = 100;
  p.lanc_emin          = -3.0;
  p.lanc_emax          =  3.0;
  p.lanc_boltzmann_cutoff = 1e-12;

  edlib::HubbardModel<double>::ModelData bath;
  bath.hopping = {
    { 0.0, -1.0,  0.0, -1.0},
    {-1.0,  0.0, -1.0,  0.0},
    { 0.0, -1.0,  0.0, -1.0},
    {-1.0,  0.0, -1.0,  0.0}
  };
  bath.U              = {5.0, 5.0, 5.0, 5.0};
  bath.mu             = {2.5, 2.5, 2.5, 2.5};
  bath.magnetic_field = {0.01, 0.01, 0.01, 0.01};

#ifdef USE_MPI
  using HamType = edlib::SRSHubbardHamiltonian;
  HamType ham(p, bath, MPI_COMM_WORLD);
#else
  using HamType = edlib::SOCSRHubbardHamiltonian;
  HamType ham(p, bath);
#endif

  ham.diag();

  if (rank == 0) {
    std::cout << std::setprecision(10);
    std::cout << "Ground state: "
              << ham.eigenpairs().begin()->eigenvalue() << std::endl;
  }

  // Static observables.
  edlib::StaticObservables<HamType> so(p);
  auto obs = so.calculate_static_observables(ham);
  if (rank == 0) {
    std::cout << "<N> per orbital: ";
    for (double v : obs[edlib::StaticObservables<HamType>::_N_]) std::cout << v << " ";
    std::cout << std::endl;
  }

  // Single-particle Green's function (Matsubara).
  edlib::MatsubaraMesh fmesh(p.lanc_beta, p.lanc_nomega, edlib::Statistics::Fermionic);
  edlib::GreensFunction<HamType, edlib::MatsubaraMesh> gf(p, ham, fmesh);
  gf.compute();
  if (rank == 0) {
    std::cout << "G(iw_0, orb=0, spin=0) = " << gf.G()(0, 0, 0) << std::endl;
  }

  // Local Sz susceptibility (bosonic Matsubara).
  edlib::MatsubaraMesh bmesh(p.lanc_beta, p.lanc_nomega, edlib::Statistics::Bosonic);
  edlib::ChiLoc<HamType, edlib::MatsubaraMesh> chi(p, ham, bmesh);
  chi.compute<edlib::SzOperator<double>>();
  if (rank == 0) {
    std::cout << "ChiSz(iw_0, orb=0) = " << chi.G()(0, 0) << std::endl;
  }

#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
