// CUDA spin-resolved storage example: 4-site Hubbard ring at half filling.
//
// Same physics as HubbardCore.cpp, but the Hamiltonian matvec runs on the
// GPU (SpinResolvedStorageCuda + cpp-arnoldi CudaBackend). This TU is built
// by nvcc because the storage launches CUDA kernels.
//
// Build: configure EDLib with -DEDLIB_USE_CUDA=ON.

#include <edlib/GreensFunction.h>
#include <edlib/Hamiltonian.h>
#include <edlib/Mesh.h>
#include <edlib/StaticObservables.h>

#include <iomanip>
#include <iostream>

int main() {
  edlib::Parameters p;
  p.nsites         = 4;
  p.nspins         = 2;
  p.arpack_nev     = 1;
  p.lanc_beta      = 10.0;
  p.lanc_nomega    = 32;
  p.lanc_nlanc     = 100;

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

  using HamType = edlib::SRSCudaHubbardHamiltonian;
  HamType ham(p, bath);

  ham.diag();
  std::cout << std::setprecision(10);
  std::cout << "Ground state: "
            << ham.eigenpairs().begin()->eigenvalue() << std::endl;

  // Single-particle Green's function (Matsubara) — exercises the Lanczos
  // path through the GPU matvec.
  edlib::MatsubaraMesh fmesh(p.lanc_beta, p.lanc_nomega, edlib::Statistics::Fermionic);
  edlib::GreensFunction<HamType, edlib::MatsubaraMesh> gf(p, ham, fmesh);
  gf.compute();
  std::cout << "G(iw_0, orb=0, spin=0) = " << gf.G()(0, 0, 0) << std::endl;

  return 0;
}
