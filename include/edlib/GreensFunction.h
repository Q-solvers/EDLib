#ifndef EDLIB_GREENSFUNCTION_H
#define EDLIB_GREENSFUNCTION_H

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "edlib/Dyson.h"
#include "edlib/EigenPair.h"
#include "edlib/ExecutionStatistic.h"
#include "edlib/Gf.h"
#include "edlib/Lanczos.h"
#include "edlib/Parameters.h"

namespace edlib {

  /**
   * Single-particle Green's function via Lanczos continued fractions.
   *
   * GF3 shape is [n_omega, n_orbital, n_spin] for the local GF, and
   * [n_omega, n_orbital * n_orbital, n_spin] for the non-local G_ij block.
   *
   * Orbital pairs to compute may be supplied at construction. If empty, all
   * diagonal local Green's functions are computed. Pairs (i, j) with i != j
   * additionally drive the non-local G_ij block.
   */
  template <class Hamiltonian, class Mesh>
  class GreensFunction : public Lanczos<Hamiltonian, Mesh> {
    using Base = Lanczos<Hamiltonian, Mesh>;
    using Base::hamiltonian;
    using Base::lanczos;
    using Base::kernel;
    using Base::compute_continued_fraction;
    using Base::suffix;
    using typename Base::precision;
    using typename Base::KVector;

  public:
    using Base::beta;
    using Base::omega;
    using ModelType = typename Hamiltonian::ModelType;
    using GF_TYPE   = GF3;

    GreensFunction(const Parameters& p, Hamiltonian& h, Mesh omega_mesh,
                   const std::vector<std::array<int, 2>>& orbital_pairs = {})
        : Base(p, h, std::move(omega_mesh)),
          _model(h.model()),
          _G_g  ({omega().extent(), h.model().interacting_orbitals(),
                  p.nspins}),
          _G_l  ({omega().extent(), h.model().interacting_orbitals(),
                  p.nspins}),
          _G    ({omega().extent(), h.model().interacting_orbitals(),
                  p.nspins}),
          _G_g_ij({omega().extent(),
                   h.model().interacting_orbitals() * h.model().interacting_orbitals(),
                   p.nspins}),
          _G_l_ij({omega().extent(),
                   h.model().interacting_orbitals() * h.model().interacting_orbitals(),
                   p.nspins}),
          _G_ij  ({omega().extent(),
                   h.model().interacting_orbitals() * h.model().interacting_orbitals(),
                   p.nspins}),
          _cutoff(static_cast<precision>(p.lanc_boltzmann_cutoff)),
          _nspins(p.nspins),
          _nsites(p.nsites) {
      if (p.eigenvalues_only) {
        throw std::logic_error(
            "GreensFunction: eigenvectors were not computed (eigenvalues_only=true).");
      }
      const int n_orb = h.model().interacting_orbitals();
      if (orbital_pairs.empty()) {
        for (int i = 0; i < n_orb; ++i) _g_orbs.push_back(i);
      } else {
        for (const auto& pr : orbital_pairs) {
          if (pr[0] == pr[1]) _g_orbs.push_back(pr[0]);
          else                _g_ij_orb_pairs.push_back({pr[0], pr[1]});
        }
        std::sort(_g_orbs.begin(), _g_orbs.end());
        _g_orbs.erase(std::unique(_g_orbs.begin(), _g_orbs.end()), _g_orbs.end());
        // Ensure each non-local pair has its two local GFs available.
        for (const auto& pr : _g_ij_orb_pairs) {
          for (int e : pr) {
            if (std::find(_g_orbs.begin(), _g_orbs.end(), e) == _g_orbs.end()) {
              _g_orbs.push_back(e);
            }
          }
        }
      }
    }

    /// Run the Lanczos continued-fraction sum over all eigenpairs.
    void compute() {
      _G_g *= std::complex<double>(0);
      _G_l *= std::complex<double>(0);
      _Z = precision(0);
      if (hamiltonian().eigenpairs().empty()) return;
#ifdef USE_MPI
      int rank;
      MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
#endif
      const auto& groundstate = *hamiltonian().eigenpairs().begin();
      for (const auto& pair : hamiltonian().eigenpairs()) {
        _Z += std::exp(-(pair.eigenvalue() - groundstate.eigenvalue()) * beta());
      }
      for (const auto& pair : hamiltonian().eigenpairs()) {
        precision boltzmann_f =
            std::exp(-(pair.eigenvalue() - groundstate.eigenvalue()) * beta());
        if (std::abs(_cutoff - boltzmann_f) > std::numeric_limits<precision>::epsilon()
            && boltzmann_f < _cutoff) {
          continue;
        }
#ifdef USE_MPI
        if (rank == 0)
#endif
          std::cout << "Compute Green's function contribution for eigenvalue E="
                    << pair.eigenvalue() << " with Boltzmann factor = "
                    << boltzmann_f << "; for sector" << pair.sector() << std::endl;
        local_contribution   (pair, groundstate);
        nonlocal_contribution(pair, groundstate);
      }
#ifdef USE_MPI
      if (rank == 0) {
#endif
        _G_g    /= std::complex<double>(_Z);
        _G_l    /= std::complex<double>(_Z);
        _G_g_ij /= std::complex<double>(_Z);
        _G_l_ij /= std::complex<double>(_Z);
#ifdef USE_MPI
      }
#endif
      non_local_gf();
      _G = _G_g + _G_l;
      edlib::statistics.updateEvent("Greens function");
    }

    /**
     * Compute the self-energy by Dyson's equation, using the model's bare
     * Green's function for the same mesh. Returns the GF3-shaped sigma.
     */
    GF_TYPE compute_selfenergy() {
      GF_TYPE bare ({omega().extent(), _nsites * _nsites, _nspins});
      GF_TYPE sigma({omega().extent(), _nsites * _nsites, _nspins});
      _model.bare_greens_function(bare, omega(), beta());
      solve_dyson(bare, _G_ij, sigma, _nsites);
      return sigma;
    }

    const GF_TYPE& G_g()    const { return _G_g; }
    const GF_TYPE& G_l()    const { return _G_l; }
    const GF_TYPE& G()      const { return _G; }
    const GF_TYPE& G_g_ij() const { return _G_g_ij; }
    const GF_TYPE& G_l_ij() const { return _G_l_ij; }
    const GF_TYPE& G_ij()   const { return _G_ij; }

  private:
    void local_contribution(const typename Hamiltonian::EigenPairType& pair,
                            const typename Hamiltonian::EigenPairType& groundstate) {
#ifdef USE_MPI
      int rank; MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
#endif
      for (std::size_t io = 0; io < _g_orbs.size(); ++io) {
        for (int ispin = 0; ispin < _model.spins(); ++ispin) {
          int orb = _g_orbs[io];
          KVector outvec = kernel().make_vector(1);
          precision expectation_value = 0;
          _model.symmetry().set_sector(pair.sector());
          if (create_particles(std::array<int, 1>{{orb}}, ispin,
                               pair.eigenvector(), outvec, expectation_value)) {
            int nlanc = lanczos(outvec);
#ifdef USE_MPI
            if (rank == 0)
#endif
            {
              std::cout << "orbital: " << orb << "   spin: "
                        << (ispin == 0 ? "up" : "down")
                        << " <n|aa*|n>=" << expectation_value
                        << " nlanc:" << nlanc << std::endl;
              compute_continued_fraction(expectation_value, pair.eigenvalue(),
                                         groundstate.eigenvalue(),
                                         nlanc, 1, _G_g, orb, ispin);
            }
          }
          _model.symmetry().set_sector(pair.sector());
          if (annihilate_particles(std::array<int, 1>{{orb}}, ispin,
                                   pair.eigenvector(), outvec, expectation_value)) {
            int nlanc = lanczos(outvec);
#ifdef USE_MPI
            if (rank == 0)
#endif
            {
              std::cout << "orbital: " << orb << "   spin: "
                        << (ispin == 0 ? "up" : "down")
                        << " <n|a*a|n>=" << expectation_value
                        << " nlanc:" << nlanc << std::endl;
              compute_continued_fraction(expectation_value, pair.eigenvalue(),
                                         groundstate.eigenvalue(),
                                         nlanc, -1, _G_l, orb, ispin);
            }
          }
        }
      }
    }

    void nonlocal_contribution(const typename Hamiltonian::EigenPairType& pair,
                               const typename Hamiltonian::EigenPairType& groundstate) {
#ifdef USE_MPI
      int rank; MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
#endif
      for (std::size_t io = 0; io < _g_ij_orb_pairs.size(); ++io) {
        for (int ispin = 0; ispin < _model.spins(); ++ispin) {
          auto orbs = _g_ij_orb_pairs[io];
          KVector outvec = kernel().make_vector(1);
          precision expectation_value = 0;
          _model.symmetry().set_sector(pair.sector());
          if (create_particles(std::array<int, 2>{{orbs[0], orbs[1]}}, ispin,
                               pair.eigenvector(), outvec, expectation_value)) {
            int nlanc = lanczos(outvec);
#ifdef USE_MPI
            if (rank == 0)
#endif
            {
              std::cout << "orbitals: " << orbs[0] << ", " << orbs[1]
                        << "   spin: " << (ispin == 0 ? "up" : "down")
                        << " <n|aa*|n>=" << expectation_value
                        << " nlanc:" << nlanc << std::endl;
              compute_continued_fraction(expectation_value, pair.eigenvalue(),
                                         groundstate.eigenvalue(), nlanc, 1, _G_g_ij,
                                         _model.interacting_orbitals() * orbs[0] + orbs[1],
                                         ispin);
            }
          }
          _model.symmetry().set_sector(pair.sector());
          if (annihilate_particles(std::array<int, 2>{{orbs[0], orbs[1]}}, ispin,
                                   pair.eigenvector(), outvec, expectation_value)) {
            int nlanc = lanczos(outvec);
#ifdef USE_MPI
            if (rank == 0)
#endif
            {
              std::cout << "orbitals: " << orbs[0] << ", " << orbs[1]
                        << "   spin: " << (ispin == 0 ? "up" : "down")
                        << " <n|a*a|n>=" << expectation_value
                        << " nlanc:" << nlanc << std::endl;
              compute_continued_fraction(expectation_value, pair.eigenvalue(),
                                         groundstate.eigenvalue(), nlanc, -1, _G_l_ij,
                                         _model.interacting_orbitals() * orbs[0] + orbs[1],
                                         ispin);
            }
          }
        }
      }
    }

    void non_local_gf() {
      const int n_orb = _model.interacting_orbitals();
      for (int iomega = 0; iomega < omega().extent(); ++iomega) {
        for (std::size_t io = 0; io < _g_ij_orb_pairs.size(); ++io) {
          for (int ispin = 0; ispin < _model.spins(); ++ispin) {
            auto orbs = _g_ij_orb_pairs[io];
            int ij = n_orb * orbs[0] + orbs[1];
            for (int jj = 0; jj < 2; ++jj) {
              _G_g_ij(iomega, ij, ispin) -= _G_g(iomega, orbs[jj], ispin);
              _G_l_ij(iomega, ij, ispin) -= _G_l(iomega, orbs[jj], ispin);
            }
            _G_g_ij(iomega, ij, ispin) *= std::complex<double>(0.5, 0.0);
            _G_l_ij(iomega, ij, ispin) *= std::complex<double>(0.5, 0.0);
          }
        }
        for (std::size_t io = 0; io < _g_orbs.size(); ++io) {
          int orb = _g_orbs[io];
          for (int ispin = 0; ispin < _model.spins(); ++ispin) {
            _G_g_ij(iomega, n_orb * orb + orb, ispin) = _G_g(iomega, orb, ispin);
            _G_l_ij(iomega, n_orb * orb + orb, ispin) = _G_l(iomega, orb, ispin);
          }
        }
      }
      _G_ij = _G_g_ij + _G_l_ij;
    }

    template <std::size_t N>
    bool create_particles(std::array<int, N> orbitals, int spin,
                          const KVector& invec,
                          KVector& outvec,
                          precision& expectation_value) {
      if (!_model.symmetry().can_create_particle(spin)) return false;
      hamiltonian().storage().reset();
      auto next_sec = _model.symmetry().create_particle(spin);
      outvec = kernel().make_vector(hamiltonian().storage().vector_size(next_sec));
      edlib::statistics.registerEvent("adag");
      for (int orb : orbitals) {
        hamiltonian().storage().init();
        KVector tmp = kernel().make_vector(kernel().size(outvec));
        kernel().a_adag(orb + spin * _model.orbitals(),
                        invec, tmp, next_sec, /*annihilate=*/false);
        kernel().add(outvec, tmp);
      }
      edlib::statistics.updateEvent("adag");
      double norm = kernel().dot(outvec, outvec
#ifdef USE_MPI
       , hamiltonian().comm()
#endif
      );
      kernel().scale(outvec, precision(1) / std::sqrt(norm));
      _model.symmetry().set_sector(next_sec);
      expectation_value = static_cast<precision>(norm);
      return std::abs(norm) > 1e-10;
    }

    template <std::size_t N>
    bool annihilate_particles(std::array<int, N> orbitals, int spin,
                              const KVector& invec,
                              KVector& outvec,
                              precision& expectation_value) {
      if (!_model.symmetry().can_destroy_particle(spin)) return false;
      hamiltonian().storage().reset();
      auto next_sec = _model.symmetry().destroy_particle(spin);
      outvec = kernel().make_vector(hamiltonian().storage().vector_size(next_sec));
      edlib::statistics.registerEvent("a");
      for (int orb : orbitals) {
        hamiltonian().storage().init();
        KVector tmp = kernel().make_vector(kernel().size(outvec));
        kernel().a_adag(orb + spin * _model.orbitals(),
                        invec, tmp, next_sec, /*annihilate=*/true);
        kernel().add(outvec, tmp);
      }
      edlib::statistics.updateEvent("a");
      double norm = kernel().dot(outvec, outvec
#ifdef USE_MPI
       , hamiltonian().comm()
#endif
      );
      kernel().scale(outvec, precision(1) / std::sqrt(norm));
      _model.symmetry().set_sector(next_sec);
      expectation_value = static_cast<precision>(norm);
      return std::abs(norm) > 1e-10;
    }

    ModelType&            _model;
    GF_TYPE               _G_g, _G_g_ij;
    GF_TYPE               _G_l, _G_l_ij;
    GF_TYPE               _G,   _G_ij;
    precision             _cutoff;
    precision             _Z = precision(0);
    int                   _nspins;
    int                   _nsites;
    std::vector<int>      _g_orbs;
    std::vector<std::array<int, 2>> _g_ij_orb_pairs;
  };

}

#endif
