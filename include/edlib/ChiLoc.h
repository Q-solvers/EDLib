#ifndef EDLIB_CHILOC_H
#define EDLIB_CHILOC_H

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "edlib/EigenPair.h"
#include "edlib/Gf.h"
#include "edlib/Lanczos.h"
#include "edlib/Mesh.h"
#include "edlib/Parameters.h"

namespace edlib {

  /**
   * Base for bosonic operators used in local susceptibility evaluation.
   */
  template <class Precision>
  class BosonicOperator {
  public:
    explicit BosonicOperator(Precision avg) : _avg(avg) {}
    Precision average() const { return _avg; }
  private:
    Precision _avg;
  };

  /// Sz operator: 0.5 * (n_up - n_down) on a site.
  template <class Precision>
  class SzOperator : public BosonicOperator<Precision> {
  public:
    explicit SzOperator(Precision avg = Precision(0))
        : BosonicOperator<Precision>(avg) {}

    template <class Model>
    Precision action(long long state, int site, const Model& model) const {
      return Precision(0.5) *
             (model.checkState(state, site,                     model.max_total_electrons()) -
              model.checkState(state, site + model.orbitals(),  model.max_total_electrons()));
    }
    std::string name() const { return "Sz"; }
  };

  /// Charge operator: n_up + n_down on a site.
  template <class Precision>
  class NOperator : public BosonicOperator<Precision> {
  public:
    explicit NOperator(Precision avg = Precision(1))
        : BosonicOperator<Precision>(avg) {}

    template <class Model>
    Precision action(long long state, int site, const Model& model) const {
      return model.checkState(state, site,                    model.max_total_electrons()) +
             model.checkState(state, site + model.orbitals(), model.max_total_electrons());
    }
    std::string name() const { return "N"; }
  };

  /**
   * Local susceptibility \chi_{ii}(\omega) and its non-local extension via
   * Lanczos continued fractions. Stores the result in GF2 = Gf<complex,2>:
   * shape [n_omega, n_orbital] for the diagonal and
   * [n_omega, n_orbital*n_orbital] for the non-local block.
   */
  template <class Hamiltonian, class Mesh>
  class ChiLoc : public Lanczos<Hamiltonian, Mesh> {
    using Base = Lanczos<Hamiltonian, Mesh>;
    using Base::zero_freq;
    using Base::lanczos;
    using Base::hamiltonian;
    using Base::compute_sym_continued_fraction;
    using typename Base::precision;
    using Sector  = typename Hamiltonian::ModelType::Sector;

  public:
    using Base::beta;
    using Base::omega;
    using GF_TYPE = GF2;

    ChiLoc(const Parameters& p, Hamiltonian& h, Mesh omega_mesh,
           const std::vector<std::array<int, 2>>& orbital_pairs = {})
        : Base(p, h, std::move(omega_mesh)),
          _model(h.model()),
          gf   ({omega().extent(), h.model().interacting_orbitals()}),
          gf_ij({omega().extent(),
                 h.model().interacting_orbitals() * h.model().interacting_orbitals()}),
          _cutoff(static_cast<precision>(p.lanc_boltzmann_cutoff)),
          _type("Sz") {
      if (p.eigenvalues_only) {
        throw std::logic_error(
            "ChiLoc: eigenvectors were not computed (eigenvalues_only=true).");
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
        for (const auto& pr : _g_ij_orb_pairs) {
          for (int e : pr) {
            if (std::find(_g_orbs.begin(), _g_orbs.end(), e) == _g_orbs.end()) {
              _g_orbs.push_back(e);
            }
          }
        }
      }
    }

    const GF_TYPE& G()    const { return gf; }
    const GF_TYPE& G_ij() const { return gf_ij; }

    template <class Op = SzOperator<precision>>
    void compute(const double* avg_ptr = nullptr) {
      static_assert(std::is_base_of<BosonicOperator<precision>, Op>::value,
                    "ChiLoc::compute: Op must derive from BosonicOperator");
      gf    *= std::complex<double>(0);
      gf_ij *= std::complex<double>(0);
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
      const Op op = (avg_ptr == nullptr ? Op() : Op(precision(*avg_ptr)));
      _type = op.name();

      for (const auto& pair : hamiltonian().eigenpairs()) {
        precision boltzmann_f =
            std::exp(-(pair.eigenvalue() - groundstate.eigenvalue()) * beta());
        if (std::abs(_cutoff - boltzmann_f) > std::numeric_limits<precision>::epsilon()
            && boltzmann_f < _cutoff) continue;
#ifdef USE_MPI
        if (rank == 0)
#endif
          std::cout << "Compute Green's function contribution for eigenvalue E="
                    << pair.eigenvalue() << " with Boltzmann factor = "
                    << boltzmann_f << "; for sector" << pair.sector() << std::endl;
        local_contribution   (groundstate, op, pair);
        nonlocal_contribution(groundstate, op, pair);
      }
#ifdef USE_MPI
      if (rank == 0) {
#endif
        gf    /= std::complex<double>(_Z);
        gf_ij /= std::complex<double>(_Z);
        local_correction(op);
        non_local_correction(op);
#ifdef USE_MPI
      }
#endif
    }

    template <class O>
    void local_correction(const O& op) {
      for (int orb : _g_orbs) zero_freq_contribution(op, gf, orb);
    }

    template <class O>
    void non_local_correction(const O& op) {
      const int n_orb = _model.interacting_orbitals();
      for (const auto& pr : _g_ij_orb_pairs) {
        zero_freq_contribution(op, gf_ij, pr[0] * n_orb + pr[1]);
      }
      for (int iomega = 0; iomega < omega().extent(); ++iomega) {
        for (const auto& pr : _g_ij_orb_pairs) {
          int ij = pr[0] * n_orb + pr[1];
          for (int jj = 0; jj < 2; ++jj) gf_ij(iomega, ij) -= gf(iomega, pr[jj]);
          gf_ij(iomega, ij) *= std::complex<double>(0.5);
        }
        for (int orb : _g_orbs) gf_ij(iomega, orb * n_orb + orb) = gf(iomega, orb);
      }
    }

    template <class Op = SzOperator<precision>>
    void local_contribution(const EigenPair<precision, Sector>& groundstate, const Op op,
                            const EigenPair<precision, Sector>& pair) {
#ifdef USE_MPI
      int rank; MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
#endif
      for (int orb : _g_orbs) {
        std::vector<precision> outvec(1, precision(0));
        precision expectation_value = 0;
        _model.symmetry().set_sector(pair.sector());
        if (operation(orb, pair.eigenvector(), outvec, expectation_value, op)) {
          int nlanc = lanczos(outvec);
#ifdef USE_MPI
          if (rank == 0)
#endif
          {
            std::cout << "orbital: " << orb << " <n|" << op.name() << op.name()
                      << "|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
            compute_sym_continued_fraction(expectation_value, pair.eigenvalue(),
                                           groundstate.eigenvalue(),
                                           nlanc, 1, gf, orb);
          }
        }
      }
    }

    template <class Op>
    void nonlocal_contribution(const EigenPair<precision, Sector>& groundstate, const Op& op,
                               const EigenPair<precision, Sector>& pair) {
#ifdef USE_MPI
      int rank; MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
#endif
      const int n_orb = _model.interacting_orbitals();
      for (const auto& pr : _g_ij_orb_pairs) {
        _model.symmetry().set_sector(pair.sector());
        std::vector<precision> outvec(1, precision(0));
        precision expectation_value = 0;
        if (operation(pr[0], pr[1], pair.eigenvector(), outvec, expectation_value, op)) {
          int nlanc = lanczos(outvec);
#ifdef USE_MPI
          if (rank == 0)
#endif
          {
            std::cout << "orbitals: " << pr[0] << ", " << pr[1] << " <n|"
                      << op.name() << op.name() << "|n>=" << expectation_value
                      << " nlanc:" << nlanc << std::endl;
            compute_sym_continued_fraction(expectation_value, pair.eigenvalue(),
                                           groundstate.eigenvalue(),
                                           nlanc, 1, gf_ij,
                                           pr[0] * n_orb + pr[1]);
          }
        }
      }
    }

  private:
    // Matsubara zero-frequency correction: removes the analytical tail of
    // Chi(W_n) ~ c2/W^2 + c4/W^4 and applies the sum rule.
    template <class O>
    void zero_freq_contribution(const O& op, GF_TYPE& G, int i) {
      if constexpr (std::is_same_v<Mesh, MatsubaraMesh>) {
        double chiSum = 0.0;
        double c2, c4, tail;
        double c2_2, c4_2, tail2;
        get_tail(i, omega().extent() - 1, G, c2,   c4,   tail);
        get_tail(i, omega().extent() - 2, G, c2_2, c4_2, tail2);
        if (std::abs(tail - tail2) / std::abs(tail) > 1e-4) {
          std::cerr << "Not enough frequencies to compute high frequency tail. "
                       "Please increase number of frequencies. Diff: "
                    << std::abs(tail - tail2) / std::abs(tail) << std::endl;
        }
        for (int iomega = 1; iomega < omega().extent(); ++iomega) {
          double om = omega().points()[iomega];
          chiSum += G(iomega, i).real() - c2 / (om * om) - c4 / (om * om * om * om);
        }
        G(0, i) -= std::complex<double>(2 * chiSum + 2 * tail
                                         - op.average() * op.average() * beta(),
                                         0.0);
      } else {
        (void)op; (void)G; (void)i;
      }
    }

    void get_tail(int i, int freq, const GF_TYPE& G,
                  double& c2, double& c4, double& tail) const {
      double om1 = omega().points()[freq];
      double om2 = omega().points()[freq - 1];
      double om1_2 = om1 * om1;
      double om2_2 = om2 * om2;
      double g1 = G(freq,     i).real();
      double g2 = G(freq - 1, i).real();
      c2   = -(g2 * om2_2 * om2_2 - g1 * om1_2 * om1_2) / (om1_2 - om2_2);
      c4   = -(g1 * om1_2 * om1_2 * om2_2 - g2 * om2_2 * om2_2 * om1_2) / (om1_2 - om2_2);
      tail =  c2 * beta() * beta() / 24.0
            + c4 * beta() * beta() * beta() * beta() / 1440.0;
    }

    template <class Op>
    bool operation(int orbital, const std::vector<precision>& invec,
                   std::vector<precision>& outvec, precision& expectation_value,
                   const Op& o) {
      hamiltonian().storage().reset();
      outvec.assign(hamiltonian().storage().vector_size(_model.symmetry().sector()),
                    precision(0));
      for (std::size_t i = 0; i < invec.size(); ++i) {
        _model.symmetry().next_state();
        long long nst = _model.symmetry().state();
        outvec[i] = o.action(nst, orbital, _model) * invec[i];
      }
      double norm = hamiltonian().storage().vv(outvec, outvec
#ifdef USE_MPI
          , hamiltonian().comm()
#endif
      );
      for (auto& v : outvec) v /= std::sqrt(norm);
      _model.symmetry().init();
      expectation_value = static_cast<precision>(norm);
      return expectation_value > precision(1e-9);
    }

    template <class Op>
    bool operation(int mu, int nu, const std::vector<precision>& invec,
                   std::vector<precision>& outvec, precision& expectation_value,
                   const Op& o) {
      hamiltonian().storage().reset();
      outvec.assign(hamiltonian().storage().vector_size(_model.symmetry().sector()),
                    precision(0));
      for (std::size_t i = 0; i < invec.size(); ++i) {
        _model.symmetry().next_state();
        long long nst = _model.symmetry().state();
        outvec[i] = (o.action(nst, mu, _model) + o.action(nst, nu, _model)) * invec[i];
      }
      double norm = hamiltonian().storage().vv(outvec, outvec
#ifdef USE_MPI
          , hamiltonian().comm()
#endif
      );
      for (auto& v : outvec) v /= std::sqrt(norm);
      _model.symmetry().init();
      expectation_value = static_cast<precision>(norm);
      return true;
    }

    typename Hamiltonian::ModelType& _model;
    GF_TYPE                          gf;
    GF_TYPE                          gf_ij;
    precision                        _cutoff;
    precision                        _Z = precision(0);
    std::string                      _type;
    std::vector<int>                 _g_orbs;
    std::vector<std::array<int, 2>>  _g_ij_orb_pairs;
  };

}

#endif
