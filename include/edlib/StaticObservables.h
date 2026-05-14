#ifndef EDLIB_STATICOBSERVABLES_H
#define EDLIB_STATICOBSERVABLES_H

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "edlib/EigenPair.h"
#include "edlib/MpiTypes.h"
#include "edlib/Parameters.h"

namespace edlib {

  /**
   * Compute thermal-averaged static observables (N, N_up, N_dn, M, D_occ,
   * M_i M_j, N_eff, E) from a diagonalized Hamiltonian.
   */
  template <class Hamiltonian>
  class StaticObservables {
  public:
    using precision = typename Hamiltonian::ModelType::precision;
    using sector    = typename Hamiltonian::ModelType::Sector;

    static constexpr const char* _N_     = "N";
    static constexpr const char* _N_UP_  = "N_up";
    static constexpr const char* _N_DN_  = "N_dn";
    static constexpr const char* _M_     = "M";
    static constexpr const char* _D_OCC_ = "D_occ";
    static constexpr const char* _N_EFF_ = "N_eff";
    static constexpr const char* _MI_MJ_ = "M_i M_j";
    static constexpr const char* _E_     = "E";

    explicit StaticObservables(const Parameters& p)
        : _beta(static_cast<precision>(p.lanc_beta)),
          _cutoff(static_cast<precision>(p.lanc_boltzmann_cutoff)) {
      if (p.eigenvalues_only) {
        throw std::logic_error(
            "StaticObservables: eigenvectors were not computed (eigenvalues_only=true).");
      }
#ifdef USE_MPI
      const int nitems = 2;
      int          blocklengths[nitems] = {1, 1};
      MPI_Datatype types[nitems] = { mpi_type<std::size_t>(), mpi_type<double>() };
      MPI_Aint     offsets[nitems];
      offsets[0] = offsetof(Element, ind);
      offsets[1] = offsetof(Element, val);
      MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_Element);
      MPI_Type_commit(&mpi_Element);
#endif
    }

    void print_static_observables(Hamiltonian& ham, std::ostream& out) {
      auto obs = calculate_static_observables(ham);
#ifdef USE_MPI
      int myid; MPI_Comm_rank(ham.comm(), &myid);
      if (myid != 0) return;
#endif
      for (const auto& kv : obs) {
        out << "<" << kv.first << "> = {";
        for (std::size_t i = 0; i < kv.second.size(); ++i) {
          if (i) out << ", ";
          out << kv.second[i];
        }
        out << "}\n";
      }
    }

    std::map<std::string, std::vector<precision>>
    calculate_static_observables(Hamiltonian& ham) {
#ifdef USE_MPI
      int myid; MPI_Comm_rank(ham.storage().comm(), &myid);
#endif
      const int n_orb = ham.model().interacting_orbitals();
      std::map<std::string, std::vector<precision>> avg = {
          {_N_,     std::vector<precision>(n_orb, precision(0))},
          {_N_UP_,  std::vector<precision>(n_orb, precision(0))},
          {_N_DN_,  std::vector<precision>(n_orb, precision(0))},
          {_M_,     std::vector<precision>(n_orb, precision(0))},
          {_D_OCC_, std::vector<precision>(n_orb, precision(0))},
          {_MI_MJ_, std::vector<precision>(n_orb * n_orb, precision(0))},
          {_N_EFF_, std::vector<precision>(1, precision(0))},
          {_E_,     std::vector<precision>(1, precision(0))},
      };
      precision sum = precision(0);
      const auto& groundstate = *ham.eigenpairs().begin();
      for (const auto& pair : ham.eigenpairs()) {
        precision boltzmann_f =
            std::exp(-(pair.eigenvalue() - groundstate.eigenvalue()) * _beta);
        if (std::abs(_cutoff - boltzmann_f) > std::numeric_limits<precision>::epsilon()
            && boltzmann_f < _cutoff) continue;
#ifdef USE_MPI
        if (myid == 0)
#endif
          std::cout << "Compute static observables contribution for eigenvalue E="
                    << pair.eigenvalue() << " with Boltzmann factor = "
                    << boltzmann_f << "; for sector" << pair.sector() << std::endl;
        auto contrib = calculate_static_observables_eigenvector(ham, pair);
        for (auto& kv : contrib) {
          for (std::size_t i = 0; i < kv.second.size(); ++i) {
            avg[kv.first][i] += kv.second[i] * boltzmann_f;
          }
        }
        sum += boltzmann_f;
      }
#ifdef USE_MPI
      if (myid == 0)
#endif
        std::cout << "Statsum: " << sum << std::endl;
      for (auto& kv : avg) for (auto& v : kv.second) v /= sum;
      return avg;
    }

    std::vector<std::pair<long long, precision>>
    find_largest_coefficients(Hamiltonian& ham,
                              const EigenPair<precision, sector>& pair,
                              std::size_t nmax, precision /*trivial*/) {
      ham.model().symmetry().set_sector(pair.sector());
      ham.storage().reset();
      int count = static_cast<int>(std::min(nmax, pair.eigenvector().size()));
      std::vector<std::size_t> largest(pair.eigenvector().size());
#ifdef USE_MPI
      int myid, nprocs;
      MPI_Comm_rank(ham.comm(), &myid);
      MPI_Comm_size(ham.comm(), &nprocs);
      std::vector<int> counts(nprocs);
      std::vector<int> displs(nprocs + 1);
      MPI_Gather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, ham.comm());
      if (myid == 0) {
        displs[0] = 0;
        for (int i = 0; i < nprocs; ++i) displs[i + 1] = displs[i] + counts[i];
      }
#endif
      for (std::size_t i = 0; i < largest.size(); ++i) largest[i] = i;
      std::partial_sort(largest.begin(), largest.begin() + count, largest.end(),
                        [&pair](std::size_t a, std::size_t b) {
                          return std::abs(pair.eigenvector()[a]) > std::abs(pair.eigenvector()[b])
                              || (std::abs(pair.eigenvector()[a]) == std::abs(pair.eigenvector()[b])
                                  && a < b);
                        });
#ifdef USE_MPI
      std::vector<Element> send(count);
      for (int i = 0; i < count; ++i) {
        send[i] = Element(largest[i] + ham.storage().offset(),
                          pair.eigenvector()[largest[i]]);
      }
      std::vector<Element> all(displs[nprocs]);
      MPI_Gatherv(send.data(), count, mpi_Element,
                  all.data(), counts.data(), displs.data(), mpi_Element,
                  0, ham.comm());
      if (myid == 0) {
        nmax = std::min(nmax, all.size());
        std::partial_sort(all.begin(), all.begin() + nmax, all.end(),
                          [](const Element& a, const Element& b) { return a > b; });
        std::vector<std::pair<long long, precision>> ret(nmax);
        for (std::size_t i = 0; i < nmax; ++i) {
          long long nst = ham.model().symmetry().state_by_index(all[i].ind);
          ret[i] = {nst, static_cast<precision>(all[i].val)};
        }
        return ret;
      }
      return {};
#else
      std::vector<std::pair<long long, precision>> ret(count);
      for (int i = 0; i < count; ++i) {
        long long nst = ham.model().symmetry().state_by_index(largest[i]);
        ret[i] = {nst, pair.eigenvector()[largest[i]]};
      }
      return ret;
#endif
    }

    std::vector<std::pair<std::size_t, precision>>
    calculate_class_contrib(Hamiltonian& ham, const EigenPair<precision, sector>& pair,
                            std::size_t nmax, precision trivial, bool cumulative) {
      auto coeffs = find_largest_coefficients(ham, pair, nmax, trivial);
      std::vector<std::pair<std::size_t, precision>> contribs;
#ifdef USE_MPI
      int myid; MPI_Comm_rank(ham.comm(), &myid);
      if (myid != 0) return contribs;
#endif
      for (std::size_t i = 0; i < coeffs.size(); ++i) {
        if (i == 0 || std::abs(coeffs[i - 1].second - coeffs[i].second) > trivial) {
          contribs.push_back({i, coeffs[i].second * coeffs[i].second});
          if (i && cumulative) contribs.back().second += contribs[contribs.size() - 2].second;
        } else {
          contribs.back().second += coeffs[i].second * coeffs[i].second;
        }
      }
      if (!contribs.empty()) contribs.pop_back();
      return contribs;
    }

    void print_major_electronic_configuration(Hamiltonian& ham,
                                              const EigenPair<precision, sector>& pair,
                                              std::size_t nmax, precision trivial,
                                              std::ostream& out) {
      auto coeffs = find_largest_coefficients(ham, pair, nmax, trivial);
#ifdef USE_MPI
      int myid; MPI_Comm_rank(ham.comm(), &myid);
      if (myid != 0) return;
#endif
      out << "Eigenvector components for eigenvalue " << pair.eigenvalue() << " ";
      pair.sector().print(out);
      out << std::endl;
      for (const auto& c : coeffs) {
        out << c.second << " * |";
        std::string spin_down = std::bitset<64>(c.first).to_string().substr(
            64 - ham.model().orbitals(), ham.model().orbitals());
        std::string spin_up = std::bitset<64>(c.first).to_string().substr(
            64 - 2 * ham.model().orbitals(), ham.model().orbitals());
        out << spin_up << "|" << spin_down << ">\n";
      }
    }

    void print_class_contrib(Hamiltonian& ham, const EigenPair<precision, sector>& pair,
                             std::size_t nmax, precision trivial, bool cumulative,
                             std::ostream& out) {
      auto contribs = calculate_class_contrib(ham, pair, nmax, trivial, cumulative);
#ifdef USE_MPI
      int myid; MPI_Comm_rank(ham.comm(), &myid);
      if (myid != 0) return;
#endif
      out << "Contributions of eigenvector component classes for eigenvalue "
          << pair.eigenvalue() << " ";
      pair.sector().print(out);
      out << std::endl;
      for (const auto& c : contribs) out << c.first << "\t" << c.second << "\n";
    }

#ifdef USE_MPI
    struct Element {
      Element() = default;
      Element(std::size_t i, double v) : ind(i), val(v) {}
      std::size_t ind;
      double      val;
      bool operator>(const Element& e) const {
        return std::abs(val) > std::abs(e.val)
            || (std::abs(val) == std::abs(e.val) && ind < e.ind);
      }
      bool operator<(const Element& e) const {
        return std::abs(val) < std::abs(e.val)
            || (std::abs(val) == std::abs(e.val) && ind > e.ind);
      }
    };
    MPI_Datatype mpi_Element;
#endif

  private:
    std::map<std::string, std::vector<precision>>
    calculate_static_observables_eigenvector(Hamiltonian& ham,
                                             const EigenPair<precision, sector>& pair) {
      const int n_orb = ham.model().interacting_orbitals();
      std::vector<precision> n     (n_orb, precision(0));
      std::vector<precision> n_up  (n_orb, precision(0));
      std::vector<precision> n_down(n_orb, precision(0));
      std::vector<precision> m     (n_orb, precision(0));
      std::vector<precision> mimj  (n_orb * n_orb, precision(0));
      std::vector<precision> d_occ (n_orb, precision(0));
      precision inverse_N_eff = precision(0);

      ham.model().symmetry().set_sector(pair.sector());
      ham.storage().reset();
      for (std::size_t i = 0; i < pair.eigenvector().size(); ++i) {
        precision weight = pair.eigenvector()[i] * pair.eigenvector()[i];
        ham.model().symmetry().next_state();
        long long nst = ham.model().symmetry().state();
        for (int orb = 0; orb < n_orb; ++orb) {
          int el_up   = ham.model().checkState(nst, orb,                       ham.model().max_total_electrons());
          int el_down = ham.model().checkState(nst, orb + ham.model().orbitals(),
                                               ham.model().max_total_electrons());
          n     [orb] += (el_up + el_down) * weight;
          n_up  [orb] += el_up   * weight;
          n_down[orb] += el_down * weight;
          m     [orb] += (el_up - el_down) * weight;
          for (int orb2 = 0; orb2 < n_orb; ++orb2) {
            int el_up2   = ham.model().checkState(nst, orb2, ham.model().max_total_electrons());
            int el_down2 = ham.model().checkState(nst, orb2 + n_orb,
                                                  ham.model().max_total_electrons());
            mimj[n_orb * orb + orb2] += (el_up - el_down) * (el_up2 - el_down2) * weight;
          }
          d_occ[orb] += el_up * el_down * weight;
        }
        inverse_N_eff += weight * weight;
      }

      std::map<std::string, std::vector<precision>> result = {
          {_N_,     std::vector<precision>(n_orb, precision(0))},
          {_N_UP_,  std::vector<precision>(n_orb, precision(0))},
          {_N_DN_,  std::vector<precision>(n_orb, precision(0))},
          {_M_,     std::vector<precision>(n_orb, precision(0))},
          {_D_OCC_, std::vector<precision>(n_orb, precision(0))},
          {_MI_MJ_, std::vector<precision>(n_orb * n_orb, precision(0))},
          {_N_EFF_, std::vector<precision>(1, precision(0))},
          {_E_,     std::vector<precision>(1, precision(0))},
      };
#ifdef USE_MPI
      MPI_Reduce(n     .data(), result[_N_]    .data(), n     .size(), mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(n_up  .data(), result[_N_UP_] .data(), n_up  .size(), mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(n_down.data(), result[_N_DN_] .data(), n_down.size(), mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(m     .data(), result[_M_]    .data(), m     .size(), mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(d_occ .data(), result[_D_OCC_].data(), d_occ .size(), mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(mimj  .data(), result[_MI_MJ_].data(), mimj  .size(), mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(&inverse_N_eff, &result[_N_EFF_][0], 1, mpi_type<precision>(), MPI_SUM, 0, ham.comm());
#else
      result[_N_]     = n;
      result[_N_UP_]  = n_up;
      result[_N_DN_]  = n_down;
      result[_M_]     = m;
      result[_D_OCC_] = d_occ;
      result[_MI_MJ_] = mimj;
      result[_N_EFF_][0] = inverse_N_eff;
#endif
      result[_E_][0] = pair.eigenvalue();
      result[_N_EFF_][0] = precision(1) / result[_N_EFF_][0];
      return result;
    }

    precision _beta;
    precision _cutoff;
  };

}

#endif
