#ifndef EDLIB_LANCZOS_H
#define EDLIB_LANCZOS_H

#include <cmath>
#include <complex>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "edlib/Gf.h"
#include "edlib/Mesh.h"
#include "edlib/Parameters.h"

namespace edlib {

  /**
   * Lanczos basis builder + continued fraction evaluator for spectral
   * functions. Templated directly on the frequency mesh type (no MeshFactory
   * layer): callers pass a mesh instance into the constructor.
   */
  template <class Hamiltonian, class Mesh>
  class Lanczos {
  protected:
    using precision = typename Hamiltonian::ModelType::precision;

  public:
    Lanczos(const Parameters& p, Hamiltonian& h, Mesh omega)
        : ham(h),
          _omega(std::move(omega)),
          _beta(static_cast<precision>(p.lanc_beta)),
          _Nl(p.lanc_nlanc),
          alfalanc(p.lanc_nlanc, precision(0)),
          betalanc(p.lanc_nlanc + 1, precision(0)),
          det(p.lanc_nlanc, std::complex<double>(0)),
          dl (p.lanc_nlanc, std::complex<double>(0)) {}

    const Mesh& omega() const { return _omega; }

  protected:
    int lanczos(std::vector<precision>& v) {
      int nlanc = 0;
      const std::size_t size = v.size();
      std::vector<precision> w(size, precision(0));
      precision alf = 0, bet = 0;
      ham.fill();
      if (size != 0) {
        ham.storage().prepare_work_arrays(v.data());
        for (int iter = 1; iter <= _Nl; ++iter) {
          ++nlanc;
          if (iter != 1) {
            for (std::size_t j = 0; j < size; ++j) {
              precision dummy = v[j];
              v[j] = w[j] / bet;
              w[j] = -bet * dummy;
            }
          }
          alf = 0;
          bet = 0;
          ham.storage().av(v.data(), w.data(), size, false);
          alf = ham.storage().vv(v, w);
          alfalanc[iter - 1] = alf;
          for (std::size_t j = 0; j < size; ++j) w[j] -= alf * v[j];
          bet = ham.storage().vv(w, w);
          bet = std::sqrt(bet);
          if (iter != _Nl) betalanc[iter] = bet;
          if (std::abs(bet) < precision(1e-10)) break;
        }
        hamiltonian().storage().finalize(0, false);
      } else {
        hamiltonian().storage().finalize(0, false);
      }
#ifdef USE_MPI
      MPI_Barrier(hamiltonian().comm());
#endif
      return nlanc;
    }

    /// Continued-fraction evaluation into a 3-index GF (frequency × orbital × spin).
    void compute_continued_fraction(double expectation_value, double excited_state,
                                    double groundstate, int nlanc, int isign,
                                    GF3& gf, int site, int spin) {
      double expb = (_beta * (excited_state - groundstate) > 25)
                  ? 0.0
                  : std::exp(-_beta * (excited_state - groundstate));
      for (int iomega = 0; iomega < _omega.extent(); ++iomega) {
        std::complex<double> ener  = freq_point(iomega) + excited_state * isign;
        double               shift = 1.0;
        std::complex<double> swp   = get_frac_point(expectation_value, nlanc, isign, expb, shift, ener);
        gf(iomega, site, spin) += swp;
      }
    }

    /// Continued-fraction evaluation into a 2-index GF (frequency × orbital), symmetric variant.
    void compute_sym_continued_fraction(double expectation_value, double excited_state,
                                        double groundstate, int nlanc, int isign,
                                        GF2& gf, int site) {
      double expb = (_beta * (excited_state - groundstate) > 25)
                  ? 0.0
                  : std::exp(-_beta * (excited_state - groundstate));
      update_static(gf, site, expectation_value, expb);
      for (int iomega = zero_freq(); iomega < _omega.extent(); ++iomega) {
        std::complex<double> ener  =  freq_point(iomega) + excited_state * isign;
        std::complex<double> ener2 = -freq_point(iomega) + excited_state * isign;
        double               shift = 1.0;
        std::complex<double> swp   = get_frac_point(expectation_value, nlanc, isign, expb, shift, ener);
        swp                       += get_frac_point(expectation_value, nlanc, isign, expb, shift, ener2);
        gf(iomega, site) += swp;
      }
    }

          Hamiltonian& hamiltonian()       { return ham; }
    const Hamiltonian& hamiltonian() const { return ham; }

  public:
    /// Inverse temperature used by the spectral evaluation.
    precision beta() const { return _beta; }

  protected:

    std::complex<double> freq_point(int index) const {
      if constexpr (std::is_same_v<Mesh, RealFreqMesh>) {
        return std::complex<double>(_omega.points()[index], M_PI / _beta);
      } else {
        return std::complex<double>(0.0, _omega.points()[index]);
      }
    }

    std::string suffix() const {
      if constexpr (std::is_same_v<Mesh, RealFreqMesh>) return "_r";
      else                                              return "";
    }

    int zero_freq() const {
      if constexpr (std::is_same_v<Mesh, MatsubaraMesh>) return 1;
      else                                               return 0;
    }

    void update_static(GF2& gf, int site, double expectation_value, double expb) {
      if constexpr (std::is_same_v<Mesh, MatsubaraMesh>) {
        gf(0, site) -= std::complex<double>(expectation_value * _beta * expb, 0.0);
      }
    }

  private:
    Hamiltonian& ham;
    Mesh         _omega;
    precision    _beta;
    int          _Nl;

    std::vector<precision> alfalanc;
    std::vector<precision> betalanc;

    std::vector<std::complex<double>> det;
    std::vector<std::complex<double>> dl;

    std::complex<double> get_frac_point(double expectation_value, int nlanc, int isign,
                                        double expb, double shift,
                                        const std::complex<double>& ener) {
      std::complex<double> swp(0, 0);
      det.assign(nlanc, std::complex<double>(0));
      for (int i = 0; i < nlanc; ++i) {
        dl[i] = ener - (double(alfalanc[i]) * isign);
      }
      if (nlanc == 1) {
        det[0] = dl[0];
        swp += expectation_value * expb / det[0];
      } else {
        det[nlanc - 1] = dl[nlanc - 1];
        det[nlanc - 2] = dl[nlanc - 2] * dl[nlanc - 1] - std::pow(betalanc[nlanc - 1], 2);
        for (int i = nlanc - 3; i >= 0; --i) {
          det[i] = dl[i] * det[i + 1] - std::pow(betalanc[i + 1], 2) * det[i + 2];
          // Avoid overflow when intermediate determinants get huge.
          if (std::abs(det[i]) > (std::numeric_limits<float>::max() / 2.0) && i != 0) {
            shift = 1.0 / (std::numeric_limits<float>::max() / 1000.0);
            det[i]     *= shift;
            det[i + 1] *= shift;
          }
        }
        swp += expectation_value * expb * det[1] / det[0];
      }
      return swp;
    }
  };

}

#endif
