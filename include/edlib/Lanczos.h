//
// Created by iskakoff on 01/08/16.
//

#ifndef HUBBARD_LANCZOS_H
#define HUBBARD_LANCZOS_H

#include <alps/gf/mesh.hpp>
#include <alps/gf/gf.hpp>
#include <alps/params.hpp>

#include <cmath>

namespace EDLib {
  namespace gf {
    template<class Hamiltonian>
    class Lanczos {
    protected:
      typedef typename Hamiltonian::ModelType::precision precision;
    public:
      Lanczos(alps::params &p, Hamiltonian &h, const alps::gf::statistics::statistics_type& type = alps::gf::statistics::statistics_type::FERMIONIC) :
        ham(h), _omega(p["lanc.BETA"].as<double>(), p["lanc.NOMEGA"].as<int>(), type),_Nl(p["lanc.NLANC"]),
        alfalanc(p["lanc.NLANC"], 0.0), betalanc(int(p["lanc.NLANC"]) + 1, 0.0), det(p["lanc.NLANC"], 0), dl(p["lanc.NLANC"], 0.0) {}

      const alps::gf::matsubara_positive_mesh &omega() const {
        return _omega;
      }
    protected:
      int lanczos(std::vector < precision > &v) {
        int nlanc = 0;
        unsigned long size = v.size();
        std::vector < precision > w(size, precision(0.0));
        precision alf = 0, bet = 0;
        ham.fill();
        if(v.size()!=0) {
          ham.storage().prepare_work_arrays(v.data());
          for (int iter = 1; iter <= _Nl; ++iter) {
            ++nlanc;
            if (iter != 1) {
              for (int j = 0; j < size; ++j) {
                precision dummy = v[j];
                v[j] = w[j] / bet;
                w[j] = -bet * dummy;
              }
            }
            alf = 0.0;
            bet = 0.0;
            ham.storage().av(v.data(), w.data(), size, false);
            alf = ham.storage().vv(v, w);
            alfalanc[iter - 1] = alf;
            for (int j = 0; j < size; ++j) {
              w[j] -= alf * v[j];
            }
            bet = ham.storage().vv(w, w);
            bet = std::sqrt(bet);

            if (iter != _Nl) betalanc[iter] = bet;
            if (std::abs(bet) < 1e-10 /*|| iter >= (2 * ham.model().symmetry().sector().size())*/) {
              break;
            }
          }
          hamiltonian().storage().finalize(0, false);
        }
#ifdef USE_MPI
        MPI_Barrier(hamiltonian().comm());
#endif
        return nlanc;
      }

      /**
       * Compute lanczos continues fraction
       */
      template<typename GF_TYPE>
      void compute_continues_fraction(double expectation_value, double excited_state, double groundstate, int nlanc, int isign, GF_TYPE &gf,
                                      const alps::gf::index_mesh::index_type & site, const alps::gf::index_mesh::index_type & spin) {
        double expb = 0;
        double shift;
        if (_omega.beta() * (excited_state - groundstate) > 25)
          expb = 0;
        else
          expb = exp(-_omega.beta() * (excited_state - groundstate));
        const std::vector < double > &freqs = _omega.points();
        for (int iomega = 0; iomega < _omega.extent(); ++iomega) {
          shift = 1.0;
          std::complex<double> swp = 0.0;
          std::complex < double > ener = std::complex < double >(0.0, freqs[iomega]) + (excited_state) * isign;
          swp = get_frac_point(expectation_value, nlanc, isign, expb, shift, ener);

          gf(alps::gf::matsubara_positive_mesh::index_type(iomega), site, spin) += swp;
        }
      }

      /**
       * Compute symmetrized lanczos continues fraction
       */
      template<typename GF_TYPE>
      void compute_sym_continues_fraction(double expectation_value, double excited_state, double groundstate, int nlanc, int isign, GF_TYPE &gf,
                                          const alps::gf::index_mesh::index_type & site) {
        double expb = 0;
        double shift;
        if (_omega.beta() * (excited_state - groundstate) > 25)
          expb = 0;
        else
          expb = exp(-_omega.beta() * (excited_state - groundstate));
        const std::vector < double > &freqs = _omega.points();
        gf(alps::gf::matsubara_positive_mesh::index_type(0), site) -= expectation_value*_omega.beta()*expb;
        for (int iomega = 1; iomega < _omega.extent(); ++iomega) {
          shift = 1.0;
          std::complex < double > ener = std::complex < double >(0.0,   freqs[iomega]) + (excited_state) * isign;
          std::complex < double > ener2 = std::complex < double >(0.0, -freqs[iomega]) + (excited_state) * isign;
          std::complex<double> swp = get_frac_point(expectation_value, nlanc, isign, expb, shift, ener);
          swp += get_frac_point(expectation_value, nlanc, isign, expb, shift, ener2);
          gf(alps::gf::matsubara_positive_mesh::index_type(iomega), site) += swp;
        }
      }

      const Hamiltonian &hamiltonian() const {
        return ham;
      };

      Hamiltonian &hamiltonian() {
        return ham;
      };

    private:
      alps::gf::matsubara_positive_mesh _omega;

      int _Nl;
      Hamiltonian &ham;

      std::vector < precision > alfalanc;
      std::vector < precision > betalanc;

      std::vector < std::complex < double > > det;
      std::vector < std::complex < double > > dl;

      /**
       *
       * @param expectation_value
       * @param nlanc - number of Lanczos operation have been performed
       * @param isign - sign
       * @param expb - Boltzman exponent value
       * @param shift - overflow avoiding shift
       * @param ener - first denominator in fraction (w - E)
       * @return GF value for specific energy point.
       */
      std::complex < double > get_frac_point(double expectation_value, int nlanc, int isign, double expb, double shift, const std::complex < double > &ener) {
        std::complex < double > swp = 0.0;
        det.assign(nlanc, 0.0);
        for (int i = 0; i < nlanc; ++i) {
          dl[i] = ener - ((double) (alfalanc[i]) * isign);
        }
        if (nlanc == 1) {
          det[0] = dl[0];
          swp += expectation_value * expb / det[0];
        } else {
          det[nlanc - 1] = dl[nlanc - 1];
          det[nlanc - 2] = dl[nlanc - 2] * dl[nlanc - 1] - std::pow(betalanc[nlanc - 1], 2);
          for (int i = nlanc - 3; i >= 0; --i) {
            det[i] = (dl[i] * det[i + 1] - std::pow(betalanc[i + 1], 2) * det[i + 2]);
            if(abs(det[i]) > (std::numeric_limits<float>::max() / 2.0) && i != 0) {
              shift = 1.0/(std::numeric_limits<float>::max() / 1000.0);
              det[i] *=shift;
              det[i + 1] *=shift;
            }
          }
          swp += expectation_value * expb * det[1] / det[0];
        }
        return swp;
      }
    };

  }
}
#endif //HUBBARD_LANCZOS_H
