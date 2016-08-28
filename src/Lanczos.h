//
// Created by iskakoff on 01/08/16.
//

#ifndef HUBBARD_LANCZOS_H
#define HUBBARD_LANCZOS_H

#include <alps/gf/mesh.hpp>
#include <alps/gf/gf.hpp>

namespace EDLib {
  namespace gf {
    template<typename precision, class Hamiltonian>
    class Lanczos {
    public:
      Lanczos(EDParams &p, Hamiltonian &h) : _omega(p["lanc.BETA"], p["lanc.NOMEGA"]), _Nl(p["lanc.NLANC"]), ham(h), alfalanc(_Nl), betalanc(_Nl + 1), det(_Nl), dl(_Nl) {

      }

      const alps::gf::matsubara_positive_mesh &omega() const {
        return _omega;
      }

      int lanczos(std::vector < precision > &v) {
        int nlanc = 0;
        unsigned long size = v.size();
        std::vector < precision > w(size, precision(0.0));
        precision alf = 0, bet = 0;
        ham.fill();
        for (int iter = 1; iter <= _Nl; iter++) {
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
          for (int k = 0; k < v.size(); ++k) {
            alf += w[k] * v[k];
          }
          alfalanc[iter - 1] = alf;
          for (int j = 0; j < size; ++j) {
            w[j] -= alf * v[j];
          }
          for (int k = 0; k < v.size(); ++k) {
            bet += w[k] * w[k];
          }
          bet = std::sqrt(bet);

          if (iter != _Nl) betalanc[iter] = bet;
          if (std::abs(bet) < 1e-10 /*|| iter >= (2 * ham.model().symmetry().sector().size())*/) {
            break;
          }
        }
        return nlanc;
      }

      /**
       * Compute lanczos continues fraction
       */
      void computefrac(double expectation_value, double excited_state, double groundstate, int nlanc, int isign, alps::gf::omega_gf &gf) {
        double expb = 0;
        if (_omega.beta() * (excited_state - groundstate) > 25)
          expb = 0;
        else
          expb = exp(-_omega.beta() * (excited_state - groundstate));
        const std::vector < double > &freqs = _omega.points();
        for (int iomega = 0; iomega < _omega.extent(); ++iomega) {
          std::complex < double > ener = std::complex < double >(0.0, freqs[iomega]) + (excited_state) * isign;
          det.assign(nlanc, 0.0);
          for (int i = 0; i < nlanc; ++i) {
            dl[i] = ener - ((double) (alfalanc[i]) * isign);
          }
          if (nlanc == 1) {
            det[0] = dl[0];
            gf(alps::gf::matsubara_positive_mesh::index_type(iomega)) += expectation_value * expb / det[0];
          } else {
            det[nlanc - 1] = dl[nlanc - 1];
            det[nlanc - 2] = dl[nlanc - 2] * dl[nlanc - 1] - std::pow(betalanc[nlanc - 1], 2);
            for (int i = nlanc - 3; i >= 0; --i) {
              det[i] = dl[i] * det[i + 1] - std::pow(betalanc[i + 1], 2) * det[i + 2];
            }
            gf(alps::gf::matsubara_positive_mesh::index_type(iomega)) += expectation_value * expb * det[1] / det[0];
          }
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
    };

  }
}
#endif //HUBBARD_LANCZOS_H
