//
// Created by iskakoff on 01/08/16.
//

#ifndef HUBBARD_LANCZOS_H
#define HUBBARD_LANCZOS_H

#include <alps/gf/mesh.hpp>
#include <alps/gf/gf.hpp>
#include <alps/params.hpp>

#include <cmath>

#include "MeshFactory.h"

namespace EDLib {
  namespace gf {
    template<class Hamiltonian, class Mesh=alps::gf::matsubara_positive_mesh, typename ... Args>
    class Lanczos {
    protected:
      typedef typename Hamiltonian::ModelType::precision precision;
      typedef typename Mesh::index_type mesh_index;
    public:
      Lanczos(alps::params &p, Hamiltonian &h, Args...args) :
        ham(h), _omega(MeshFactory<Mesh, Args...>::createMesh(p, args...)),_Nl(p["lanc.NLANC"]),
        alfalanc(p["lanc.NLANC"], 0.0), betalanc(int(p["lanc.NLANC"]) + 1, 0.0), det(p["lanc.NLANC"], 0), dl(p["lanc.NLANC"], 0.0), _beta(p["lanc.BETA"].as<precision>()) {}

      const Mesh &omega() const {
        return _omega;
      }
    protected:
      /**
       * Lanczos basis construction
       * @param v - initial vector
       * @return number of Lanczos iteration
       */
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
      void compute_continued_fraction(double expectation_value, double excited_state, double groundstate, int nlanc, int isign, GF_TYPE &gf,
                                      const alps::gf::index_mesh::index_type &site, const alps::gf::index_mesh::index_type &spin) {
        double expb = 0;
        double shift;
        if (_beta * (excited_state - groundstate) > 25)
          expb = 0;
        else
          expb = exp(-_beta * (excited_state - groundstate));
        const std::vector < double > &freqs = _omega.points();
        for (int iomega = 0; iomega < _omega.extent(); ++iomega) {
          shift = 1.0;
          std::complex<double> swp = 0.0;
          std::complex < double > ener = freq_point(iomega) + (excited_state) * isign;
          swp = get_frac_point(expectation_value, nlanc, isign, expb, shift, ener);

          gf(mesh_index(iomega), site, spin) += swp;
        }
      }

      /**
       * Compute symmetrized lanczos continues fraction
       */
      template<typename GF_TYPE>
      void compute_sym_continued_fraction(double expectation_value, double excited_state, double groundstate, int nlanc, int isign, GF_TYPE &gf,
                                          const alps::gf::index_mesh::index_type &site) {
        double expb = 0;
        double shift;
        if (_beta * (excited_state - groundstate) > 25)
          expb = 0;
        else
          expb = exp(-_beta * (excited_state - groundstate));
        const std::vector < double > &freqs = _omega.points();
        update_static(gf, site, expectation_value, expb);
        for (int iomega = zero_freq(); iomega < _omega.extent(); ++iomega) {
          shift = 1.0;
          std::complex < double > ener =   freq_point(iomega) + (excited_state) * isign;
          std::complex < double > ener2 = -freq_point(iomega) + (excited_state) * isign;
          std::complex<double> swp = get_frac_point(expectation_value, nlanc, isign, expb, shift, ener);
          swp += get_frac_point(expectation_value, nlanc, isign, expb, shift, ener2);
          gf(mesh_index(iomega), site) += swp;
        }
      }

      const Hamiltonian &hamiltonian() const {
        return ham;
      };

      Hamiltonian &hamiltonian() {
        return ham;
      };

      precision beta() const {
        return _beta;
      }

      /**
       * Computes complex value for frequency.
       * For Matsubara frequecy z = i*omega_n. For real frequency add small imaginary temperature dependent broadering z = omega_n + i\delta
       *
       * @tparam M - mesh type
       * @param index frequency index
       * @return proper complex representation for current frequency
       */
      template<typename M=Mesh>
      typename std::enable_if<std::is_base_of<alps::gf::matsubara_positive_mesh, M>::value, std::complex<double>>::type
      freq_point(int index) {
        return std::complex<double>(0.0, _omega.points()[index]);
      };

      template<typename M=Mesh>
      typename std::enable_if<std::is_base_of<alps::gf::real_frequency_mesh, M>::value, std::complex<double>>::type
      freq_point(int index) {
        return std::complex<double>(_omega.points()[index], M_PI/_beta);
      };
      
      template<typename M=Mesh>
      typename std::enable_if<std::is_base_of<alps::gf::matsubara_positive_mesh, M>::value, std::string>::type
      suffix() {
        return "";
      };

      template<typename M=Mesh>
      typename std::enable_if<std::is_base_of<alps::gf::real_frequency_mesh, M>::value, std::string>::type
      suffix() {
        return "_r";
      };

      /**
       * Compute smallest index for frequency. Since Lanczos continued fraction can not compute zero Matsubara frequency bosonic Green's function
       * we should compute it from first non-zero Matsubara and treat zero frequency separately.
       *
       * @tparam M - mesh type
       * @return valid smallest index for frequency
       */
      template<typename M=Mesh>
      typename std::enable_if<std::is_base_of<alps::gf::matsubara_positive_mesh, M>::value, int>::type
      zero_freq() {
        return 1;
      };

      template<typename M=Mesh>
      typename std::enable_if<std::is_base_of<alps::gf::real_frequency_mesh, M>::value, int>::type
      zero_freq() {
        return 0;
      };

      template<typename GF_TYPE, typename M=Mesh>
      typename std::enable_if<std::is_base_of<alps::gf::matsubara_positive_mesh, M>::value, void>::type
      update_static(GF_TYPE& gf, const alps::gf::index_mesh::index_type &site, double expectation_value, double expb) {
        gf(mesh_index(0), site) -= expectation_value*_beta*expb;
      };

      template<typename GF_TYPE, typename M=Mesh>
      typename std::enable_if<std::is_base_of<alps::gf::real_frequency_mesh, M>::value, void>::type
      update_static(GF_TYPE& gf, const alps::gf::index_mesh::index_type &site, double expectation_value, double expb) {

      };

    private:
      /// frequency mesh
      Mesh _omega;
      /// inverse temperature
      precision _beta;

      /// maximum number of Lanczos iteration
      int _Nl;
      /// Hamiltonain object
      Hamiltonian &ham;

      /// Lanczos tridiagonal matrix
      std::vector < precision > alfalanc;
      std::vector < precision > betalanc;

      /// continued fraction arrays
      std::vector < std::complex < double > > det;
      std::vector < std::complex < double > > dl;

      /**
       * Computes continued fraction for specific frequency
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
