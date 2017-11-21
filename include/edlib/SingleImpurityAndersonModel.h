//
// Created by iskakoff on 23/08/16.
//

#ifndef HUBBARD_SINGLEIMPURITYANDERSONMODEL_H
#define HUBBARD_SINGLEIMPURITYANDERSONMODEL_H

#include <alps/gf/mesh.hpp>
#include <alps/gf/gf.hpp>
#include "FermionicModel.h"
#include "CommonUtils.h"

namespace EDLib {
  namespace Model {
    namespace SingleImpurityAnderson {
      template<typename prec>
      class InnerState {
      public:
        virtual int valid(long long, int) const {return 0;};
        virtual void set(long long,long long&, int&, int) const {};
        int inline checkState(long long nst, int im, int Ns) const {
          return (int)((nst & (1ll << (2*Ns - 1 - im))) >> (2*Ns - 1 - im));
        }

        /**
         * Anihilate particle
         * \param i [in] - site to anihilate particle
         * \param jold [in] - current state
         * \param k [out] - resulting state
         * \param isign [out] - fermionic sign
         * \param Ip [in] - number of fermionic sites
         */
        void inline a(int i, long long jold, long long &k, int &isign, int Ip) const {
          long long sign = 0;
          for (int ll = 0; ll < i; ++ll) {
            sign += ((jold & (1ll << (Ip - ll - 1))) != 0) ? 1 : 0;
          }
          isign = (sign % 2) == 0 ? 1 : -1;
          k = jold - (1ll << (Ip - i - 1));
        }

        /**
         * Create particle
         * \param i [in] - site to create particle
         * \param jold [in] - current state
         * \param k [out] - resulting state
         * \param isign [out] - fermionic sign
         * \param Ip [in] - number of fermionic sites
         */
        void inline adag(int i, long long jold, long long &k, int &isign, int Ip) const {
          long long sign = 0;
          for (int ll = 0; ll < i; ++ll) {
            sign += ((jold & (1ll << (Ip - ll - 1))) != 0) ? 1 : 0;
          }
          isign = (sign % 2) == 0 ? 1 : -1;
          k = jold + (1ll << (Ip - i - 1));
        }

        virtual inline prec value() const { return 0.0; }
      };
      template<typename prec>
      class InnerHybridizationState : public InnerState<prec> {
        using InnerState<prec>::checkState;
        using InnerState<prec>::a;
        using InnerState<prec>::adag;
      public:
        InnerHybridizationState(int ii, int jj, int spin, prec val) : _indicies(ii, jj), _spin(spin), _value(val) {};

        const inline std::pair < int, int > &indicies() const { return _indicies; }

        virtual inline prec value() const { return _value; }

        inline int spin() const { return _spin; }
        virtual int valid(long long nst, int Ns) const {
          return (checkState(nst, _indicies.first + _spin * Ns, Ns) * (1 - checkState(nst, _indicies.second + _spin * Ns, Ns)));
        }
        virtual void set(long long nst,long long&k, int&sign, int Ns) const {
          long long k1, k2;
          int isign1, isign2;
          a(_indicies.first + _spin * Ns, nst, k1, isign1, 2*Ns);
          adag(_indicies.second + _spin * Ns, k1, k2, isign2, 2*Ns);
          k = k2;
          sign = isign1 * isign2;
        }

      private:
        std::pair < int, int > _indicies;
        int _spin;
        prec _value;
      };
      template<typename prec>
      class InnerInteractionState : public InnerState<prec> {
        using InnerState<prec>::checkState;
        using InnerState<prec>::a;
        using InnerState<prec>::adag;
      public:
        InnerInteractionState(int i, int j, int k, int l, int sigma, int sigmaprime, prec U) :
          _i(i), _j(j), _k(k), _l(l), _sigma(sigma), _sigmaprime(sigmaprime), _U(U) {}

        int i() const {
          return _i;
        }

        int j() const {
          return _j;
        }

        int k() const {
          return _k;
        }

        int l() const {
          return _l;
        }

        prec U() const {
          return _U;
        }
        /**
         * @brief Check the possible transition
         *
         * Evaluate the following four operators product:
         * a^*_i a^*_j a_l a_k | nst>
         *
         * @param nst - current state
         * @param Ns - number of fermionic sites
         * @return One if transition is possible
         */
        virtual int valid(long long nst, int Ns) const {
          int Ip = 2*Ns;
          if(checkState(nst, _k + _sigma * Ns, Ns) != 0) {
            long long k3 = nst - (1ll << (Ip - 1 - _k - _sigma * Ns));
            if (checkState(k3, _l + _sigmaprime * Ns, Ns) != 0) {
              long long k4 = k3 - (1ll << (Ip - 1 - _l - _sigmaprime * Ns));
              if (checkState(k4, _j + _sigmaprime * Ns, Ns) == 0) {
                long long k2 = k4 | (1ll << (Ip - 1 - _j - _sigmaprime * Ns));
                return (1-checkState(k2, _i + _sigma * Ns, Ns));
              }
            }
          }
          return 0;
        }
        /**
         * @brief Computes the new state for inter-orbital Coulomb transition
         *
         * |k> = a^*_i a^*_j a_l a_k | nst>
         *
         * @param nst - current state
         * @param k - next state
         * @param sign - sign of transition
         * @param Ns - number of fermionic sites
         */
        virtual void set(long long nst,long long&k, int&sign, int Ns) const {
          long long k1, k2, k3, k4;
          int isign1, isign2, isign3, isign4;
          a(_k + _sigma * Ns, nst, k3, isign1, 2*Ns);
          a(_l + _sigmaprime * Ns, k3, k4, isign2, 2*Ns);
          adag(_j + _sigmaprime * Ns, k4, k2, isign3, 2*Ns);
          adag(_i + _sigma * Ns, k2, k1, isign4, 2*Ns);
          k = k1;
          sign = isign1 * isign2*isign3*isign4;
        }

        /**
         * @brief Return the interaction strength for current spin-orbital combination
         * @return U_{ijkl}
         */
        virtual inline prec value() const { return 0.5*_U; }

      private:
        int _i;
        int _j;
        int _k;
        int _l;
        int _sigma;
        int _sigmaprime;
        prec _U;
      };
    }
    /**
     * Single multi-orbital Impurity Anderson Model class
     *
     * @tparam prec - floating point precision
     */
    template<typename prec>
    class SingleImpurityAndersonModel : public FermionicModel {
    public:
      typedef prec precision;
      typedef typename Symmetry::SzSymmetry SYMMETRY;
      typedef typename SingleImpurityAnderson::InnerState<precision> St;
      typedef typename SingleImpurityAnderson::InnerHybridizationState<precision> HSt;
      typedef typename SingleImpurityAnderson::InnerInteractionState<precision> USt;
      typedef typename Symmetry::SzSymmetry::Sector Sector;

      SingleImpurityAndersonModel(alps::params &p): FermionicModel(p), _symmetry(p), _ml(p["siam.NORBITALS"]),
                                                    _Epsk(p["siam.NORBITALS"], std::vector<std::vector<double> >()),
                                                    _Vk(p["siam.NORBITALS"], std::vector<std::vector<double> >()),
                                                    _t0(p["siam.NORBITALS"], std::vector<std::vector<double> >(p["siam.NORBITALS"], std::vector<double>(_ms, 0.0))),
                                                    _bath_ind(p["siam.NORBITALS"], 0) {
        std::string input = p["INPUT_FILE"];
        alps::hdf5::archive input_data(input.c_str(), "r");
        if (_ml > _Ns) {
          throw std::invalid_argument("Incorrect values for the total number of sites and the number of orbitals. Please check input file.");
        }
        if(_ms != 2) {
          throw std::invalid_argument("Incorrect values for the number of spins. Please check input file.");
        }
        _Ip = _ms * _Ns;
        for (int im = 0; im < _ml; ++im) {
          std::stringstream s;
          s<<"Bath/Vk_"<<im<<"/values";
          input_data >> alps::make_pvp(s.str().c_str(), _Vk[im]);
          s.str("");
          s<<"Bath/Epsk_"<<im<<"/values";
          input_data >> alps::make_pvp(s.str().c_str(), _Epsk[im]);
          s.str("");
          s<<"t0_"<<im<<"/values";
          input_data >> alps::make_pvp(s.str().c_str(), _t0[im]);
        }
        input_data >> alps::make_pvp("Eps0/values", _Eps0);
        input_data >> alps::make_pvp("mu", _xmu);
        input_data >> alps::make_pvp("interaction/values", _U);
        input_data.close();
        if(_U.size() != _ml) {
          throw std::invalid_argument("Incorrect number of orbitals. Please check input file.");
        }
        int b_ind = 0;
        for(int im = 0; im< _ml; ++im ){
          _bath_ind[im] = b_ind;
          b_ind += _Vk[im].size();
        }
        for(int im = 0; im< _ml; ++im ){
          if(_t0[im].size()>_ml) {
            throw std::invalid_argument("Inter orbital hoppings array dimension are bigger than number of impurity orbitals");
          }
        }
        if(b_ind != _Ns - _ml) {
          throw std::invalid_argument("Total number of state does not equal to sum of the total number of bath levels and the number of impurity orbitals");
        }
        // interorbital hoppings
        for (int im = 0; im < _ml; ++im) {
          for (int jm = 0; jm < im; ++jm) {
            for (int is = 0; is < _ms; ++is) {
              if (std::abs(_t0[im][jm][is]) > 1e-10) {
                _T_states.push_back(HSt(im, jm, is, _t0[im][jm][is]));
                _T_states.push_back(HSt(jm, im, is, _t0[im][jm][is]));
              }
            }
          }
        }
        // fill hybridization part
        for (int im = 0; im < _ml; ++im) {
          for (int ik = 0; ik < _Vk[im].size(); ++ik) {
            for (int is = 0; is < _ms; ++is) {
              if (std::abs(_Vk[im][ik][is]) > 1e-10) {
                int imk = ik + _bath_ind[im] + _ml;
                _T_states.push_back(HSt(im, imk, is, _Vk[im][ik][is]));
                _T_states.push_back(HSt(imk, im, is, _Vk[im][ik][is]));
              }
            }
          }
        }
        // fill off-diagonal interaction term
        for (int is1 = 0; is1 < _ms; ++is1) {
          for (int is2 = 0; is2 < _ms; ++is2) {
            for (int i = 0; i < _ml; ++i) {
              for (int j = 0; j < _ml; ++j) {
                for (int k = 0; k < _ml; ++k) {
                  for (int l = 0; l < _ml; ++l) {
                    // skip density-density contribution
                    if ( ( (i == l) && (j == k) && (is1==is2) ) || ( (i == k) && (j == l) ) ) {
                      continue;
                    }
                    if(std::abs(_U[i][j][k][l]) != 0.0) {
                      _V_states.push_back(USt(i, j, k, l, is1, is2, _U[i][j][k][l]));
                    }
                  }
                }
              }
            }
          }
        }
      }

      /**
       * computes diagonal contribution for the specific occupation basis state
       * @param state - occupation basis state
       * @return <state | H | state>
       */
      inline const precision diagonal(long long state) const {
        precision xtemp = 0.0;
        for (int im = 0; im < _ml; ++im) {
          for (int is = 0; is < _ms; ++is) {
            for (int ik = 0; ik < _Epsk[im].size(); ++ik) {
              int ikm = ik + _bath_ind[im] + _ml;
              xtemp+=(_Epsk[im][ik][is] * checkState(state, ikm + is * _Ns, _Ip));
            }
            xtemp += (_Eps0[im][is] - _xmu) * checkState(state, im + is * _Ns, _Ip);
          }
          xtemp += _U[im][im][im][im] * checkState(state, im, _Ip) * checkState(state, im + _Ns, _Ip);
          for(int jm = 0; jm < _ml; ++jm) {
            for(int is = 0; is< _ms; ++is)
            if(im!=jm) {
              xtemp += 0.5 * (_U[im][jm][im][jm] - _U[im][jm][jm][im]) * checkState(state, im + is*_Ns, _Ip) * checkState(state, jm + is*_Ns, _Ip);
              xtemp += 0.5 * (_U[im][jm][im][jm]) * checkState(state, im + is*_Ns, _Ip) * checkState(state, jm + (1-is)*_Ns, _Ip);
            }
          }
        }
        return xtemp;
      }

      /**
       * @deprecated
       */
      inline long long interacting_states(long long nst) {
        long long up = 0;
        for (int is = 0; is < _ms; ++is) {
          up = nst >> (_Ip - _ml);
        }
        long long down = (nst & ((1ll<<_Ns) - 1))>>(_Ns-_ml);
        return (up<<_ml) + down;
      }

      /**
       * Check that state describes valid transition for basis vector |nst>
       * @param state - transition state
       * @param nst - occupation basis vector
       * @return 1 or 0 wheater the transition is possible or not respectively
       */
      inline int valid(const St &state, long long nst) {
        return state.valid(nst, _Ns);
      }

      /**
       * Perform transition "state" from state |nst> to |k>
       * @param state
       * @param nst
       * @param k
       * @param sign
       * @return contribution to off-diagonal Hamilonian element
       */
      inline precision set(const St &state, long long nst, long long &k, int &sign) {
        state.set(nst, k, sign, _Ns);
        return state.value();
      }

      /**
       * Model symmetry type
       */
      SYMMETRY &symmetry() {
        return _symmetry;
      }

      /**
       * @return Hybridization transitions
       */
      inline const std::vector<HSt>& T_states() const {
        return _T_states;
      }
      /**
       * @return Off-diagonal Coulomb transitions
       */
      inline const std::vector<USt>& V_states() const {
        return _V_states;
      }

      /**
       * Only impurity orbitals have Coulomb interaction. Bath is non-interacting.
       * @return number of impurity orbitals
       */
      int interacting_orbitals() const {
        return _ml;
      }

      /**
       * Computes bare Green's function
       * @tparam Mesh - Green's function frequency mesh
       * @param bare_gf - Bare Green's function container
       * @param beta - inverse temperature
       */
      template<typename Mesh>
      void bare_greens_function(alps::gf::three_index_gf<std::complex<double>, Mesh, alps::gf::index_mesh, alps::gf::index_mesh >& bare_gf, double beta) {
        for(int iw = 0; iw< bare_gf.mesh1().points().size(); ++iw) {
          typename Mesh::index_type w(iw);
          for (int im: bare_gf.mesh2().points()) {
            for (int is : bare_gf.mesh3().points()) {
              std::complex<double> delta = 0;
              for(int ik = 0; ik< _Epsk[im].size(); ++ik) {
                delta += _Vk[im][ik][is]*_Vk[im][ik][is]/(common::freq_point(iw, bare_gf.mesh1(), beta) - _Epsk[im][ik][is]);
              }
              bare_gf(w, alps::gf::index_mesh::index_type(im), alps::gf::index_mesh::index_type(is)) = 1.0/(common::freq_point(iw, bare_gf.mesh1(), beta) - _Eps0[im][is] - delta);
            }
          }
        }
      }

    private:
      /// model symmetry
      SYMMETRY _symmetry;
      /// number of impurity orbitals
      int _ml;
      /// Coulomb interaction matrix
      std::vector < std::vector < std::vector < std::vector < precision > > > > _U;
      /// chemical potential
      precision _xmu;
      /// impurity orbitals energy
      std::vector < std::vector < precision > > _Eps0;
      /// number of bath states
      int _Nk;
      /// Hoppings between orbitals
      std::vector < std::vector < std::vector < precision > > > _t0;
      /// Hybridization with bath
      std::vector < std::vector < std::vector < precision > > > _Vk;
      /// Bath energy levels
      std::vector < std::vector < std::vector < precision > > > _Epsk;
      /// indices for bath
      std::vector<int> _bath_ind;

      /// Kinetic part of the off-diagonal Hamiltonian elements
      std::vector < HSt > _T_states;
      /// Interaction part of the off-diagonal Hamiltonian elements
      std::vector < USt > _V_states;
    };

  }
}
#endif //HUBBARD_SINGLEIMPURITYANDERSONMODEL_H
