//
// Created by iskakoff on 23/08/16.
//

#ifndef HUBBARD_SINGLEIMPURITYANDERSONMODEL_H
#define HUBBARD_SINGLEIMPURITYANDERSONMODEL_H

#include "FermionicModel.h"

namespace EDLib {
  namespace Model {
    namespace SingleImpurityAnderson {
      class InnerState {
      public:
        virtual int valid(long long, int) const {return 0;};
        virtual void set(long long,long long&, int&, int) const {};
        int inline checkState(long long nst, int im, int Ns) const {
          return (int)((1ll << (2*Ns - 1 - im)) >> (2*Ns - 1 - im));
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
      };
      template<typename prec>
      class InnerHybridizationState : public InnerState {
      public:
        InnerHybridizationState(int ii, int jj, int spin, prec val) : _indicies(ii, jj), _spin(spin), _value(val) {};

        const inline std::pair < int, int > &indicies() const { return _indicies; }

        const inline prec &value() const { return _value; }

        const inline int spin() const { return _spin; }
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
      class InnerInteractionState : public InnerState {
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
          a(_l + _sigmaprime * Ns, k3, k4, isign1, 2*Ns);
          adag(_j + _sigmaprime * Ns, k4, k2, isign2, 2*Ns);
          adag(_i + _sigma * Ns, k2, k1, isign2, 2*Ns);
          k = k1;
          sign = isign1 * isign2;
        }

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
    template<typename precision>
    class SingleImpurityAndersonModel : public FermionicModel {
    public:
      typedef typename Symmetry::ImpuritySzSymmetry SYMMETRY;
      typedef typename SingleImpurityAnderson::InnerState St;
      typedef typename SingleImpurityAnderson::InnerHybridizationState<precision> HSt;
      typedef typename SingleImpurityAnderson::InnerInteractionState<precision> USt;
      typedef typename Symmetry::ImpuritySzSymmetry::Sector Sector;

      SingleImpurityAndersonModel(EDParams &p): _symmetry(p, p["NSITES"], p["siam.NORBITALS"]), _Ns(p["NSITES"]), _ml(p["siam.NORBITALS"]) {
        std::string input = p["INPUT_FILE"];
        alps::hdf5::archive input_data(input.c_str(), "r");
        input_data >> alps::make_pvp("NSPINS", _ms);
        _symmetry = SYMMETRY(p, _Ns, _ml);
        if ((_Ns - _ml) % _ml != 0 || _ml > _Ns) {
          throw std::logic_error("Incorrect values for the total number of sites and the number of orbitals. Please check input file.");
        }
        if(_ms != 2) {
          throw std::logic_error("Incorrect values for the number of spins. Please check input file.");
        }
        _Ip = _ms * _Ns;
        _Nk = (_Ns - _ml) / _ml;
        _Eps.assign(_ml, std::vector < precision >(_ms, precision(0.0)));
        _Epsk.assign(_ml, std::vector < std::vector < precision > >(_Nk, std::vector < precision >(_ms, precision(0.0)))),
          _Vk.assign(_ml, std::vector < std::vector < precision > >(_Nk, std::vector < precision >(_ms, precision(0.0)))),
          _U.assign(_ml, std::vector < std::vector < std::vector < precision > > >
            (_ml, std::vector < std::vector < precision > >
              (_ml, std::vector < precision >(_ml, precision(0.0)))));
        input_data >> alps::make_pvp("hopping/values", _Vk);
        input_data >> alps::make_pvp("hopping/values", _Epsk);
        input_data >> alps::make_pvp("hopping/values", _Eps);
        input_data >> alps::make_pvp("XMU", _xmu);
        input_data >> alps::make_pvp("interaction/values", _U);
        input_data.close();
        for (int im = 0; im < _ml; ++im) {
          for (int ik = 0; ik < _Nk; ++ik) {
            for (int is = 0; is < _ms; ++is) {
              if (std::abs(_Vk[im][ik][is]) > 1e-10) {
                _T_states.push_back(HSt(im, ik, is, _Vk[im][ik][is]));
              }
            }
          }
        }
        for (int is1 = 0; is1 < _ms; ++is1) {
          for (int is2 = 0; is2 < _ms; ++is2) {
            for (int i = 0; i < _ml; ++i) {
              for (int j = 0; j < _ml; ++j) {
                for (int k = 0; k < _ml; ++k) {
                  for (int l = 0; l < _ml; ++l) {
                    if (i == j && i == k && i == l) {
                      // skip diagonal contribution
                      continue;
                    }
                    _V_states.push_back(USt(i, j, k, l, is1, is2, _U[i][j][k][l]));
                  }
                }
              }
            }
          }
        }
      }

      inline const precision diagonal(long long state) const {
        precision xtemp = 0.0;
        for (int im = 0; im < _ml; ++im) {
          for (int is = 0; is < _ms; ++is) {
            for (int ik = 0; ik < _Nk; ++ik) {
              int ikm = ik + (im) * _Nk + _ml;
              xtemp+=(_Epsk[im][ik][is] * checkState(state, ikm + is * _Ns, _Ip));
            }
            xtemp += (_Eps[im][is] - _xmu) * checkState(state, im + is * _Ns, _Ip);
          }
          xtemp += _U[im][im][im][im] * checkState(state, im, _Ip) * checkState(state, im + _Ns, _Ip);
        }
        return xtemp;
      }

      inline long long interacting_states(long long nst) {
        long long up = 0;
        for (int is = 0; is < _ms; ++is) {
          up = nst >> (_Ip - _ml);
        }
        long long down = (nst & ((1ll<<_Ns) - 1))>>(_Ns-_ml);
        return (up<<_ml) + down;
      }

      inline int valid(const St &state, long long nst) {
        return state.valid(nst, _Ns);
      }

      inline void set(const St &state, long long nst, long long &k, int &sign) {
        state.set(nst, k, sign, _Ns);
      }

      SYMMETRY &symmetry() {
        return _symmetry;
      }

      inline const std::vector<St>& T_states() const {
        return _T_states;
      }

      inline const std::vector < St > &V_states() const {
        return _V_states;
      }

      const int interacting_orbitals() const {
        return _ml;
      }

    private:
      SYMMETRY _symmetry;
      int _Ns;
      int _ms;
      int _ml;
      int _Nk;
      int _Ip;
      precision _xmu;
      std::vector < std::vector < precision > > _Eps;
      std::vector < std::vector < std::vector < precision > > > _Vk;
      std::vector < std::vector < std::vector < precision > > > _Epsk;
      std::vector < std::vector < std::vector < std::vector < precision > > > > _U;

      // Kinetic part of the off-diagonal Hamiltonian elements
      std::vector < St > _T_states;
      // Interaction part of the off-diagonal Hamiltonian elements
      std::vector < St > _V_states;
    };

  }
}
#endif //HUBBARD_SINGLEIMPURITYANDERSONMODEL_H
