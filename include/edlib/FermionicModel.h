//
// Created by iskakoff on 28/08/16.
//

#ifndef HUBBARD_FERMIONICMODEL_H
#define HUBBARD_FERMIONICMODEL_H

#include <alps/params.hpp>

namespace EDLib {
  namespace Model {
/**
 * @brief FermionicModel base class
 *
 * Define common fermionic routines for binary represented state
 *
 * @author iskakoff
 */
    class FermionicModel {
    public:
      FermionicModel(alps::params &p) : _Ns(p["NSITES"]), _ms(p["NSPINS"]), _Ip(int(p["NSPINS"]) * int(p["NSITES"])) {
      }

      /**
       * @brief Check that im state is occupated
       *
       * @param nst - current state
       * @param im - state to check
       * @param Ip - total number of fermionic spins for all sites
       *
       * @return 0 if state is empty, 1 - otherwise
       */
      int inline checkState(long long nst, const int &im, int Ip) const {
        return (int) ((nst & (1ll << (Ip - 1 - im))) >> (Ip - 1 - im));
      }
      /**
       * @brief Anihilate particle
       * @param i [in] - site to anihilate particle
       * @param jold [in] - current state
       * @param k [out] - resulting state
       * @param isign [out] - fermionic sign
       */
      void inline a(int i, long long jold, long long &k, int &isign) {
        long long sign = 0;
        for (int ll = 0; ll < i; ++ll) {
          sign += ((jold & (1ll << (_Ip - ll - 1))) != 0) ? 1 : 0;
        }
        isign = (sign % 2) == 0 ? 1 : -1;
        k = jold - (1ll << (_Ip - i - 1));
      }

      /**
       * @brief Create particle
       * \param i [in] - site to create particle
       * \param jold [in] - current state
       * \param k [out] - resulting state
       * \param isign [out] - fermionic sign
       */
      void inline adag(int i, long long jold, long long &k, int &isign) {
        long long sign = 0;
        for (int ll = 0; ll < i; ++ll) {
          sign += ((jold & (1ll << (_Ip - ll - 1))) != 0) ? 1 : 0;
        }
        isign = (sign % 2) == 0 ? 1 : -1;
        k = jold + (1ll << (_Ip - i - 1));
      }


      const int orbitals() const {
        return _Ns;
      }

      const int max_total_electrons() const {
        return _Ip;
      }

      const int spins() const {
        return _ms;
      }

    protected:
      /**
       * _Ns - number of lattice sites
       * _ms - number of electron spins
       * _Ip - maximum number of electrons
       */
      int _Ns;
      int _ms;
      int _Ip;
    };
  }
}

#endif //HUBBARD_FERMIONICMODEL_H
