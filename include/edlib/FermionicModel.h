#ifndef EDLIB_FERMIONICMODEL_H
#define EDLIB_FERMIONICMODEL_H

#include "edlib/Parameters.h"

namespace edlib {

  /**
   * Base class for fermionic models with binary-state representation.
   *
   * solve_dyson and bare_greens_function (present in the legacy
   * EDLib::Model::FermionicModel) are deliberately not here — they belong to
   * derived models that know about Green's-function machinery (HubbardModel
   * etc.) and are ported in Phase 3 alongside Eigen-based linear algebra.
   */
  class FermionicModel {
  public:
    explicit FermionicModel(const Parameters& p)
        : _Ns(p.nsites), _ms(p.nspins), _Ip(p.nspins * p.nsites) {}

    int orbitals()             const { return _Ns; }
    int spins()                const { return _ms; }
    int max_total_electrons()  const { return _Ip; }

    inline int checkState(long long nst, int im, int Ip) const {
      return static_cast<int>((nst & (1ll << (Ip - 1 - im))) >> (Ip - 1 - im));
    }

    inline void a(int i, long long jold, long long& k, int& isign) const {
      long long sign = 0;
      for (int ll = 0; ll < i; ++ll) {
        sign += ((jold & (1ll << (_Ip - ll - 1))) != 0) ? 1 : 0;
      }
      isign = (sign % 2) == 0 ? 1 : -1;
      k = jold - (1ll << (_Ip - i - 1));
    }

    inline void adag(int i, long long jold, long long& k, int& isign) const {
      long long sign = 0;
      for (int ll = 0; ll < i; ++ll) {
        sign += ((jold & (1ll << (_Ip - ll - 1))) != 0) ? 1 : 0;
      }
      isign = (sign % 2) == 0 ? 1 : -1;
      k = jold + (1ll << (_Ip - i - 1));
    }

  protected:
    int _Ns;   ///< number of lattice sites
    int _ms;   ///< number of electron spins
    int _Ip;   ///< maximum number of electrons (= _Ns * _ms)
  };

}

#endif
