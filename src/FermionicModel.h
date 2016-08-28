//
// Created by iskakoff on 28/08/16.
//

#ifndef HUBBARD_FERMIONICMODEL_H
#define HUBBARD_FERMIONICMODEL_H

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
    };
  }
}

#endif //HUBBARD_FERMIONICMODEL_H
