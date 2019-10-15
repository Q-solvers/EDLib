//
// Created by iskakoff on 01/06/17.
//

#ifndef EDLIB_HOLSTEINANDERSONPARAMETER_H
#define EDLIB_HOLSTEINANDERSONPARAMETER_H

namespace EDLib {
  namespace Ext {
    /**
     * Defines additional parameters for Holstein Anderson impurity model
     * @param params - ALPSCore parameters container
     */
    void define_parameters(alps::params &params) {
      params.define < int >("NBBITS", 3, "Total maximum number of bits per bosonic level in the system.");
      params.define < int >("NBLEVEL", 1, "Total maximum number of bosonic levels in the system.");
    }
  }
}

#endif //EDLIB_HOLSTEINANDERSONPARAMETER_H
