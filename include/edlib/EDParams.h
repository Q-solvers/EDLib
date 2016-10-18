//
// Created by iskakoff on 18/08/16.
//

#ifndef HUBBARD_EDPARAMS_H
#define HUBBARD_EDPARAMS_H


#include <alps/params.hpp>
namespace EDLib {
  class EDParams : public alps::params {
  public:
    EDParams() : params() {
      define_parameters();
    }

    EDParams(int argc, const char **argv) : alps::params(argc, argv) {
      define_parameters();
    }

  private:
    /**
   * Define parameters used in programm
   */
    void define_parameters() {

      define < int >("NSITES_BOSE", 2, "Number of bosonic sites. Should be 2^K - 1");
    }

  };

}
#endif //HUBBARD_EDPARAMS_H
