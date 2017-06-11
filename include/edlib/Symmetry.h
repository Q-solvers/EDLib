//
// Created by iskakoff on 19/07/16.
//

#ifndef HUBBARD_SYMMETRY_H
#define HUBBARD_SYMMETRY_H

#include <alps/params.hpp>

namespace EDLib {
  namespace Symmetry {
/**
 * Base class for symmetries
 */
    class Symmetry {
    public:
      Symmetry() : _state(0) {};

      virtual ~Symmetry() {};

      /**
       * Check and, if possible, change basis state to the next state
       * @return true if there was next basis state
       */
      virtual bool next_state() = 0;

      /**
       * @param spin
       * @return true if we can create/destroy particle for specific spin
       */
      virtual bool can_create_particle(int spin) = 0;
      virtual bool can_destroy_particle(int spin) = 0;

      /**
       * @return current basis state
       */
      long long state() const {
        return _state;
      };
      long long &state() {
        return _state;
      };

      /**
       * @param state - basis state
       * @return index in the ordered basis for specific basis state
       */
      virtual int index(long long state) = 0;

      /**
       * Reset stateful object ot the initial state
       */
      virtual void reset() = 0;

      /**
       * init symmetry stateful object
       */
      virtual void init() = 0;

      /**
       * check that there is the next simmetry sector. if exist set current sector to the next available
       * @return true if there is next state
       */
      virtual bool next_sector() = 0;

    private:
      /// basis state
      long long _state;
    };
  }
}
#endif //HUBBARD_SYMMETRY_H
