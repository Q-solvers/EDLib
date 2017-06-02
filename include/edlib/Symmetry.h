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

      virtual bool next_state() = 0;

      virtual bool can_create_particle(int spin) = 0;
      virtual bool can_destroy_particle(int spin) = 0;

      long long state() const {
        return _state;
      };

      long long &state() {
        return _state;
      };

      virtual int index(long long combination) = 0;

      virtual void reset() = 0;

      virtual void init() = 0;

      // check that there is the next simmetry sector. if exist set current sector to the next available
      virtual bool next_sector() = 0;

    private:
      long long _state;
    };
  }
}
#endif //HUBBARD_SYMMETRY_H
