//
// Created by iskakoff on 19/07/16.
//

#ifndef HUBBARD_SYMMETRY_H
#define HUBBARD_SYMMETRY_H

#include "EDParams.h"

/**
 *  TODO: decide do we really need this abstract class or should we remove it and just make specification for Symmetry
 *  The profit of this class is clear specification of methods.
 *  But compiler will never see it because we specify symmetry as a template parameter
 *  Therefore it gives us a virtual methods overhead.
 */
/**
 * Base class for symmetries
 */
class Symmetry {
public:
  Symmetry(): _state(0){};
  virtual ~Symmetry(){};
  virtual bool next_state() = 0;
  const long long state() const {
    return _state;
  };
  long long & state() {
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


#endif //HUBBARD_SYMMETRY_H
