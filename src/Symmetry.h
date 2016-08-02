//
// Created by iskakoff on 19/07/16.
//

#ifndef HUBBARD_COMBINATION_H
#define HUBBARD_COMBINATION_H

#include <alps/params.hpp>

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
  Symmetry(alps::params& p){};
  virtual ~Symmetry(){};
  virtual bool next_state() = 0;
  virtual long long state() = 0;
//  virtual int index(long long combination) = 0;
  virtual long long state(int index) = 0;
  virtual void reset() = 0;
  virtual void init() = 0;
  // check that there is the next simmetry sector. if exist set current sector to the next available
  virtual bool next_sector() = 0;
};


#endif //HUBBARD_COMBINATION_H
