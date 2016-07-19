//
// Created by iskakoff on 19/07/16.
//

#ifndef HUBBARD_COMBINATION_H
#define HUBBARD_COMBINATION_H

#include <alps/params.hpp>

/**
 * Base class for combinatorics
 */
class Combination {
public:
  Combination(alps::params& p){};
  virtual ~Combination(){};
  virtual bool next_combination() = 0;
  virtual long long combination() = 0;
  virtual int index(long long combination) = 0;
  virtual long long combination(int index) = 0;
  virtual void reset() = 0;
};


#endif //HUBBARD_COMBINATION_H
