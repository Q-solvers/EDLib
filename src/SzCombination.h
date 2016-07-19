//
// Created by iskakoff on 19/07/16.
//

#ifndef HUBBARD_SZCOMBINATION_H
#define HUBBARD_SZCOMBINATION_H

#include "Combination.h"
#include <alps/params.hpp>

class SzCombination: public Combination {
public:

  SzCombination(alps::params& p): Combination(p), _state(0) {};
  virtual ~SzCombination() {};

  virtual bool next_combination() override {
    return false;
  }

  virtual long long combination() override {
    return _state;
  }

  virtual int index(long long combination) override {
    return 0;
  }

  virtual long long combination(int index) override {
    return 0;
  }

  virtual void reset() override {_state = 0ll;}

private:



  long long _state;
};


#endif //HUBBARD_SZCOMBINATION_H
