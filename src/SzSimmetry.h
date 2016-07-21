//
// Created by iskakoff on 19/07/16.
//

#ifndef HUBBARD_SZCOMBINATION_H
#define HUBBARD_SZCOMBINATION_H

#include <queue>

#include <alps/params.hpp>

#include "Simmetry.h"

/**
 * Sz simmetry class
 */
class SzSimmetry: public Simmetry {
public:

  SzSimmetry(alps::params& p): Simmetry(p), _state(0), _current_sector(-1,-1,0), _Ns(p["NSITES"]) {
    //TODO: read sectors from parameter file
  };
  virtual ~SzSimmetry() {};

  virtual bool next_state() override {
    return false;
  }

  virtual long long state() override {
    return _state;
  }

  virtual int index(long long state) override {
    return 0;
  }

  virtual long long state(int index) override {
    return 0;
  }

  virtual void reset() override {_state = 0ll;}

  virtual void init() override {
    // TODO: Decide what we should have to init
    reset();
  };

  virtual bool next_sector() override {
    if(_sectors.empty())
      return false;
    _current_sector = _sectors.front();
    _sectors.pop();
    return true;
  }

private:
  class Sector {
    friend class SzSimmetry;
  protected:
    Sector(int up, int down, size_t size) : _nup(up), _ndown(down), _size(size)  {};
    bool operator<(Sector s) {
      return _size< s._size;
    }
  private:
    int _nup;
    int _ndown;
    size_t _size;
  };

  SzSimmetry::Sector _current_sector;
  std::queue<SzSimmetry::Sector> _sectors;
  long long _state;
  int _Ns;
};


#endif //HUBBARD_SZCOMBINATION_H
