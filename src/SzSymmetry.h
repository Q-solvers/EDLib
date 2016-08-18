//
// Created by iskakoff on 19/07/16.
//

#ifndef HUBBARD_SZCOMBINATION_H
#define HUBBARD_SZCOMBINATION_H

#include <queue>

#include <alps/params.hpp>

#include "Symmetry.h"

/**
 * Sz symmetry class
 */
class SzSymmetry: public Symmetry {
public:
  class Sector {
  public:
    friend class SzSymmetry;
    friend std::ostream& operator << (std::ostream& o, const SzSymmetry::Sector& c) { return o << " (nup: " <<c._nup <<" ndown: "<<c._ndown<<") size: "<<c._size;  }
    Sector(int up, int down, size_t size) : _nup(up), _ndown(down), _size(size)  {};
  protected:
    bool operator<(Sector s) {
      return _size< s._size;
    }

  public:
    int nup() const {return _nup;}

    int ndown() const {return _ndown;}

    size_t size() const {return _size;}

    void print() const {
      std::cout<<_nup<<" "<<_ndown;
    }

  private:
    int _nup;
    int _ndown;
    size_t _size;
  };

  SzSymmetry(alps::params& p): Symmetry(p), _state(0), _current_sector(-1,-1,0), _Ns(p["NSITES"]), upstate(_Ns+1), dostate(_Ns+1),
                               _c_n_k(_Ns+1, std::vector<int>(_Ns+1, 0)), basis(_Ns+1), ninv(_Ns+1, std::vector<int>(1<<_Ns, 0)),
                               _first(true){
    //TODO: read sectors from parameter file
    _Ip = 2*_Ns;
    _ind = 0;
    for(int i = 0; i<=_Ns;++i) {
      int cnk = C_n_k_i(_Ns,i);
      basis[i].reserve(cnk);
      for (int k = 0; k < cnk; ++k) {
        basis[i][k] = next_basis(_Ns, i, upstate, k==0);
        ninv[i][basis[i][k]] = k;
      }
      for(int j = 0; j<= _Ns;++j) {
        _c_n_k[i][j] = C_n_k_i(i, j);
      }
    }
    if(p.exists("arpack.SECTOR") && bool(p["arpack.SECTOR"])) {
      std::vector<std::vector<int> > sectors;
      std::string input = p["INPUT_FILE"];
      alps::hdf5::archive input_file(input, "r");
      input_file>>alps::make_pvp("sectors/values", sectors);
      input_file.close();
      for(auto& sector : sectors) {
        _sectors.push(SzSymmetry::Sector(sector[0], sector[1], (size_t) (C_n_k_i(_Ns, sector[0]) * C_n_k_i(_Ns, sector[1]))));
      }
    } else {
      for(int i = 0; i<=_Ns;++i) {
        for(int j = 0; j<= _Ns;++j) {
          _sectors.push(SzSymmetry::Sector(i, j, (size_t) (C_n_k_i(_Ns, i) * C_n_k_i(_Ns, j))));
        }
      }
    }
  };
  virtual ~SzSymmetry() {};

  virtual bool next_state() override {
    long long res = 0;
    if(_first) {
      _first = false;
    }
    if(_ind>=_current_sector.size()) {
      return false;
    }
    int u = 0, d = 0;
    u = _ind / _c_n_k[_Ns][_current_sector.ndown()];
    d = _ind % _c_n_k[_Ns][_current_sector.ndown()];
    res=basis[_current_sector.nup()][u];
    res<<=_Ns;
    res+=basis[_current_sector.ndown()][d];
    _ind++;
    _state = res;
    return true;
  }

  virtual long long state() override {
    return _state;
  }

  int index(long long state, const SzSymmetry::Sector & sector) {
    long long up = state>>_Ns;
    long long down = state & ((1ll<<_Ns) - 1);
    int cup= _c_n_k[_Ns][sector.nup()];
    int cdo= _c_n_k[_Ns][sector.ndown()];
    return ninv[sector.nup()][(int)up]*(cdo) + ninv[sector.ndown()][(int)down];
  }

  virtual int index(long long state) override {
    return index(state, _current_sector);
  }

  virtual long long state(int index) override {
    return 0;
  }

  virtual void reset() override {
    _state = 0ll;
    _first = true;
    _ind = 0;
  }

  virtual void init() override {
    // TODO: Decide what we should have to init
    reset();
    initState(_current_sector.nup(), upstate);
    initState(_current_sector.ndown(), dostate);
  };

  virtual bool next_sector() override {
    if(_sectors.empty())
      return false;
    _current_sector = _sectors.front();
    _sectors.pop();
    std::cout<<"Diagonalizating Sz-symmetry sector with nup: "<<_current_sector.nup()<<" ndown: "<<_current_sector.ndown()<<" size: "<<_current_sector.size()<<std::endl;
    return true;
  }

  void set_sector(const SzSymmetry::Sector& sector) {
    _current_sector = sector;
    init();
  }

  const SzSymmetry::Sector& sector() const {
    return _current_sector;
  }

  /**
   * Calculate n2!/n1!
   */
  int variation(int n1, int n2) {
    int result = 1;
    for(int i=n1; i<=n2; i++) {
      result *= i;
    }
    return result;
  }

  /**
   * Calculate number of combinations:
   * C_k^n = n!/(k!*(n-k)!)
   */
  int C_n_k_i(int n, int k) {
    if((n - k) > k) {
      return variation(n - k + 1, n) / variation(1, k);
    }
    return variation(k + 1, n) / variation(1, n - k);
  }

  int c_n_k(int n, int k) {
    return _c_n_k[n][k];
  }

private:

  /**
   *
   */
  void initState(int ik, std::vector<int>& vec) {
    for(int i=0;i<ik;i++) {
      vec[i] = i;
    }
  }

  /**
   *
   */
  bool next_combination(int n, int k, std::vector < int > &old) {
    for (int i = k - 1; i >= 0; i--) {
      if (old[i] < (n - 1 - k + (i + 1))) {
        old[i] += 1;
        for (int j = i + 1; j < k; j++) {
          old[j] = old[j - 1] + 1;
        }
        return true;
      }
    }
    return false;
  }

  int next_basis(int n, int k, std::vector < int > &old, bool start){
    int res = 0;
    if(start) {
      initState(k, old);
      //pass
    }
    else {
      next_combination(n, k, old);
    }
    for (int i = 0; i < k; i++) {
      res += (1 << old[i]);
    }
    return res;
  }

  SzSymmetry::Sector _current_sector;
  std::queue<SzSymmetry::Sector> _sectors;
  long long _state;
  int _Ns;
  int _Ip;
  int _ind;
  std::vector<int> upstate;
  std::vector<int> dostate;
  std::vector<std::vector<int> > _c_n_k;
  std::vector<std::vector<int> > basis;
  std::vector<std::vector<int> > ninv;
  bool _first;
};


#endif //HUBBARD_SZCOMBINATION_H
