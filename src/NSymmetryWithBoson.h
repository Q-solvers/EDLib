//
// Created by iskakoff on 21/08/16.
//

#ifndef HUBBARD_NSYMMETRYWITHBOSON_H
#define HUBBARD_NSYMMETRYWITHBOSON_H


#include <queue>
#include "Symmetry.h"
#include "Combination.h"

class NSymmetryWithBoson: public Symmetry {
public:
  class Sector {
  public:
    friend class NSymmetryWithBoson;
    friend std::ostream& operator << (std::ostream& o, const NSymmetryWithBoson::Sector& c) { return o << " (nup+ndown: " <<c._n<<") size: "<<c._size;  }
    Sector(int n, size_t size) : _n(n), _size(size)  {};
  protected:
    bool operator<(Sector s) {
      return _size< s._size;
    }

  public:
    int n() const {return _n;}

    size_t size() const {return _size;}

    void print() const {
      std::cout<<_n;
    }

  private:
    int _n;
    size_t _size;
  };

  NSymmetryWithBoson(EDParams &p) : Symmetry(p), _Nf(2*int(p["NSITES"])), _Nb(p["NSITES_BOSE"]), _totstate(_Nf, 0.0), _current_sector(-1, 0),
                           _comb(_Nf) {
    if(p.exists("arpack.SECTOR") && bool(p["arpack.SECTOR"])) {
      std::vector<std::vector<int> > sectors;
      std::string input = p["INPUT_FILE"];
      alps::hdf5::archive input_file(input, "r");
      input_file>>alps::make_pvp("sectors/values", sectors);
      input_file.close();
      for(auto& sector : sectors) {
        _sectors.push(NSymmetryWithBoson::Sector(sector[0], (size_t) (_comb.c_n_k(_Nf, sector[0]) * (1<<_Nb))));
      }
    } else {
      for(int i = 0; i<=_Nf;++i) {
        _sectors.push(NSymmetryWithBoson::Sector(i, (size_t) (_comb.c_n_k(_Nf, i) * (1<<_Nb))));
      }
    }
  }

  virtual bool next_state() override {
    long long res = 0;
    if(_ind>=_current_sector.size()) {
      return false;
    }
    int f = 0, b = 0;
    int nb2 = 1 << _Nb;
    // get current bosonic state
    b = _ind % ((1<<_Nb));
    // get current fermionic state
    f = (_ind/((1<<_Nb))) % _comb.c_n_k(_Nf, _current_sector.n());
    res=next_basis(_Nf, _current_sector.n(), _totstate, f == 0, b == 0);
    res<<=_Nb;
    res+=b;
    state() = res;
    _ind++;
    return true;
  }

  virtual int index(long long st) override {
    long long f = st>>(_Nb);
    int bos = st & ((1<<_Nb)-1);
    int cup = _comb.c_n_k(_Nf, _current_sector.n())*(1<<_Nb);
    return cup - ((num(f, _Nf, _current_sector.n())<<_Nb) - bos) - ((1<<_Nb) - 1) - 1;
  }

  virtual void reset() override {
    state() = 0ll;
    _first = true;
    _ind = 0;
  }

  virtual void init() override {
    reset();
    initState(_current_sector.n(), _totstate);
  }

  virtual bool next_sector() override {
    if(_sectors.empty())
      return false;
    _current_sector = _sectors.front();
    _sectors.pop();
    std::cout<<"Diagonalizating N-symmetry sector with total number of particles: "<<_current_sector.n()<<" size: "<<_current_sector.size()<<std::endl;
    return true;
  }

private:

  int _Nf;
  int _Nb;
  bool _first;
  int _ind;
  NSymmetryWithBoson::Sector _current_sector;
  std::queue<NSymmetryWithBoson::Sector> _sectors;
  std::vector<int> _totstate;
  Combination _comb;

  int next_basis(int n, int k, std::vector < int > &old, bool u, bool b){
    int res = 0;
    if(u) {
      initState(k, old);
      //pass
    }
    else if(b) {
      next_combination(n, k, old);
    }
    for (int i = 0; i < k; i++) {
      res += (1 << old[i]);
    }
    return res;
  }
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

  inline const int num(long long b, int n, int m) const {
    int res = 0;
    if (((b & (1ll << (_Nf - n))) == 0) and ((n - 1) > 0) and (m > 0) and (m < n))
      res = num(b, n - 1, m);
    else if (((n - 1) > 0) and (m > 0) and (m < n))
      res = _comb.c_n_k(n - 1, m) + num(b, n - 1, m - 1);
    return res;
  }
};


#endif //HUBBARD_NSYMMETRYWITHBOSON_H
