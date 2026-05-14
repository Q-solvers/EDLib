#ifndef EDLIB_SZSYMMETRY_H
#define EDLIB_SZSYMMETRY_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <queue>
#include <vector>

#include "edlib/Combination.h"
#include "edlib/Parameters.h"
#include "edlib/Symmetry.h"

namespace edlib {

  /**
   * Sz-conserving symmetry: sectors indexed by (n_up, n_down).
   *
   * Sector restrictions, if any, are passed as an array of {n_up, n_down}
   * pairs at construction. An empty list (default) means "all sectors".
   */
  class SzSymmetry : public Symmetry {
  public:
    class Sector {
    public:
      friend class SzSymmetry;

      friend std::ostream& operator<<(std::ostream& o, const Sector& c) {
        return o << " (nup: " << c._nup << " ndown: " << c._ndown
                 << ") size: " << c._size;
      }

      Sector(int up, int down, std::size_t size)
          : _nup(up), _ndown(down), _size(size) {}

      int         nup()   const { return _nup; }
      int         ndown() const { return _ndown; }
      std::size_t size()  const { return _size; }

      void print(std::ostream& out) const { out << _nup << " " << _ndown; }
      void print() const { print(std::cout); }

      bool operator<(const Sector& s) const {
        return _size < s._size
            || (_size == s._size && _nup <  s._nup && _ndown <  s._ndown)
            || (_size == s._size && _nup == s._nup && _ndown <  s._ndown);
      }
      bool operator>(const Sector& s) const { return s < *this; }

    private:
      int         _nup;
      int         _ndown;
      std::size_t _size;
    };

    explicit SzSymmetry(int N,
                        const std::vector<std::array<int, 2>>& sector_list = {})
        : Symmetry(), _current_sector(-1, -1, 0), _Ns(N),
          upstate(N + 1), dostate(N + 1),
          _comb(N), basis(N + 1),
          ninv(N + 1, std::vector<int>(1 << N, 0)),
          _first(true) {
      initial_fill();
      populate_sectors(sector_list);
    }

    SzSymmetry(const Parameters& p,
               const std::vector<std::array<int, 2>>& sector_list = {})
        : SzSymmetry(p.nsites, sector_list) {}

    bool next_state() override {
      if (_first) _first = false;
      if (_ind >= static_cast<int>(_current_sector.size())) return false;
      state() = state_by_index(_ind);
      ++_ind;
      return true;
    }

    inline long long state_by_index(int ind) {
      int u = ind / _comb.c_n_k(_Ns, _current_sector.ndown());
      int d = ind % _comb.c_n_k(_Ns, _current_sector.ndown());
      long long res = basis[_current_sector.nup()][u];
      res <<= _Ns;
      res += basis[_current_sector.ndown()][d];
      return res;
    }

    int index(long long state, const Sector& sector) {
      long long up = state >> _Ns;
      long long down = state & ((1ll << _Ns) - 1);
      int cdo = _comb.c_n_k(_Ns, sector.ndown());
      return ninv[sector.nup()][static_cast<int>(up)] * cdo
           + ninv[sector.ndown()][static_cast<int>(down)];
    }

    int index(long long state) override { return index(state, _current_sector); }

    void reset() override {
      state() = 0ll;
      _first = true;
      _ind   = 0;
    }

    void init() override {
      reset();
      _comb.init_state(_current_sector.nup(),   upstate);
      _comb.init_state(_current_sector.ndown(), dostate);
    }

    bool next_sector() override {
      if (_sectors.empty()) return false;
      _current_sector = _sectors.front();
      _sectors.pop();
      return true;
    }

    void set_sector(const Sector& s) { _current_sector = s; init(); }
    const Sector& sector() const     { return _current_sector; }

    inline const Combination& comb() const { return _comb; }

    bool can_create_particle(int spin) override {
      return spin == 0 ? _current_sector.nup()   < _Ns
                       : _current_sector.ndown() < _Ns;
    }
    bool can_destroy_particle(int spin) override {
      return spin == 0 ? _current_sector.nup()   > 0
                       : _current_sector.ndown() > 0;
    }

    Sector destroy_particle(int spin) {
      return Sector(_current_sector.nup()   - (1 - spin),
                    _current_sector.ndown() - spin,
                    _comb.c_n_k(_Ns, _current_sector.nup()   - (1 - spin))
                  * _comb.c_n_k(_Ns, _current_sector.ndown() - spin));
    }
    Sector create_particle(int spin) {
      return Sector(_current_sector.nup()   + (1 - spin),
                    _current_sector.ndown() + spin,
                    _comb.c_n_k(_Ns, _current_sector.nup()   + (1 - spin))
                  * _comb.c_n_k(_Ns, _current_sector.ndown() + spin));
    }

    std::queue<Sector>& sectors() { return _sectors; }

#ifdef USE_MPI
    void set_offset(std::size_t offset) { _ind += offset; }
#endif

  private:
    void populate_sectors(const std::vector<std::array<int, 2>>& sector_list) {
      std::vector<Sector> sectors;
      if (!sector_list.empty()) {
        for (const auto& s : sector_list) {
          sectors.emplace_back(s[0], s[1],
              static_cast<std::size_t>(_comb.c_n_k(_Ns, s[0]) * _comb.c_n_k(_Ns, s[1])));
        }
      } else {
        for (int i = 0; i <= _Ns; ++i) {
          for (int j = 0; j <= _Ns; ++j) {
            sectors.emplace_back(i, j,
                static_cast<std::size_t>(_comb.c_n_k(_Ns, i) * _comb.c_n_k(_Ns, j)));
          }
        }
      }
      std::sort(sectors.begin(), sectors.end(), std::less<Sector>());
      for (const auto& e : sectors) _sectors.push(e);
    }

    void initial_fill() {
      _Ip = 2 * _Ns;
      _ind = 0;
      for (int i = 0; i <= _Ns; ++i) {
        int cnk = _comb.c_n_k(_Ns, i);
        basis[i].resize(cnk);
        for (int k = 0; k < cnk; ++k) {
          basis[i][k] = next_basis(_Ns, i, upstate, k == 0);
          ninv[i][basis[i][k]] = k;
        }
      }
    }

    int next_basis(int n, int k, std::vector<int>& old, bool start) {
      int res = 0;
      if (start) _comb.init_state(k, old);
      else       _comb.next_combination(n, k, old);
      for (int i = 0; i < k; ++i) res += (1 << old[i]);
      return res;
    }

    Sector              _current_sector;
    std::queue<Sector>  _sectors;
    int                 _Ns;
    int                 _Ip = 0;
    int                 _ind = 0;
    std::vector<int>    upstate;
    std::vector<int>    dostate;
    std::vector<std::vector<int>> basis;
    std::vector<std::vector<int>> ninv;
    Combination         _comb;
    bool                _first;
  };

}

#endif
