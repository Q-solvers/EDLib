#ifndef EDLIB_NSYMMETRY_H
#define EDLIB_NSYMMETRY_H

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <queue>
#include <vector>

#include "edlib/Combination.h"
#include "edlib/Parameters.h"
#include "edlib/Symmetry.h"

namespace edlib {

  /**
   * Total-N conserving symmetry: sectors indexed by a single integer n.
   *
   * sector_list (default empty = "all sectors") restricts iteration to the
   * supplied n values. _N is 2 * nsites (spinful), matching legacy semantics.
   */
  class NSymmetry : public Symmetry {
  public:
    class Sector {
    public:
      friend class NSymmetry;

      friend std::ostream& operator<<(std::ostream& o, const Sector& c) {
        return o << " (nup+ndown: " << c._n << ") size: " << c._size;
      }

      Sector(int n, std::size_t size) : _n(n), _size(size) {}

      int         n()    const { return _n; }
      std::size_t size() const { return _size; }

      void print(std::ostream& out) const { out << _n; }
      void print() const { print(std::cout); }

      bool operator<(const Sector& s) const {
        return _size < s._size || (_size == s._size && _n < s._n);
      }
      bool operator>(const Sector& s) const { return s < *this; }

    private:
      int         _n;
      std::size_t _size;
    };

    explicit NSymmetry(int N, const std::vector<int>& sector_list = {})
        : Symmetry(), _N(N), _totstate(N, 0), _current_sector(-1, 0), _comb(N) {
      populate_sectors(sector_list);
    }

    NSymmetry(const Parameters& p, const std::vector<int>& sector_list = {})
        : NSymmetry(2 * p.nsites, sector_list) {}

    void set_sector(const Sector& s) { _current_sector = s; init(); }
    const Sector& sector() const     { return _current_sector; }

    bool next_state() override {
      if (_ind >= static_cast<int>(_current_sector.size())) return false;
      long long res = next_basis(_N, _current_sector.n(), _totstate);
      ++_ind;
      state() = res;
      return true;
    }

    int index(long long st) override {
      return _comb.c_n_k(_N, _current_sector.n()) - num(st, _N, _current_sector.n()) - 1;
    }

    void reset() override {
      state() = 0ll;
      _first  = true;
      _ind    = 0;
    }

    void init() override {
      reset();
      _comb.init_state(_current_sector.n(), _totstate);
    }

    bool next_sector() override {
      if (_sectors.empty()) return false;
      _current_sector = _sectors.front();
      _sectors.pop();
      return true;
    }

    const Combination& comb() const { return _comb; }
    std::queue<Sector>& sectors()   { return _sectors; }

    bool can_create_particle(int /*spin*/) override {
      return _current_sector.n() < _N - 1;
    }
    bool can_destroy_particle(int /*spin*/) override {
      return _current_sector.n() > 0;
    }

  protected:
    int              _ind = 0;
    std::vector<int> _totstate;
    int              _N;
    bool             _first = true;
    Combination      _comb;

    int next_basis(int n, int k, std::vector<int>& old) {
      int res = 0;
      if (_first) { _comb.init_state(k, old); _first = false; }
      else        { _comb.next_combination(n, k, old); }
      for (int i = 0; i < k; ++i) res += (1 << old[i]);
      return res;
    }

    inline int num(long long b, int n, int m) const {
      int res = 0;
      if (((b & (1ll << (_N - n))) == 0) && ((n - 1) > 0) && (m > 0) && (m < n))
        res = num(b, n - 1, m);
      else if (((n - 1) > 0) && (m > 0) && (m < n))
        res = _comb.c_n_k(n - 1, m) + num(b, n - 1, m - 1);
      return res;
    }

  private:
    void populate_sectors(const std::vector<int>& sector_list) {
      std::vector<Sector> sectors;
      if (!sector_list.empty()) {
        for (int n : sector_list) {
          sectors.emplace_back(n, static_cast<std::size_t>(_comb.c_n_k(_N, n)));
        }
      } else {
        for (int i = 0; i <= _N; ++i) {
          sectors.emplace_back(i, static_cast<std::size_t>(_comb.c_n_k(_N, i)));
        }
      }
      std::sort(sectors.begin(), sectors.end(), std::less<Sector>());
      for (const auto& e : sectors) _sectors.push(e);
    }

    Sector             _current_sector;
    std::queue<Sector> _sectors;
  };

}

#endif
