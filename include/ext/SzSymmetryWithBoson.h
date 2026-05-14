#ifndef EDLIB_EXT_SZSYMMETRYWITHBOSON_H
#define EDLIB_EXT_SZSYMMETRYWITHBOSON_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <queue>
#include <vector>

#include "edlib/Combination.h"
#include "edlib/Parameters.h"
#include "edlib/Symmetry.h"

namespace edlib { namespace ext {

  /**
   * Sz-conserving fermionic symmetry augmented with a bosonic state encoded
   * in the low `boson_bits_count` bits of each basis state. Sector
   * dimensions are size_fermion(nup, ndown) * 2^boson_bits_count.
   */
  class SzSymmetryWithBoson : public Symmetry {
  public:
    class Sector {
    public:
      friend class SzSymmetryWithBoson;

      friend std::ostream& operator<<(std::ostream& o, const Sector& c) {
        return o << " (nup: " << c._nup << " ndown: " << c._ndown
                 << " bosons cutoff: " << ((1 << c._bosons) - 1)
                 << ") size: " << c._size;
      }

      Sector(int up, int down, int bosonic_bits, std::size_t size)
          : _nup(up), _ndown(down), _size(size), _bosons(bosonic_bits) {}

      int         nup()   const { return _nup; }
      int         ndown() const { return _ndown; }
      std::size_t size()  const { return _size; }
      std::size_t fsize() const { return _size / (1 << _bosons); }

      void print() const { std::cout << _nup << " " << _ndown; }

      bool operator<(const Sector& s) const { return _size < s._size; }

    private:
      int         _nup;
      int         _ndown;
      std::size_t _size;
      int         _bosons;
    };

    SzSymmetryWithBoson(int N, int boson_bits_count,
                        const std::vector<std::array<int, 2>>& sector_list = {})
        : Symmetry(),
          _current_sector(-1, -1, boson_bits_count, 0),
          _Ns(N),
          _boson_bits_count(boson_bits_count),
          _maximum_bosons((1 << boson_bits_count) - 1),
          upstate(N + 1), dostate(N + 1),
          _comb(N), basis(N + 1),
          ninv(N + 1, std::vector<int>(1 << N, 0)),
          _first(true) {
      initial_fill();
      populate_sectors(sector_list);
    }

    /**
     * @param p Core parameters (uses p.nsites).
     * @param boson_bits_count total bits reserved for the bosonic state
     *        (legacy: NBBITS * NBLEVEL).
     * @param sector_list optional sector restriction; empty = all sectors.
     */
    SzSymmetryWithBoson(const Parameters& p,
                        int boson_bits_count,
                        const std::vector<std::array<int, 2>>& sector_list = {})
        : SzSymmetryWithBoson(p.nsites, boson_bits_count, sector_list) {}

    bool next_state() override {
      int u = (_ind % _current_sector.fsize()) / _comb.c_n_k(_Ns, _current_sector.ndown());
      int d = (_ind % _current_sector.fsize()) % _comb.c_n_k(_Ns, _current_sector.ndown());
      long long st = basis[_current_sector.nup()][u];
      st <<= _Ns;
      st += basis[_current_sector.ndown()][d];
      _bosons = _ind / _current_sector.fsize();
      st <<= _boson_bits_count;
      st += _bosons;
      state() = st;
      ++_ind;
      return _ind <= static_cast<int>(_current_sector.size());
    }

    int index(long long st) override { return index(st, _current_sector); }

    int index(long long state, const Sector& sector) {
      long long fnst = state >> _boson_bits_count;
      long long up   = fnst >> _Ns;
      long long down = fnst & ((1ll << _Ns) - 1);
      int cdo = _comb.c_n_k(_Ns, sector.ndown());
      int fermions_index = ninv[sector.nup()][static_cast<int>(up)] * cdo
                         + ninv[sector.ndown()][static_cast<int>(down)];
      int bosons = state & _maximum_bosons;
      return sector.fsize() * bosons + fermions_index;
    }

    void reset() override {
      state() = 0ll;
      _ind    = 0;
      _bosons = 0;
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

    int maximum_bosons() const { return _maximum_bosons; }

    bool can_create_particle(int spin) override {
      return spin == 0 ? _current_sector.nup()   < _Ns - 1
                       : _current_sector.ndown() < _Ns - 1;
    }
    bool can_destroy_particle(int spin) override {
      return spin == 0 ? _current_sector.nup()   > 0
                       : _current_sector.ndown() > 0;
    }

    Sector destroy_particle(int spin) {
      return Sector(_current_sector.nup()   - (1 - spin),
                    _current_sector.ndown() - spin,
                    _boson_bits_count,
                    _comb.c_n_k(_Ns, _current_sector.nup()   - (1 - spin))
                  * _comb.c_n_k(_Ns, _current_sector.ndown() - spin)
                  * (_maximum_bosons + 1));
    }
    Sector create_particle(int spin) {
      return Sector(_current_sector.nup()   + (1 - spin),
                    _current_sector.ndown() + spin,
                    _boson_bits_count,
                    _comb.c_n_k(_Ns, _current_sector.nup()   + (1 - spin))
                  * _comb.c_n_k(_Ns, _current_sector.ndown() + spin)
                  * (_maximum_bosons + 1));
    }

    const Combination&   comb()    const { return _comb; }
    std::queue<Sector>&  sectors()       { return _sectors; }

#ifdef USE_MPI
    void set_offset(std::size_t offset) { _ind += offset; }
#endif

  private:
    void populate_sectors(const std::vector<std::array<int, 2>>& sector_list) {
      if (!sector_list.empty()) {
        for (const auto& s : sector_list) {
          _sectors.push(Sector(s[0], s[1], _boson_bits_count,
              static_cast<std::size_t>(_comb.c_n_k(_Ns, s[0])
                                     * _comb.c_n_k(_Ns, s[1])
                                     * (_maximum_bosons + 1))));
        }
      } else {
        for (int i = 0; i <= _Ns; ++i) {
          for (int j = 0; j <= _Ns; ++j) {
            _sectors.push(Sector(i, j, _boson_bits_count,
                static_cast<std::size_t>(_comb.c_n_k(_Ns, i)
                                       * _comb.c_n_k(_Ns, j)
                                       * (_maximum_bosons + 1))));
          }
        }
      }
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

    Sector                            _current_sector;
    std::queue<Sector>                _sectors;
    int                               _Ns;
    int                               _Ip = 0;
    int                               _ind = 0;
    int                               _bosons = 0;
    int                               _boson_bits_count;
    int                               _maximum_bosons;

    std::vector<int>                  upstate;
    std::vector<int>                  dostate;
    std::vector<std::vector<int>>     basis;
    std::vector<std::vector<int>>     ninv;
    Combination                       _comb;
    bool                              _first;
  };

}}  // namespace edlib::ext

#endif
