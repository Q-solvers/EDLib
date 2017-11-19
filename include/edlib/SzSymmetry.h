//
// Created by iskakoff on 19/07/16.
//

#ifndef HUBBARD_SZCOMBINATION_H
#define HUBBARD_SZCOMBINATION_H

#include <queue>

#include "Symmetry.h"
#include "Combination.h"
#include "NSymmetry.h"

namespace EDLib {
  namespace Symmetry {
/**
 * Sz symmetry class
 */
    class SzSymmetry : public Symmetry {
    public:
      class Sector {
      public:
        friend class SzSymmetry;

        friend std::ostream &operator<<(std::ostream &o, const SzSymmetry::Sector &c) { return o << " (nup: " << c._nup << " ndown: " << c._ndown << ") size: " << c._size; }

        Sector(int up, int down, size_t size) : _nup(up), _ndown(down), _size(size) {};

        /**
         * @return number of electrons with spin-up
         */
        int nup() const { return _nup; }
        /**
         * @return number of electrons with spin-down
         */
        int ndown() const { return _ndown; }

        /**
         * @return sector dimension
         */
        size_t size() const { return _size; }

        void print() const {
          std::cout << _nup << " " << _ndown;
        }
      protected:
        bool operator<(Sector s) {
          return _size < s._size;
        }

      private:
        int _nup;
        int _ndown;
        size_t _size;
      };

      SzSymmetry(int N) :Symmetry(), _current_sector(-1, -1, 0), _Ns(N), upstate(N + 1), dostate(N + 1),
                         _comb(N), basis(N + 1), ninv(N + 1, std::vector < int >(1 << N, 0)),
                         _first(true) {
        initial_fill();
      };

      SzSymmetry(alps::params &p) : Symmetry(), _current_sector(-1, -1, 0), _Ns(p["NSITES"]), upstate(_Ns + 1), dostate(_Ns + 1),
                                _comb(_Ns), basis(_Ns + 1), ninv(_Ns + 1, std::vector < int >(1 << _Ns, 0)),
                                _first(true) {
        initial_fill();
        if (p.exists("arpack.SECTOR") && bool(p["arpack.SECTOR"])) {
          std::vector < std::vector < int > > sectors;
          std::string input = p["INPUT_FILE"];
          alps::hdf5::archive input_file(input, "r");
          input_file >> alps::make_pvp("sectors/values", sectors);
          input_file.close();
          for (int i = 0; i < sectors.size(); ++i) {
            _sectors.push(SzSymmetry::Sector(sectors[i][0], sectors[i][1], (size_t) (_comb.c_n_k(_Ns, sectors[i][0]) * _comb.c_n_k(_Ns, sectors[i][1]))));
          }
        } else {
          for (int i = 0; i <= _Ns; ++i) {
            for (int j = 0; j <= _Ns; ++j) {
              _sectors.push(SzSymmetry::Sector(i, j, (size_t) (_comb.c_n_k(_Ns, i) * _comb.c_n_k(_Ns, j))));
            }
          }
        }
      }

      virtual ~SzSymmetry() {};

      virtual bool next_state() override {
        if (_first) {
          _first = false;
        }
        if (_ind >= _current_sector.size()) {
          return false;
        }
        state() = state_by_index(_ind);
        _ind++;
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

      virtual int index(long long state, const SzSymmetry::Sector &sector) {
        long long up = state >> _Ns;
        long long down = state & ((1ll << _Ns) - 1);
        int cup = _comb.c_n_k(_Ns, sector.nup());
        int cdo = _comb.c_n_k(_Ns, sector.ndown());
        return ninv[sector.nup()][(int) up] * (cdo) + ninv[sector.ndown()][(int) down];
      }

      virtual int index(long long state)  override{
        return index(state, _current_sector);
      }

      virtual void reset()  override{
        state() = 0ll;
        _first = true;
        _ind = 0;
      }

      virtual void init()  override {
        // TODO: Decide what we should have to init
        reset();
        _comb.init_state(_current_sector.nup(), upstate);
        _comb.init_state(_current_sector.ndown(), dostate);
      };

      virtual bool next_sector() override {
        if (_sectors.empty())
          return false;
        _current_sector = _sectors.front();
        _sectors.pop();
        return true;
      }

      void set_sector(const SzSymmetry::Sector &sector) {
        _current_sector = sector;
        init();
      }

      const SzSymmetry::Sector &sector() const {
        return _current_sector;
      }

      inline const Combination &comb() const {
        return _comb;
      }

      bool can_create_particle(int spin) override {
        return spin == 0 ? _current_sector.nup() < _Ns : _current_sector.ndown() < _Ns;
      }

      bool can_destroy_particle(int spin) override {
        return spin == 0 ? _current_sector.nup() > 0 : _current_sector.ndown() > 0;
      }

      SzSymmetry::Sector destroy_partice(int spin) {
        return Sector(_current_sector.nup() - (1 - spin), _current_sector.ndown() - spin,
                      _comb.c_n_k(_Ns, _current_sector.nup() - (1 - spin)) * _comb.c_n_k(_Ns, _current_sector.ndown() - spin));
      }

      SzSymmetry::Sector create_partice(int spin) {
        return Sector(_current_sector.nup() + (1 - spin), _current_sector.ndown() + spin,
                      _comb.c_n_k(_Ns, _current_sector.nup() + (1 - spin)) * _comb.c_n_k(_Ns, _current_sector.ndown() + spin));
      }

#ifdef USE_MPI
      void set_offset(size_t offset) {_ind += offset;}
#endif
    private:
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
      };

      int next_basis(int n, int k, std::vector < int > &old, bool start) {
        int res = 0;
        if (start) {
          _comb.init_state(k, old);
          //pass
        } else {
          _comb.next_combination(n, k, old);
        }
        for (int i = 0; i < k; i++) {
          res += (1 << old[i]);
        }
        return res;
      }

      SzSymmetry::Sector _current_sector;
      std::queue < SzSymmetry::Sector > _sectors;
      int _Ns;
      int _Ip;
      int _ind;
      std::vector < int > upstate;
      std::vector < int > dostate;
      std::vector < std::vector < int > > basis;
      std::vector < std::vector < int > > ninv;
      Combination _comb;
      bool _first;
    protected:
    public:
      std::queue<SzSymmetry::Sector> &sectors() {
        return _sectors;
      }
    };
  }
}

#endif //HUBBARD_SZCOMBINATION_H
