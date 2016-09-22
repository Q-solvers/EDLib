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
      protected:
        bool operator<(Sector s) {
          return _size < s._size;
        }

      public:
        int nup() const { return _nup; }

        int ndown() const { return _ndown; }

        size_t size() const { return _size; }

        void print() const {
          std::cout << _nup << " " << _ndown;
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

      SzSymmetry(EDParams &p) : Symmetry(), _current_sector(-1, -1, 0), _Ns(p["NSITES"]), upstate(_Ns + 1), dostate(_Ns + 1),
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

      virtual bool next_state() {
        long long res = 0;
        if (_first) {
          _first = false;
        }
        if (_ind >= _current_sector.size()) {
          return false;
        }
        int u = 0, d = 0;
        u = _ind / _comb.c_n_k(_Ns, _current_sector.ndown());
        d = _ind % _comb.c_n_k(_Ns, _current_sector.ndown());
        res = basis[_current_sector.nup()][u];
        res <<= _Ns;
        res += basis[_current_sector.ndown()][d];
        _ind++;
        state() = res;
        return true;
      }

      virtual long long state_indexed(int ind) {
        long long res = 0;
/* FIXME ?
        if (ind >= _current_sector.size()) {
          throw;
        }
*/
        int u = 0, d = 0;
        u = ind / _comb.c_n_k(_Ns, _current_sector.ndown());
        d = ind % _comb.c_n_k(_Ns, _current_sector.ndown());
        res = basis[_current_sector.nup()][u];
        res <<= _Ns;
        res += basis[_current_sector.ndown()][d];
        return res;
      }

      int index(long long state, const SzSymmetry::Sector &sector) {
        long long up = state >> _Ns;
        long long down = state & ((1ll << _Ns) - 1);
        int cup = _comb.c_n_k(_Ns, sector.nup());
        int cdo = _comb.c_n_k(_Ns, sector.ndown());
        return ninv[sector.nup()][(int) up] * (cdo) + ninv[sector.ndown()][(int) down];
      }

      virtual int index(long long state) {
        return index(state, _current_sector);
      }

      virtual void reset() {
        state() = 0ll;
        _first = true;
        _ind = 0;
      }

      virtual void init() {
        // TODO: Decide what we should have to init
        reset();
        _comb.init_state(_current_sector.nup(), upstate);
        _comb.init_state(_current_sector.ndown(), dostate);
      };

      virtual bool next_sector() {
        if (_sectors.empty())
          return false;
        _current_sector = _sectors.front();
        _sectors.pop();
        std::cout << "Diagonalizating Sz-symmetry sector with nup: " << _current_sector.nup() << " ndown: " << _current_sector.ndown() << " size: " << _current_sector.size()
                  << std::endl;
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
#ifdef ALPS_HAVE_MPI
      void set_offset(size_t offset) {_ind += offset;}
#endif
    private:
      void initial_fill() {
        _Ip = 2 * _Ns;
        _ind = 0;
        for (int i = 0; i <= _Ns; ++i) {
          int cnk = _comb.c_n_k(_Ns, i);
          basis[i].reserve(cnk);
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
    class ImpuritySzSymmetry : public SzSymmetry {
    public:
      ImpuritySzSymmetry(EDParams &p, int Ns, int ml) : SzSymmetry(p), _up_symmetry(Ns-ml), _down_symmetry(Ns-ml), _imp_symmetry(ml), _ml(ml), _Nk(Ns-ml) {
      };
      virtual bool next_state() {
        long long res = 0;
        if(_imp_max_down==_imp_down && _imp_max_up==_imp_up && !_imp_symmetry.next_state())
          return false;
        if(!_imp_symmetry.next_state()) {

        }
        state() = res;
        return true;
      }

      virtual void reset() {
        state() = 0ll;
        _up_symmetry.reset();
        _down_symmetry.reset();
        _imp_symmetry.reset();
//        _first = true;
//        _ind = 0;
      }

      virtual void init() {
        reset();
        int ndown = sector().ndown();
        int nup = sector().nup();
        _up_symmetry.sectors();
        _imp_max_down = std::min(_ml, ndown);
        _imp_max_up = std::min(_ml, nup);
        _up_symmetry.set_sector(NSymmetry::Sector(std::min(nup, _Nk), _up_symmetry.comb().c_n_k(_Nk, std::min(nup, _Nk))));
        _down_symmetry.set_sector(NSymmetry::Sector(std::min(ndown, _Nk), _down_symmetry.comb().c_n_k(_Nk, std::min(ndown, _Nk))));
        _imp_symmetry.set_sector(SzSymmetry::Sector(std::max(0, nup - _Nk), std::max(0, ndown-_Nk),
                                                    _imp_symmetry.comb().c_n_k(_Nk, std::min(ndown, _Nk))));
      };
      
      

    private:
      NSymmetry _up_symmetry;
      NSymmetry _down_symmetry;
      SzSymmetry _imp_symmetry;
      int _ml;
      int _Nk;
      int _imp_up;
      int _imp_down;
      int _imp_max_up;
      int _imp_max_down;
    };
  }
}

#endif //HUBBARD_SZCOMBINATION_H
