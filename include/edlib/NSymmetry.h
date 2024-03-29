//
// Created by iskakoff on 21/08/16.
//

#ifndef EDLIB_NSYMMETRY_H
#define EDLIB_NSYMMETRY_H


#include <queue>
#include <alps/hdf5.hpp>
#include "Symmetry.h"
#include "Combination.h"
namespace EDLib {
  namespace Symmetry {
    class NSymmetry : public Symmetry {
    public:
      class Sector {
      public:
        friend class NSymmetry;

        friend std::ostream &operator<<(std::ostream &o, const NSymmetry::Sector &c) { return o << " (nup+ndown: " << c._n << ") size: " << c._size; }

        Sector(int n, size_t size) : _n(n), _size(size) {};

        bool operator<(const Sector & s) const{
          return _size < s._size || (_size == s._size && _n < s._n);
        }

        bool operator>(const Sector & s) const{
          return s < *this;
        }

      public:
        int n() const { return _n; }

        size_t size() const { return _size; }

        void print(std::ostream & out) const {
          out << _n;
        }
        
        void print() const {print(std::cout);}

      private:
        int _n;
        size_t _size;
      };

      NSymmetry(int N) : Symmetry(), _N(N), _totstate(N, 0.0), _current_sector(-1, 0), _comb(N) {}

      NSymmetry(alps::params &p) : Symmetry(), _N(2 * int(p["NSITES"])), _totstate(2 * p["NSITES"].as<int>(), 0.0), _current_sector(-1, 0),
                               _comb(2 * p["NSITES"].as<int>()) {
        std::vector<Sector> sectors;
        if (p.exists("arpack.SECTOR") && bool(p["arpack.SECTOR"])) {
          std::vector < std::vector < int > > sectors_list;
          std::string input = p["INPUT_FILE"];
          alps::hdf5::archive input_file(input, "r");
          input_file >> alps::make_pvp("sectors/values", sectors_list);
          input_file.close();
          for (int kkk = 0; kkk < sectors_list.size(); ++kkk) {
            sectors.push_back(NSymmetry::Sector(sectors_list[kkk][0], (size_t) (_comb.c_n_k(_N, sectors_list[kkk][0]))));
          }
        } else {
          for (int i = 0; i <= _N; ++i) {
            sectors.push_back(NSymmetry::Sector(i, (size_t) (_comb.c_n_k(_N, i))));
          }
        }
        std::sort(sectors.begin(), sectors.end(), std::less<Sector>());
        for(auto const& e : sectors) {
          _sectors.push(e);
        }
      }

      void set_sector(const NSymmetry::Sector &sector) {
        _current_sector = sector;
        init();
      }

      const Sector &sector() const {
        return _current_sector;
      }

      virtual bool next_state() override {
        long long res = 0;
        if (_ind >= _current_sector.size()) {
          return false;
        }
        res = next_basis(_N, _current_sector.n(), _totstate);
        _ind++;
        state() = res;
        return true;
      }

      virtual int index(long long st) override  {
        return _comb.c_n_k(_N, _current_sector.n()) - num(st, _N, _current_sector.n()) - 1;
      }

      virtual void reset() override {
        state() = 0ll;
        _first = true;
        _ind = 0;
      }

      virtual void init() override {
        reset();
        _comb.init_state(_current_sector.n(), _totstate);
      }

      virtual bool next_sector() override {
        if (_sectors.empty())
          return false;
        _current_sector = _sectors.front();
        _sectors.pop();
        return true;
      }

      const Combination &comb() const {
        return _comb;
      }

      std::queue<NSymmetry::Sector> &sectors() {
        return _sectors;
      }

      bool can_create_particle(int spin) override {
        return _current_sector.n() < _N - 1;
      }

      bool can_destroy_particle(int spin) override {
        return _current_sector.n() > 0;
      }

    protected:
      int _ind;
      std::vector < int > _totstate;
      int _N;
      bool _first;
      Combination _comb;

      int next_basis(int n, int k, std::vector < int > &old) {
        int res = 0;
        if (_first) {
          _comb.init_state(k, old);
          _first = false;
          //pass
        } else {
          _comb.next_combination(n, k, old);
        }
        for (int i = 0; i < k; i++) {
          res += (1 << old[i]);
        }
        return res;
      }

      inline int num(long long b, int n, int m) const {
        int res = 0;
        if (((b & (1ll << (_N - n))) == 0) and ((n - 1) > 0) and (m > 0) and (m < n))
          res = num(b, n - 1, m);
        else if (((n - 1) > 0) and (m > 0) and (m < n))
          res = _comb.c_n_k(n - 1, m) + num(b, n - 1, m - 1);
        return res;
      }

    private:
      NSymmetry::Sector _current_sector;
      std::queue < NSymmetry::Sector > _sectors;
    };
  }
}

#endif //EDLIB_NSYMMETRY_H
