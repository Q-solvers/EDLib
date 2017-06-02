//
// Created by iskakoff on 20/04/17.
//

#ifndef HUBBARD_NSYMMETRYWITHBOSON_H
#define HUBBARD_NSYMMETRYWITHBOSON_H

#include <edlib/Symmetry.h>
#include <edlib/Combination.h>
#include <queue>
#include <edlib/HDF5Utils.h>

namespace EDLib {
  namespace Symmetry {

    /**
     * @brief SzSymmetryWithBoson class
     *
     * @author iskakoff
     */
    class SzSymmetryWithBoson : public Symmetry {
    public:
      class Sector {
      public:
        friend class SzSymmetryWithBoson;

        friend std::ostream &operator<<(std::ostream &o, const SzSymmetryWithBoson::Sector &c) { return o << " (nup: " << c._nup << " ndown: " << c._ndown<<" bosons cutoff: "<< (1<<c._bosons)-1 << ") size: " << c._size; }

        Sector(int up, int down, int bosonic_bits, size_t size) : _nup(up), _ndown(down), _size(size), _bosons(bosonic_bits) {};

        int nup() const { return _nup; }

        int ndown() const { return _ndown; }

        size_t size() const { return _size; }
        size_t fsize() const { return _size/(1<<_bosons); }

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
        int _bosons;
      };


      SzSymmetryWithBoson(alps::params &p) : _boson_bits_count(p["NBBITS"].as<int>() * p["NBLEVEL"].as<int>()),
                                             _maximum_bosons((1<<(p["NBBITS"].as<int>() * p["NBLEVEL"].as<int>())) - 1),
                                             _current_sector(-1, -1, 0, 0), _Ns(p["NSITES"]), upstate(_Ns + 1), dostate(_Ns + 1),
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
            _sectors.push(SzSymmetryWithBoson::Sector(sectors[i][0], sectors[i][1], _boson_bits_count, (size_t) (_comb.c_n_k(_Ns, sectors[i][0]) * _comb.c_n_k(_Ns, sectors[i][1]) * (_maximum_bosons+1))));
          }
        } else {
          for (int i = 0; i <= _Ns; ++i) {
            for (int j = 0; j <= _Ns; ++j) {
              _sectors.push(SzSymmetryWithBoson::Sector(i, j, _boson_bits_count, (size_t) (_comb.c_n_k(_Ns, i) * _comb.c_n_k(_Ns, j) * (_maximum_bosons + 1))));
            }
          }
        }
      }
      SzSymmetryWithBoson(int N, int Nb) : _current_sector(-1, -1, 0, 0), _Ns(N), upstate(N + 1), dostate(N + 1),
                                           _comb(N), basis(N + 1), ninv(N + 1, std::vector < int >(1 << N, 0)),
                                           _first(true) {
        _boson_bits_count = Nb;
        _maximum_bosons = (1<<Nb) - 1;
        initial_fill();
      }

      virtual bool next_state()  override {
        int u = (_ind % (sector().fsize())) / _comb.c_n_k(_Ns, _current_sector.ndown());
        int d = (_ind % (sector().fsize())) % _comb.c_n_k(_Ns, _current_sector.ndown());
        long long st = basis[_current_sector.nup()][u];
        st <<= _Ns;
        st += basis[_current_sector.ndown()][d];
        _bosons = _ind / sector().fsize();
        st<<=_boson_bits_count;
        st+= _bosons;
        state() = st;
        _ind++;
        return _ind <= sector().size();
      }

      virtual int index(long long st) override {
        return index(st, sector());
      }

      virtual int index(long long state, const SzSymmetryWithBoson::Sector &sector) {
        long long fnst = state >> _boson_bits_count;
        long long up = fnst >> _Ns;
        long long down = fnst & ((1ll << _Ns) - 1);
        int cup = _comb.c_n_k(_Ns, sector.nup());
        int cdo = _comb.c_n_k(_Ns, sector.ndown());
        int fermions_index = ninv[sector.nup()][(int) up] * (cdo) + ninv[sector.ndown()][(int) down];
        int bosons = state & (_maximum_bosons);
        return sector.fsize() * bosons + fermions_index;
      }

      virtual void reset() override {
        state() = 0ll;
        _ind = 0;
        _bosons = 0;
      }

      virtual void init() override {
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

      void set_sector(const SzSymmetryWithBoson::Sector &sector) {
        _current_sector = sector;
        init();
      }

      const SzSymmetryWithBoson::Sector &sector() const {
        return _current_sector;
      }

      inline const Combination &comb() const {
        return _comb;
      }
#ifdef USE_MPI
      void set_offset(size_t offset) {_ind += offset;}
#endif

      int maximum_bosons() const {
        return _maximum_bosons;
      }

      bool can_create_particle(int spin) override {
        return spin == 0 ? _current_sector.nup() < _Ns - 1 : _current_sector.ndown() < _Ns - 1;
      }

      bool can_destroy_particle(int spin) override {
        return spin == 0 ? _current_sector.nup() > 0 : _current_sector.ndown() > 0;
      }

      SzSymmetryWithBoson::Sector destroy_partice(int spin) {
        return Sector(_current_sector.nup() - (1 - spin), _current_sector.ndown() - spin, _boson_bits_count,
                      _comb.c_n_k(_Ns, _current_sector.nup() - (1 - spin)) * _comb.c_n_k(_Ns, _current_sector.ndown() - spin) * (_maximum_bosons  + 1));
      }

      SzSymmetryWithBoson::Sector create_partice(int spin) {
        return Sector(_current_sector.nup() + (1 - spin), _current_sector.ndown() + spin, _boson_bits_count,
                      _comb.c_n_k(_Ns, _current_sector.nup() + (1 - spin)) * _comb.c_n_k(_Ns, _current_sector.ndown() + spin)* (_maximum_bosons  + 1));
      }

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

      SzSymmetryWithBoson::Sector _current_sector;
      std::queue < SzSymmetryWithBoson::Sector > _sectors;
      int _Ns;
      int _Ip;
      int _ind;
      int _bosons;
      int _boson_bits_count;
      int _maximum_bosons;

      std::vector < int > upstate;
      std::vector < int > dostate;
      std::vector < std::vector < int > > basis;
      std::vector < std::vector < int > > ninv;
      Combination _comb;
      bool _first;
    protected:
    public:
      std::queue<SzSymmetryWithBoson::Sector> &sectors() {
        return _sectors;
      }
    };

  }

  namespace hdf5 {
    template<>
    void EDLib::hdf5::HDF5Utils<typename EDLib::Symmetry::SzSymmetryWithBoson::Sector>::save(const typename EDLib::Symmetry::SzSymmetryWithBoson::Sector& s, alps::hdf5::archive & ar, const std::string& path) {
      ar[path + "/nup"]<<s.nup();
      ar[path + "/ndown"]<<s.ndown();
      ar[path + "/size"]<<s.size();
    }
  }
}


#endif //HUBBARD_NSYMMETRYWITHBOSON_H
