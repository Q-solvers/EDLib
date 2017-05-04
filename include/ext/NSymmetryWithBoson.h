//
// Created by iskakoff on 20/04/17.
//

#ifndef HUBBARD_NSYMMETRYWITHBOSON_H
#define HUBBARD_NSYMMETRYWITHBOSON_H

#include <edlib/NSymmetry.h>

namespace EDLib {
  namespace Symmetry {

    /**
     * @brief NSymmetryWithBoson class
     *
     * @author iskakoff
     */
    class NSymmetryWithBoson : public NSymmetry {
    public:
      NSymmetryWithBoson(alps::params &p) : NSymmetry(p) {}
      NSymmetryWithBoson(int N, int Nb) : NSymmetry(N) {
        _boson_bits_count = Nb;
        _maximum_bosons = 1<<Nb - 1;
      }

      virtual bool next_state() {
        long long res = 0;
        if (_ind >= sector().size() && _bosons>=_maximum_bosons) {
          return false;
        } else if (_ind >= sector().size() && _bosons < _maximum_bosons) {
          _comb.init_state(sector().n(), _totstate);
          ++_bosons;
        }
        res = next_basis(_N, sector().n(), _totstate);
        _ind++;
        state() = res<<;
        return true;
      }

      virtual int index(long long st) {
        return _comb.c_n_k(_N, sector().n()) - num(st, _N, sector().n()) - 1;
      }

      virtual void reset() {
        state() = 0ll;
        _first = true;
        _ind = 0;
        _bosons = 0;
      }

      virtual void init() {
        reset();
        _comb.init_state(sector().n(), _totstate);
      }

      int maximum_bosons() const {
        return _maximum_bosons;
      }

    private:
      int _bosons;
      int _boson_bits_count;
      int _maximum_bosons;
    };
  }
}


#endif //HUBBARD_NSYMMETRYWITHBOSON_H
