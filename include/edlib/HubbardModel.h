//
// Created by iskakoff on 29/07/16.
//

#ifndef HUBBARD_HUBBARDMODEL_H
#define HUBBARD_HUBBARDMODEL_H

#include <vector>

#include <alps/params.hpp>
#include "SzSymmetry.h"
#include "FermionicModel.h"

namespace EDLib {
  namespace Model {
    namespace Hubbard {
      template<typename prec>
      class InnerState {
      public:
        InnerState(int ii, int jj, int spin, prec val) : _indicies(ii, jj), _spin(spin), _value(val) {};

        const inline std::pair < int, int > &indicies() const { return _indicies; }

        const inline prec &value() const { return _value; }

        inline int spin() const { return _spin; }

      private:
        std::pair < int, int > _indicies;
        int _spin;
        prec _value;
      };
    }

    template<typename prcsn>
    class HubbardModel: public FermionicModel {
    public:
      typedef prcsn precision;
      typedef typename Symmetry::SzSymmetry SYMMETRY;
      typedef typename Hubbard::InnerState < precision > St;
      typedef typename Symmetry::SzSymmetry::Sector Sector;

      HubbardModel(alps::params &p) : FermionicModel(p), _symmetry(p) {
        Eps.assign(p["NSITES"], std::vector < precision >(p["NSPINS"], precision(0.0)));
        t.assign(p["NSITES"], std::vector < precision >(p["NSITES"], precision(0.0)));
        U.assign(p["NSITES"], precision(0.0));
        _xmu.assign(p["NSITES"], precision(0.0));
        _Hmag = 0.0;
        std::string input = p["INPUT_FILE"];
        alps::hdf5::archive input_data(input.c_str(), "r");
        input_data >> alps::make_pvp("BETA", _beta);
        if(input_data.is_data("magnetic_field")) {
          input_data >> alps::make_pvp("magnetic_field", _Hmag);
        } else{
          _Hmag = 0.0;
        }

        input_data >> alps::make_pvp("hopping/values", t);
        input_data >> alps::make_pvp("interaction/values", U);
        input_data >> alps::make_pvp("chemical_potential/values", _xmu);
        input_data.close();
        for (int ii = 0; ii < _Ns; ++ii) {
          for (int jj = 0; jj < _Ns; ++jj) {
            if (std::abs(t[ii][jj]) > 1e-10) {
              for (int is = 0; is < _ms; ++is) {
                _states.push_back(St(ii, jj, is, t[ii][jj]));
              }
            }
          }
        }
      };

      inline int valid(const St &state, long long nst) {
        return (checkState(nst, state.indicies().first + state.spin() * _Ns, _Ip) * (1 - checkState(nst, state.indicies().second + state.spin() * _Ns, _Ip)));
      }

      inline void set(const St &state, long long nst, long long &k, int &sign) {
        long long k1, k2;
        int isign1, isign2;
        a(state.indicies().first + state.spin() * _Ns, nst, k1, isign1);
        adag(state.indicies().second + state.spin() * _Ns, k1, k2, isign2);
        k = k2;
        // -t c^+ c
        sign = -isign1 * isign2;
      }

      inline precision diagonal(long long state) const {
        precision xtemp = 0.0;
        for (int im = 0; im < _Ns; ++im) {
          for (int is = 0; is < _ms; ++is) {
            xtemp += (Eps[im][is] - _xmu[is]) * checkState(state, im + is * _Ns, _Ip);
          }
          xtemp += U[im] * checkState(state, im, _Ip) * checkState(state, im + _Ns, _Ip);
          xtemp += _Hmag * (checkState(state, im + _Ns, _Ip) - checkState(state, im, _Ip));
        }
        return xtemp;
      }

      inline long long interacting_states(long long nst) {
        return nst;
      }

      const std::vector < St > &T_states() const { return _states; };
      // We have only diagonal interaction
      const std::vector < St > &V_states() const { return _V_states; };

      int interacting_orbitals() const {
        return _Ns;
      }

      inline const Symmetry::SzSymmetry &symmetry() const {
        return _symmetry;
      }

      inline Symmetry::SzSymmetry &symmetry() {
        return _symmetry;
      }

    private:
      // Symmetry
      Symmetry::SzSymmetry _symmetry;
      // Hopping
      std::vector < std::vector < precision > > t;
      // Interaction
      std::vector < precision > U;
      // Chemical potential
      std::vector < precision > _xmu;
      // Magnetic field
      precision _Hmag;
      // site energy shift
      std::vector < std::vector < precision > > Eps;
      // Inverse temperature
      double _beta;

      // Non-diagonal states iterator
      std::vector < St > _states;
      std::vector < St > _V_states;
    };

  }
}
#endif //HUBBARD_HUBBARDMODEL_H
