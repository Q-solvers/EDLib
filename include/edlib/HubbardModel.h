//
// Created by iskakoff on 29/07/16.
//

#ifndef HUBBARD_HUBBARDMODEL_H
#define HUBBARD_HUBBARDMODEL_H

#include <vector>

#include <alps/params.hpp>
#include <alps/gf/mesh.hpp>
#include <alps/gf/gf.hpp>
#include "SzSymmetry.h"
#include "FermionicModel.h"
#include "CommonUtils.h"

#include <Eigen/Core>

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
        _Eps.assign(p["NSITES"], std::vector < precision >(p["NSPINS"], precision(0.0)));
        t.assign(p["NSITES"], std::vector < precision >(p["NSITES"], precision(0.0)));
        U.assign(p["NSITES"], precision(0.0));
        J.assign(p["NSITES"], std::vector < precision >(p["NSITES"], precision(0.0)));
        _xmu.assign(p["NSITES"], precision(0.0));
        _Hmag.assign(p["NSITES"], precision(0.0));
        std::string input = p["INPUT_FILE"];
        alps::hdf5::archive input_data(input.c_str(), "r");
        if(input_data.is_data("magnetic_field/values")) {
          input_data >> alps::make_pvp("magnetic_field/values", _Hmag);
        }

        input_data >> alps::make_pvp("hopping/values", t);
        input_data >> alps::make_pvp("interaction/values", U);
        if(input_data.is_data("exchange/values")) {
          input_data >> alps::make_pvp("exchange/values", J);
        }
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

      /**
       * check that current basis vector get non-zero contribution
       *
       * @param state - electron state combination of spin and site indices
       * @param nst - current basis state
       * @return 1 if there is nonzero contribution, otherwise 0
       */
      inline int valid(const St &state, long long nst) {
        return (checkState(nst, state.indicies().first + state.spin() * _Ns, _Ip) * (1 - checkState(nst, state.indicies().second + state.spin() * _Ns, _Ip)));
      }

      /**
       * Compute off-diagonal term for transition from nst-state to k-state
       *
       * @param state - transition state
       * @param nst - initial basis state
       * @param k - resulting basis state
       * @param sign - fermionic sign for transition
       * @return contribution to off-diagonal element for transition state
       */
      inline precision set(const St &state, long long nst, long long &k, int &sign) {
        long long k1, k2;
        int isign1, isign2;
        a(state.indicies().first + state.spin() * _Ns, nst, k1, isign1);
        adag(state.indicies().second + state.spin() * _Ns, k1, k2, isign2);
        k = k2;
        // -t c^+ c
        sign = -isign1 * isign2;
        return state.value();
      }

      /**
       * Computes diagonal contribution for state s: <s| H |s>
       *
       * @param state - current basis state
       * @return value of <s| H |s>
       */
      inline precision diagonal(long long state) const {
        precision xtemp = 0.0;
        for (int im = 0; im < _Ns; ++im) {
          for (int is = 0; is < _ms; ++is) {
            xtemp += (_Eps[im][is] - _xmu[is]) * checkState(state, im + is * _Ns, _Ip);
          }
          xtemp += U[im] * checkState(state, im, _Ip) * checkState(state, im + _Ns, _Ip);
          xtemp += _Hmag[im] * (checkState(state, im + _Ns, _Ip) - checkState(state, im, _Ip));
          for (int im2 = 0; im2 < _Ns; ++im2) {
            xtemp +=
              J[im][im2] *
              (checkState(state, im, _Ip) - checkState(state, im + _Ns, _Ip)) *
              (checkState(state, im2, _Ip) - checkState(state, im2 + _Ns, _Ip));
          }
        }
        return xtemp;
      }

      /**
       * @deprecated
       */
      inline long long interacting_states(long long nst) {
        return nst;
      }


      /**
       * @return hopping transitions
       */
      const std::vector < St > &T_states() const { return _states; };
      // We have only diagonal interaction
      /**
       * @return off-diagonal interaction transitions
       */
      const std::vector < St > &V_states() const { return _V_states; };

      /**
       * For Hubbard model all orbitals are interacting
       *
       * @return total number of sites
       */
      int interacting_orbitals() const {
        return _Ns;
      }

      /**
       * For Hubbard model Hamiltonian comutes with spin-operator ([Sz, H] = 0)
       * @return symmetry object
       */
      inline const Symmetry::SzSymmetry &symmetry() const {
        return _symmetry;
      }

      inline Symmetry::SzSymmetry &symmetry() {
        return _symmetry;
      }

      /**
       * Compute bare Green's function for specific mesh
       *
       * @tparam Mesh - mesh-type
       * @param bare_gf - bare Green's function container
       * @param beta - inverse temperature
       */
      template<typename Mesh>
      void bare_greens_function(alps::gf::three_index_gf<std::complex<double>, Mesh, alps::gf::index_mesh, alps::gf::index_mesh >& bare_gf, double beta) {
        for(int iw = 0; iw< bare_gf.mesh1().points().size(); ++iw) {
          typename Mesh::index_type w(iw);
          for (int is : bare_gf.mesh3().points()) {
            Eigen::MatrixXcd G_inv = Eigen::MatrixXcd::Zero(_Ns, _Ns);
            for(int I = 0; I<_Ns; ++I) {
              G_inv(I, I) = (common::freq_point(iw, bare_gf.mesh1(), beta) + _xmu[I] - _Eps[I][is]);
              for(int J = 0; J<_Ns; ++J){
                 int im = I*_Ns + J;
                 G_inv(I, J) += t[I][J];
              }
            }
            G_inv = G_inv.inverse().eval();
            for (int im: bare_gf.mesh2().points()) {
              int I = im / _Ns;
              int J = im % _Ns;
              bare_gf(w, alps::gf::index_mesh::index_type(im), alps::gf::index_mesh::index_type(is)) = G_inv(I, J);
            }
          }
        }
      }

    private:
      /// Symmetry
      Symmetry::SzSymmetry _symmetry;
      /// Hopping
      std::vector < std::vector < precision > > t;
      /// Interaction
      std::vector < precision > U;
      /// Exchange
      std::vector < std::vector < precision > > J;
      /// Chemical potential
      std::vector < precision > _xmu;
      /// Magnetic field
      std::vector < precision > _Hmag;
      /// site energy shift
      std::vector < std::vector < precision > > _Eps;

      /// Non-diagonal states iterators
      std::vector < St > _states;
      std::vector < St > _V_states;
    };

  }
}
#endif //HUBBARD_HUBBARDMODEL_H
