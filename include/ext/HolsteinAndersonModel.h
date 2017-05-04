//
// Created by iskakoff on 20/04/17.
//

#ifndef HUBBARD_HOLSTEINANDERSONMODEL_H
#define HUBBARD_HOLSTEINANDERSONMODEL_H


#include <edlib/FermionicModel.h>
#include <edlib/CommonUtils.h>
#include <alps/gf/gf.hpp>
#include <edlib/SingleImpurityAndersonModel.h>
#include "NSymmetryWithBoson.h"

namespace EDLib {
  namespace Model {
    namespace HolsteinAnderson {
      template<typename prec>
      class InnerState : public SingleImpurityAnderson::InnerState<prec> {
      public:
        virtual int valid(long long, int, int) const {return 0;};
        virtual void set(long long,long long&, int&, int, int) const {};
      };

      template<typename prec>
      class HybridisationInnerState : public InnerState<prec> {

      };
      template<typename prec>
      class BosonInnerState : public InnerState<prec> {

      };
    }
    /**
     * @brief HolsteinAndersonModel class
     *
     * @author iskakoff
     */
    template<typename prec>
    class HolsteinAndersonModel : public FermionicModel {
    public:
      typedef prec precision;
      typedef typename Symmetry::NSymmetryWithBoson SYMMETRY;
      typedef typename HolsteinAnderson::InnerState < precision > St;
      typedef typename HolsteinAnderson::HybridisationInnerState < precision > HSt;
      typedef typename HolsteinAnderson::BosonInnerState < precision > BSt;
      typedef typename Symmetry::NSymmetryWithBoson::Sector Sector;


      HolsteinAndersonModel(alps::params &p) : FermionicModel(p), _symmetry(p) {}

      inline const precision diagonal(long long full_state) const {
        precision xtemp = 0.0;
        long long state = full_state>>_Nb;
        long long bosons = full_state & ((1<<_Nb) - 1);

        for (int is = 0; is < _ms; ++is) {
          for (int ik = 0; ik < _Epsk.size(); ++ik) {
            int ikm = ik + 1;
            xtemp+=(_Epsk[ik][is] * checkState(state, ikm + is * _Ns, _Ip));
          }
          xtemp += (_Eps[is] - _xmu) * checkState(state, is * _Ns, _Ip);
        }
        xtemp += _U * checkState(state, 0, _Ip) * checkState(state, _Ns, _Ip);
        return xtemp;
      }

      inline int valid(const St &state, long long nst) {
        return state.valid(nst, _Ns, _Nb);
      }

      inline void set(const St &state, long long nst, long long &k, int &sign) {
        state.set(nst, k, sign, _Ns, _Nb);
      }

      SYMMETRY &symmetry() {
        return _symmetry;
      }

      inline const std::vector<HSt>& T_states() const {
        return _F_states;
      }

      inline const std::vector<BSt>& V_states() const {
        return _B_states;
      }

      int interacting_orbitals() const {
        return 1;
      }

      template<typename Mesh>
      void bare_greens_function(alps::gf::three_index_gf<std::complex<double>, Mesh, alps::gf::index_mesh, alps::gf::index_mesh >& bare_gf, double beta) {
      }

      inline const Symmetry::NSymmetryWithBoson& symmetry() const {
        return _symmetry;
      }

      inline Symmetry::NSymmetryWithBoson& symmetry() {
        return _symmetry;
      }
    private:
      Symmetry::NSymmetryWithBoson _symmetry;
      int _Nk;
      int _bd;
      int _Nb;
      precision _xmu;
      std::vector < prec > _Eps;
      std::vector < std::vector < precision > > _Vk;
      std::vector < std::vector < precision > > _Epsk;
      precision _U;
      // bosonic part
      std::vector < prec > _w0;
      std::vector < prec > _W;


      std::vector<HSt> _F_states;
      std::vector<BSt> _B_states;
    };

  }
}

#endif //HUBBARD_HOLSTEINANDERSONMODEL_H
