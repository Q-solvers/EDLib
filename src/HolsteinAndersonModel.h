//
// Created by iskakoff on 22/08/16.
//

#ifndef HUBBARD_HOLSTEINANDERSONMODEL_H
#define HUBBARD_HOLSTEINANDERSONMODEL_H

#include <utility>
#include "NSymmetryWithBoson.h"

namespace HolsteinAnderson {
  template<typename prec>
  class InnerState {
  public:
    InnerState(int ii, int spin, prec val) : _index(ii), _spin(spin), _value(val) {};
    const inline int& index() const {return _index;}
    const inline prec& value() const {return _value;}
    const inline int spin() const {return _spin;}
  private:
    int _index;
    int _spin;
    prec _value;
  };
}

template<typename precision>
class HolsteinAndersonModel {
public:

  typedef typename HolsteinAnderson::InnerState<precision> St;
  typedef typename NSymmetryWithBoson::Sector Sector;

  HolsteinAndersonModel(EDParams& p) : _symmetry(p), _Nf(2*int(p["NSITES"])), _Ns(p["NSITES"]), _Nb(p["NSITES_BOSE"]), _ms(p["NSPINS"]) {
    std::string input = p["INPUT_FILE"];
    alps::hdf5::archive input_data(input.c_str(), "r");
    input_data>>alps::make_pvp("BETA", _beta);
    input_data>>alps::make_pvp("U", _U);
    input_data>>alps::make_pvp("Eps", _Eps);
    input_data>>alps::make_pvp("fermion_bath/hybridization/values", _Vk);
    input_data>>alps::make_pvp("fermion_bath/level/values", _Epsk);
    input_data>>alps::make_pvp("bosonic_bath/energy", _W);
    input_data>>alps::make_pvp("bosonic_bath/coupling", _g);
    input_data.close();
    for(int ii = 0; ii< _Nf-1; ++ii) {
      if (std::abs(_Vk[ii]) > 1e-10) {
        for(int is = 0; is < _ms; ++is) {
          _states.push_back(St(ii, is, _Vk[ii]));
        }
      }
    }
  }


  inline const precision diagonal(long long state) const {
    precision xtemp = 0.0;
    for(int is = 0 ; is< _ms; ++is) {
      xtemp += _Eps * checkState(state, is* _Ns);
    }
    xtemp += _U * checkState(state, 0)* checkState(state, _Ns);
    xtemp += _W*number_of_bosons(state);
    return xtemp;
  }

private:
  // Symmetry
  NSymmetryWithBoson _symmetry;
  // Onsite Coulomb interaction
  precision _U;
  // Impurity level
  precision _Eps;
  // Fermionic bath energy levels
  std::vector<precision> _Epsk;
  // Energy of bosonic level
  precision _W;
  // Hybridization with fermionic bath
  std::vector<precision> _Vk;
  // Fermion-bosonic coupling
  precision _g;
  // Inverse temperature
  precision _beta;

  /**
   * _Nf - number of fermionic sites for both spins
   * _Ns - number of fermionic sites for single spin
   * _Nb - number of bosonic sites
   * _ms - number of spins
   */
  int _Nf;
  int _Ns;
  int _Nb;
  int _ms;

  // Non-diagonal states iterator
  std::vector<St> _states;

  /**
   * Check that im state is occupated
   *
   * \param nst - current state
   * \param im - state to check
   *
   * \return 0 if state is empty, 1 - otherwise
   */
  int inline checkState(const long long& nst, const int& im) const {
    return (int) ( ( (nst>>_Nb) & (1ll << (_Nf - 1 - im) ) ) >> (_Nf - 1 - im));
  }

  /**
   * \param nst - current state
   * \return number of bosons for current state
   */
  int inline number_of_bosons(const long long& nst) const{
    return nst & ((1ll<<_Nb) - 1);
  }
};


#endif //HUBBARD_HOLSTEINANDERSONMODEL_H
