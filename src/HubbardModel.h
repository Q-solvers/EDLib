//
// Created by iskakoff on 29/07/16.
//

#ifndef HUBBARD_HUBBARDMODEL_H
#define HUBBARD_HUBBARDMODEL_H

#include <vector>
#include <alps/params.hpp>


namespace Hubbard {
  template<typename prec>
  class InnerState {
  public:
    InnerState(int ii, int jj, int spin, prec val) : _indicies(ii, jj), _spin(spin), _value(val) {};
    const inline std::pair<int, int>& indicies() const {return _indicies;}
    const inline prec& value() const {return _value;}
    const inline int spin() const {return _spin;}
  private:
    std::pair<int, int> _indicies;
    int _spin;
    prec _value;
  };
}

template<typename precision>
class HubbardModel {
public:
  typedef typename Hubbard::InnerState<precision> St;
  typedef typename SzSymmetry::Sector Sector;

  HubbardModel(alps::params &p) : _symmetry(p), _Ns(p["NSITES"]), _ms(p["NSPINS"]), _Ip(_ms*_Ns),
                                  Eps(p["NSITES"], std::vector<precision>(p["NSPINS"], precision(0.0))),
                                  t(p["NSITES"], std::vector<precision>(p["NSITES"], precision(0.0))),
                                  U(p["NSITES"], precision(0.0)) {
    std::string input = p["INPUT_FILE"];
    alps::hdf5::archive input_data(input.c_str(), "r");
    input_data>>alps::make_pvp("BETA", _beta);
    input_data>>alps::make_pvp("hopping/values", t);
    input_data>>alps::make_pvp("interaction/values", U);
    input_data.close();
    for(int i = 0; i< _Ns; ++i ) {
      // HARDCODED Half-filling paramagnetic shift
      Eps[i][0] = Eps[i][1] = -U[i]/2.0;
    }
    for(int ii = 0; ii< _Ns; ++ii) {
      for (int jj = 0; jj < _Ns; ++jj) {
        if (std::abs(t[ii][jj]) > 1e-10) {
          for(int is = 0; is < _ms; ++is) {
            _states.push_back(St(ii, jj, is, t[ii][jj]));
          }
        }
      }
    }
    // TODO: move to input file and make site-dependent
    _xmu = precision(0.0);
  };

  inline int valid(const St & state, const long long &nst) {
    return (checkState(nst, state.indicies().first + state.spin() * _Ns)*  (1 - checkState(nst, state.indicies().second + state.spin() * _Ns)));
  }

  inline void set(const St & state, const long long &nst, long long &k, int &sign) {
    long long k1, k2;
    int isign1, isign2;
    a(state.indicies().first + state.spin() * _Ns, nst, k1, isign1);
    adag(state.indicies().second + state.spin() * _Ns, k1, k2, isign2);
    k = k2;
    // -t c^+ c
    sign = -isign1*isign2;
  }

  inline const precision diagonal(long long state) const {
    precision xtemp = 0.0;
    for(int im = 0;im < _Ns;++im){
      for(int is = 0 ; is< _ms; ++is) {
        xtemp += (Eps[im][is] - _xmu) * checkState(state, im + is* _Ns);
      }
      xtemp += U[im]* checkState(state, im)* checkState(state, im + _Ns);
    }
    return xtemp;
  }
  /**
   * Check that im state is occupated
   *
   * \param nst - current state
   * \param im - state to check
   *
   * \return 0 if state is empty, 1 - otherwise
   */
  int inline checkState(const long long& nst, const int& im) const {
    return (int) ((nst & (1ll << (_Ip - 1 - im))) >> (_Ip - 1 - im));
  }
  /**
   * Anihilate particle
   * \param i [in] - site to anihilate particle
   * \param jold [in] - current state
   * \param k [out] - resulting state
   * \param isign [out] - fermionic sign
   */
  void inline a(const int& i,const long long& jold,long long &k,int &isign) {
    long long sign=0;
    for(int ll=0; ll<i; ++ll) {
      sign+= ((jold&(1ll<<(_Ip-ll-1)))!=0) ? 1 : 0;
    }
    isign = (sign % 2) == 0 ? 1 : -1;
    k=jold-(1ll<<(_Ip-i-1));
  }

  /**
   * Create particle
   * \param i [in] - site to create particle
   * \param jold [in] - current state
   * \param k [out] - resulting state
   * \param isign [out] - fermionic sign
   */
  void inline adag(const int& i, const long long &jold, long long& k, int& isign) {
    long long sign=0;
    for(int ll=0; ll<i; ++ll) {
      sign+= ((jold&(1ll<<(_Ip-ll-1)))!=0) ? 1 : 0;
    }
    isign = (sign % 2) == 0 ? 1 : -1;
    k = jold + (1ll << (_Ip - i -1));
  }

  const std::vector<St>& states() const {return _states;};

  const int orbitals() const {
    return _Ns;
  }

  const int spins() const {
    return _ms;
  }

  inline const SzSymmetry &symmetry() const {
    return _symmetry;
  }

  inline SzSymmetry &symmetry() {
    return _symmetry;
  }

  /**
 * Actions on v>
 */
  bool create_particle(int orbital, int spin, const std::shared_ptr<precision>& invec, std::vector<precision>& outvec, double & expectation_value) {
    // check that the particle can be annihilated
    if((spin == 0 && _symmetry.sector().nup()==_Ns) || (spin == 1 && _symmetry.sector().ndown()==_Ns)) {
      return false;
    }
    _symmetry.init();
    long long k = 0;
    int sign = 0;
    int nup_new = _symmetry.sector().nup() + (1-spin);
    int ndn_new = _symmetry.sector().ndown() + spin;
    Sector next_sec(nup_new, ndn_new, _symmetry.c_n_k(_Ns, nup_new) * _symmetry.c_n_k(_Ns,ndn_new));
    outvec.assign(next_sec.size(), 0.0);
    double norm = 0.0;
    int i = 0;
    while(_symmetry.next_state()){
      long long nst = _symmetry.state();
      if(checkState(nst, orbital + spin*_Ns) == 0) {
        adag(orbital + spin*_Ns, nst, k, sign);
        int i1 = _symmetry.index(k, next_sec);
        outvec[i1] = sign * invec.get()[i];
        norm += std::norm(outvec[i1]);
      }
      ++i;
    };
//    norm = std::sqrt(norm);
    for (int j = 0; j < next_sec.size(); ++j) {
      outvec[j] /= std::sqrt(norm);
    }
    _symmetry.set_sector(next_sec);
    expectation_value = norm;
    return true;
  };
  bool annihilate_particle(int orbital, int spin, const std::shared_ptr<precision>& invec, std::vector<precision>& outvec, double & expectation_value) {
    // check that the particle can be annihilated
    if((spin == 0 && _symmetry.sector().nup()==0) || (spin == 1 && _symmetry.sector().ndown()==0)) {
      return false;
    }
    _symmetry.init();
    long long k = 0;
    int sign = 0;
    int nup_new = _symmetry.sector().nup() - (1-spin);
    int ndn_new = _symmetry.sector().ndown() - spin;
    Sector next_sec(nup_new, ndn_new, _symmetry.c_n_k(_Ns, nup_new) * _symmetry.c_n_k(_Ns, ndn_new));
    outvec.assign(next_sec.size(), precision(0.0));
    double norm = 0.0;
    int i = 0;
    while(_symmetry.next_state()){
      long long nst = _symmetry.state();
      if(checkState(nst, orbital + spin*_Ns)) {
        a(orbital + spin*_Ns, nst, k, sign);
        int i1 = _symmetry.index(k, next_sec);
        outvec[i1] = sign * invec.get()[i];
        // v_i * v_i^{\star}
        norm += std::norm(outvec[i1]);
      }
      ++i;
    };
    for (int j = 0; j < next_sec.size(); ++j) {
      outvec[j] /= std::sqrt(norm);
    }
    _symmetry.set_sector(next_sec);
    // <v|a^{\star}a|v>
    expectation_value = norm;
    return true;
  };

private:
  // Symmetry
  SzSymmetry _symmetry;
  // Hopping
  std::vector<std::vector<precision> > t;
  // Interaction
  std::vector<precision> U;
  // Chemical potential
  precision _xmu;
  // site energy shift
  std::vector<std::vector<precision> > Eps;
  // Inverse temperature
  double _beta;

  // Non-diagonal states iterator
  std::vector<St> _states;

  /**
   * _Ns - number of lattice sites
   * _ms - number of electron spins
   * _Ip - maximum number of electrons
   */
  int _Ns;
  int _ms;
  int _Ip;
};


#endif //HUBBARD_HUBBARDMODEL_H
