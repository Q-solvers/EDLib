//
// Created by iskakoff on 01/08/16.
//

#ifndef HUBBARD_LANCZOS_H
#define HUBBARD_LANCZOS_H

#include <alps/gf/mesh.hpp>
#include <alps/gf/gf.hpp>
#include <alps/params.hpp>

template<typename precision, class Hamiltonian>
class Lanczos {
public:
  Lanczos(alps::params& p, Hamiltonian &h): _omega(p["lanc.BETA"], p["lanc.NOMEGA"]), _Nl(p["lanc.NLANC"]), ham(h), alfalanc(_Nl), betalanc(_Nl+1) {

  }

  const alps::gf::matsubara_positive_mesh & omega() const {
    return _omega;
  }

  int lanczos(std::vector<precision> &v) {
    int nlanc = 0;
    unsigned long size = v.size();
    std::vector<precision> w(size, precision(0.0));
    precision alf = 0, bet = 0;
    ham.fill();
    for (int iter = 1; iter <= _Nl; iter++) {
      ++nlanc;
      if(iter != 1) {
        for (int j = 0; j < size; ++j) {
          precision dummy = v[j];
          v[j] = w[j]/bet;
          w[j] = -bet*dummy;
        }
      }
      ham.storage().av(v.data(), w.data(), size, false);
      for (int k = 0; k < v.size(); ++k) {
        alf += w[k] * v[k];
      }
      alfalanc[iter - 1] = alf;
      for (int j = 0; j < size; ++j) {
        w[j] -= alf*v[j];
      }
      for (int k = 0; k < v.size(); ++k) {
        bet += w[k] * w[k];
      }
      bet = std::sqrt(bet);

      if (iter != _Nl) betalanc[iter] = bet;
      if (std::abs(bet) < 1e-10 && iter >= (2 * ham.symmetry().sector().size())) {
        break;
      }
    }
    return nlanc;
  }

  const Hamiltonian &hamiltonian() const{
    return ham;
  };

  Hamiltonian &hamiltonian() {
    return ham;
  };

private:
  alps::gf::matsubara_positive_mesh _omega;

  int _Nl;
  Hamiltonian &ham;

  std::vector<precision> alfalanc;
  std::vector<precision> betalanc;
};


#endif //HUBBARD_LANCZOS_H
