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

  int lanczos() {
    int nlanc = 0;
    MPI_Win v_win;
    precision alf = 0, bet = 0;
    for (int iter = 1; iter <= _Nl; iter++) {
      nlanc = nlanc + 1;
      alfalanc[iter - 1] = alf;
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

private:
  alps::gf::matsubara_positive_mesh _omega;

  int _Nl;
  Hamiltonian &ham;

  std::vector<precision> alfalanc;
  std::vector<precision> betalanc;
};


#endif //HUBBARD_LANCZOS_H
