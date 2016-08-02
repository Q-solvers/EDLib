//
// Created by iskakoff on 01/08/16.
//

#ifndef HUBBARD_GREENSFUNCTION_H
#define HUBBARD_GREENSFUNCTION_H


#include "Lanczos.h"

template<typename precision, class Hamiltonian, class Model>
class GreensFunction: public Lanczos<precision, Hamiltonian> {
  using Lanczos<precision, Hamiltonian>::hamiltonian;
public:
  GreensFunction(alps::params& p, Hamiltonian& h) : Lanczos<precision, Hamiltonian>(p, h), gf(Lanczos<precision, Hamiltonian>::omega()), _model(p) {
  }

  void compute() {
    for(auto & eigenpair : hamiltonian().eigenpairs()) {
      hamiltonian().symmetry().set_sector(eigenpair.sector());
      for(int i = 0; i<_model.orbitals(); ++i) {
        for(int is = 0; is< _model.spins())
      }
    }
  }
private:
  alps::gf::omega_gf gf;
  Model _model;
};


#endif //HUBBARD_GREENSFUNCTION_H
