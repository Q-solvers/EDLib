//
// Created by iskakoff on 01/08/16.
//

#ifndef HUBBARD_GREENSFUNCTION_H
#define HUBBARD_GREENSFUNCTION_H


#include "Lanczos.h"

template<typename precision, class Hamiltonian, class Model>
class GreensFunction: public Lanczos<precision, Hamiltonian> {
  using Lanczos<precision, Hamiltonian>::hamiltonian;
  using Lanczos<precision, Hamiltonian>::lanczos;
public:
  GreensFunction(alps::params& p, Hamiltonian& h) : Lanczos<precision, Hamiltonian>(p, h), gf(Lanczos<precision, Hamiltonian>::omega()), _model(p) {
  }

  void compute() {
    for(auto & eigenpair : hamiltonian().eigenpairs()) {
      hamiltonian().symmetry().set_sector(eigenpair.sector());
      for(int i = 0; i<1/*_model.orbitals()*/; ++i) {
        for(int is = 0; is< _model.spins() ; ++is) {
          std::vector<precision> outvec(eigenpair.sector().size(), precision(0.0));
          if(hamiltonian().symmetry().template create_particle<precision, Model>(i, is, eigenpair.eigenvector(), outvec, _model)){
            lanczos(outvec);
          }
          if(hamiltonian().symmetry().template annihilate_particle<precision, Model>(i, is, eigenpair.eigenvector(), outvec, _model)){
            lanczos(outvec);
          }
        }
      }
    }
  }
private:
  alps::gf::omega_gf gf;
  Model _model;
};


#endif //HUBBARD_GREENSFUNCTION_H
