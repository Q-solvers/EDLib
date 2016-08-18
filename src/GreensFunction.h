//
// Created by iskakoff on 01/08/16.
//

#ifndef HUBBARD_GREENSFUNCTION_H
#define HUBBARD_GREENSFUNCTION_H


#include <iomanip>
#include "Lanczos.h"

template<typename precision, class Hamiltonian, class Model>
class GreensFunction: public Lanczos<precision, Hamiltonian> {
  using Lanczos<precision, Hamiltonian>::hamiltonian;
  using Lanczos<precision, Hamiltonian>::lanczos;
  using Lanczos<precision, Hamiltonian>::omega;
  using Lanczos<precision, Hamiltonian>::computefrac;
public:
  GreensFunction(alps::params& p, Hamiltonian& h) : Lanczos<precision, Hamiltonian>(p, h), _model(p),
                                                    gf(p["NSPINS"], std::vector<alps::gf::omega_gf>(p["NSITES"], alps::gf::omega_gf(Lanczos<precision, Hamiltonian>::omega() ) ) ),
                                                    _cutoff(p["lanc.BOLTZMANN_CUTOFF"]){
  }

  void compute() {
    auto & ground_state = hamiltonian().eigenpairs()[0];
    for(auto & eigenpair : hamiltonian().eigenpairs()) {
      std::cout<<"Compute Green's function contribution for eigenvalue E="<<eigenpair.eigenvalue()<<" for sector"<<eigenpair.sector()<<std::endl;
      if(std::exp(-(eigenpair.eigenvalue() - ground_state.eigenvalue())*omega().beta()) < _cutoff) {
        std::cout<<"Skipped by Boltzmann factor."<<std::endl;
        continue;
      }
      for(int i = 0; i<1/*_model.orbitals()*/; ++i) {
        for(int is = 0; is< _model.spins() ; ++is) {
          std::vector<precision> outvec(eigenpair.sector().size(), precision(0.0));
          double expectation_value = 0;
          hamiltonian().symmetry().set_sector(eigenpair.sector());
          if(hamiltonian().symmetry().template create_particle<precision, Model>(i, is, eigenpair.eigenvector(), outvec, _model, expectation_value)){
            int nlanc = lanczos(outvec);
            std::cout<<"orbital: "<<i<<"   spin: "<<is<<" "<<expectation_value<<std::endl;
            computefrac(expectation_value, eigenpair.eigenvalue(), ground_state.eigenvalue(), nlanc, 1, gf[is][i]);
          }
          hamiltonian().symmetry().set_sector(eigenpair.sector());
          if(hamiltonian().symmetry().template annihilate_particle<precision, Model>(i, is, eigenpair.eigenvector(), outvec, _model, expectation_value)){
            int nlanc = lanczos(outvec);
            std::cout<<"orbital: "<<i<<"   spin: "<<is<<" "<<expectation_value<<std::endl;
            computefrac(expectation_value, eigenpair.eigenvalue(), ground_state.eigenvalue(), nlanc, -1, gf[is][i]);
          }
        }
      }
    }
    for(int i = 0; i<1/*_model.orbitals()*/; ++i) {
      for (int is = 0; is < _model.spins(); ++is) {
        std::ostringstream Gomega_name;
        Gomega_name<<"G_omega"<<"_"<<(is ? "up": "down")<<"_"<<i;
        std::ofstream G_omega_file(Gomega_name.str().c_str());
        G_omega_file << std::setprecision(14) << gf[is][i];
        Gomega_name<<".h5";
        alps::hdf5::archive ar(Gomega_name.str().c_str(), alps::hdf5::archive::WRITE);
        gf[is][i].save(ar, "/G_omega");
        G_omega_file.close();
        ar.close();
      }
    }
  }
private:
  std::vector< std::vector<alps::gf::omega_gf> > gf;
  Model _model;
  double _cutoff;
};


#endif //HUBBARD_GREENSFUNCTION_H
