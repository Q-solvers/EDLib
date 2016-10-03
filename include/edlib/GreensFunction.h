//
// Created by iskakoff on 01/08/16.
//

#ifndef HUBBARD_GREENSFUNCTION_H
#define HUBBARD_GREENSFUNCTION_H


#include <iomanip>
#include "Lanczos.h"

namespace EDLib {
  namespace gf {
    template<typename precision, class Hamiltonian>
    class GreensFunction : public Lanczos < precision, Hamiltonian > {
      using Lanczos < precision, Hamiltonian >::hamiltonian;
      using Lanczos < precision, Hamiltonian >::lanczos;
      using Lanczos < precision, Hamiltonian >::omega;
      using Lanczos < precision, Hamiltonian >::compute_continues_fraction;
    public:
      GreensFunction(EDParams &p, Hamiltonian &h) : Lanczos < precision, Hamiltonian >(p, h), _model(h.model()),
                                                    gf(p["NSPINS"],
                                                       std::vector < alps::gf::omega_gf >(p["NSITES"], alps::gf::omega_gf(Lanczos < precision, Hamiltonian >::omega()))),
                                                    _cutoff(p["lanc.BOLTZMANN_CUTOFF"]) {
        if(p["storage.EIGENVALUES_ONLY"] == 1) {
          throw std::logic_error("Eigenvectors have not been computed. Green's function can not be evaluated.");
        }
      }

      void compute() {
        for (int kkk = 0; kkk < hamiltonian().eigenpairs().size(); ++kkk) {
          _Z += std::exp(-(hamiltonian().eigenpairs()[kkk].eigenvalue() - hamiltonian().eigenpairs()[0].eigenvalue()) * omega().beta());
        }
        for (int kkk = 0; kkk < hamiltonian().eigenpairs().size(); ++kkk) {
          precision boltzmann_f = std::exp(-(hamiltonian().eigenpairs()[kkk].eigenvalue() - hamiltonian().eigenpairs()[0].eigenvalue()) * omega().beta());
          if (boltzmann_f < _cutoff) {
//        std::cout<<"Skipped by Boltzmann factor."<<std::endl;
            continue;
          }
          std::cout << "Compute Green's function contribution for eigenvalue E=" << hamiltonian().eigenpairs()[kkk].eigenvalue() << " with Boltzmann factor = "
                    << boltzmann_f << "; for sector" << hamiltonian().eigenpairs()[kkk].sector() << std::endl;
          for (int i = 0; i < 1/*_model.orbitals()*/; ++i) {
            for (int is = 0; is < _model.spins(); ++is) {
              std::vector < precision > outvec(hamiltonian().eigenpairs()[kkk].sector().size(), precision(0.0));
              precision expectation_value = 0;
              _model.symmetry().set_sector(hamiltonian().eigenpairs()[kkk].sector());
              if (create_particle(i, is, hamiltonian().eigenpairs()[kkk].eigenvector(), outvec, expectation_value)) {
                int nlanc = lanczos(outvec);
                std::cout << "orbital: " << i << "   spin: " << (is == 0 ? "up" : "down") << " <n|aa*|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
                compute_continues_fraction(expectation_value, hamiltonian().eigenpairs()[kkk].eigenvalue(), hamiltonian().eigenpairs()[0].eigenvalue(), nlanc, 1, gf[is][i]);
              }
              _model.symmetry().set_sector(hamiltonian().eigenpairs()[kkk].sector());
              if (annihilate_particle(i, is, hamiltonian().eigenpairs()[kkk].eigenvector(), outvec, expectation_value)) {
                int nlanc = lanczos(outvec);
                std::cout << "orbital: " << i << "   spin: " << (is == 0 ? "up" : "down") << " <n|a*a|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
                compute_continues_fraction(expectation_value, hamiltonian().eigenpairs()[kkk].eigenvalue(), hamiltonian().eigenpairs()[0].eigenvalue(), nlanc, -1, gf[is][i]);
              }
            }
          }
        }
        for (int i = 0; i < 1/*_model.orbitals()*/; ++i) {
          for (int is = 0; is < _model.spins(); ++is) {
            gf[is][i] /= _Z;
            std::ostringstream Gomega_name;
            Gomega_name << "G_omega" << "_" << (is ? "down" : "up") << "_" << i;
            std::ofstream G_omega_file(Gomega_name.str().c_str());
            G_omega_file << std::setprecision(14) << gf[is][i];
            Gomega_name << ".h5";
            alps::hdf5::archive ar(Gomega_name.str().c_str(), alps::hdf5::archive::WRITE);
            gf[is][i].save(ar, "/G_omega");
            G_omega_file.close();
            ar.close();
          }
        }
        std::cout << "Statsum: " << _Z << std::endl;
      }

    private:
      std::vector < std::vector < alps::gf::omega_gf > > gf;
      typename Hamiltonian::ModelType &_model;
      precision _cutoff;
      precision _Z;

      /**
       * @brief Perform the create operator action to the eigenstate
       *
       * @param orbital - the orbital to create a particle
       * @param spin - the spin of a particle to create
       * @param invec - current eigenstate
       * @param outvec - Op-vec product
       * @param expectation_value - expectation value of aa*
       * @return true if the particle has been created
       */
      bool create_particle(int orbital, int spin, const std::shared_ptr < precision > &invec, std::vector < precision > &outvec, double &expectation_value) {
        // check that the particle can be annihilated
        if ((spin == 0 && _model.symmetry().sector().nup() == _model.orbitals()) || (spin == 1 && _model.symmetry().sector().ndown() == _model.orbitals())) {
          return false;
        }
        _model.symmetry().init();
        long long k = 0;
        int sign = 0;
        int nup_new = _model.symmetry().sector().nup() + (1 - spin);
        int ndn_new = _model.symmetry().sector().ndown() + spin;
        typename Hamiltonian::ModelType::Sector next_sec(nup_new, ndn_new, _model.symmetry().comb().c_n_k(_model.orbitals(), nup_new) * _model.symmetry().comb().c_n_k(_model.orbitals(), ndn_new));
        outvec.assign(next_sec.size(), 0.0);
        double norm = 0.0;
        int i = 0;
        while (_model.symmetry().next_state()) {
          long long nst = _model.symmetry().state();
          if (_model.checkState(nst, orbital + spin * _model.orbitals(), _model.max_total_electrons()) == 0) {
            _model.adag(orbital + spin * _model.orbitals(), nst, k, sign);
            int i1 = _model.symmetry().index(k, next_sec);
            outvec[i1] = sign * invec.get()[i];
            norm += std::norm(outvec[i1]);
          }
          ++i;
        };
        for (int j = 0; j < next_sec.size(); ++j) {
          outvec[j] /= std::sqrt(norm);
        }
        _model.symmetry().set_sector(next_sec);
        expectation_value = norm;
        return true;
      };

      /**
       * @brief Perform the annihilator operator action to the eigenstate
       *
       * @param orbital - the orbital to destroy a particle
       * @param spin - the spin of a particle to destroy
       * @param invec - current eigenstate
       * @param outvec - Op-vec product
       * @param expectation_value - expectation value of a*a
       * @return true if the particle has been destroyed
       */
      bool annihilate_particle(int orbital, int spin, const std::shared_ptr < precision > &invec, std::vector < precision > &outvec, double &expectation_value) {
        // check that the particle can be annihilated
        if ((spin == 0 && _model.symmetry().sector().nup() == 0) || (spin == 1 && _model.symmetry().sector().ndown() == 0)) {
          return false;
        }
        _model.symmetry().init();
        long long k = 0;
        int sign = 0;
        int nup_new = _model.symmetry().sector().nup() - (1 - spin);
        int ndn_new = _model.symmetry().sector().ndown() - spin;
        typename Hamiltonian::ModelType::Sector next_sec(nup_new, ndn_new, _model.symmetry().comb().c_n_k(_model.orbitals(), nup_new) * _model.symmetry().comb().c_n_k(_model.orbitals(), ndn_new));
        outvec.assign(next_sec.size(), precision(0.0));
        double norm = 0.0;
        int i = 0;
        while (_model.symmetry().next_state()) {
          long long nst = _model.symmetry().state();
          if (_model.checkState(nst, orbital + spin * _model.orbitals(), _model.max_total_electrons())) {
            _model.a(orbital + spin * _model.orbitals(), nst, k, sign);
            int i1 = _model.symmetry().index(k, next_sec);
            outvec[i1] = sign * invec.get()[i];
            // v_i * v_i^{\star}
            norm += std::norm(outvec[i1]);
          }
          ++i;
        };
        for (int j = 0; j < next_sec.size(); ++j) {
          outvec[j] /= std::sqrt(norm);
        }
        _model.symmetry().set_sector(next_sec);
        // <v|a^{\star}a|v>
        expectation_value = norm;
        return true;
      };
    };
  }
}

#endif //HUBBARD_GREENSFUNCTION_H
