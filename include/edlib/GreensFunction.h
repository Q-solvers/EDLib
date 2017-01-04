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
      GreensFunction(alps::params &p, Hamiltonian &h) : Lanczos < precision, Hamiltonian >(p, h), _model(h.model()),
                                                        gf(Lanczos < precision, Hamiltonian >::omega(), alps::gf::index_mesh(p["NSITES"].as<int>()), alps::gf::index_mesh(p["NSPINS"].as<int>())),
                                                        _cutoff(p["lanc.BOLTZMANN_CUTOFF"]) {
        if(p["storage.EIGENVALUES_ONLY"] == 1) {
          throw std::logic_error("Eigenvectors have not been computed. Green's function can not be evaluated.");
        }
      }

      void compute() {
        if(hamiltonian().eigenpairs().empty())
          return;
        const EigenPair<precision, typename Hamiltonian::ModelType::Sector> &groundstate =  *hamiltonian().eigenpairs().begin();
        for (typename std::set<EigenPair<precision, typename Hamiltonian::ModelType::Sector> >::iterator kkk = hamiltonian().eigenpairs().begin(); kkk != hamiltonian().eigenpairs().end(); kkk++) {
          const EigenPair<precision, typename Hamiltonian::ModelType::Sector> &eigenpair = *kkk;
          _Z += std::exp(-(eigenpair.eigenvalue() - groundstate.eigenvalue()) * omega().beta());
        }
        for (typename std::set<EigenPair<precision, typename Hamiltonian::ModelType::Sector> >::iterator kkk = hamiltonian().eigenpairs().begin(); kkk != hamiltonian().eigenpairs().end(); kkk++) {
          const EigenPair<precision, typename Hamiltonian::ModelType::Sector>& pair = *kkk;
          precision boltzmann_f = std::exp(-(pair.eigenvalue() - groundstate.eigenvalue()) * omega().beta());
          if (boltzmann_f < _cutoff) {
//        std::cout<<"Skipped by Boltzmann factor."<<std::endl;
            continue;
          }
          std::cout << "Compute Green's function contribution for eigenvalue E=" << pair.eigenvalue() << " with Boltzmann factor = "
                    << boltzmann_f << "; for sector" << pair.sector() << std::endl;
          for (int i = 0; i < _model.orbitals(); ++i) {
            for (int is = 0; is < _model.spins(); ++is) {
              std::vector < precision > outvec(pair.sector().size(), precision(0.0));
              precision expectation_value = 0;
              _model.symmetry().set_sector(pair.sector());
              if (create_particle(i, is, pair.eigenvector(), outvec, expectation_value)) {
                hamiltonian().storage().prepare_work_arrays(outvec.data());
                int nlanc = lanczos(outvec);
                hamiltonian().storage().finalize();
                std::cout << "orbital: " << i << "   spin: " << (is == 0 ? "up" : "down") << " <n|aa*|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
                compute_continues_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, 1, gf, alps::gf::index_mesh::index_type(i), alps::gf::index_mesh::index_type(is));
              }
              _model.symmetry().set_sector(pair.sector());
              if (annihilate_particle(i, is, pair.eigenvector(), outvec, expectation_value)) {
                hamiltonian().storage().prepare_work_arrays(outvec.data());
                int nlanc = lanczos(outvec);
                hamiltonian().storage().finalize();
                std::cout << "orbital: " << i << "   spin: " << (is == 0 ? "up" : "down") << " <n|a*a|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
                compute_continues_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, -1, gf, alps::gf::index_mesh::index_type(i), alps::gf::index_mesh::index_type(is));
              }
            }
          }
        }
#ifdef USE_MPI
        int rank;
        MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
        if(rank == 0) {
#endif
        gf /= _Z;
        std::ostringstream Gomega_name;
        Gomega_name << "G_omega";
        std::ofstream G_omega_file(Gomega_name.str().c_str());
        G_omega_file << std::setprecision(14) << gf;
        Gomega_name << ".h5";
        alps::hdf5::archive ar(Gomega_name.str().c_str(), alps::hdf5::archive::WRITE);
        gf.save(ar, "/G_omega");
        G_omega_file.close();
        ar.close();
        std::cout << "Statsum: " << _Z << std::endl;
#ifdef USE_MPI
        }
#endif
      }

    private:
      typedef alps::gf::three_index_gf<std::complex<double>, alps::gf::matsubara_mesh<alps::gf::mesh::POSITIVE_ONLY>, alps::gf::index_mesh, alps::gf::index_mesh >  GF_TYPE;
      GF_TYPE gf;
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
      bool create_particle(int orbital, int spin, const std::vector < precision > &invec, std::vector < precision > &outvec, double &expectation_value) {
        // check that the particle can be annihilated
        if ((spin == 0 && _model.symmetry().sector().nup() == _model.orbitals()) || (spin == 1 && _model.symmetry().sector().ndown() == _model.orbitals())) {
          return false;
        }
        hamiltonian().storage().reset();
        long long k = 0;
        int sign = 0;
        int nup_new = _model.symmetry().sector().nup() + (1 - spin);
        int ndn_new = _model.symmetry().sector().ndown() + spin;
        typename Hamiltonian::ModelType::Sector next_sec(nup_new, ndn_new, _model.symmetry().comb().c_n_k(_model.orbitals(), nup_new) * _model.symmetry().comb().c_n_k(_model.orbitals(), ndn_new));
        outvec.assign(hamiltonian().storage().vector_size(next_sec), 0.0);
        int i = 0;
        hamiltonian().storage().a_adag(orbital + spin * _model.orbitals(), invec, outvec, next_sec, false);
        double norm = hamiltonian().storage().vv(outvec, outvec);
        for (int j = 0; j < outvec.size(); ++j) {
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
      bool annihilate_particle(int orbital, int spin, const std::vector < precision > &invec, std::vector < precision > &outvec, double &expectation_value) {
        // check that the particle can be annihilated
        if ((spin == 0 && _model.symmetry().sector().nup() == 0) || (spin == 1 && _model.symmetry().sector().ndown() == 0)) {
          return false;
        }
        hamiltonian().storage().reset();
        long long k = 0;
        int sign = 0;
        int nup_new = _model.symmetry().sector().nup() - (1 - spin);
        int ndn_new = _model.symmetry().sector().ndown() - spin;
        typename Hamiltonian::ModelType::Sector next_sec(nup_new, ndn_new, _model.symmetry().comb().c_n_k(_model.orbitals(), nup_new) * _model.symmetry().comb().c_n_k(_model.orbitals(), ndn_new));
        outvec.assign(hamiltonian().storage().vector_size(next_sec), precision(0.0));
        int i = 0;
        hamiltonian().storage().a_adag(orbital + spin * _model.orbitals(), invec, outvec, next_sec, true);
        double norm = hamiltonian().storage().vv(outvec, outvec);
        for (int j = 0; j < outvec.size(); ++j) {
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
