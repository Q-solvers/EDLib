//
// Created by iskakoff on 03/01/17.
//

#ifndef HUBBARD_CHILOC_H
#define HUBBARD_CHILOC_H



#include <iomanip>
#include "Lanczos.h"
#include "EigenPair.h"


namespace EDLib {
  namespace gf {
    class BosonicOperator{};

    template<typename precision, class ModelType>
    class SzOperator: public BosonicOperator {
    public:
      SzOperator(ModelType &m) : _model(m){};
      precision action(long long state, int iii) const {
        return 0.5*(_model.checkState(state, iii, _model.max_total_electrons()) - _model.checkState(state, iii - _model.orbitals(), _model.max_total_electrons()));
      }
    private:
      const ModelType &_model;
    };

    template<typename precision, class Hamiltonian>
    class ChiLoc : public Lanczos < precision, Hamiltonian > {
      using Lanczos < precision, Hamiltonian >::hamiltonian;
      using Lanczos < precision, Hamiltonian >::lanczos;
      using Lanczos < precision, Hamiltonian >::omega;
      using Lanczos < precision, Hamiltonian >::compute_sym_continues_fraction;
    public:
      ChiLoc(alps::params &p, Hamiltonian &h) : Lanczos < precision, Hamiltonian >(p, h, alps::gf::statistics::statistics_type::BOSONIC), _model(h.model()),
                                                        gf(Lanczos < precision, Hamiltonian >::omega(), alps::gf::index_mesh(p["NSITES"].as<int>())),
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
        SzOperator<precision, typename Hamiltonian::ModelType> op(_model);
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
              std::vector < precision > outvec(pair.sector().size(), precision(0.0));
              precision expectation_value = 0;
              _model.symmetry().set_sector(pair.sector());
              if (operation(i, pair.eigenvector(), outvec, expectation_value, op)) {
                hamiltonian().storage().prepare_work_arrays(outvec.data());
                int nlanc = lanczos(outvec);
                hamiltonian().storage().finalize();
                std::cout << "orbital: " << i << " <n|O|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
                compute_sym_continues_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, 1, gf, alps::gf::index_mesh::index_type(i));
              }
          }
        }
#ifdef USE_MPI
        int rank;
        MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
        if(rank == 0) {
#endif
        gf/= _Z;
        /// Compute static susceptibility
        for (int i = 0; i < _model.orbitals(); ++i) {
          double chiSum = 0.0;
          for (int iomega = 1; iomega < omega().extent(); ++iomega) {
            chiSum = chiSum + gf(alps::gf::matsubara_positive_mesh::index_type(iomega), alps::gf::index_mesh::index_type(i)).real();
          }
          gf(alps::gf::matsubara_positive_mesh::index_type(0), alps::gf::index_mesh::index_type(i)) -= 2 * chiSum;
        }
        std::ostringstream Gomega_name;
        Gomega_name << "Chi_omega";
        std::ofstream G_omega_file(Gomega_name.str().c_str());
        G_omega_file << std::setprecision(14) << gf;
        Gomega_name << ".h5";
        alps::hdf5::archive ar(Gomega_name.str().c_str(), alps::hdf5::archive::WRITE);
        gf.save(ar, "/Chi_omega");
        G_omega_file.close();
        ar.close();
        std::cout << "Statsum: " << _Z << std::endl;
#ifdef USE_MPI
        }
#endif
      }

    private:
      typedef alps::gf::two_index_gf<std::complex<double>, alps::gf::matsubara_mesh<alps::gf::mesh::POSITIVE_ONLY>, alps::gf::index_mesh>  GF_TYPE;
      GF_TYPE gf;
      typename Hamiltonian::ModelType &_model;
      precision _cutoff;
      precision _Z;

      /**
       * @brief Perform the operation to the eigenstate
       *
       * @param orbital - the orbital to action
       * @param invec - current eigenstate
       * @param outvec - Op-vec product
       * @param expectation_value - expectation value of aa*
       * @return true if the particle has been created
       */
      template<typename Op>
      bool operation(int orbital, const std::vector < precision > &invec, std::vector < precision > &outvec, double &expectation_value, const Op& o) {
        hamiltonian().storage().reset();
        long long k = 0;
        int sign = 0;
        outvec.assign(hamiltonian().storage().vector_size(_model.symmetry().sector()), 0.0);
        int i = 0;
        while (_model.symmetry().next_state()) {
          long long nst = _model.symmetry().state();
          outvec[i] = o.action(nst, orbital) * invec[i];
          ++i;
        };
        double norm = hamiltonian().storage().vv(outvec, outvec);
        for (int j = 0; j < outvec.size(); ++j) {
          outvec[j] /= std::sqrt(norm);
        }
        _model.symmetry().init();
        expectation_value = norm;
        return true;
      };
    };
  }
}


#endif //HUBBARD_CHILOC_H
