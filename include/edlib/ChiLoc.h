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

    template<typename precision>
    class SzOperator: public BosonicOperator {
    public:
      SzOperator() {};
      template<class ModelType>
      precision action(long long state, int iii, const ModelType & model) const {
        /// 0.5 * (n_up - n_down)
        return 0.5*(model.checkState(state, iii, model.max_total_electrons()) -
                    model.checkState(state, iii + model.orbitals(), model.max_total_electrons()));
      }
      /// Current implementation restricted to paramagnetic solution
      precision average () const {
        return 0.0;
      }

      std::string name() const {return "Sz";};
    };

    template<typename precision>
    class NOperator: public BosonicOperator {
    public:
      NOperator() {};
      template<class ModelType>
      precision action(long long state, int iii, const ModelType & model) const {
        /// (n_up + n_down)
        return (model.checkState(state, iii, model.max_total_electrons()) +
                model.checkState(state, iii + model.orbitals(), model.max_total_electrons()));
      }
      /// Current implementation restricted to half-filled
      precision average() const {
        return 1.0;
      }
      std::string name() const {return "N";};
    };

    template<class Hamiltonian>
    class ChiLoc : public Lanczos < Hamiltonian > {
      using Lanczos < Hamiltonian >::hamiltonian;
      using Lanczos < Hamiltonian >::lanczos;
      using Lanczos < Hamiltonian >::omega;
      using Lanczos < Hamiltonian >::compute_sym_continues_fraction;
      using typename Lanczos < Hamiltonian >::precision;
    public:
      ChiLoc(alps::params &p, Hamiltonian &h) : Lanczos < Hamiltonian >(p, h, alps::gf::statistics::statistics_type::BOSONIC), _model(h.model()),
                                                        gf(Lanczos < Hamiltonian >::omega(), alps::gf::index_mesh(h.model().interacting_orbitals())),
                                                        _cutoff(p["lanc.BOLTZMANN_CUTOFF"]) {
        if(p["storage.EIGENVALUES_ONLY"] == 1) {
          throw std::logic_error("Eigenvectors have not been computed. Green's function can not be evaluated.");
        }
      }

      /**
       * Compute two-particle Green's function for specific operator Op.
       *
       * @tparam Op type of bosonic operator (should be either SzOperator or NOperator)
       */
      template<typename Op = SzOperator<precision> >
      void compute() {
        static_assert(std::is_base_of<BosonicOperator, Op>::value, "Wrong bosonic operator.");
        // cleanup
        gf *= 0.0;
        _Z = 0.0;
        if(hamiltonian().eigenpairs().empty())
          return;
#ifdef USE_MPI
        int rank;
        MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
#endif
        const EigenPair<precision, typename Hamiltonian::ModelType::Sector> &groundstate =  *hamiltonian().eigenpairs().begin();
        for (auto kkk = hamiltonian().eigenpairs().begin(); kkk != hamiltonian().eigenpairs().end(); kkk++) {
          const EigenPair<precision, typename Hamiltonian::ModelType::Sector> &eigenpair = *kkk;
          _Z += std::exp(-(eigenpair.eigenvalue() - groundstate.eigenvalue()) * omega().beta());
        }
        const Op op;
        for (auto kkk = hamiltonian().eigenpairs().begin(); kkk != hamiltonian().eigenpairs().end(); kkk++) {
          const EigenPair<precision, typename Hamiltonian::ModelType::Sector>& pair = *kkk;
          precision boltzmann_f = std::exp(-(pair.eigenvalue() - groundstate.eigenvalue()) * omega().beta());
          if (boltzmann_f < _cutoff) {
//        std::cout<<"Skipped by Boltzmann factor."<<std::endl;
            continue;
          }
          std::cout << "Compute Green's function contribution for eigenvalue E=" << pair.eigenvalue() << " with Boltzmann factor = "
                    << boltzmann_f << "; for sector" << pair.sector() << std::endl;
          for (int i = 0; i < _model.interacting_orbitals(); ++i) {
            std::vector < precision > outvec(1, precision(0.0));
            precision expectation_value = 0;
            _model.symmetry().set_sector(pair.sector());
            if (operation(i, pair.eigenvector(), outvec, expectation_value, op)) {
              int nlanc = lanczos(outvec);
#ifdef USE_MPI
              if(rank==0){
#endif
              std::cout << "orbital: " << i << " <n|O|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
              compute_sym_continues_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, 1, gf, alps::gf::index_mesh::index_type(i));
#ifdef USE_MPI
            }
#endif
            }
          }
        }
#ifdef USE_MPI
        if(rank == 0) {
#endif
        gf/= _Z;
        /// Compute static susceptibility
        for (int i = 0; i < _model.orbitals(); ++i) {
          double chiSum = 0.0;
          for (int iomega = 1; iomega < omega().extent(); ++iomega) {
            chiSum = chiSum + gf(alps::gf::matsubara_positive_mesh::index_type(iomega), alps::gf::index_mesh::index_type(i)).real();
          }
          gf(alps::gf::matsubara_positive_mesh::index_type(0), alps::gf::index_mesh::index_type(i)) -= 2 * chiSum + op.average()*omega().beta();
        }
        std::ostringstream Gomega_name;
        Gomega_name << "Chi"<<op.name()<<"_omega";
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
          outvec[i] = o.action(nst, orbital, _model) * invec[i];
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
