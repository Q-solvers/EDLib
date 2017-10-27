//
// Created by iskakoff on 03/01/17.
//

#ifndef HUBBARD_CHILOC_H
#define HUBBARD_CHILOC_H


#include <iostream>
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
      /// Current implementation restricted to half-filled case
      precision average() const {
        return 1.0;
      }
      std::string name() const {return "N";};
    };

    template<class Hamiltonian, class Mesh, typename ... Args>
    class ChiLoc : public Lanczos < Hamiltonian, Mesh, Args... > {
      using Lanczos < Hamiltonian, Mesh, Args... >::hamiltonian;
      using Lanczos < Hamiltonian, Mesh, Args... >::lanczos;
      using Lanczos < Hamiltonian, Mesh, Args... >::omega;
      using Lanczos < Hamiltonian, Mesh, Args... >::beta;
      using Lanczos < Hamiltonian, Mesh, Args... >::compute_sym_continued_fraction;
      using typename Lanczos < Hamiltonian, Mesh, Args... >::precision;
    public:
      ChiLoc(alps::params &p, Hamiltonian &h, Args... args) : Lanczos < Hamiltonian, Mesh, Args... >(p, h, args...), _model(h.model()),
                                                        gf(Lanczos < Hamiltonian, Mesh, Args... >::omega(), alps::gf::index_mesh(h.model().interacting_orbitals())),
                                                        _cutoff(p["lanc.BOLTZMANN_CUTOFF"]), _type("Sz") {
        if(p["storage.EIGENVALUES_ONLY"] == 1) {
          throw std::logic_error("Eigenvectors have not been computed. Green's function can not be evaluated.");
        }
        std::string input = p["INPUT_FILE"];
        alps::hdf5::archive input_file(input.c_str(), "r");
        if(input_file.is_data("ChiLoc_orbitals/values")){
          input_file >> alps::make_pvp("ChiLoc_orbitals/values", gf_orbs);
        }else{
          // Or calculate only the diagonal part.
          gf_orbs.clear();
          for(int i = 0; i < h.model().interacting_orbitals(); ++i){
            gf_orbs.push_back({i, i});
          }
        }
        input_file.close();
        // Find all unique indices for the diagonal part.
        diagonal_orbs.resize(gf_orbs.size() * 2);
        for(int i = 0; i < gf_orbs.size(); ++i){
         for(int j = 0; j < 2; ++j){
           diagonal_orbs[i + j * gf_orbs.size()] = gf_orbs[i][j];
         }
        }
        std::sort(diagonal_orbs.begin(), diagonal_orbs.end());
        diagonal_orbs.erase(std::unique(diagonal_orbs.begin(), diagonal_orbs.end()), diagonal_orbs.end());
        // Find offdiagonal indices.
        offdiagonal_orbs.clear();
        for(int i = 0; i < gf_orbs.size(); ++i){
          if(gf_orbs[i][0] != gf_orbs[i][1]){
            throw std::logic_error("Offdiagonal susceptibility is not supported yet.");
            offdiagonal_orbs.push_back(gf_orbs[i]);
          }
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
          _Z += std::exp(-(eigenpair.eigenvalue() - groundstate.eigenvalue()) * beta());
        }
        const Op op;
        _type = op.name();
        for (auto kkk = hamiltonian().eigenpairs().begin(); kkk != hamiltonian().eigenpairs().end(); kkk++) {
          const EigenPair<precision, typename Hamiltonian::ModelType::Sector>& pair = *kkk;
          precision boltzmann_f = std::exp(-(pair.eigenvalue() - groundstate.eigenvalue()) * beta());
          if (boltzmann_f < _cutoff) {
            // Skipped by Boltzmann factor.
            continue;
          }
#ifdef USE_MPI
          if(rank == 0)
#endif
          std::cout << "Compute Green's function contribution for eigenvalue E=" << pair.eigenvalue() << " with Boltzmann factor = "
                    << boltzmann_f << "; for sector" << pair.sector() << std::endl;
          for (int iorb = 0; iorb < diagonal_orbs.size(); ++iorb) {
            std::vector < precision > outvec(1, precision(0.0));
            precision expectation_value = 0;
            _model.symmetry().set_sector(pair.sector());
            if (operation(diagonal_orbs[iorb], pair.eigenvector(), outvec, expectation_value, op)) {
              int nlanc = lanczos(outvec);
#ifdef USE_MPI
              if(rank==0){
#endif
              std::cout << "orbital: " << diagonal_orbs[iorb] << " <n|" + op.name() + op.name() + "|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
                compute_sym_continued_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, 1, gf, alps::gf::index_mesh::index_type(diagonal_orbs[iorb]));
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
        zero_freq_contribution(op);
#ifdef USE_MPI
        }
#endif
      }

      template<typename O, typename M= Mesh>
      typename std::enable_if<std::is_base_of<alps::gf::matsubara_positive_mesh, M>::value, void>::type zero_freq_contribution(const O& op) {
        /// Compute static susceptibility
        for (int i = 0; i < _model.orbitals(); ++i) {
          double chiSum = 0.0;
          /// Susceptibility decays as c2/w^2
          /// compute c2 and c4 from two largest freq points
          /// and for next pair as well to check convergence
          double c2, c4;
          double c2_2, c4_2;
          double tail, tail2;
          get_tail(i, omega().extent()-1, c2, c4, tail);
          get_tail(i, omega().extent()-2, c2_2, c4_2, tail2);
          if(std::abs(tail-tail2)/std::abs(tail) > 1e-4) {
            std::cerr<<"Not enough frequencies to compute high frequency tail. Please increase number of frequencies. Diff: "<<std::abs(tail-tail2)/std::abs(tail)<<std::endl;
          }
          for (int iomega = 1; iomega < omega().extent(); ++iomega) {
            double om = omega().points()[iomega];
            chiSum = chiSum + gf(typename Mesh::index_type(iomega), alps::gf::index_mesh::index_type(i)).real() - c2/(om*om) - c4/(om*om*om*om);
          }
          gf(typename Mesh::index_type(0), alps::gf::index_mesh::index_type(i)) -= 2 * chiSum + 2* tail - op.average()*beta();
        }
      };

      template<typename O, typename M= Mesh>
      typename std::enable_if<std::is_base_of<alps::gf::real_frequency_mesh, M>::value, void>::type zero_freq_contribution(const O& op) {

      };

      /**
       * High frequency behaviour of Chi(W_n) = c_2/(W_n)^2 + c_4/(W_n)^4
       *
       * @param i -- current orbital
       * @param freq -- frequency
       * @param c2 -- c_2 coefficient
       * @param c4 -- c_4 coefficient
       * @param tail -- analytical value of tail
       */
      void get_tail(int i, int freq, double &c2, double& c4, double &tail) const {
        double om1 = omega().points()[freq];
        double om2 = omega().points()[freq-1];
        double om1_2 = om1*om1;
        double om2_2 = om2*om2;
        double g1 = gf(typename Mesh::index_type(freq), alps::gf::index_mesh::index_type(i)).real();
        double g2 = gf(typename Mesh::index_type(freq - 1), alps::gf::index_mesh::index_type(i)).real();
        c2 = - (g2*om2_2*om2_2 - g1*om1_2*om1_2)/(om1_2-om2_2);
        c4 = - (g1*om1_2*om1_2*om2_2 - g2*om2_2*om2_2*om1_2)/(om1_2-om2_2);
        tail= c2 * beta() * beta() / 24.0 + c4 * beta() * beta() * beta() * beta() / 1440.0;
      }

      void save(alps::hdf5::archive& ar, const std::string & path) {
#ifdef USE_MPI
        int rank;
        MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
        if(rank == 0) {
#endif
          gf.save(ar, path + "/Chi" + _type +"_omega");
          std::ostringstream Gomega_name;
          Gomega_name << "Chi"<<_type<<"_omega";
          std::ofstream G_omega_file(Gomega_name.str().c_str());
          G_omega_file << std::setprecision(14) << gf;
          G_omega_file.close();
          std::cout << "Statsum: " << _Z << std::endl;
          ar[path + "/@Statsum"] << _Z;
#ifdef USE_MPI
        }
#endif
      }

    private:
      typedef alps::gf::two_index_gf<std::complex<double>, Mesh, alps::gf::index_mesh>  GF_TYPE;
      GF_TYPE gf;
      typename Hamiltonian::ModelType &_model;
      precision _cutoff;
      precision _Z;
      /// type of recently computed susceptibility
      std::string _type;
      /// Orbital pairs to calculate the Green's function.
      std::vector<std::vector<int>> gf_orbs;
      /// Orbitals to calculate the diagonal Green's function.
      std::vector<int> diagonal_orbs;
      /// Orbital pairs to calculate the offdiagonal Green's function.
      std::vector<std::vector<int>> offdiagonal_orbs;

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
        for(int i = 0; i< invec.size(); ++i) {
          _model.symmetry().next_state();
          long long nst = _model.symmetry().state();
          outvec[i] = o.action(nst, orbital, _model) * invec[i];
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
