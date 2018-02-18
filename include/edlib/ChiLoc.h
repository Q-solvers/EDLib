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
    /**
     * Base class for susceptibility bosonic operator
     * @tparam precision - floating point precision
     */
    template<typename precision>
    class BosonicOperator{
    public:
      /**
       * @param avg - average value for zero Matsubara frequency correction
       */
      BosonicOperator(precision avg) : _avg(avg) {}
      precision average () const {
        return _avg;
      }
    private:
      /// operator average value
      precision _avg;
    };

    /**
     * Spin operator
     * @tparam precision -- floating point precision
     */
    template<typename precision>
    class SzOperator: public BosonicOperator<precision> {
    public:
      /**
       * @param avg -- average value of spin operator. By default we consider paramagnetic system and
       * average spin is 0.
       */
      SzOperator(precision avg = 0.0) : BosonicOperator<precision>(avg) {};
      /**
       *
       * @tparam ModelType -- specific model type
       * @param state - occupation basis state
       * @param iii - electonic orbital
       * @param model - model instance
       * @return value of the total spin for site 'iii' in the basis state 'state'
       */
      template<class ModelType>
      precision action(long long state, int iii, const ModelType & model) const {
        // 0.5 * (n_up - n_down)
        return 0.5*(model.checkState(state, iii, model.max_total_electrons()) -
                    model.checkState(state, iii + model.orbitals(), model.max_total_electrons()));
      }
      /**
       * @return operator name
       */
      std::string name() const {return "Sz";};
    };

    /**
     * Charge operator
     * @tparam precision - floating point precision
     */
    template<typename precision>
    class NOperator: public BosonicOperator<precision> {
    public:
      NOperator(precision avg = 1.0) : BosonicOperator<precision>(avg) {};
      /**
       *
       * @tparam ModelType -- specific model type
       * @param state - occupation basis state
       * @param iii - electonic orbital
       * @param model - model instance
       * @return value of the total charge for site 'iii' in the basis state 'state'
       */
      template<class ModelType>
      precision action(long long state, int iii, const ModelType & model) const {
        // (n_up + n_down)
        return (model.checkState(state, iii, model.max_total_electrons()) +
                model.checkState(state, iii + model.orbitals(), model.max_total_electrons()));
      }
      std::string name() const {return "N";};
    };

    /**
     * General class for local susceptibilities calculation.
     *
     * @tparam Hamiltonian - Hamiltonian instance type
     * @tparam Mesh - type of frequency mesh. can be either alps::gf::real_frequency_mesh or alps::gf::matsubara_positive_only
     * @tparam Args - additional parameters for Mesh. For matsubara mesh should be alps::gf::statistics::statistics_type
     */
    template<class Hamiltonian, class Mesh, typename ... Args>
    class ChiLoc : public Lanczos < Hamiltonian, Mesh, Args... > {
      using Lanczos < Hamiltonian, Mesh, Args... >::zero_freq;
      using Lanczos < Hamiltonian, Mesh, Args... >::omega;
      using Lanczos < Hamiltonian, Mesh, Args... >::lanczos;
      using Lanczos < Hamiltonian, Mesh, Args... >::hamiltonian;
      using Lanczos < Hamiltonian, Mesh, Args... >::beta;
      using Lanczos < Hamiltonian, Mesh, Args... >::compute_sym_continued_fraction;
      using typename Lanczos < Hamiltonian, Mesh, Args... >::precision;
      using Sector = typename Hamiltonian::ModelType::Sector;
      /// Green's function conatainer type
      typedef alps::gf::two_index_gf<std::complex<double>, Mesh, alps::gf::index_mesh>  GF_TYPE;
    public:
      /**
       *
       * @param p - ALPSCore parameters
       * @param h - Hamiltonian instance
       * @param args - additional parameters for mesh. For example for matsubara mesh should it be alps::gf::statistics::statistics_type::BOSONIC
       */
      ChiLoc(alps::params &p, Hamiltonian &h, Args... args) : Lanczos < Hamiltonian, Mesh, Args... >(p, h, args...), _model(h.model()),
                                                        gf(Lanczos < Hamiltonian, Mesh, Args... >::omega(), alps::gf::index_mesh(h.model().interacting_orbitals())),
                                                        gf_ij(Lanczos < Hamiltonian, Mesh, Args... >::omega(), alps::gf::index_mesh(h.model().interacting_orbitals()*h.model().interacting_orbitals())),
                                                        _cutoff(p["lanc.BOLTZMANN_CUTOFF"]), _type("Sz") {
        // we can not evaluate Green's function if eigenvectors have not been computed.
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
            offdiagonal_orbs.push_back(gf_orbs[i]);
          }
        }
      }

      /**
       * Compute two-particle Green's function for specific operator Op.
       *
       * @tparam Op type of bosonic operator (should be either SzOperator or NOperator)
       * @param avg_ptr - pointer to the operator average value. If null -- the default average vlues will be used for specific operator.
       */
      template<typename Op = SzOperator<precision> >
      void compute(const double * avg_ptr = nullptr) {
        static_assert(std::is_base_of<BosonicOperator<precision>, Op>::value, "Wrong bosonic operator.");
        // reset to 0
        gf *= 0.0;
        gf_ij *= 0.0;
        _Z = 0.0;
        // check that at least one eigen-value have computed
        if(hamiltonian().eigenpairs().empty())
          return;
#ifdef USE_MPI
        int rank;
        MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
#endif
        // get groundstate
        const EigenPair<precision, typename Hamiltonian::ModelType::Sector> &groundstate =  *hamiltonian().eigenpairs().begin();
        // compute partition function
        for (auto kkk = hamiltonian().eigenpairs().begin(); kkk != hamiltonian().eigenpairs().end(); kkk++) {
          const EigenPair<precision, typename Hamiltonian::ModelType::Sector> &eigenpair = *kkk;
          _Z += std::exp(-(eigenpair.eigenvalue() - groundstate.eigenvalue()) * beta());
        }
        // init bosonic operator
        const Op op = (avg_ptr == nullptr ? Op() : Op(*avg_ptr));
        _type = op.name();
        // loop over all eigen-pairs
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
          local_contribution(groundstate, op, pair);
          nonlocal_contribution(groundstate, op, pair);
        }
#ifdef USE_MPI
        if(rank == 0) {
#endif
        // normalize Green's function
        gf /= _Z;
        gf_ij /= _Z;
        // Compute zero-frequency contribution
        local_correction(op);
        // Compute non-local correction
        non_local_correction(op);
#ifdef USE_MPI
        }
#endif
      }

      /**
       * Evaluates static correction at zero matsubara frequency
       *
       * @tparam O -- bosonic operator type
       * @param op -- bosonic operator
       */
      template<typename O>
      void local_correction(const O&op) {
        for (int i = 0; i < _model.interacting_orbitals(); ++i) {
          zero_freq_contribution(op, gf, i);
        }
      }

      /**
       * Computes non-local Green's function from symmetrized averaged: G_ij = 0.5( <(op_i + op_j)(op_i + op_j)> - <op_i op_i> - <op_j op_j> ),
       * and evaluates static correction at zero matsubara frequency
       *
       *
       * @tparam O -- bosonic operator type
       * @param op -- bosonic operator
       */
      template<typename O>
      void non_local_correction(const O& op) {
        for (int i = 0; i < _model.interacting_orbitals(); ++i) {
          for (int j = 0; j < _model.interacting_orbitals(); ++j) {
            zero_freq_contribution(op, gf_ij, i*_model.interacting_orbitals() + j);
            for (int iomega = 0; iomega < omega().extent(); ++iomega) {
              gf_ij(typename Mesh::index_type(iomega), alps::gf::index_mesh::index_type(i*_model.interacting_orbitals() + j)) -=
                 gf(typename Mesh::index_type(iomega), alps::gf::index_mesh::index_type(i)) + gf(typename Mesh::index_type(iomega), alps::gf::index_mesh::index_type(j));
            }
          }
        }
        gf_ij *= 0.5;
      };

      /**
       * Compute local Green's function G_ii
       *
       * @tparam Op -- type of bosonic operator
       * @param groundstate -- system groundstate
       * @param op -- bosonic operator
       * @param pair -- current Eigen-Pair
       */
      template<typename Op = SzOperator<precision> >
      void local_contribution(const EigenPair <precision, Sector> &groundstate, const Op op,
                              const EigenPair <precision, Sector> &pair) {
#ifdef USE_MPI
        int rank;
        MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
#endif
        for (int iorb = 0; iorb < diagonal_orbs.size(); ++iorb) {
          int i = diagonal_orbs[iorb];
          std::vector < precision > outvec(1, precision(0.0));
          precision expectation_value = 0;
          _model.symmetry().set_sector(pair.sector());
          if (operation(i, pair.eigenvector(), outvec, expectation_value, op)) {
            int nlanc = lanczos(outvec);
#ifdef USE_MPI
            if(rank==0){
#endif
              std::cout << "orbital: " << i << " <n|" + op.name() + op.name() + "|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
              // compute symmetrized Lanczos continued fraction
              compute_sym_continued_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, 1,
                                             gf, alps::gf::index_mesh::index_type(i));
#ifdef USE_MPI
            }
#endif
          }
        }
      }
      /**
       * Computes symmetrized non-local Green's function for operator op: G_ij = <(op_i + op_j)(op_i + op_j)>
       *
       * @tparam Op -- type of bosonic operator
       * @param groundstate -- system groundstate
       * @param op -- bosonic operator
       * @param pair -- current Eigen-Pair
       */
      template<typename Op>
      void nonlocal_contribution(const EigenPair <precision, Sector> &groundstate, const Op &op,
                                 const EigenPair <precision, Sector> &pair) {
#ifdef USE_MPI
        int rank;
        MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
#endif
        for (int i = 0; i < _model.interacting_orbitals(); ++i) {
          for (int j = 0; j < _model.interacting_orbitals(); ++j) {
            _model.symmetry().set_sector(pair.sector());
            std::vector < precision > outvec(1, precision(0.0));
            precision expectation_value = 0;
            if(operation(i, j, pair.eigenvector(), outvec, expectation_value, op)) {
              int nlanc = lanczos(outvec);
#ifdef USE_MPI
        if(rank == 0) {
#endif
              std::cout << "orbital: " << i << " <n|" + op.name() + op.name() + "|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
              // compute symmetrized Lanczos continued fraction
              compute_sym_continued_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, 1,
                                             gf_ij, alps::gf::index_mesh::index_type(i * _model.interacting_orbitals() + j));
#ifdef USE_MPI
              }
#endif
            }
          }
        }
      }

      /**
       * zero Matsubara frequency contribution. Since for bosonic Green's function Lanczos continued fraction can not
       * properly evaluate zero frequency contribution we compute it from the following sum-rule:
       * \Chi(w_0) = \beta \Chi(tau = 0) - 2 \sum_{n!=0} \Chi(w_n)
       *
       * @tparam O - type of operator
       * @tparam M - Mesh type
       * @param op - operator instance
       */
      template<typename O, typename M= Mesh>
      typename std::enable_if<std::is_base_of<alps::gf::matsubara_positive_mesh, M>::value, void>::type zero_freq_contribution(const O& op, GF_TYPE& G, int i) {
        // Compute static susceptibility for each orbital

        double chiSum = 0.0;
        // Susceptibility decays as c2/w^2 + c4/w^4
        // compute c2 and c4 from two largest freq points
        // and for next pair as well to check convergence
        double c2, c4;
        double c2_2, c4_2;
        double tail, tail2;
        // check that we have enough Matsubara frequencies and we already sit in the high-frequency tail regime
        get_tail(i, omega().extent()-1, G, c2, c4, tail);
        get_tail(i, omega().extent()-2, G, c2_2, c4_2, tail2);
        if(std::abs(tail-tail2)/std::abs(tail) > 1e-4) {
          std::cerr<<"Not enough frequencies to compute high frequency tail. Please increase number of frequencies. Diff: "<<std::abs(tail-tail2)/std::abs(tail)<<std::endl;
        }
        // loop over non-zero Matsubara frequencies
        for (int iomega = 1; iomega < omega().extent(); ++iomega) {
          double om = omega().points()[iomega];
          chiSum = chiSum + G(typename Mesh::index_type(iomega), alps::gf::index_mesh::index_type(i)).real() - c2/(om*om) - c4/(om*om*om*om);
        }
        // computes zero-frequency contribution
        G(typename Mesh::index_type(0), alps::gf::index_mesh::index_type(i)) -= 2 * chiSum + 2* tail - op.average()*beta();
      };

      /**
       * For real frequency mesh there is nothing to do. Zero frequency has been already correctly computed.
       */
      template<typename O, typename M= Mesh>
      typename std::enable_if<std::is_base_of<alps::gf::real_frequency_mesh, M>::value, void>::type zero_freq_contribution(const O& op, GF_TYPE& G, int i) {

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
      void get_tail(int i, int freq, GF_TYPE& G, double &c2, double& c4, double &tail) const {
        double om1 = omega().points()[freq];
        double om2 = omega().points()[freq-1];
        double om1_2 = om1*om1;
        double om2_2 = om2*om2;
        double g1 = G(typename Mesh::index_type(freq), alps::gf::index_mesh::index_type(i)).real();
        double g2 = G(typename Mesh::index_type(freq - 1), alps::gf::index_mesh::index_type(i)).real();
        c2 = - (g2*om2_2*om2_2 - g1*om1_2*om1_2)/(om1_2-om2_2);
        c4 = - (g1*om1_2*om1_2*om2_2 - g2*om2_2*om2_2*om1_2)/(om1_2-om2_2);
        tail= c2 * beta() * beta() / 24.0 + c4 * beta() * beta() * beta() * beta() / 1440.0;
      }

      /**
       * Save susceptibility to HDF5 archive and to the text-file
       * @param ar - hdf5 archive file
       * @param path - root path in hdf5 archive
       */
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

          std::ostringstream Gomega_name2;
          Gomega_name2 << "Chi_ij_"<<_type<<"_omega";
          std::ofstream G_omega_file2(Gomega_name2.str().c_str());
          G_omega_file2<< std::setprecision(14) << gf_ij;
          G_omega_file2.close();
#ifdef USE_MPI
        }
#endif
      }

    private:
      /// Green's function conatainer
      GF_TYPE gf;
      /// Nonlocal
      GF_TYPE gf_ij;
      /// Specific model
      typename Hamiltonian::ModelType &_model;
      /// Boltzmann-factor cut-off
      precision _cutoff;
      /// Partition function
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
        double norm = hamiltonian().storage().vv(outvec, outvec
#ifdef USE_MPI
            , hamiltonian().comm()
#endif
        );
        for (int j = 0; j < outvec.size(); ++j) {
          outvec[j] /= std::sqrt(norm);
        }
        _model.symmetry().init();
        expectation_value = norm;
        return expectation_value > 1e-9;
      };

      /**
       * @brief Perform the symmetrized operation on the eigenstate
       *
       * @param mu - the first orbital to action
       * @param mu - the second orbital to action
       * @param invec - current eigenstate
       * @param outvec - Op-vec product
       * @param expectation_value - expectation value of aa*
       * @return true if the particle has been created
       */
      template<typename Op>
      bool operation(int mu, int nu, const std::vector < precision > &invec, std::vector < precision > &outvec, double &expectation_value, const Op& o) {
        hamiltonian().storage().reset();
        long long k = 0;
        int sign = 0;
        outvec.assign(hamiltonian().storage().vector_size(_model.symmetry().sector()), 0.0);
        for(int i = 0; i< invec.size(); ++i) {
          _model.symmetry().next_state();
          long long nst = _model.symmetry().state();
          outvec[i] = (o.action(nst, mu, _model) + o.action(nst, nu, _model) )* invec[i];
        };
        double norm = hamiltonian().storage().vv(outvec, outvec
#ifdef USE_MPI
            , hamiltonian().comm()
#endif
        );
        for (int i = 0; i < outvec.size(); ++i) {
          outvec[i] /= std::sqrt(norm);
        }
        _model.symmetry().init();
        expectation_value = norm;
        return true;//expectation_value > 1e-9;
      };
    };
  }
}


#endif //HUBBARD_CHILOC_H
