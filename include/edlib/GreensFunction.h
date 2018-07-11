//
// Created by iskakoff on 01/08/16.
//

#ifndef HUBBARD_GREENSFUNCTION_H
#define HUBBARD_GREENSFUNCTION_H


#include <iomanip>
#include "Lanczos.h"
#include "EigenPair.h"
#include "ExecutionStatistic.h"

namespace EDLib {
  namespace gf {
    /**
     * Class for evaluation of the single-particle Green's function.
     *
     *
     *
     * @tparam Hamiltonian - type of Hamiltonian object
     * @tparam Mesh - type of frequency mesh. can be either alps::gf::real_frequency_mesh or alps::gf::matsubara_positive_only
     * @tparam Args - additional Mesh parametrization for alps::gf::matsubara_positive_only
     */
    template<class Hamiltonian, typename Mesh, typename... Args>
    class GreensFunction : public Lanczos < Hamiltonian, Mesh, Args...> {
      using Lanczos < Hamiltonian, Mesh, Args... >::hamiltonian;
      using Lanczos < Hamiltonian, Mesh, Args... >::lanczos;
      using Lanczos < Hamiltonian, Mesh, Args... >::omega;
      using Lanczos < Hamiltonian, Mesh, Args... >::beta;
      using Lanczos < Hamiltonian, Mesh, Args... >::compute_continued_fraction;
      using Lanczos < Hamiltonian, Mesh, Args... >::suffix;
      using typename Lanczos < Hamiltonian, Mesh, Args... >::precision;
    public:
      /**
       * Construct Green's function class for calculating
       * G(t-t') = -<T[c(t')c^{*}(t)]>, where T is time ordering operator
       *
       * By default, all local Green functions are computed. An array of orbital
       * pairs can be supplied in the input file as the GreensFunction_orbitals
       * group.
       *
       * @param p - AlpsCore parameter object
       * @param h - Hamiltonain instance
       * @param args - additional parameters for Mesh. For example for Matsubara mesh should it be alps::gf::statistics::statistics_type::FERMIONIC
       */
      GreensFunction(alps::params &p, Hamiltonian &h, Args ... args) : Lanczos < Hamiltonian, Mesh, Args... >(p, h, args...), _model(h.model()),
                                                        _G_g(Lanczos < Hamiltonian, Mesh, Args... >::omega(), alps::gf::index_mesh(h.model().interacting_orbitals()), alps::gf::index_mesh(p["NSPINS"].as<int>())),
                                                        _G_l(Lanczos < Hamiltonian, Mesh, Args... >::omega(), alps::gf::index_mesh(h.model().interacting_orbitals()), alps::gf::index_mesh(p["NSPINS"].as<int>())),
                                                        _G(Lanczos < Hamiltonian, Mesh, Args... >::omega(), alps::gf::index_mesh(h.model().interacting_orbitals()), alps::gf::index_mesh(p["NSPINS"].as<int>())),
                                                        _G_g_ij(Lanczos < Hamiltonian, Mesh, Args... >::omega(), alps::gf::index_mesh(h.model().interacting_orbitals() * h.model().interacting_orbitals()), alps::gf::index_mesh(p["NSPINS"].as<int>())),
                                                        _G_l_ij(Lanczos < Hamiltonian, Mesh, Args... >::omega(), alps::gf::index_mesh(h.model().interacting_orbitals() * h.model().interacting_orbitals()), alps::gf::index_mesh(p["NSPINS"].as<int>())),
                                                        _G_ij(Lanczos < Hamiltonian, Mesh, Args... >::omega(), alps::gf::index_mesh(h.model().interacting_orbitals() * h.model().interacting_orbitals()), alps::gf::index_mesh(p["NSPINS"].as<int>())),
                                                        _cutoff(p["lanc.BOLTZMANN_CUTOFF"]) {
        if(p["storage.EIGENVALUES_ONLY"] == 1) {
          throw std::logic_error("Eigenvectors have not been computed. Green's function can not be evaluated.");
        }
        std::string input = p["INPUT_FILE"];
        alps::hdf5::archive input_file(input.c_str(), "r");
        std::vector<std::vector<int>> gf_orbs;
        if(input_file.is_data("GreensFunction_orbitals/values")){
          input_file >> alps::make_pvp("GreensFunction_orbitals/values", gf_orbs);
        }else{
          // Or calculate only the diagonal part.
          gf_orbs.clear();
          for(int i = 0; i < h.model().interacting_orbitals(); ++i){
            gf_orbs.push_back({i, i});
          }
        }
        input_file.close();
        for(size_t ii = 0; ii < gf_orbs.size(); ++ii){
          if(gf_orbs[ii][0] == gf_orbs[ii][1]){
           _g_orbs.push_back(gf_orbs[ii][0]);
          }else{
           _g_ij_orb_pairs.push_back(std::array<size_t, 2>{size_t(gf_orbs[ii][0]), size_t(gf_orbs[ii][1])});
          }
        }
        // Find all unique indices for the diagonal part.
        std::sort(_g_orbs.begin(), _g_orbs.end());
        _g_orbs.erase(std::unique(_g_orbs.begin(), _g_orbs.end()), _g_orbs.end());
        // Check that we will have the two local GFs required by nonlocal GF.
        for(size_t ii = 0; ii < _g_ij_orb_pairs.size(); ++ii){
          for(size_t jj = 0; jj < 2; ++jj){
            bool found = false;
            for(size_t kk = 0; kk < _g_orbs.size(); ++kk){
              if(_g_orbs[kk] == _g_ij_orb_pairs[ii][jj]){
                found = true;
                break;
              }
            }
            if(!found){
              std::ostringstream msg;
              throw std::logic_error("Nonlocal GF requires local GFs for both orbitals.");
            }
          }
        }
      }

      /**
       * Evaluates Green's function by Lanczos continued fraction method.
       * For each eigenvalue checks Boltsman factor cut-off
       */
      void compute() {
        // init Green's function with 0
        _G_g *= 0.0;
        _G_l *= 0.0;
        // init Partition function
        _Z = 0.0;
        // check that at least one eigen-value have been computed
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
        // iterate over eigen-pairs
        for (auto kkk = hamiltonian().eigenpairs().begin(); kkk != hamiltonian().eigenpairs().end(); kkk++) {
          const EigenPair<precision, typename Hamiltonian::ModelType::Sector>& pair = *kkk;
          // compute Boltzmann-factor
          precision boltzmann_f = std::exp(-(pair.eigenvalue() - groundstate.eigenvalue()) * beta());
          // Skip all eigenvalues with Boltzmann-factor smaller than cutoff
          if (boltzmann_f < _cutoff) {
            continue;
          }
#ifdef USE_MPI
          if(rank == 0)
#endif
          std::cout << "Compute Green's function contribution for eigenvalue E=" << pair.eigenvalue() << " with Boltzmann factor = " << boltzmann_f << "; for sector" << pair.sector() << std::endl;
          local_contribution(pair, groundstate);
          nonlocal_contribution(pair, groundstate);
        }
#ifdef USE_MPI
        if(rank == 0) {
#endif
        // normalize Green's function
        _G_g /= _Z;
        _G_l /= _Z;
        _G_g_ij /= _Z;
        _G_l_ij /= _Z;
#ifdef USE_MPI
        }
#endif
        non_local_gf();
        _G = _G_g + _G_l;
        common::statistics.updateEvent("Greens function");
      }

      /**
       * Save Green's function into the hdf5 archive and in plain text file
       *
       * @param ar -- hdf5 archive to save Green's function
       * @param path -- root path in hdf5 archive
       */
      void save(alps::hdf5::archive& ar, const std::string & path) {
#ifdef USE_MPI
        int rank;
        MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
        if(rank == 0) {
#endif
          if(_g_orbs.size()){
            _G.save(ar, path + "/G_omega" + suffix());
            std::ostringstream Gomega_name;
            Gomega_name << "G_omega"<<suffix();
            std::ofstream G_omega_file(Gomega_name.str().c_str());
            G_omega_file << std::setprecision(14) << _G;
            G_omega_file.close();
          }
          std::cout << "Statsum: " << _Z << std::endl;
          ar[path + "/@Statsum"] << _Z;
          if(_g_ij_orb_pairs.size()){
            std::ostringstream Gomega_name2;
            Gomega_name2 << "G_ij_omega"<<suffix();
            std::ofstream G_omega_file2(Gomega_name2.str().c_str());
            G_omega_file2<< std::setprecision(14) << _G_ij;
            G_omega_file2.close();
          }
#ifdef USE_MPI
        }
#endif
      }

      /**
       * Comutes self-energy from Dyson equation based on model specific bare Green's function
       *
       * @param ar -- hdf5 archive to save self-energy
       * @param path -- root path in hdf5 archive
       */
      void compute_selfenergy(alps::hdf5::archive &ar, const std::string &path){
        // Bare Green's function and Self-energy should be defined on the same grid
        GF_TYPE bare(_G.mesh1(), _G.mesh2(), _G.mesh3());
        GF_TYPE sigma(_G.mesh1(), _G.mesh2(), _G.mesh3());
        // obtain model-specific bare Green's function
        _model.bare_greens_function(bare, beta());
        bare.save(ar, path + "/G0_omega");
        std::ostringstream Gomega_name;
        Gomega_name << "G0_omega";
        std::ofstream G_omega_file(Gomega_name.str().c_str());
        G_omega_file << std::setprecision(14) << bare;
        G_omega_file.close();
        // solve Dyson equation
        for(int iw = 0; iw< bare.mesh1().points().size(); ++iw) {
          typename Mesh::index_type w(iw);
          for (int im: bare.mesh2().points()) {
            for (int is : bare.mesh3().points()) {
              sigma(w, alps::gf::index_mesh::index_type(im), alps::gf::index_mesh::index_type(is)) =
                1.0/bare(w, alps::gf::index_mesh::index_type(im), alps::gf::index_mesh::index_type(is)) - 1.0/_G(w, alps::gf::index_mesh::index_type(im), alps::gf::index_mesh::index_type(is));
            }
          }
        }
        // store to file
        sigma.save(ar, path + "/Sigma_omega");
        Gomega_name.str("");
        Gomega_name << "Sigma_omega";
        G_omega_file.open(Gomega_name.str().c_str());
        G_omega_file << std::setprecision(14) << sigma;
        G_omega_file.close();
      }

      /// Green's function type
      typedef alps::gf::three_index_gf<std::complex<double>, Mesh, alps::gf::index_mesh, alps::gf::index_mesh >  GF_TYPE;
      typedef typename alps::gf::index_mesh::index_type index_mesh_index;
      typedef typename Mesh::index_type frequency_mesh_index;

      const GF_TYPE &G_g() const {return _G_g;}
      const GF_TYPE &G_l() const {return _G_l;}
      const GF_TYPE &G() const {return _G;}
      const GF_TYPE &G_g_ij() const {return _G_g_ij;}
      const GF_TYPE &G_l_ij() const {return _G_l_ij;}
      const GF_TYPE &G_ij() const {return _G_ij;}

    private:

      /**
       * Compute local Green's function G_ii
       *
       * @param groundstate -- system groundstate
       * @param pair -- current Eigen-Pair
       */
      void local_contribution(const EigenPair<precision, typename Hamiltonian::ModelType::Sector>& pair, const EigenPair<precision, typename Hamiltonian::ModelType::Sector>& groundstate) {
#ifdef USE_MPI
        int rank;
        MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
#endif
        /// iterate over orbitals for the diagonal Green's function
        for (int iorb = 0; iorb < _g_orbs.size(); ++iorb) {
          /// iterate over spins
          for (int ispin = 0; ispin < _model.spins(); ++ispin) {
            int orb = _g_orbs[iorb];
            std::vector < precision > outvec(1, precision(0.0));
            precision expectation_value = 0.0;
            _model.symmetry().set_sector(pair.sector());
            /// first we are going to compute G_g = <c(t)c^+(0)>, for this we need to create particle in current eigenstate
            if (create_particles(std::array<size_t, 1>{{size_t(orb)}}, ispin, pair.eigenvector(), outvec, expectation_value)) {
              /// Perform Lanczos factorization for starting vector |outvec>
              int nlanc = lanczos(outvec);
#ifdef USE_MPI
              if(!rank)
#endif
              {
                std::cout << "orbital: " << orb << "   spin: " << (ispin == 0 ? "up" : "down") << " <n|aa*|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
                /// Using computed Lanczos factorization compute approximation for \frac{1}{z - H} by calculation of a continued fraction
                compute_continued_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, 1, _G_g, index_mesh_index(orb), index_mesh_index(ispin));
              }
            }
            /// restore symmetry sector
            _model.symmetry().set_sector(pair.sector());
            /// then compute G_l = <c^+(t)c(0)>, for this we need to destroy particle
            if (annihilate_particles(std::array<size_t, 1>{{size_t(orb)}}, ispin, pair.eigenvector(), outvec, expectation_value)) {
              int nlanc = lanczos(outvec);
#ifdef USE_MPI
              if(!rank)
#endif
              {
                std::cout << "orbital: " << orb << "   spin: " << (ispin == 0 ? "up" : "down") << " <n|a*a|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
                compute_continued_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, -1, _G_l, index_mesh_index(orb), index_mesh_index(ispin));
              }
            }
          }
        }
      }

      /**
       * Computes the part "<(c_i + c_j)(c_i^+ + c_j^+)>" of non-local Green's function G_ij
       *
       * @param groundstate -- system groundstate
       * @param pair -- current Eigen-Pair
       */
      void nonlocal_contribution(const EigenPair<precision, typename Hamiltonian::ModelType::Sector>& pair, const EigenPair<precision, typename Hamiltonian::ModelType::Sector>& groundstate) {
#ifdef USE_MPI
        int rank;
        MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
#endif
        for (int iorb = 0; iorb < _g_ij_orb_pairs.size(); ++iorb) {
          for (int ispin = 0; ispin < _model.spins(); ++ispin) {
            auto orbs = _g_ij_orb_pairs[iorb];
            std::vector < precision > outvec(1, precision(0.0));
            bool found[2];
            precision expectation_value = 0.0;
            /// create particle on two different orbs, sum the resulting vectors and compute contribution to Green's function
            _model.symmetry().set_sector(pair.sector());
            if(create_particles(orbs, ispin, pair.eigenvector(), outvec, expectation_value)) {
              int nlanc = lanczos(outvec);
#ifdef USE_MPI
              if(!rank)
#endif
              {
                std::cout << "orbitals: " << orbs[0] << ", " << orbs[1] << "   spin: " << (ispin == 0 ? "up" : "down") << " <n|aa*|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
                compute_continued_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, 1, _G_g_ij, index_mesh_index(_model.interacting_orbitals() * orbs[0] + orbs[1]), index_mesh_index(ispin));
              }
            }
            /// perform the same for destroying of a particle
            /// restore symmetry sector
            _model.symmetry().set_sector(pair.sector());
            if(annihilate_particles(orbs, ispin, pair.eigenvector(), outvec, expectation_value) ) {
              int nlanc = lanczos(outvec);
#ifdef USE_MPI
              if(!rank)
#endif
              {
                std::cout << "orbitals: " << orbs[0] << ", " << orbs[1] << "   spin: " << (ispin == 0 ? "up" : "down") << " <n|a*a|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
                compute_continued_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, -1, _G_l_ij, index_mesh_index(_model.interacting_orbitals() * orbs[0] + orbs[1]), index_mesh_index(ispin));
              }
            }
          }
        }
      }

      /**
       * Computes non-local Green's function G_ij = 0.5( <(c_i + c_j)(c_i^+ + c_j^+)> - <c_i c_i^+> - <c_j c_j^+> ).
       *
       * @tparam O -- bosonic operator type
       * @param op -- bosonic operator
       */
      void non_local_gf() {
        for (int iomega = 0; iomega < omega().extent(); ++iomega) {
          for (int iorb = 0; iorb < _g_ij_orb_pairs.size(); ++iorb) {
            for (int ispin = 0; ispin < _model.spins(); ++ispin) {
              auto orbs = _g_ij_orb_pairs[iorb];
              for (int jj = 0; jj < 2; ++jj) {
                _G_g_ij(frequency_mesh_index(iomega), index_mesh_index(_model.interacting_orbitals() * orbs[0] + orbs[1]), index_mesh_index(ispin)) -= _G_g(frequency_mesh_index(iomega), index_mesh_index(orbs[jj]), index_mesh_index(ispin));
                _G_l_ij(frequency_mesh_index(iomega), index_mesh_index(_model.interacting_orbitals() * orbs[0] + orbs[1]), index_mesh_index(ispin)) -= _G_l(frequency_mesh_index(iomega), index_mesh_index(orbs[jj]), index_mesh_index(ispin));
              }
              _G_g_ij(frequency_mesh_index(iomega), index_mesh_index(_model.interacting_orbitals() * orbs[0] + orbs[1]), index_mesh_index(ispin)) *= 0.5;
              _G_l_ij(frequency_mesh_index(iomega), index_mesh_index(_model.interacting_orbitals() * orbs[0] + orbs[1]), index_mesh_index(ispin)) *= 0.5;
            }
          }
          // and copy diagonal G for completeness
          for (int iorb = 0; iorb < _g_orbs.size(); ++iorb) {
            size_t orb = _g_orbs[iorb];
            for (int ispin = 0; ispin < _model.spins(); ++ispin) {
             _G_g_ij(frequency_mesh_index(iomega), index_mesh_index(_model.interacting_orbitals() * orb + orb), index_mesh_index(ispin)) = _G_g(frequency_mesh_index(iomega), index_mesh_index(orb), index_mesh_index(ispin));
             _G_l_ij(frequency_mesh_index(iomega), index_mesh_index(_model.interacting_orbitals() * orb + orb), index_mesh_index(ispin)) = _G_l(frequency_mesh_index(iomega), index_mesh_index(orb), index_mesh_index(ispin));
            }
          }
        }
        _G_ij = _G_g_ij + _G_l_ij;
      }

      /// Green's function container object
      GF_TYPE _G_g;
      GF_TYPE _G_g_ij;
      GF_TYPE _G_l;
      GF_TYPE _G_l_ij;
      GF_TYPE _G;
      GF_TYPE _G_ij;
      /// Model we are solving
      typename Hamiltonian::ModelType &_model;
      /// Boltzmann-factor cut-off
      precision _cutoff;
      /// Partition function
      precision _Z;
      /// Orbitals used for the diagonal Green's function calculation
      std::vector<size_t> _g_orbs;
      /// Orbital pairs used for the offdiagonal Green's function calculation
      std::vector<std::array<size_t, 2> > _g_ij_orb_pairs;

      /**
       * @brief Perform the create operator action to the eigenstate
       *
       * @param orbitals - the orbital to create a particle
       * @param spin - the spin of a particle to create
       * @param invec - current eigenstate
       * @param outvec - Op-vec product
       * @param expectation_value - expectation value of aa*
       * @return true if the particle has been created
       */
      template<size_t N>
      bool create_particles(std::array<size_t, N> orbitals, int spin, const std::vector<precision> &invec, std::vector<precision> &outvec, double &expectation_value) {
        // check that the particle can be annihilated
        if (!_model.symmetry().can_create_particle(spin)) {
          return false;
        }
        hamiltonian().storage().reset();
        int nup_new = _model.symmetry().sector().nup() + (1 - spin);
        int ndn_new = _model.symmetry().sector().ndown() + spin;
        typename Hamiltonian::ModelType::Sector next_sec = _model.symmetry().create_partice(spin);
        outvec.assign(hamiltonian().storage().vector_size(next_sec), 0.0);
        common::statistics.registerEvent("adag");
        for(auto orb : orbitals) {
          hamiltonian().storage().init();
          std::vector<precision> tmpout(outvec.size());
          hamiltonian().storage().a_adag(orb + spin * _model.orbitals(), invec, tmpout, next_sec, false);
          std::transform(tmpout.begin(), tmpout.end(), outvec.begin(), outvec.begin(), std::plus<precision>());
        }
        common::statistics.updateEvent("adag");
        double norm = hamiltonian().storage().vv(outvec, outvec
#ifdef USE_MPI
            , hamiltonian().comm()
#endif
        );
        for (int j = 0; j < outvec.size(); ++j) {
          outvec[j] /= std::sqrt(norm);
        }
        _model.symmetry().set_sector(next_sec);
        expectation_value = norm;
        return std::abs(norm) > 1.e-10;
      };

      /**
       * @brief Perform the annihilator operator action to the eigenstate
       *
       * @param orbitals - the orbitals list to destroy a particle
       * @param spin - the spin of a particle to destroy
       * @param invec - current eigenstate
       * @param outvec - Op-vec product
       * @param expectation_value - expectation value of a*a
       * @return true if the particle has been destroyed
       */
      template<size_t N>
      bool annihilate_particles(std::array<size_t, N> orbitals, int spin, const std::vector<precision> &invec, std::vector<precision> &outvec, double &expectation_value) {
        // check that the particle can be annihilated
        if (!_model.symmetry().can_destroy_particle(spin)) {
          return false;
        }
        hamiltonian().storage().reset();
        int nup_new = _model.symmetry().sector().nup() - (1 - spin);
        int ndn_new = _model.symmetry().sector().ndown() - spin;
        typename Hamiltonian::ModelType::Sector next_sec = _model.symmetry().destroy_partice(spin);
        outvec.assign(hamiltonian().storage().vector_size(next_sec), precision(0.0));
        common::statistics.registerEvent("a");
        for(auto orb : orbitals) {
          hamiltonian().storage().init();
          std::vector<precision> tmpout(outvec.size());
          hamiltonian().storage().a_adag(orb + spin * _model.orbitals(), invec, tmpout, next_sec, true);
          std::transform(tmpout.begin(), tmpout.end(), outvec.begin(), outvec.begin(), std::plus<precision>());
        }
        common::statistics.updateEvent("a");
        double norm = hamiltonian().storage().vv(outvec, outvec
#ifdef USE_MPI
            , hamiltonian().comm()
#endif
            );
        for (int j = 0; j < outvec.size(); ++j) {
          outvec[j] /= std::sqrt(norm);
        }
        _model.symmetry().set_sector(next_sec);
        // <v|a^{\star}a|v>
        expectation_value = norm;
        return std::abs(norm) > 1.e-10;
      };
    };
  }
}

#endif //HUBBARD_GREENSFUNCTION_H
