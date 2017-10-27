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
    template<class Hamiltonian, typename Mesh, typename... Args>
    class GreensFunction : public Lanczos < Hamiltonian, Mesh, Args...> {
      using Lanczos < Hamiltonian, Mesh, Args... >::hamiltonian;
      using Lanczos < Hamiltonian, Mesh, Args... >::lanczos;
      using Lanczos < Hamiltonian, Mesh, Args... >::omega;
      using Lanczos < Hamiltonian, Mesh, Args... >::beta;
      using Lanczos < Hamiltonian, Mesh, Args... >::compute_continued_fraction;
      using typename Lanczos < Hamiltonian, Mesh, Args... >::precision;
    public:
      GreensFunction(alps::params &p, Hamiltonian &h, Args ... args) : Lanczos < Hamiltonian, Mesh, Args... >(p, h, args...), _model(h.model()),
                                                        gf(Lanczos < Hamiltonian, Mesh, Args... >::omega(), alps::gf::index_mesh(h.model().interacting_orbitals()), alps::gf::index_mesh(p["NSPINS"].as<int>())),
                                                        gf_ij(Lanczos < Hamiltonian, Mesh, Args... >::omega(), alps::gf::index_mesh(h.model().interacting_orbitals() * h.model().interacting_orbitals()), alps::gf::index_mesh(p["NSPINS"].as<int>())),
                                                        _cutoff(p["lanc.BOLTZMANN_CUTOFF"]) {
        if(p["storage.EIGENVALUES_ONLY"] == 1) {
          throw std::logic_error("Eigenvectors have not been computed. Green's function can not be evaluated.");
        }
        std::string input = p["INPUT_FILE"];
        alps::hdf5::archive input_file(input.c_str(), "r");
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

      void compute() {
        gf *= 0.0;
        gf_ij *= 0.0;
        _Z = 0.0;
        if(hamiltonian().eigenpairs().empty())
          return;
#ifdef USE_MPI
        int rank;
        MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
#endif
        /// get groundstate
        const EigenPair<precision, typename Hamiltonian::ModelType::Sector> &groundstate =  *hamiltonian().eigenpairs().begin();
        /// compute statsum
        for (auto kkk = hamiltonian().eigenpairs().begin(); kkk != hamiltonian().eigenpairs().end(); kkk++) {
          const EigenPair<precision, typename Hamiltonian::ModelType::Sector> &eigenpair = *kkk;
          _Z += std::exp(-(eigenpair.eigenvalue() - groundstate.eigenvalue()) * beta());
        }
        common::statistics.registerEvent("Greens function");
        /// iterate over eigen-pairs
        for (auto kkk = hamiltonian().eigenpairs().begin(); kkk != hamiltonian().eigenpairs().end(); kkk++) {
          const EigenPair<precision, typename Hamiltonian::ModelType::Sector>& pair = *kkk;
          /// compute Boltzmann-factor
          precision boltzmann_f = std::exp(-(pair.eigenvalue() - groundstate.eigenvalue()) * beta());
          /// Skip all eigenvalues with Boltzmann-factor smaller than cutoff
          if (boltzmann_f < _cutoff) {
//        std::cout<<"Skipped by Boltzmann factor."<<std::endl;
            continue;
          }
#ifdef USE_MPI
          if(rank==0)
#endif
          std::cout << "Compute Green's function contribution for eigenvalue E=" << pair.eigenvalue() << " with Boltzmann factor = " << boltzmann_f << "; for sector" << pair.sector() << std::endl;
          local_contribution(pair, groundstate);
          nonlocal_contribution(pair, groundstate);
        }
#ifdef USE_MPI
        if(rank == 0) {
#endif
        /// normalize Green's function by statsum Z.
        gf /= _Z;
        gf_ij /= _Z;
#ifdef USE_MPI
        }
#endif
        non_local_gf();
        common::statistics.updateEvent("Greens function");
      }

      /**
       * Save Green's function in the hdf5 archive and in plain text file
       * @param ar -- hdf5 archive to save Green's function
       * @param path -- root path in hdf5 archive
       */
      void save(alps::hdf5::archive& ar, const std::string & path) {
#ifdef USE_MPI
        int rank;
        MPI_Comm_rank(hamiltonian().storage().comm(), &rank);
        if(rank == 0) {
#endif
          if(diagonal_orbs.size()){
            gf.save(ar, path + "/G_omega");
            std::ostringstream Gomega_name;
            Gomega_name << "G_omega";
            std::ofstream G_omega_file(Gomega_name.str().c_str());
            G_omega_file << std::setprecision(14) << gf;
            G_omega_file.close();
          }
          std::cout << "Statsum: " << _Z << std::endl;
          ar[path + "/@Statsum"] << _Z;
          if(offdiagonal_orbs.size()){
            std::ostringstream Gomega_name2;
            Gomega_name2 << "G_ij_omega";
            std::ofstream G_omega_file2(Gomega_name2.str().c_str());
            G_omega_file2<< std::setprecision(14) << gf_ij;
            G_omega_file2.close();
          }
#ifdef USE_MPI
        }
#endif
      }

      void compute_selfenergy(alps::hdf5::archive &ar, const std::string &path){
        alps::gf::three_index_gf<std::complex<double>, Mesh, alps::gf::index_mesh, alps::gf::index_mesh> bare(gf.mesh1(), gf.mesh2(), gf.mesh3());
        alps::gf::three_index_gf<std::complex<double>, Mesh, alps::gf::index_mesh, alps::gf::index_mesh> sigma(gf.mesh1(), gf.mesh2(), gf.mesh3());
        _model.bare_greens_function(bare, beta());
        bare.save(ar, path + "/G0_omega");
        std::ostringstream Gomega_name;
        Gomega_name << "G0_omega";
        std::ofstream G_omega_file(Gomega_name.str().c_str());
        G_omega_file << std::setprecision(14) << bare;
        G_omega_file.close();
        for(int iw = 0; iw< bare.mesh1().points().size(); ++iw) {
          typename Mesh::index_type w(iw);
          for (int im: bare.mesh2().points()) {
            for (int is : bare.mesh3().points()) {
              sigma(w, alps::gf::index_mesh::index_type(im), alps::gf::index_mesh::index_type(is)) =
                1.0/bare(w, alps::gf::index_mesh::index_type(im), alps::gf::index_mesh::index_type(is)) - 1.0/gf(w, alps::gf::index_mesh::index_type(im), alps::gf::index_mesh::index_type(is));
            }
          }
        }
        sigma.save(ar, path + "/Sigma_omega");
        Gomega_name.str("");
        Gomega_name << "Sigma_omega";
        G_omega_file.open(Gomega_name.str().c_str());
        G_omega_file << std::setprecision(14) << sigma;
        G_omega_file.close();
      }

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
        for (int iorb = 0; iorb < diagonal_orbs.size(); ++iorb) {
          /// iterate over spins
          for (int ispin = 0; ispin < _model.spins(); ++ispin) {
            int orb = diagonal_orbs[iorb];
            std::vector < precision > outvec(1, precision(0.0));
            precision expectation_value = 0.0;
            _model.symmetry().set_sector(pair.sector());
            /// first we are going to create particle and compute contribution to Green's function
            if (create_particle(orb, ispin, pair.eigenvector(), outvec, expectation_value)) {
              /// Perform Lanczos factorization for starting vector |outvec>
              int nlanc = lanczos(outvec);
#ifdef USE_MPI
              if(!rank)
#endif
              {
                std::cout << "orbital: " << orb << "   spin: " << (ispin == 0 ? "up" : "down") << " <n|aa*|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
                /// Using computed Lanczos factorization compute approximation for \frac{1}{z - H} by calculation of a continued fraction
                compute_continued_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, 1, gf, index_mesh_index(orb), index_mesh_index(ispin));
              }
            }
            /// restore symmetry sector
            _model.symmetry().set_sector(pair.sector());
            /// perform the same for destroying of a particle
            if (annihilate_particle(orb, ispin, pair.eigenvector(), outvec, expectation_value)) {
              int nlanc = lanczos(outvec);
#ifdef USE_MPI
              if(!rank)
#endif
              {
                std::cout << "orbital: " << orb << "   spin: " << (ispin == 0 ? "up" : "down") << " <n|a*a|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
                compute_continued_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, -1, gf, index_mesh_index(orb), index_mesh_index(ispin));
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
        for (int iorb = 0; iorb < offdiagonal_orbs.size(); ++iorb) {
          for (int ispin = 0; ispin < _model.spins(); ++ispin) {
            std::vector<int> orbs = offdiagonal_orbs[iorb];
            std::vector<std::vector<precision>> outvec(2, std::vector<precision>(1, precision(0.0)));
            std::vector < precision > sumvec(1, precision(0.0));
            bool found[2];
            precision expectation_value = 0.0;
            /// create particle on two different orbs, sum the resulting vectors and compute contribution to Green's function
            for(int i = 0; i < 2; ++i){
              _model.symmetry().set_sector(pair.sector());
              found[i] = create_particle(orbs[i], ispin, pair.eigenvector(), outvec[i], expectation_value);
            }
            if(found[0] || found[1]){
              if(found[0] && found[1]){
               sumvec = outvec[0];
               for(size_t i = 0; i < sumvec.size(); ++i){
                 sumvec[i] += outvec[1][i];
               }
              }else{
                sumvec = (found[0] ? outvec[0] : outvec[1]);
              }
              int nlanc = lanczos(sumvec);
#ifdef USE_MPI
              if(!rank)
#endif
              {
                std::cout << "orbitals: " << orbs[0] << ", " << orbs[1] << "   spin: " << (ispin == 0 ? "up" : "down") << " <n|aa*|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
                compute_continued_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, 1, gf_ij, index_mesh_index(_model.interacting_orbitals() * orbs[0] + orbs[1]), index_mesh_index(ispin));
              }
            }
            /// perform the same for destroying of a particle
            for(int i = 0; i < 2; ++i){
              /// restore symmetry sector
              _model.symmetry().set_sector(pair.sector());
              found[i] = annihilate_particle(orbs[i], ispin, pair.eigenvector(), outvec[i], expectation_value);
            }
            if(found[0] || found[1]){
              if(found[0] && found[1]){
               sumvec = outvec[0];
               for(size_t i = 0; i < sumvec.size(); ++i){
                 sumvec[i] += outvec[1][i];
               }
              }else{
                sumvec = (found[0] ? outvec[0] : outvec[1]);
              }
              int nlanc = lanczos(sumvec);
#ifdef USE_MPI
              if(!rank)
#endif
              {
                std::cout << "orbitals: " << orbs[0] << ", " << orbs[1] << "   spin: " << (ispin == 0 ? "up" : "down") << " <n|a*a|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
                compute_continued_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, -1, gf_ij, index_mesh_index(_model.interacting_orbitals() * orbs[0] + orbs[1]), index_mesh_index(ispin));
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
          for (int iorb = 0; iorb < offdiagonal_orbs.size(); ++iorb) {
            for (int ispin = 0; ispin < _model.spins(); ++ispin) {
              std::vector<int> orbs = offdiagonal_orbs[iorb];
              for (int j = 0; j < 2; ++j) {
                gf_ij(frequency_mesh_index(iomega), index_mesh_index(_model.interacting_orbitals() * orbs[0] + orbs[1]), index_mesh_index(ispin)) -= gf(frequency_mesh_index(iomega), index_mesh_index(orbs[j]), index_mesh_index(ispin));
              }
              gf_ij(frequency_mesh_index(iomega), index_mesh_index(_model.interacting_orbitals() * orbs[0] + orbs[1]), index_mesh_index(ispin)) *= 0.5;
            }
          }
        }
      }

      /// Green's function type
      typedef alps::gf::three_index_gf<std::complex<double>, Mesh, alps::gf::index_mesh, alps::gf::index_mesh >  GF_TYPE;
      typedef typename alps::gf::index_mesh::index_type index_mesh_index;
      typedef typename Mesh::index_type frequency_mesh_index;
      /// Green's function container object
      GF_TYPE gf;
      GF_TYPE gf_ij;
      /// Model we are solving
      typename Hamiltonian::ModelType &_model;
      /// Boltzmann-factor cutoff
      precision _cutoff;
      /// Statsum
      precision _Z;
      /// Orbital pairs to calculate the Green's function.
      std::vector<std::vector<int>> gf_orbs;
      /// Orbitals to calculate the diagonal Green's function.
      std::vector<int> diagonal_orbs;
      /// Orbital pairs to calculate the offdiagonal Green's function.
      std::vector<std::vector<int>> offdiagonal_orbs;

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
        common::statistics.registerEvent("adag");
        hamiltonian().storage().a_adag(orbital + spin * _model.orbitals(), invec, outvec, next_sec, false);
        common::statistics.updateEvent("adag");
        std::cout<<"adag in "<<common::statistics.event("adag").first<<" s \n";
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
        common::statistics.registerEvent("a");
        hamiltonian().storage().a_adag(orbital + spin * _model.orbitals(), invec, outvec, next_sec, true);
        common::statistics.updateEvent("a");
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
