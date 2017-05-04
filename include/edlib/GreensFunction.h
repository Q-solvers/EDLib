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
                                                        _cutoff(p["lanc.BOLTZMANN_CUTOFF"]) {
        if(p["storage.EIGENVALUES_ONLY"] == 1) {
          throw std::logic_error("Eigenvectors have not been computed. Green's function can not be evaluated.");
        }
      }

      void compute() {
        gf *= 0.0;
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
          /// iterate over interacting orbitals
          for (int i = 0; i < _model.interacting_orbitals(); ++i) {
            /// iterate over spins
            for (int is = 0; is < _model.spins(); ++is) {
              std::vector < precision > outvec(1, precision(0.0));
              precision expectation_value = 0;
              _model.symmetry().set_sector(pair.sector());
              /// first we are going to create particle and compute contribution to Green's function
              if (create_particle(i, is, pair.eigenvector(), outvec, expectation_value)) {
                /// Perform Lanczos factorization for starting vector |outvec>
                int nlanc = lanczos(outvec);
#ifdef USE_MPI
               if(rank==0){
#endif
                std::cout << "orbital: " << i << "   spin: " << (is == 0 ? "up" : "down") << " <n|aa*|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
                /// Using computed Lanczos factorization compute approximation for \frac{1}{z - H} by calculation of a continued fraction
                 compute_continued_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, 1, gf, alps::gf::index_mesh::index_type(i),
                                            alps::gf::index_mesh::index_type(is));
#ifdef USE_MPI
               }
#endif
              }
              /// restore symmetry sector
              _model.symmetry().set_sector(pair.sector());
              /// perform the same for destroying of a particle
              if (annihilate_particle(i, is, pair.eigenvector(), outvec, expectation_value)) {
                int nlanc = lanczos(outvec);
#ifdef USE_MPI
                if(rank==0){
#endif
                std::cout << "orbital: " << i << "   spin: " << (is == 0 ? "up" : "down") << " <n|a*a|n>=" << expectation_value << " nlanc:" << nlanc << std::endl;
                  compute_continued_fraction(expectation_value, pair.eigenvalue(), groundstate.eigenvalue(), nlanc, -1, gf, alps::gf::index_mesh::index_type(i),
                                             alps::gf::index_mesh::index_type(is));
#ifdef USE_MPI
                }
#endif
              }
            }
          }
        }
#ifdef USE_MPI
        if(rank == 0) {
#endif
        /// normalize Green's function by statsum Z.
        gf /= _Z;
#ifdef USE_MPI
        }
#endif
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
          gf.save(ar, path + "/G_omega");
          std::ostringstream Gomega_name;
          Gomega_name << "G_omega";
          std::ofstream G_omega_file(Gomega_name.str().c_str());
          G_omega_file << std::setprecision(14) << gf;
          G_omega_file.close();
          std::cout << "Statsum: " << _Z << std::endl;
          ar[path + "/@Statsum"] << _Z;
#ifdef USE_MPI
        }
#endif
      }

      void compute_selfenergy(alps::hdf5::archive &ar, const std::string &path){
        GF_TYPE bare(gf.mesh1(), gf.mesh2(), gf.mesh3());
        GF_TYPE sigma(gf.mesh1(), gf.mesh2(), gf.mesh3());
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
      /// Green's function type
      typedef alps::gf::three_index_gf<std::complex<double>, Mesh, alps::gf::index_mesh, alps::gf::index_mesh >  GF_TYPE;
      /// Green's function container object
      GF_TYPE gf;
      /// Model we are solving
      typename Hamiltonian::ModelType &_model;
      /// Boltzmann-factor cutoff
      precision _cutoff;
      /// Statsum
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
