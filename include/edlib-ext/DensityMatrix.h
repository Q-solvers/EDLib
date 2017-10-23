#ifndef EDLIBEXT_DENSITYMATRIX_H
#define EDLIBEXT_DENSITYMATRIX_H

#include "edlib/EigenPair.h"

#include <alps/params.hpp>

namespace EDLib {
  namespace ext {
    template<class Hamiltonian>
    class DensityMatrix {
    protected:
      typedef typename Hamiltonian::ModelType::precision precision;
      typedef typename Hamiltonian::ModelType::Sector sector;
      typedef typename Hamiltonian::ModelType::SYMMETRY symmetry;

    public:

      DensityMatrix(alps::params &p, Hamiltonian& ham__) :
        Ns(int(p["NSITES"])),
        Ip(int(p["NSPINS"]) * int(p["NSITES"])),
        Nspins(int(p["NSPINS"])),
        beta(p["lanc.BETA"].as<precision>()),
        cutoff(p["lanc.BOLTZMANN_CUTOFF"]),
        ham_(ham__)
      {
        std::string input = p["INPUT_FILE"];
        alps::hdf5::archive input_file(input, "r");
        if(input_file.is_data("densitymatrix_orbs/values") && p["storage.EIGENVALUES_ONLY"] == 0){
          input_file >> alps::make_pvp("densitymatrix_orbs/values", orbsA);
          std::sort(orbsA.begin(), orbsA.end());
          orbsA.erase(std::unique(orbsA.begin(), orbsA.end()), orbsA.end());
          symA = std::vector<symmetry>(2, symmetry(orbsA.size()));
          symB = std::vector<symmetry>(1, symmetry(Ns - orbsA.size()));
          Ns_A = orbsA.size();
          Ns_B = Ns - Ns_A;
          for(int iorb = 0; iorb < Ns; ++iorb){
            orbsB.push_back(iorb);
          }
          for(int iorb = orbsA.size() - 1; iorb >= 0; --iorb){
            orbsB.erase(orbsB.begin() + orbsA[iorb]);
          }
          if(p.exists("densitymatrix.SECTOR") && bool(p["densitymatrix.SECTOR"])){
            std::vector<std::vector<int>> _sectors;
            input_file >> alps::make_pvp("densitymatrix_sectors/values", _sectors);
            for (int i = 0; i < _sectors.size(); ++i) {
             secA.push_back(sector(_sectors[i][0], _sectors[i][1], (size_t)(
                 symA[0].comb().c_n_k(Ns_A, _sectors[i][0]) *
                 symA[0].comb().c_n_k(Ns_A, _sectors[i][1])
             )));
            }
          }else{
            for (int i = 0; i <= Ns_A; ++i) {
              for (int j = 0; j <= Ns_A; ++j) {
                secA.push_back(sector(i, j, (size_t)(
                  symA[0].comb().c_n_k(Ns_A, i) *
                  symA[0].comb().c_n_k(Ns_A, j)
                )));
              }
            }
          }
          input_file.close();
          for(size_t isect = 0; isect < secA.size(); ++isect){
            rho.insert(
              std::pair<size_t, std::vector<std::vector<precision>>>(
                isect,
                std::vector<std::vector<precision>>(
                  secA[isect].size(),
                  std::vector<precision>(secA[isect].size(), 0.0)
                )
              )
            );
          }
        }else{
         Ns_A = 0;
         std::cout << "Density matrix can not be calculated. " << std::endl;
        }
      }

      /**
       * Compute reduced density matrix.
       *
       * @return sectors of the reduced density matrix
       */
      const std::map<size_t, std::vector<std::vector<precision>>> &compute() {
        for(size_t isect = 0; isect < secA.size(); ++isect){
          for(size_t jj = 0; jj < secA[isect].size(); ++jj){
            for(size_t kk = 0; kk < secA[isect].size(); ++kk){
             rho[isect][jj][kk] = 0.0;
            }
          }
        }
        precision sum = 0.0;
        const EigenPair<precision, sector> &groundstate = *ham_.eigenpairs().begin();
        // Loop over all eigenpairs.
        for(auto ipair = ham_.eigenpairs().begin(); ipair != ham_.eigenpairs().end(); ++ipair){
          const EigenPair<precision, sector>& pair = *ipair;
          // Calculate Boltzmann factor, skip the states with trivial contribution.
          precision boltzmann_f = std::exp(
           -(pair.eigenvalue() - groundstate.eigenvalue()) * beta
          );
          if(boltzmann_f < cutoff){
            continue;
          }
          // Sum the contributions.
          compute_eigenvector(pair, boltzmann_f);
          sum += boltzmann_f;
        }
        for(size_t isect = 0; isect < secA.size(); ++isect){
          for(size_t jj = 0; jj < secA[isect].size(); ++jj){
            for(size_t kk = 0; kk < secA[isect].size(); ++kk){
             rho[isect][jj][kk] /= sum;
            }
          }
        }
        return rho;
      }

      /**
       * Print the reduced density matrix.
       */
      void print() {
  #ifdef USE_MPI
        int myid;
        MPI_Comm_rank(ham_.comm(), &myid);
        if(!myid)
  #endif
        for(size_t isect = 0; isect < secA.size(); ++isect){
          std::cout << "Density matrix sector " << secA[isect].nup() << " " << secA[isect].ndown() << std::endl;
          for(size_t jj = 0; jj < secA[isect].size(); ++jj){
            for(size_t kk = 0; kk < secA[isect].size(); ++kk){
              if(kk){
                std::cout << "\t";
              }
              std::cout << rho[isect][jj][kk];
            }
            std::cout << std::endl;
          }
        }
      }

      /**
       * Combine the sectors of reduced density matrix.
       *
       * @return the whole reduced density matrix
       */
      const std::vector<std::vector<precision>> full() const {
        size_t fullsize = std::pow(2, Ns_A * Nspins);
        std::vector<std::vector<precision>> result(
          fullsize, std::vector<precision>(fullsize, 0.0)
        );
        size_t shift = 0;
        for(size_t isect = 0; isect < secA.size(); ++isect){
          for(size_t ii = 0; ii < secA[isect].size(); ++ii){
            for(size_t jj = 0; jj < secA[isect].size(); ++jj){
              result[shift + ii][shift + jj] = rho.at(isect)[ii][jj];
            }
          }
          shift += secA[isect].size();
        }
        return result;
      }

    inline const std::vector<sector> &sectors()
    const {
     return secA;
    }

    inline const std::map<size_t, std::vector<std::vector<precision>>> &matrix()
    const {
     return rho;
    }

    inline const Hamiltonian &ham()
    const {
     return ham_;
    }

    private:

      /**
       * Compute the contribution of an eigenvector to the density matrix.
       *
       * @param pair   the eigenpair
       * @param weight additional multiplier (Boltzmann factor)
       */
      void compute_eigenvector(const EigenPair<precision, sector>& pair, precision weight) {
  #ifdef USE_MPI
        int myid;
        int nprocs;
        MPI_Comm_rank(ham_.comm(), &myid);
        MPI_Comm_size(ham_.comm(), &nprocs);
        std::vector<int> counts(nprocs);
        std::vector<int> displs(nprocs + 1);
        std::vector<precision> evec(pair.sector().size());
        int size = pair.eigenvector().size();
        MPI_Allgather(&size, 1,
                      alps::mpi::detail::mpi_type<int>(),
                      counts.data(), 1,
                      alps::mpi::detail::mpi_type<int>(),
                      ham_.comm()
        );
        displs[0] = 0;
        for(size_t i = 0; i < nprocs; ++i){
         displs[i + 1] = displs[i] + counts[i];
        }
        MPI_Allgatherv(pair.eigenvector().data(), pair.eigenvector().size(),
                       alps::mpi::detail::mpi_type<precision>(),
                       evec.data(), counts.data(), displs.data(),
                       alps::mpi::detail::mpi_type<precision>(),
                       ham_.comm()
        );
  #endif
        ham_.model().symmetry().set_sector(pair.sector());
        for(size_t isect = 0; isect < secA.size(); ++isect){
          symA[0].set_sector(secA[isect]);
          int nupB = pair.sector().nup() - secA[isect].nup();
          int ndownB = pair.sector().ndown() - secA[isect].ndown();
          if(
           (nupB < 0) ||
           (ndownB < 0) ||
           (nupB > Ns_B) ||
           (ndownB > Ns_B)
          ){
           continue;
          }
          sector secB = sector(nupB, ndownB,
            symB[0].comb().c_n_k(Ns_B, nupB) * symB[0].comb().c_n_k(Ns_B, ndownB)
          );
          symB[0].set_sector(secB);
          for(size_t ii = 0; ii < secB.size(); ++ii){
            symB[0].next_state();
            long long stateB = symB[0].state();
            symA[0].set_sector(secA[isect]);
            for(size_t jj = 0; jj < secA[isect].size(); ++jj){
              symA[0].next_state();
              long long state0 = mergestate(symA[0].state(), stateB);
              symA[1].set_sector(secA[isect]);
              for(size_t kk = 0; kk < secA[isect].size(); ++kk){
                symA[1].next_state();
                long long state1 = mergestate(symA[1].state(), stateB);
                rho[isect][jj][kk] += weight *
  #ifdef USE_MPI
                  evec[ham_.model().symmetry().index(state0)] *
                  evec[ham_.model().symmetry().index(state1)];
  #else
                  pair.eigenvector()[ham_.model().symmetry().index(state0)] *
                  pair.eigenvector()[ham_.model().symmetry().index(state1)];
  #endif
              }
            }
          }
        }
      }

      /**
       * Combine basis vectors.
       *
       * @param  stateA basis vector of subsystem A
       * @param  stateB basis vector of subsystem B
       * @return        basis vector of the whole system
       */
      long long mergestate(long long stateA, long long stateB){
        long long state = 0;
        long long newstate = 0;
        int isign;
        for(size_t ispin = 0; ispin < Nspins; ++ispin){
          for(size_t iorb = 0; iorb < orbsA.size(); ++iorb){
            if(ham_.model().checkState(stateA, iorb + ispin * Ns_A, Nspins * Ns_A)){
              ham_.model().adag(orbsA[iorb]  + ispin * Ns, state, newstate, isign);
              state = newstate;
            }
          }
          for(size_t iorb = 0; iorb < orbsB.size(); ++iorb){
            if(ham_.model().checkState(stateB, iorb + ispin * Ns_B, Nspins * Ns_B)){
              ham_.model().adag(orbsB[iorb]  + ispin * Ns, state, newstate, isign);
              state = newstate;
            }
          }
        }
        return state;
      }

      /// The hamiltonian of the system A+B
      Hamiltonian& ham_;
      /// Reduced density matrix for the system A
      std::map<size_t, std::vector<std::vector<precision>>> rho;

      /// Symmetries of the system A.
      std::vector<symmetry> symA;
      /// Symmetry of the system B.
      // FIXME Have to use a vector even for one symmetry, else the constructor won't compile.
      std::vector<symmetry> symB;
      /// Symmetry sectors of the system A.
      std::vector<sector> secA;
      /// The orbitals for which the density matrix is calculated.
      std::vector<int> orbsA;
      std::vector<int> orbsB;
      /// The number of sites.
      int Ns, Ns_A, Ns_B;
      int Ip;
      int Nspins;
      /// Inverse temperature.
      precision beta;
      /// Minimal Boltzmann-factor.
      precision cutoff;

    };

  }
}

#endif
