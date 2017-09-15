#ifndef HUBBARD_DENSITYMATRIX_H
#define HUBBARD_DENSITYMATRIX_H

#include "EigenPair.h"

#include <alps/params.hpp>

namespace EDLib {
  template<class Hamiltonian>
  class DensityMatrix {
  protected:
    typedef typename Hamiltonian::ModelType::precision precision;
    typedef typename Hamiltonian::ModelType::Sector sector;
    typedef typename Hamiltonian::ModelType::SYMMETRY symmetry;

  public:

    DensityMatrix(alps::params &p) :
      Ns(int(p["NSITES"])),
      Ip(int(p["NSPINS"]) * int(p["NSITES"])),
      Nspins(int(p["NSPINS"])),
      beta(p["lanc.BETA"].as<precision>()),
      cutoff(p["lanc.BOLTZMANN_CUTOFF"])
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
      }else{
       Ns_A = 0;
       std::cout << "Density matrix can not be calculated. " << std::endl;
      }
    }

    std::map<size_t, std::vector<std::vector<precision>>> compute(Hamiltonian& ham) {
      std::map<size_t, std::vector<std::vector<precision>>> rho;
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
      precision sum = 0.0;
      const EigenPair<precision, sector> &groundstate = *ham.eigenpairs().begin();
      // Loop over all eigenpairs.
      for(auto ipair = ham.eigenpairs().begin(); ipair != ham.eigenpairs().end(); ++ipair){
        const EigenPair<precision, sector>& pair = *ipair;
        // Calculate Boltzmann factor, skip the states with trivial contribution.
        precision boltzmann_f = std::exp(
         -(pair.eigenvalue() - groundstate.eigenvalue()) * beta
        );
        if(boltzmann_f < cutoff){
          continue;
        }
        // Sum the contributions.
        compute_eigenvector(rho, ham, pair, boltzmann_f);
        sum += boltzmann_f;
      }
      for(size_t isect = 0; isect < secA.size(); ++isect){
        for(size_t jj = 0; jj < secA[isect].size(); ++jj){
          for(size_t kk = 0; kk < secA[isect].size(); ++kk){
           rho[isect][jj][kk] /= sum;
          }
        }
      }
#ifdef USE_MPI
      int myid;
      MPI_Comm_rank(ham.comm(), &myid);
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
      return rho;
    }

  private:

    void compute_eigenvector(std::map<size_t, std::vector<std::vector<precision>>>& rho, Hamiltonian& ham, const EigenPair<precision, sector>& pair, precision multiplier) {
#ifdef USE_MPI
      int myid;
      int nprocs;
      MPI_Comm_rank(ham.comm(), &myid);
      MPI_Comm_size(ham.comm(), &nprocs);
      std::vector<int> counts(nprocs);
      std::vector<int> displs(nprocs + 1);
      std::vector<precision> evec(pair.sector().size());
      int size = pair.eigenvector().size();
      MPI_Allgather(&size, 1,
                    alps::mpi::detail::mpi_type<int>(),
                    counts.data(), 1,
                    alps::mpi::detail::mpi_type<int>(),
                    ham.comm()
      );
      displs[0] = 0;
      for(size_t i = 0; i < nprocs; ++i){
       displs[i + 1] = displs[i] + counts[i];
      }
      MPI_Allgatherv(pair.eigenvector().data(), pair.eigenvector().size(),
                     alps::mpi::detail::mpi_type<precision>(),
                     evec.data(), counts.data(), displs.data(),
                     alps::mpi::detail::mpi_type<precision>(),
                     ham.comm()
      );
#endif
      ham.model().symmetry().set_sector(pair.sector());
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
          long long stateA[2];
          long long state[2];
          for(size_t jj = 0; jj < secA[isect].size(); ++jj){
            symA[0].next_state();
            stateA[0] = symA[0].state();
            state[0] = mergestate(ham, stateA[0], stateB);
            symA[1].set_sector(secA[isect]);
            for(size_t kk = 0; kk < secA[isect].size(); ++kk){
              symA[1].next_state();
              stateA[1] = symA[1].state();
              state[1] = mergestate(ham, stateA[1], stateB);
              rho[isect][jj][kk] += multiplier *
#ifdef USE_MPI
                evec[ham.model().symmetry().index(state[0])] *
                evec[ham.model().symmetry().index(state[1])];
#else
                pair.eigenvector()[ham.model().symmetry().index(state[0])] *
                pair.eigenvector()[ham.model().symmetry().index(state[1])];
#endif
            }
          }
        }
      }
    }


    long long mergestate(Hamiltonian& ham, long long stateA, long long stateB){
      long long state = 0;
      long long newstate = 0;
      int isign;
      for(size_t ispin = 0; ispin < Nspins; ++ispin){
        for(size_t iorb = 0; iorb < orbsA.size(); ++iorb){
          if(ham.model().checkState(stateA, iorb + ispin * Ns_A, Nspins * Ns_A)){
            ham.model().adag(orbsA[iorb]  + ispin * Ns, state, newstate, isign);
            state = newstate;
          }
        }
        for(size_t iorb = 0; iorb < orbsB.size(); ++iorb){
          if(ham.model().checkState(stateB, iorb + ispin * Ns_B, Nspins * Ns_B)){
            ham.model().adag(orbsB[iorb]  + ispin * Ns, state, newstate, isign);
            state = newstate;
          }
        }
      }
      return state;
    }

    std::vector<symmetry> symA;
    std::vector<symmetry> symB; // FIXME What's the fucking difference from just symmetry symB -- it won't compile!
    std::vector<sector> secA;
    /// The orbitals for which the density matrix is calculated
    std::vector<int> orbsA;
    std::vector<int> orbsB;
    // The number of sites
    int Ns, Ns_A, Ns_B;
    int Ip;
    int Nspins;
    /// Inverse temperature
    precision beta;
    /// Boltzmann-factor cutoff
    precision cutoff;

  };

}

#endif
