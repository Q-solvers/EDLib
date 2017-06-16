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
      for(int iorb = 0; iorb < orbsA.size(); ++iorb){
        orbsB.erase(orbsB.begin() + orbsA[iorb]);
      }
      if(p.exists("densitymatrix.SECTOR") && bool(p["densitymatrix.SECTOR"])){
        std::vector<std::vector<int>> _sectors;
        input_file >> alps::make_pvp("densitymatrix_sectors/values", _sectors);
        for (int i = 0; i < _sectors.size(); ++i) {
         secA.push_back(sector(_sectors[i][0], _sectors[i][1], (size_t)(symA[0].comb().c_n_k(Ns_A, _sectors[i][0]) * symA[0].comb().c_n_k(Ns_A, _sectors[i][1]))));
        }
      }
      input_file.close();
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
//          std::cout<<"Skipped by Boltzmann factor."<<std::endl;
          continue;
        }
        // Sum the contributions.
        compute_eigenvector(rho, ham, pair, boltzmann_f);
        sum += boltzmann_f;
      }
      return rho;
    }

  private:

    void compute_eigenvector(std::map<size_t, std::vector<std::vector<precision>>>& rho, Hamiltonian& ham, const EigenPair<precision, sector>& pair, precision multiplier) {
      for(size_t isect = 0; isect < secA.size(); ++isect){
        symA[0].set_sector(secA[isect]);
        sector secB = sector(pair.sector().nup() - secA[isect].nup(), pair.sector().ndown() - secA[isect].ndown(), pair.sector().size() / secA[isect].size());
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
              rho[isect][ii][jj] += multiplier * pair.eigenvector()[state[0]] * pair.eigenvector()[state[1]];
            }
          }
        }
      }
    }


    long long mergestate(Hamiltonian& ham, long long stateA, long long stateB){
      long long state = 0;
      for(size_t ispin = 0; ispin < Nspins; ++ispin){
        for(size_t iorb = 0; iorb < orbsA.size(); ++iorb){
          state += ham.model().checkState(stateA, iorb + ispin * Ns_A, Ip) * std::pow(2, orbsA[iorb]  + ispin * Ns_A);
        }
        for(size_t iorb = 0; iorb < orbsB.size(); ++iorb){
          state += ham.model().checkState(stateB, iorb + ispin * Ns_B, Ip) * std::pow(2, orbsB[iorb]  + ispin * Ns_B);
        }
      }
      return state;
    }

/* FIXME Copypasta from FermionicModel.h */
    int inline checkState(long long nst, const int &im, int Ip) const {
      return (int) ((nst & (1ll << (Ip - 1 - im))) >> (Ip - 1 - im));
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
