#ifndef HUBBARD_STATEDESCRIPTION_H
#define HUBBARD_STATEDESCRIPTION_H

#include "EigenPair.h"
#include <vector>
#include <bitset>
#include <string>

namespace EDLib {
  template<class Hamiltonian>
  class StateDescription {
  protected:
    typedef typename Hamiltonian::ModelType::precision precision;
  public:

    StateDescription(alps::params &p) :
      _beta(p["lanc.BETA"].as<precision>()),
      _cutoff(p["lanc.BOLTZMANN_CUTOFF"])
    {
      if(p["storage.EIGENVALUES_ONLY"] == 1) {
        throw std::logic_error("Eigenvectors have not been computed. StateDescription can not continue.");
      }
#ifdef USE_MPI
      const int nitems=2;
      int blocklengths[nitems] = {1, 1};
      MPI_Datatype types[nitems] = {alps::mpi::detail::mpi_type<size_t>(), alps::mpi::detail::mpi_type<precision>()};
      MPI_Aint offsets[nitems];
      offsets[0] = offsetof(Element, ind);
      offsets[1] = offsetof(Element, val);
      MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_Element);
      MPI_Type_commit(&mpi_Element);
#endif
    };

    /**
     * @brief Compute average number of electrons
     *
     * @param orb - the site for calculation
     * @param diff - 0: average number or electrons, 1: average spin
     * @param _ham - the Hamiltonian
     */
    precision avgn(Hamiltonian& _ham, int diff, int orb){
      precision avg = 0.0;
      precision sum = 0.0;
      const EigenPair<precision, typename Hamiltonian::ModelType::Sector> &groundstate =  *_ham.eigenpairs().begin();
      // Loop over all eigenpairs.
      for(auto ipair = _ham.eigenpairs().begin(); ipair != _ham.eigenpairs().end(); ++ipair){
        const EigenPair<precision, typename Hamiltonian::ModelType::Sector>& pair = *ipair;
        // Calculate Boltzmann factor, skip the states with trivial contribution.
        precision boltzmann_f = std::exp(
         -(pair.eigenvalue() - groundstate.eigenvalue()) * _beta
        );
        if(boltzmann_f < _cutoff){
//          std::cout<<"Skipped by Boltzmann factor."<<std::endl;
          continue;
        }
        // Sum the contributions.
        avg += avgn1(_ham, pair, diff, orb) * boltzmann_f;
        sum += boltzmann_f;
      }
      return avg / sum;
    }

    std::vector<std::pair<long long, precision>> find(Hamiltonian& _ham, const EigenPair<typename Hamiltonian::ModelType::precision, typename Hamiltonian::ModelType::Sector>& pair, size_t nmax, precision trivial){
      _ham.model().symmetry().set_sector(pair.sector());
      _ham.storage().reset();
      int count = std::min(nmax, pair.eigenvector().size());
      std::vector<size_t> largest = std::vector<size_t>(pair.eigenvector().size());
#ifdef USE_MPI
      int myid;
      int nprocs;
      MPI_Comm_rank(_ham.comm(), &myid);
      MPI_Comm_size(_ham.comm(), &nprocs);
      std::vector<int> counts(nprocs);
      std::vector<int> displs(nprocs + 1);
      MPI_Gather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, _ham.comm());
      if(!myid){
        displs[0] = 0;
        for(size_t i = 0; i < nprocs; i++){
         displs[i + 1] = displs[i] + counts[i];
        }
      }
#endif
      for(size_t i = 0; i < largest.size(); ++i){
        largest[i] = i;
      }
      std::partial_sort(largest.begin(), largest.begin()+count, largest.end(), [&pair] (int a, int b) -> bool {
        return ((std::abs(pair.eigenvector()[a]) > std::abs(pair.eigenvector()[b]))
                || ((std::abs(pair.eigenvector()[a]) == std::abs(pair.eigenvector()[b])) && (a < b))
               );
      });
#ifdef USE_MPI
      std::vector<Element> send(count);
      for(size_t i = 0; i < count; ++i){
       send[i] = Element(largest[i] + _ham.storage().offset(), pair.eigenvector()[largest[i]]);
      }
      std::vector<Element> all(displs[nprocs]);
      MPI_Gatherv(send.data(), count, mpi_Element, all.data(), counts.data(), displs.data(), mpi_Element, 0, _ham.comm());
      if (!myid) {
        nmax = std::min(nmax, all.size());
        std::partial_sort(all.begin(), all.begin()+nmax, all.end(), [] (Element a, Element b) -> bool {return (a > b);});
        std::vector<std::pair<long long, precision>> ret(nmax);
        for(size_t i = 0; i < nmax; ++i){
          long long nst = _ham.model().symmetry().state_by_index(all[i].ind);
          ret[i] = std::pair<long long, precision>(nst, all[i].val);
        }
        return ret;
      } else {
        std::vector<std::pair<long long, precision>> ret(0);
        return ret;
      }
#else
      std::vector<std::pair<long long, precision>> ret(count);
      for(size_t i = 0; i < count; ++i){
        long long nst = _ham.model().symmetry().state_by_index(largest[i]);
        ret[i] = std::pair<long long, precision>(nst, pair.eigenvector()[largest[i]]);
      }
      return ret;
#endif
    }

    void print(Hamiltonian& _ham, const EigenPair<typename Hamiltonian::ModelType::precision, typename Hamiltonian::ModelType::Sector>& pair, size_t nmax, precision trivial){
      std::vector<std::pair<long long, precision>> contribs = find(_ham, pair, nmax, trivial);
#ifdef USE_MPI
      int myid;
      MPI_Comm_rank(_ham.comm(), &myid);
      if(!myid)
#endif
      {
        std::cout << "Eigenvector components for eigenvalue " << pair.eigenvalue() << " ";
        pair.sector().print();
        std::cout << std::endl;
        for(size_t i = 0; i < contribs.size(); ++i){
          std::cout << contribs[i].second << " * |";
          std::string spin_down = std::bitset< 64 >( contribs[i].first ).to_string().substr(64-  _ham.model().orbitals(), _ham.model().orbitals());
          std::string spin_up   = std::bitset< 64 >( contribs[i].first ).to_string().substr(64-2*_ham.model().orbitals(), _ham.model().orbitals());
          std::cout<<spin_up<< "|"<<spin_down;
          std::cout << ">" << std::endl;
        }
      }
    }

#ifdef USE_MPI
    struct Element{
      Element() {};
      Element(size_t _ind, double _val) : ind(_ind), val(_val) {}
      size_t ind;
      double val;

      bool operator>(const Element &el) const {
        return ((std::abs(val) > std::abs(el.val))
                || ((std::abs(val) == std::abs(el.val)) && (ind < el.ind))
               );
      };

      bool operator<(const Element &el) const {
        return ((std::abs(val) < std::abs(el.val))
                || ((std::abs(val) == std::abs(el.val)) && (ind > el.ind))
               );
      };
    };

     MPI_Datatype mpi_Element;
#endif

  private:

    /**
     * @brief Compute average number of electrons or average spin on the site in one eigenvector.
     *
     * @param orb - the site for calculation
     * @param diff - 0: average number or electrons, 1: average spin
     * @param _ham - the Hamiltonian
     * @param pair - the eigenpair
     */
    precision avgn1(Hamiltonian& _ham, const EigenPair<typename Hamiltonian::ModelType::precision, typename Hamiltonian::ModelType::Sector>& pair, int diff, int orb){
      _ham.model().symmetry().set_sector(pair.sector());
      _ham.storage().reset();
      precision result = 0.0;
      // Loop over basis vectors.
      for(int i = 0; i < pair.eigenvector().size(); ++i){
        // Calculate sum or difference of occupations with different spin on the site of choice.
        _ham.model().symmetry().next_state();
        long long nst = _ham.model().symmetry().state();
        int electrons = 0;
        for(int is = 0; is < _ham.model().spins(); ++is){
          electrons += (diff ? 1 - 2 * is : 1) *  _ham.model().checkState(nst, orb + is * _ham.model().interacting_orbitals(), _ham.model().max_total_electrons());
        }
        // Sum, weighted by square of eigenvector component.
        result += precision(electrons) * pair.eigenvector()[i] * pair.eigenvector()[i];
      }
#ifdef USE_MPI
      precision all;
      // Add the sums from all processes, divide by the norm of eigenvector.
      MPI_Reduce(&result, &all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, _ham.comm());
      return all / _ham.storage().vv(pair.eigenvector(), pair.eigenvector());
#else
      return result / _ham.storage().vv(pair.eigenvector(), pair.eigenvector());
#endif
    }

    /// Inverse temperature
    precision _beta;
    /// Boltzmann-factor cutoff
    precision _cutoff;

  };

}

#endif
