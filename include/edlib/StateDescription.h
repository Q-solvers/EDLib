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

    StateDescription(alps::params &p, Hamiltonian& ham) :
      _ham(ham)
    {
      if(p["storage.EIGENVALUES_ONLY"] == 1) {
        throw std::logic_error("Eigenvectors have not been computed. Green's function can not be evaluated.");
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

    std::vector<std::pair<long long, precision>> find(const EigenPair<typename Hamiltonian::ModelType::precision, typename Hamiltonian::ModelType::Sector>& pair, size_t nmax, precision trivial){
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

    void print(const EigenPair<typename Hamiltonian::ModelType::precision, typename Hamiltonian::ModelType::Sector>& pair, size_t nmax, precision trivial){
      std::vector<std::pair<long long, precision>> contribs = find(pair, nmax, trivial);
#ifdef USE_MPI
      int myid;
      MPI_Comm_rank(_ham.comm(), &myid);
      if(!myid)
#endif
      {
        std::cout << "Eigenvector components for eigenvalue " << pair.eigenvalue() << " ";
        // FIXME Depends on find() having set the sector.
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

  private:

    Hamiltonian& _ham;

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

  };

}

#endif
