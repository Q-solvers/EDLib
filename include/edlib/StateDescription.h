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

    StateDescription(Hamiltonian& ham) :
      _ham(ham)
    {
#ifdef USE_MPI
      const int nitems=2;
      int blocklengths[nitems] = {1, 1};
      MPI_Datatype types[nitems] = {MPI_INT, alps::mpi::detail::mpi_type<precision>()};
      MPI_Aint offsets[nitems];
      offsets[0] = offsetof(Element, ind);
      offsets[1] = offsetof(Element, val);
      MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_Element);
      MPI_Type_commit(&mpi_Element);
#endif
    };

    void print(const EigenPair<typename Hamiltonian::ModelType::precision, typename Hamiltonian::ModelType::Sector>& pair, int nmax, precision trivial){
      std::cout << "Eigenvector components for eigenvalue " << pair.eigenvalue() << " ";
      pair.sector().print();
      std::cout << std::endl;
      _ham.model().symmetry().set_sector(pair.sector());
      //std::vector<size_t> largest = std::vector<size_t>(std::min(nmax, pair.eigenvector().size()), 0);
      std::vector<size_t> largest(pair.eigenvector().size(), 0);
      for(size_t i = 0; i < largest.size(); ++i){
        largest[i] = i;
      }
      std::partial_sort(largest.begin(), largest.begin()+nmax, largest.end(), [&pair] (int a, int b) -> bool {return std::abs(pair.eigenvector()[a]) > std::abs(pair.eigenvector()[b]);} );
#ifdef USE_MPI
      std::vector<Element> send;
      for(size_t i = 0; i < nmax; ++i){
       send.push_back(Element(pair.eigenvector()[largest[i]], largest[i] + _ham.storage().offset()));
      }
      int myid;
      MPI_Comm_rank(_ham.comm(), &myid);
      std::vector<Element> all;
      if (myid == 0) {
        int nprocs;
        MPI_Comm_size(_ham.comm(), &nprocs);
        all.resize(nprocs * nmax);
      }
      MPI_Gather(send.data(), nmax, mpi_Element, all.data(), nmax, mpi_Element, 0, _ham.comm());
      if (myid == 0) {
        std::partial_sort(all.begin(), all.begin()+nmax, all.end());
        for(size_t i = 0; i < nmax; ++i){
          std::cout << all[i].val << " * |";
          long long nst = _ham.model().symmetry().state_by_index(all[i].ind);
          std::string spin_down = std::bitset< 64 >( nst ).to_string().substr(64-  _ham.model().orbitals(), _ham.model().orbitals());
          std::string spin_up   = std::bitset< 64 >( nst ).to_string().substr(64-2*_ham.model().orbitals(), _ham.model().orbitals());
          std::cout<<spin_up<< "|"<<spin_down;
          std::cout << ">" << std::endl;
          //std::cout << " * |" <<  std::bitset<sizeof(long long)*8>(nst) << ">" << std::endl;
        }
      }
#else
      for(size_t i = 0; i < nmax; ++i){
        std::cout << pair.eigenvector()[largest[i]] << " * ";
        long long nst = _ham.model().symmetry().state_by_index(i);
        for (int j = 1; j >= 0; --j){
          std::cout << "|";
          for (int i = (j + 1)*_ham.model().interacting_orbitals(); i >= j * _ham.model().interacting_orbitals(); --i){
            std::cout << ((nst >> i) & 1);
          }
        }
        std::cout << ">" << std::endl;
      }
      std::cout << std::endl;
#endif
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
        return (val > el.val);
      };

      bool operator<(const Element &el) const {
        return (val < el.val);
      };
    };

     MPI_Datatype mpi_Element;
#endif

  };

}

#endif
