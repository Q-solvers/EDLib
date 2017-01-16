#ifndef HUBBARD_STATEDESCRIPTION_H
#define HUBBARD_STATEDESCRIPTION_H

#include "EigenPair.h"
#include <vector>

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
      MPI_Datatype types[nitems] = {MPI_INT, MPI_DOUBLE};
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
      std::vector<size_t> largest = std::vector<size_t>(pair.eigenvector().size(), 0);
      for(size_t i = 0; i < largest.size(); ++i){
        largest[i] = i;
      }
      std::partial_sort(largest.begin(), largest.begin()+nmax, largest.end(), Comp(pair.eigenvector()) );
#ifdef USE_MPI
      std::vector<Element> send = std::vector<Element>(nmax, Element(0, 0.0));
      for(size_t i = 0; i < nmax; ++i){
       send[i].val = pair.eigenvector()[largest[i]];
       send[i].ind = largest[i] + _ham.storage().offset();
      }
      int myid;
      MPI_Comm_rank(_ham.comm(), &myid);
      std::vector<size_t> all_inds;
      std::vector<double> all_vals;
      std::vector<Element> all;
      if (myid == 0) {
        int nprocs;
        MPI_Comm_size(_ham.comm(), &nprocs);
        all = std::vector<Element>(nprocs * nmax, Element(0, 0.0));
      }
      MPI_Gather(send.data(), nmax, mpi_Element, all.data(), nmax, mpi_Element, 0, _ham.comm());
      if (myid == 0) {
        std::partial_sort(all.begin(), all.begin()+nmax, all.end(), CompElement() );
        for(size_t i = 0; i < nmax; ++i){
          std::cout << all[i].val << " * ";
          long long nst = _ham.model().symmetry().state_by_index(i);
          for (int j = 1; j >= 0; --j){
            std::cout << "|";
            for (int i = (j + 1)*_ham.model().interacting_orbitals(); i >= j * _ham.model().interacting_orbitals(); --i){
              std::cout << ((nst >> i) & 1);
            }
          }
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

    struct Comp{
      Comp( const std::vector<double> & v ) : _v(v) {}
      bool operator ()(int a, int b) { return std::abs(_v[a]) > std::abs(_v[b]);}
      const std::vector<double> & _v;
    };

#ifdef USE_MPI
    struct Element{
      Element( size_t _ind, double _val ) : ind(_ind), val(_val) {}
      size_t ind;
      double val;
    };

    struct CompElement{
      bool operator ()(Element a, Element b) { return std::abs(a.val) > std::abs(b.val);}
    };

     MPI_Datatype mpi_Element;
#endif

  };

}

#endif
