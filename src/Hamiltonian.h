//
// Created by iskakoff on 19/07/16.
//

#ifndef HUBBARD_HAMILTONIAN_H
#define HUBBARD_HAMILTONIAN_H

#include <vector>
#include <alps/params.hpp>
#include <type_traits>

#include <CRSStorage.h>
#include "Simmetry.h"

template<typename prec, class Comb=Simmetry, class Storage=CRSStorage<double> >
class Hamiltonian {
  static_assert(std::is_base_of<Simmetry, Comb>::value, "Comb should extend Combinatorics");
public:
  /*
   * Allocate space for Hamiltonian matrix
   * \param [in] max_size - the maximum size of values array
   * \param [in] max_dim - the maximum dimension of Hamiltonian matrix
   */
  Hamiltonian(size_t max_size, size_t max_dim, alps::params& p) :
    storage(max_size, max_dim, p),
    combination(p),
    Eps(p["NSITES"], std::vector<double>(p["NSPINS"], 0.0)),
    t(p["NSITES"], std::vector<double>(p["NSITES"], 0.0)),
    U(p["NSITES"], 0.0),
    Ns(p["NSITES"]),
    ms(p["NSPINS"]),
    Ip(Ns * ms)
  {};

  /**
   * fill current sector
   */
  void fill() {
    combination.init();
    storage.reset();
    int i =0;
    long long k1, k2;
    int isign1, isign2;
    prec xtemp;
    while (combination.next_state()) {
      long long nst = combination.state();
      for(int im = 0;im < Ns;++im){
        xtemp += (Eps[im][0]      - xmu) * checkState(nst, im)
               + (Eps[im][ms - 1] - xmu) * checkState(nst, im + Ns);
        xtemp += U[im]* checkState(nst, im)* checkState(nst, im + Ns);

      }
      storage.addDiagonal(i, xtemp);
      for(int ii = 0; ii< Ns; ++ii) {
        for(int jj = 0; jj< Ns; ++jj) {
          for(int spin = 0; spin< ms; ++spin) {
            if (ii != jj) {
              if (checkState(nst, ii + ms * Ns)) {
                if (!checkState(nst, jj + ms * Ns)) {
                  a(ii + 1, nst, k1, isign1);
                  adag(jj + 1, k1, k2, isign2);
                  hopping(i, nst, k2, isign1 * isign2 * t[ii][jj]);
                }
              }
            }
          }
        }
      }
    }
  }
  /**
   * perform Hamiltonian diagonalization
   * result will be stored in evals and evecs
   */
  void diag() {
    while(combination.next_sector()) {
      fill();
      /**
       * perform ARPACK call
       */
    }
  }

private:
  // CSR format Hamiltonian matrix storage
  Storage storage;
  Comb combination;

  /**
   * xmu - chemical potential
   * Eps - level shift for each spin
   * t - hoppings
   * U - onsite Coulomb interaction
   * Ns - number of orbitals
   * ms - number of spins
   * Ip - maximum total number of electrons ms*Ns
   */
  double xmu;
  std::vector<std::vector<double> > Eps;
  std::vector<std::vector<double> > t;
  std::vector<double> U;
  int Ns;
  int ms;
  int Ip;

  /**
   * Check that im state is occupated
   *
   * \param nst - current state
   * \param im - state to check
   *
   * \return 0 if state is empty, 1 - otherwise
   */
  int inline checkState(const long long& nst, const int& im) {
    return (int) ((nst & (1ll << (Ip - 1 - im))) >> (Ip - 1 - im));
  }

  /**
   * \param i - current index
   * \param nst - curent state
   * \param k - state to hop between
   * \param v - hopping value
   * \param sector - current conservation law sector
   */
  void inline hopping(const int& i, const long long& nst, const long long& k, const prec &v) {
    int k_index = combination.index(k) - 1;
    storage.addElement(i, k_index, v);
  }

  void inline a(const int& i,const long long& jold,long long &k,int &isign) {
    long long sign=0;
    for(int ll=0; ll<i-1; ll++) {
      sign+= ((jold&(1ll<<(Ip-ll-1)))!=0) ? 1 : 0;
    }
    isign = (sign % 2) == 0 ? 1 : -1;
    k=jold-(1ll<<(Ip-i));
  }

  void inline adag(const int& i, const long long &jold, long long& k, int& isign) {
    long long sign=0;
    for(int ll=0; ll<i-1; ll++) {
      sign+= ((jold&(1ll<<(Ip-ll-1)))!=0) ? 1 : 0;
    }
    isign = (sign % 2) == 0 ? 1 : -1;
    k = jold + (1ll << (Ip - i));
  }
};

#endif //HUBBARD_HAMILTONIAN_H
