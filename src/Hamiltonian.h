//
// Created by iskakoff on 19/07/16.
//

#ifndef HUBBARD_HAMILTONIAN_H
#define HUBBARD_HAMILTONIAN_H

#include <vector>
#include <alps/params.hpp>
#include <type_traits>

#include "Combination.h"
#include "Sector.h"

template<typename prec, class Comb=Combination>
class Hamiltonian {
  static_assert(std::is_base_of<Combination, Comb>::value, "Comb should extend Combinatorics");
public:
  /*
   * Allocate space for Hamiltonian matrix
   * \param [in] max_size - the maximum size of values array
   * \param [in] max_dim - the maximum dimension of Hamiltonian matrix
   */
  Hamiltonian(size_t max_size, size_t max_dim, alps::params& p) :
    values(max_size, prec(0.0)),
    col_ind(max_size, 0),
    row_ptr(max_dim+1),
    combination(p),
    Eps(p["NS"], std::vector<double>(p["SPIN"], 0.0)),
    t(p["NS"], std::vector<double>(p["NS"], 0.0)),
    U(p["NS"], 0.0),
    Ns(p["NS"]),
    ms(p["SPIN"]),
    Ip(Ns * ms)
  {};

  /**
   * fill current sector
   */
  void fill(Sector& sector) {
    combination.init(sector);
    int i =0;
    int vind = 0;
    long long k1, k2;
    int isign1, isign2;
    prec xtemp;
    while (combination.next_combination()) {
      long long nst = combination.combination();
      for(int im = 0;im < Ns;++im){
        xtemp += (Eps[im][0]      - xmu) * checkState(nst, im)
               + (Eps[im][ms - 1] - xmu) * checkState(nst, im + Ns);
        xtemp += U[im]* checkState(nst, im)* checkState(nst, im + Ns);

      }
      row_ptr[i] = vind;
      col_ind[vind] = i + 1;
      values[vind] = xtemp;
      vind++;
      for(int ii = 0; ii< Ns; ++ii) {
        for(int jj = 0; jj< Ns; ++jj) {
          if(ii!=jj) {
            if (checkState(nst, ii)) {
              if (!checkState(nst, jj)) {
                a(ii + 1, nst, k1, isign1);
                adag(jj + 1, k1, k2, isign2);
                hopping(i, nst, k2, isign1 * isign2, ii, jj, vind, sector);
              }
            }
            if (checkState(nst, ii+Ns)) {
              if (!checkState(nst, jj+Ns)) {
                a(ii +Ns + 1, nst, k1, isign1);
                adag(jj +Ns + 1, k1, k2, isign2);
                hopping(i, nst, k2, isign1 * isign2, ii, jj, vind, sector);
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
    for(auto &sector : combination.next_sector()) {
      fill(sector);
      /**
       * perform ARPACK call
       */
    }
  }

private:
  // CSR format Hamiltonian matrix storage
  std::vector<prec> values;
  std::vector<int> row_ptr;
  std::vector<int> col_ind;
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
   * \param sign - fermionic sign
   * \param ii - first site
   * \param jj - second site
   * \param vind - current CRS index
   * \param sector - current conservation law sector
   */
  void inline hopping(const int& i, const long long& nst, const long long& k, int sign, const int& ii, const int& jj, int & vind, const Sector& sector) {
    int k_index = combination.ninv_value(k, Ns, sector) - 1;
    int findedstate = 0;
    bool hasstate = false;
    // check that there is no any data on the k state
    for (int iii = row_ptr[i]; iii <= vind; iii++) {
      if (col_ind[iii] == (k_index + 1)) {
        hasstate = true;
        findedstate = iii;
        break;
      }
    }
    // if there is data add value
    if (hasstate) {
      values[findedstate] += t[ii][jj] * sign;
    } else {
    // create new element in CRS arrays
      col_ind[vind] = k_index + 1;
      values[vind] = t[ii][jj] * sign;
      vind++;
    }
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
