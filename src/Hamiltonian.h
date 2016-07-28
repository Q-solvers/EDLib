//
// Created by iskakoff on 19/07/16.
//

#ifndef HUBBARD_HAMILTONIAN_H
#define HUBBARD_HAMILTONIAN_H

#include <vector>
#include <alps/params.hpp>
#include <type_traits>

#include <CRSStorage.h>
#include <fstream>
#include "Symmetry.h"
#include "EigenPair.h"

template<typename prec, class Symmetry=Symmetry, class Storage=CRSStorage<double> >
class Hamiltonian {
  static_assert(std::is_base_of<Symmetry, Symmetry>::value, "Symmetry should extend base Symmetry class");
public:
  /*
   * Allocate space for Hamiltonian matrix
   * \param [in] max_size - the maximum size of values array
   * \param [in] max_dim - the maximum dimension of Hamiltonian matrix
   * \param [in] p - alps::parameters
   */
  Hamiltonian(alps::params& p) :
    storage(p),
    symmetry(p),
    Eps(p["NSITES"], std::vector<double>(p["NSPINS"], 0.0)),
    t(p["NSITES"], std::vector<double>(p["NSITES"], 0.0)),
    U(p["NSITES"], 0.0),
    Ns(p["NSITES"]),
    ms(p["NSPINS"]),
    Ip(Ns * ms),
    _max_dim(p["MAX_DIM"]){
    std::string input = p["INPUT_FILE"];
    alps::hdf5::archive input_data(input.c_str(), "r");
    input_data>>alps::make_pvp("BETA", _beta);
    input_data>>alps::make_pvp("hopping/values", t);
    input_data>>alps::make_pvp("interaction/values", U);
    input_data.close();
    for(int i = 0; i< Ns; ++i ) {
      // HARDCODED Half-filling
      Eps[i][0] = Eps[i][1] = -U[i]/2.0;
    }
    // TODO: move to input file and make site-dependent
    xmu = 0;
  };

  /**
   * fill current sector
   */
  void fill() {
    symmetry.init();
    storage.reset();
    int i =0;
    long long k1, k2;
    int isign1, isign2;
    prec xtemp;
    while (symmetry.next_state()) {
      xtemp = 0.0;
      long long nst = symmetry.state();
      // Compute diagonal element for current i state
      for(int im = 0;im < Ns;++im){
        xtemp += (Eps[im][0]      - xmu) * checkState(nst, im)
               + (Eps[im][ms - 1] - xmu) * checkState(nst, im + Ns);
        xtemp += U[im]* checkState(nst, im)* checkState(nst, im + Ns);

      }
      storage.addDiagonal(i, xtemp);
      // non-diagonal terms calculation
      for(int ii = 0; ii< Ns; ++ii) {
        for(int jj = 0; jj< Ns; ++jj) {
          for(int spin = 0; spin< ms; ++spin) {
            // check that ii site is not equal to jj site and that there is non-zero hopping between these lattice sites
            if (ii != jj && std::abs(t[ii][jj])>1e-10) {
              // check if we can anihilate particle with spin=spin on ii lattice site
              if (checkState(nst, ii + spin * Ns)) {
                // check if we can create particle with spin=spin on jj lattice site
                if (!checkState(nst, jj + spin * Ns)) {
                  a(ii + spin * Ns + 1, nst, k1, isign1);
                  adag(jj + spin * Ns + 1, k1, k2, isign2);
                  hopping(i, nst, k2, -isign1 * isign2 * t[ii][jj]);
                }
              }
            }
          }
        }
      }
      i++;
    }
    // additional steps after all data
    storage.endMatrix();
  }
  /**
   * perform Hamiltonian diagonalization
   * result will be stored in evals and evecs
   */
  void diag() {
    while(symmetry.next_sector()) {
      size_t sector_size = symmetry.sector().size();
      if(sector_size>_max_dim) {
        std::stringstream s;
        s<<"Current sector request more memory than allocated. Increase MAX_DIM parameter. Requested "<<sector_size<<", allocated "<<_max_dim<<".";
        throw std::runtime_error(s.str().c_str());
      }
      fill();
      /**
       * perform ARPACK call
       */
      int info = storage.diag();
      if(info != 0) {

      } else {
        const std::vector<prec>& evals = storage.eigenvalues();
        const std::vector<std::vector<prec> >& evecs = storage.eigenvectors();
        for(int i = 0; i<evals.size(); ++i) {
//          std::vector<prec> evec(evecs[i]);
          eigenpairs.push_back(EigenPair<prec, typename Symmetry::Sector>(evals[i], evecs[i], symmetry.sector()));
        }
      }
    }
    std::sort(eigenpairs.begin(), eigenpairs.end());
    std::cout<<"Here is the list of eigenvalues:"<<std::endl;
    for(auto& eigenpair : eigenpairs) {
      std::cout<<eigenpair.eigenvalue()<<" ";
      eigenpair.sector().print();
      std::cout<<std::endl;
    }
  }

private:
  // CSR format Hamiltonian matrix storage
  Storage storage;
  Symmetry symmetry;

  // Eigen-pairs
  std::vector<EigenPair<prec, typename Symmetry::Sector> > eigenpairs;

  /**
   * xmu - chemical potential
   * _beta - inverse temperature
   * Eps - level shift for each spin
   * t - hoppings
   * U - onsite Coulomb interaction
   * Ns - number of orbitals
   * ms - number of spins
   * Ip - maximum total number of electrons ms*Ns
   */
  double xmu;
  double _beta;
  std::vector<std::vector<double> > Eps;
  std::vector<std::vector<double> > t;
  std::vector<double> U;
  int Ns;
  int ms;
  int Ip;
  size_t _max_dim;

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
    int k_index = symmetry.index(k) - 1;
    storage.addElement(i, k_index, v);
  }

  /**
   * Anihilate particle
   * \param i [in] - site to anihilate particle
   * \param jold [in] - current state
   * \param k [out] - resulting state
   * \param isign [out] - fermionic sign
   */
  void inline a(const int& i,const long long& jold,long long &k,int &isign) {
    long long sign=0;
    for(int ll=0; ll<i-1; ll++) {
      sign+= ((jold&(1ll<<(Ip-ll-1)))!=0) ? 1 : 0;
    }
    isign = (sign % 2) == 0 ? 1 : -1;
    k=jold-(1ll<<(Ip-i));
  }

  /**
   * Create particle
   * \param i [in] - site to create particle
   * \param jold [in] - current state
   * \param k [out] - resulting state
   * \param isign [out] - fermionic sign
   */
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
