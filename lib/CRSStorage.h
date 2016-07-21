//
// Created by iskakoff on 20/07/16.
//

#ifndef HUBBARD_CRSSTORAGE_H
#define HUBBARD_CRSSTORAGE_H


#include <vector>
#include <alps/params.hpp>
#include "fortranbinding.h"

template<typename prec>
class CRSStorage {
public:
  CRSStorage(size_t max_size, size_t max_dim, alps::params & p): _vind(0), _max_size(max_size), _max_dim(max_dim) {
    // init what you need from parameters
  };

  void reset(){
    _vind = 0;
    row_ptr.assign(_max_dim+1, 0);
    col_ind.assign(_max_size, 0);
    values.assign(_max_size, prec(0.0));
  }

  void inline addDiagonal(const int &i, prec v) {
    row_ptr[i] = _vind;
    col_ind[_vind] = i + 1;
    values[_vind] = v;
    _vind++;
  }

  /**
   * Add off-diagonal H(i,j) element
   */
  void inline addElement(const int &i, int j, prec t) {
    int findedstate = 0;
    bool hasstate = false;
    // check that there is no any data on the k state
    for (int iii = row_ptr[i]; iii <= _vind; iii++) {
      if (col_ind[iii] == (j + 1)) {
        hasstate = true;
        findedstate = iii;
        break;
      }
    }
    // if there is data add value
    if (hasstate) {
      values[findedstate] += t;
    } else {
      // create new element in CRS arrays
      col_ind[_vind] = j + 1;
      values[_vind] = t;
      _vind++;
    }
  }

  void diag() {
    int ido = 0;
    int n = _n;
    char which[3] = "SA";
    double sigma = 0.0;
    char bmat[2] = "I";
    // TODO: move to parameters
    int nev = 1;
    int ncv = std::max(2*nev, 20);
    int lworkl = ncv*(ncv+8);
    double tol = 1e-14;
    int info = 0;
    std::vector<int> iparam(11, 0);
    std::vector<int> ipntr(11, 0);
    int ishfts = 1;
    int maxitr = 1000;
    int mode = 1;
    int ldv = n;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;
    std::vector<prec> v(n);
    std::vector<prec> resid(n);
    std::vector<prec> workd(3*n);
    std::vector<prec> workl(lworkl);
    do {
      saupd(&ido, bmat, &n, which, &nev, &tol, &resid[0],&ncv, &v[0], &ldv, &iparam[0], &ipntr[0], &workd[0], &workl[0],&lworkl, &info);
      if (ido == -1 || ido == 1) {
        typename std::vector<prec>::const_iterator in = workd.begin() + ipntr[0] - 1;
        typename std::vector<prec>::iterator out = workd.begin() + ipntr[1] - 1;
        av(in, out, n);
      }
    } while (ido!=99);
    // TODO: need to perform seupd to recover the eigen pairs
  }

  void saupd(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, prec *resid, int *ncv, prec *v, int *ldv, int *iparam, int *ipntr,
             prec *workd, prec *workl, int *lworkl, int *info);
  inline void seupd(int *rvec, char *All, int *select, prec *d,
                    prec *z, int *ldz, double *sigma,
                    char *bmat, int *n, char *which, int *nev,
                    double *tol, prec *resid, int *ncv, prec *v,
                    int *ldv, int *iparam, int *ipntr, prec *workd,
                    prec *workl, int *lworkl, int *ierr);

  void av(const typename std::vector<prec>::const_iterator & v, typename std::vector<prec>::iterator & w, int n){

  }

private:
  std::vector<prec> values;
  std::vector<int> row_ptr;
  std::vector<int> col_ind;
  size_t _max_size;
  size_t _max_dim;

  int _n;

  int _vind;
};

template<>
void CRSStorage<double>::saupd(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr,
                               double *workd, double *workl, int *lworkl, int *info) {
  dsaupd_(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
}

template<>
void CRSStorage<double>::seupd(int *rvec, char *All, int *select, double *d,
                               double *z, int *ldz, double *sigma,
                               char *bmat, int *n, char *which, int *nev,
                               double *tol, double *resid, int *ncv, double *v,
                               int *ldv, int *iparam, int *ipntr, double *workd,
                               double *workl, int *lworkl, int *ierr) {
  dseupd_(rvec, All, select, d, z, ldz, sigma, bmat, n, which, nev, tol, resid, ncv, v,
          ldv, iparam, ipntr, workd, workl, lworkl, ierr);
}

#endif //HUBBARD_CRSSTORAGE_H
