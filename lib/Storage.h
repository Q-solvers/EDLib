//
// Created by iskakoff on 21/07/16.
//

#ifndef HUBBARD_STORAGE_H
#define HUBBARD_STORAGE_H

#include "fortranbinding.h"
#include <iostream>

template<typename prec>
class Storage {
public:
  Storage(size_t max_dim) {
    v.reserve(max_dim);
    resid.reserve(max_dim);
    workd.reserve(3*max_dim);
  }
  /**
   * Diagonalize current Hamiltonian
   */
  void diag() {
    int ido = 0;
    int n = _n;
    char which[3] = "SA";
    double sigma = 0.0;
    char bmat[2] = "I";
    // TODO: move to parameters
    int nev = 2;
    int ncv = std::max(2*nev, 10);
    int lworkl = ncv*(ncv+8);
    double tol = 1e-14;
    int info = 0;
    std::vector<int> iparam(11, 0);
    std::vector<int> ipntr(11, 0);
    std::vector<int> select(ncv, 0);
    int ishfts = 1;
    int maxitr = 1000;
    int mode = 1;
    int ldv = n;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;
    v.resize(n);
    resid.resize(n);
    workd.resize(3*n);
    workl.resize(lworkl);
    do {
      saupd(&ido, bmat, &n, which, &nev, &tol, &resid[0],&ncv, &v[0], &ldv, &iparam[0], &ipntr[0], &workd[0], &workl[0],&lworkl, &info);
      if (ido == -1 || ido == 1) {
        av(&workd[ipntr[0] - 1], &workd[ipntr[1] - 1], n);
      }
    } while (ido!=99);
    if(info<0) {
      std::cout<<"' '"<<std::endl;
      std::cout<<"' Error with _saupd, info = '  "<<info<<std::endl;
      std::cout<<"' Check documentation in _saupd '  "<<iparam[4]<<std::endl;
      std::cout<<"' '"<<std::endl;
      return;
    }
    int rvec = 1;
    char howmny[2] = "A";
    evals.resize(nev);
    seupd(&rvec, howmny, &select[0], &evals[0], &v[0], &ldv, &sigma, bmat, &n, which, &nev, &tol, &resid[0], &ncv, &v[0],
          &ldv,&iparam[0], &ipntr[0], &workd[0], &workl[0], &lworkl, &info);
    // TODO: need to recover the eigenvectors from v
    if(info < 0) {
      std::cout<<"' '"<<std::endl;
      std::cout<<"' Error with _seupd, info = '  "<<info<<std::endl;
      std::cout<<"' Check the documentation of _seupd. '"<<std::endl;
      std::cout<<"' '"<<std::endl;
      return;
    }
    // TODO: save eigenvalues for current sector in local array. Merge all fouded eigen pairs together and keep only N smallest
    for(auto ev : evals) {
      std::cout<<ev<<std::endl;
    }
  }

  /**
   * Matrix-Vector product
   * Should be implemented based on storage type
   */
  virtual void av(prec* v, prec* w, int n) = 0;

  void saupd(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, prec *resid, int *ncv, prec *v, int *ldv, int *iparam, int *ipntr,
             prec *workd, prec *workl, int *lworkl, int *info){};
  inline void seupd(int *rvec, char *All, int *select, prec *d,
                    prec *z, int *ldz, double *sigma,
                    char *bmat, int *n, char *which, int *nev,
                    double *tol, prec *resid, int *ncv, prec *v,
                    int *ldv, int *iparam, int *ipntr, prec *workd,
                    prec *workl, int *lworkl, int *ierr){};

protected:
  int & n() {return _n;}
private:
  // TODO: should set _n before diagonalize
  int _n;
  std::vector<prec> v;
  std::vector<prec> resid;
  std::vector<prec> workd;
  std::vector<prec> workl;

  std::vector<prec> evals;
  std::vector<std::vector<prec> > evecs;
};


template<>
void Storage<double>::saupd(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr,
                               double *workd, double *workl, int *lworkl, int *info) {
  dsaupd_(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
}

template<>
void Storage<double>::seupd(int *rvec, char *All, int *select, double *d,
                               double *z, int *ldz, double *sigma,
                               char *bmat, int *n, char *which, int *nev,
                               double *tol, double *resid, int *ncv, double *v,
                               int *ldv, int *iparam, int *ipntr, double *workd,
                               double *workl, int *lworkl, int *ierr) {
  dseupd_(rvec, All, select, d, z, ldz, sigma, bmat, n, which, nev, tol, resid, ncv, v,
          ldv, iparam, ipntr, workd, workl, lworkl, ierr);
}

#endif //HUBBARD_STORAGE_H
