//
// Created by iskakoff on 20/07/16.
//

#ifndef HUBBARD_FORTRANBINDING_H
#define HUBBARD_FORTRANBINDING_H


#include <vector>
#include <alps/config.hpp>

#ifdef __cplusplus
extern "C" {
// TODO: add headers for double complex
void dseupd_(int *rvec, char *All, int *select, double *d, double *z, int *ldz, double *sigma, char *bmat,
        int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv,
        int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *ierr);
void dsaupd_(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv,
        double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info);
void sseupd_(int *rvec, char *All, int *select, float *d, float *z, int *ldz, float *sigma, char *bmat,
        int *n, char *which, int *nev, float *tol, float *resid, int *ncv, float *v, int *ldv,
        int *iparam, int *ipntr, float *workd, float *workl, int *lworkl, int *ierr);
void ssaupd_(int *ido, char *bmat, int *n, char *which, int *nev, float *tol, float *resid, int *ncv,
        float *v, int *ldv, int *iparam, int *ipntr, float *workd, float *workl, int *lworkl, int *info);
void dmout(int* lout, int *m, int*n, double*A, int*lda, int* idigit,char* ifmt);
void smout(int* lout, int *m, int*n, float*A, int*lda, int* idigit,char* ifmt);
#ifdef USE_MPI
void pdseupd_(int *comm, int *rvec, char *All, int *select, double *d, double *z, int *ldz, double *sigma, char *bmat,
             int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv,
             int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *ierr);
void pdsaupd_(int *comm, int *ido, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv,
             double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info);
void psseupd_(int *comm, int *rvec, char *All, int *select, float *d, float *z, int *ldz, float *sigma, char *bmat,
             int *n, char *which, int *nev, float *tol, float *resid, int *ncv, float *v, int *ldv,
             int *iparam, int *ipntr, float *workd, float *workl, int *lworkl, int *ierr);
void pssaupd_(int *comm, int *ido, char *bmat, int *n, char *which, int *nev, float *tol, float *resid, int *ncv,
             float *v, int *ldv, int *iparam, int *ipntr, float *workd, float *workl, int *lworkl, int *info);
void pdmout(int* comm, int* lout, int *m, int*n, double*A, int*lda, int* idigit,char* ifmt);
void psmout(int* comm, int* lout, int *m, int*n, float*A, int*lda, int* idigit,char* ifmt);
#endif
};
#endif

#endif //HUBBARD_FORTRANBINDING_H
