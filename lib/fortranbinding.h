//
// Created by iskakoff on 20/07/16.
//

#ifndef HUBBARD_FORTRANBINDING_H
#define HUBBARD_FORTRANBINDING_H


#include <vector>

#ifdef __cplusplus
extern "C" {
// TODO: add headers for double complex
void dseupd_(int *rvec, char *All, int *select, double *d, double *z, int *ldz, double *sigma, char *bmat,
        int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv,
        int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *ierr);
void dsaupd_(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv,
        double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info);
};
#endif

#endif //HUBBARD_FORTRANBINDING_H
