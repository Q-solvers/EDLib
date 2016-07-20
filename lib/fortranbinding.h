//
// Created by iskakoff on 20/07/16.
//

#ifndef HUBBARD_FORTRANBINDING_H
#define HUBBARD_FORTRANBINDING_H


#include <vector>

#ifdef __cplusplus
extern "C" {

void darnoldi_(int* nloc2,double* vout,double *eout, int*ncv,double* Hstate0, int*nev,int*ierr,int*info);

};

void darnoldi(int& nloc,std::vector<double>& vout,std::vector<double>& eout, int& ncv,double& Hstate0, int& nev,int& ierr,int& info) {
  darnoldi_(&nloc,&vout[0],&eout[0], &ncv,&Hstate0, &nev,&ierr,&info);
}
#endif

#endif //HUBBARD_FORTRANBINDING_H
