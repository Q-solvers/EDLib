//
// Created by iskakoff on 21/07/16.
//

#ifndef HUBBARD_STORAGE_H
#define HUBBARD_STORAGE_H

#include "fortranbinding.h"
#include <iostream>
#include <EDParams.h>

namespace EDLib {
  namespace Storage {

    template<typename prec>
    class Storage {
    public:
#ifdef ALPS_HAVE_MPI
      Storage(EDParams &p, alps::mpi::communicator& comm) : _nev(p["arpack.NEV"]), _eval_only(p["storage.EIGENVALUES_ONLY"]), _comm(comm) {
#else
      Storage(EDParams &p) : _nev(p["arpack.NEV"]), _eval_only(p["storage.EIGENVALUES_ONLY"]) {
#endif
        v.reserve(size_t(p["storage.MAX_DIM"]));
        resid.reserve(size_t(p["storage.MAX_DIM"]));
        workd.reserve(3 * size_t(p["storage.MAX_DIM"]));
        if (p.exists("arpack.NCV")) {
          _ncv = p["arpack.NCV"];
        } else {
          _ncv = 2 * _nev + 3;
        }
      }

      /**
       * For matrix size equals to 1 we do not need to perform diagonalization.
       * Eigenvalue = A(0,0)
       * Eigenvector = [1.0]
       */
      virtual void zero_eigenapair() = 0;

      /**
       * Diagonalize current Hamiltonian
       */
      int diag() {
        int ido = 0;
        int n = _n;
        if(n==0) {
#ifdef ALPS_HAVE_MPI
//          std::cout<<"Wait"<<comm().rank()<<std::endl;
          MPI_Barrier(_comm);
#endif
          return 0;
        }
        if (n*comm().size() == 1) {
          zero_eigenapair();
          MPI_Barrier(MPI_COMM_WORLD);
          return 0;
        }
        std::cout<<"diag matrix:"<<n<<std::endl;
        int ncv = std::min(_ncv, _n*comm().size());
        int nev = std::min(_nev, ncv - 1);
        char which[3] = "SA";
        prec sigma = 0.0;
        char bmat[2] = "I";
        int lworkl = ncv * (ncv + 8);
        prec tol = 1e-14;
        int info = 0;
        std::vector < int > iparam(11, 0);
        std::vector < int > ipntr(11, 0);
        std::vector < int > select(ncv, 0);
        int ishfts = 1;
        int maxitr = 1000;
        int mode = 1;
        int ldv = n;

        iparam[0] = ishfts;
        iparam[2] = maxitr;
        iparam[6] = mode;
        v.assign(size_t(n) * ncv, prec(0.0));
        resid.assign(size_t(n), prec(0.0));
        workd.assign(3 * size_t(n), prec(0.0));
        workl.assign(lworkl, prec(0.0));
        prepare_work_arrays(&workd[0]);
        do {
          saupd(&ido, bmat, &n, which, &nev, &tol, &resid[0], &ncv, &v[0], &ldv, &iparam[0], &ipntr[0], &workd[0], &workl[0], &lworkl, &info);
          if (ido == -1 || ido == 1) {
            av(&workd[ipntr[0] - 1], &workd[ipntr[1] - 1], n);
          }
        } while (ido != 99);
        if (info < 0) {
          std::cout << "' '" << std::endl;
          std::cout << "' Error with _saupd, info = '  " << info << std::endl;
          std::cout << "' Check documentation in _saupd '  " << iparam[4] << std::endl;
          std::cout << "' '" << std::endl;
          finalize();
          return info;
        }
        int rvec = 1-_eval_only;
        char howmny[2] = "A";
        evals.resize(nev);
        seupd(&rvec, howmny, &select[0], &evals[0], &v[0], &ldv, &sigma, bmat, &n, which, &nev, &tol, &resid[0], &ncv, &v[0],
              &ldv, &iparam[0], &ipntr[0], &workd[0], &workl[0], &lworkl, &info);
        // TODO: need to recover the eigenvectors from v
        if (info < 0) {
          std::cout << "' '" << std::endl;
          std::cout << "' Error with _seupd, info = '  " << info << std::endl;
          std::cout << "' Check the documentation of _seupd. '" << std::endl;
          std::cout << "' '" << std::endl;
          finalize();
          return info;
        }
        // TODO: save eigenvalues for current sector in local array. Merge all fouded eigen pairs together and keep only N smallest
        if(_eval_only == 0) {
          int nconv = iparam[4];
          evecs.assign(nconv, std::vector < prec >(n, prec(0.0)));
          for (int i = 0; i < nconv; ++i) {
            int offset = i * n;
            std::memcpy(&evecs[i][0], &v[offset], n * sizeof(prec));
          }
        } else {
          int nconv = iparam[4];
          evecs.assign(nconv, std::vector < prec >(1, prec(0.0)));
        };
        finalize();
#ifdef ALPS_HAVE_MPI
        MPI_Barrier(_comm);
#endif
        return 0;
      }

      const std::vector < prec > &eigenvalues() const {
        return evals;
      }

      const std::vector < std::vector < prec > > &eigenvectors() const {
        return evecs;
      }

      std::vector < prec > &eigenvalues() {
        return evals;
      }

      std::vector < std::vector < prec > > &eigenvectors() {
        return evecs;
      }

      /**
       * Matrix-Vector product
       * Should be implemented based on storage type
       */
      virtual void av(prec *v, prec *w, int n, bool clear = true) = 0;
      virtual void prepare_work_arrays(prec *w){};
      virtual void finalize(){};

      void saupd(int *ido, char *bmat, int *n, char *which, int *nev, prec *tol, prec *resid, int *ncv, prec *v, int *ldv, int *iparam, int *ipntr,
                 prec *workd, prec *workl, int *lworkl, int *info) {};

      inline void seupd(int *rvec, char *All, int *select, prec *d,
                        prec *z, int *ldz, prec *sigma,
                        char *bmat, int *n, char *which, int *nev,
                        prec *tol, prec *resid, int *ncv, prec *v,
                        int *ldv, int *iparam, int *ipntr, prec *workd,
                        prec *workl, int *lworkl, int *ierr) {};

    protected:
      int &n() { return _n; }

#ifdef ALPS_HAVE_MPI
      virtual alps::mpi::communicator & comm() {
        return _comm;
      }
#endif
    private:
      int _n;
      int _nev;
      int _ncv;
      int _eval_only;
      std::vector < prec > v;
      std::vector < prec > resid;
      std::vector < prec > workd;
      std::vector < prec > workl;

      std::vector < prec > evals;
      std::vector < std::vector < prec > > evecs;
#ifdef ALPS_HAVE_MPI
      alps::mpi::communicator &_comm;
#endif
    };

    template<>
    void Storage < double >::saupd(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr,
                                   double *workd, double *workl, int *lworkl, int *info) {
#ifdef ALPS_HAVE_MPI
      int scomm = PMPI_Comm_c2f(comm());
      pdsaupd_(&scomm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
#else
      dsaupd_(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
#endif
    }

    template<>
    void Storage < double >::seupd(int *rvec, char *All, int *select, double *d,
                                   double *z, int *ldz, double *sigma,
                                   char *bmat, int *n, char *which, int *nev,
                                   double *tol, double *resid, int *ncv, double *v,
                                   int *ldv, int *iparam, int *ipntr, double *workd,
                                   double *workl, int *lworkl, int *ierr) {
#ifdef ALPS_HAVE_MPI
      int scomm = PMPI_Comm_c2f(comm());
      pdseupd_(&scomm, rvec, All, select, d, z, ldz, sigma, bmat, n, which, nev, tol, resid, ncv, v,
               ldv, iparam, ipntr, workd, workl, lworkl, ierr);
#else
      dseupd_(rvec, All, select, d, z, ldz, sigma, bmat, n, which, nev, tol, resid, ncv, v,
              ldv, iparam, ipntr, workd, workl, lworkl, ierr);
#endif
    }

    template<>
    void Storage < float >::saupd(int *ido, char *bmat, int *n, char *which, int *nev, float *tol, float *resid, int *ncv, float *v, int *ldv, int *iparam, int *ipntr,
                                  float *workd, float *workl, int *lworkl, int *info) {
#ifdef ALPS_HAVE_MPI
      int scomm = PMPI_Comm_c2f(comm());
      pssaupd_(&scomm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
#else
      ssaupd_(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
#endif
    }

    template<>
    void Storage < float >::seupd(int *rvec, char *All, int *select, float *d,
                                  float *z, int *ldz, float *sigma,
                                  char *bmat, int *n, char *which, int *nev,
                                  float *tol, float *resid, int *ncv, float *v,
                                  int *ldv, int *iparam, int *ipntr, float *workd,
                                  float *workl, int *lworkl, int *ierr) {
#ifdef ALPS_HAVE_MPI
      int scomm = PMPI_Comm_c2f(comm());
      psseupd_(&scomm, rvec, All, select, d, z, ldz, sigma, bmat, n, which, nev, tol, resid, ncv, v,
               ldv, iparam, ipntr, workd, workl, lworkl, ierr);
#else
      sseupd_(rvec, All, select, d, z, ldz, sigma, bmat, n, which, nev, tol, resid, ncv, v,
              ldv, iparam, ipntr, workd, workl, lworkl, ierr);
#endif
    }
  }
}
#endif //HUBBARD_STORAGE_H
