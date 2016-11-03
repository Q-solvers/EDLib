//
// Created by iskakoff on 21/07/16.
//

#ifndef HUBBARD_STORAGE_H
#define HUBBARD_STORAGE_H

#include "fortranbinding.h"
#include <iostream>
#include <alps/params.hpp>

namespace EDLib {
  namespace Storage {

    template<typename prec>
    class Storage {
    public:
#ifdef USE_MPI
      Storage(alps::params &p, MPI_Comm comm) : _comm(comm), _nev(p["arpack.NEV"]), _eval_only(p["storage.EIGENVALUES_ONLY"]) {
#else
      Storage(alps::params &p) {
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
        if (n == 0) {
#ifdef USE_MPI
          broadcast_evals(true);
#endif
          return 0;
        }
        if (_ntot == 1) {
          zero_eigenapair();
#ifdef USE_MPI
          broadcast_evals();
#endif
          return 0;
        }
        std::cout << "diag matrix:" << n << std::endl;
        int ncv = std::min(_ncv, _ntot);
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
        prepare_work_arrays(&workd[0], size_t(2 * n));
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
        int rvec = 1 - _eval_only;
        char howmny[2] = "A";
        int nconv = iparam[4];
        evals.resize(nconv);
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
        if (_eval_only == 0) {
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
#ifdef USE_MPI
        broadcast_evals();
        int myid;
        MPI_Comm_rank(comm(), &myid);
        if (myid == 0) {
          std::cout<<"Here is eigenvalues"<<std::endl;
          for (int j = 0; j < evals.size(); ++j) {
            std::cout<<evals[j]<<std::endl<<std::flush;
          }
        }
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
      virtual void prepare_work_arrays(prec *w, size_t shift = 0){};
      virtual void finalize(){};

      void saupd(int *ido, char *bmat, int *n, char *which, int *nev, prec *tol, prec *resid, int *ncv, prec *v, int *ldv, int *iparam, int *ipntr,
                 prec *workd, prec *workl, int *lworkl, int *info) {};

      inline void seupd(int *rvec, char *All, int *select, prec *d,
                        prec *z, int *ldz, prec *sigma,
                        char *bmat, int *n, char *which, int *nev,
                        prec *tol, prec *resid, int *ncv, prec *v,
                        int *ldv, int *iparam, int *ipntr, prec *workd,
                        prec *workl, int *lworkl, int *ierr) {};
#ifdef USE_MPI
      virtual MPI_Comm comm() {
        return _comm;
      }
#endif
    protected:
      int &n() { return _n; }
      int &ntot() { return _ntot; }

#ifdef USE_MPI
      void broadcast_evals(bool empty = false) {
        MPI_Barrier(_comm);
        int nconv = evals.size();
        MPI_Bcast(&nconv, 1, MPI_INT, 0, _comm);
        int rank;
        MPI_Comm_rank(_comm, &rank);
        if(rank != 0) {
          evals.resize(nconv);
          if(empty) {
            evecs.assign(nconv, std::vector<prec>(0, prec(0.0)));
          }
        }
        MPI_Bcast(evals.data(), nconv, alps::mpi::detail::mpi_type<prec>(), 0, _comm);
      }

#endif
    private:
      int _ntot;
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
#ifdef USE_MPI
      MPI_Comm _comm;
#endif
    };

    template<>
    void Storage < double >::saupd(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr,
                                   double *workd, double *workl, int *lworkl, int *info) {
#ifdef USE_MPI
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
#ifdef USE_MPI
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
#ifdef USE_MPI
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
#ifdef USE_MPI
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
