//
// Created by iskakoff on 21/07/16.
//

#ifndef EDLIB_STORAGE_H
#define EDLIB_STORAGE_H

#include <iostream>
#include <alps/params.hpp>

#include <arnoldi/arnoldi.hpp>
#ifdef USE_MPI
#include <arnoldi/mpi.hpp>
#endif

namespace EDLib {
  namespace Storage {

    template<typename prec>
    class Storage {
    public:
#ifdef USE_MPI
      Storage(alps::params &p, MPI_Comm comm) : _comm(comm), _nev(p["arpack.NEV"]), _eval_only(p["storage.EIGENVALUES_ONLY"]) {
#else
      Storage(alps::params &p) : _nev(p["arpack.NEV"]), _eval_only(p["storage.EIGENVALUES_ONLY"]) {
#endif
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
       * Diagonalize current Hamiltonian sector using cpp-arnoldi (callback-based
       * C++ port of the ARPACK symmetric driver).
       */
      int diag() {
        if (_n == 0) {
          return finalize(0, true, true);
        }
        if (_ntot == 1) {
          zero_eigenapair();
          return finalize(0);
        }

        const int ncv = std::min(_ncv, _ntot);
        const int nev = std::min(_nev, ncv - 1);

#ifdef USE_MPI
        arnoldi::MPIComm arnoldi_comm(comm());
        arnoldi::Arnoldi<arnoldi::Kind::Sym, prec, arnoldi::MPIComm>
            solver("I", _n, "SA", nev, ncv, arnoldi_comm);
#else
        arnoldi::Arnoldi<arnoldi::Kind::Sym, prec>
            solver("I", _n, "SA", nev, ncv);
#endif
        solver.tol(prec(1e-14)).maxiter(1000).mode(1).ishift(1);

        // Open MPI_Win on the solver's workd buffer (SpinResolvedStorage),
        // or no-op for serial backends (CRSStorage, SOCRSStorage).
        // Shift 2*n matches the existing convention in prepare_work_arrays.
        prepare_work_arrays(solver.workd(), size_t(2 * _n));

        solver.solve([this](const prec* x, prec* y) {
          this->av(const_cast<prec*>(x), y, _n, /*clear=*/true);
        });

        const int info = solver.info();
        if (info < 0) {
          std::cout << "' '" << std::endl;
          std::cout << "' Error with saupd, info = '  " << info << std::endl;
          std::cout << "' '" << std::endl;
          return finalize(info);
        }

        const int nconv = solver.num_converged();
        auto r = solver.eigenpairs(/*compute_vectors=*/_eval_only == 0, prec(0));

        evals.assign(r.values.begin(), r.values.begin() + nconv);
        if (_eval_only == 0) {
          evecs.assign(nconv, std::vector<prec>(_n, prec(0.0)));
          for (int i = 0; i < nconv; ++i) {
            std::memcpy(evecs[i].data(),
                        r.vectors.data() + size_t(i) * _n,
                        _n * sizeof(prec));
          }
        } else {
          evecs.assign(nconv, std::vector<prec>(1, prec(0.0)));
        }

#ifdef USE_MPI
        int myid;
        MPI_Comm_rank(comm(), &myid);
        if (myid == 0) {
#endif
          std::cout << "Here is eigenvalues" << std::endl;
          for (int j = 0; j < (int)evals.size(); ++j) {
            std::cout << evals[j] << std::endl << std::flush;
          }
          if (info == 1) {
            std::cout << "Maximum number of iterations reached." << std::endl;
          } else if (info == 3) {
            std::cout << " No shifts could be applied during implicit Arnoldi update, try increasing NCV." << std::endl;
          }
          std::cout << " ========================= " << std::endl;
          std::cout << " Size of the matrix is " << _ntot << std::endl;
          std::cout << " The number of Ritz values requested is: " << nev << std::endl;
          std::cout << " The number of Arnoldi vectors generated: " << ncv << std::endl;
          std::cout << " What portion of the spectrum: SA" << std::endl;
          std::cout << " The number of converged Ritz values is:  " << nconv << std::endl;
          std::cout << " The number of Implicit Arnoldi update iterations taken is: " << solver.num_iterations() << std::endl;
          std::cout << " The number of OP*x is: " << solver.num_op_applies() << std::endl;
          std::cout << " The convergence criterion is:  " << prec(1e-14) << std::endl;
          std::cout << " ========================= " << std::endl;
#ifdef USE_MPI
        }
#endif
        return finalize(info);
      }

      /**
       * @return eigen-values
       */
      const std::vector < prec > &eigenvalues() const {
        return evals;
      }

      std::vector < prec > &eigenvalues() {
        return evals;
      }

      /**
       * @return eigen-vectors
       */
      const std::vector < std::vector < prec > > &eigenvectors() const {
        return evecs;
      }

      std::vector < std::vector < prec > > &eigenvectors() {
        return evecs;
      }

      /**
       * Matrix-Vector product.
       * Must be implemented by each concrete storage type.
       */
      virtual void av(prec *v, prec *w, int n, bool clear = true) = 0;

      /**
       * Perform additional setup on the solver's working array before the solve
       * loop begins. Used by SpinResolvedStorage to open the MPI_Win.
       * @param w     - pointer to the solver's workd buffer
       * @param shift - offset within w (convention: 2*n)
       */
      virtual void prepare_work_arrays(prec *w, size_t shift = 0){};

      /**
       * Finalize diagonalization for the current Hamiltonian matrix sector.
       * @param info  - result code
       * @param bcast - broadcast eigenvalues across MPI ranks
       * @param empty - arrays are empty on this CPU
       * @return result code
       */
      virtual int finalize(int info, bool bcast = true, bool empty = false){return info;};

#ifdef USE_MPI
      virtual MPI_Comm comm() {
        return _comm;
      }
#endif
    protected:
      /// local CPU dimension
      int &n() { return _n; }
      /// total matrix dimension
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
      /// number of eigenvalues to be computed
      int _nev;
      /// number of Arnoldi vectors (Krylov subspace size)
      int _ncv;
      /// compute only eigenvalues to reduce memory requirements
      int _eval_only;

      std::vector < prec > evals;
      std::vector < std::vector < prec > > evecs;
#ifdef USE_MPI
      MPI_Comm _comm;
#endif
    };

  }
}
#endif //EDLIB_STORAGE_H
