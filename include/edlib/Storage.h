#ifndef EDLIB_STORAGE_H
#define EDLIB_STORAGE_H

#include <algorithm>
#include <cstring>
#include <iostream>
#include <vector>

#include <arnoldi/arnoldi.hpp>
#ifdef USE_MPI
#include <arnoldi/mpi.hpp>
#endif

#include "edlib/MpiTypes.h"
#include "edlib/Parameters.h"

namespace edlib {

  /**
   * Abstract base storage. Concrete derived storages implement matrix-vector
   * product (av) and sector fill (fill); diag() drives the cpp-arnoldi
   * symmetric solver against av().
   */
  template <class Prec>
  class Storage {
  public:
#ifdef USE_MPI
    Storage(const Parameters& p, MPI_Comm comm)
        : _nev(p.arpack_nev),
          _ncv(p.arpack_ncv > 0 ? p.arpack_ncv : 2 * p.arpack_nev + 3),
          _eval_only(p.eigenvalues_only ? 1 : 0),
          _comm(comm) {}
#else
    explicit Storage(const Parameters& p)
        : _nev(p.arpack_nev),
          _ncv(p.arpack_ncv > 0 ? p.arpack_ncv : 2 * p.arpack_nev + 3),
          _eval_only(p.eigenvalues_only ? 1 : 0) {}
#endif

    virtual ~Storage() = default;

    /// For a 1x1 matrix the eigenpair is trivial; concrete storages set it.
    virtual void zero_eigenapair() = 0;

    /// Matrix-vector product. Implemented by concrete storage.
    virtual void av(Prec* v, Prec* w, int n, bool clear = true) = 0;

    /// Optional hook before the solve loop (used by SpinResolvedStorage for
    /// MPI_Win setup). Default is no-op.
    virtual void prepare_work_arrays(Prec* /*w*/, std::size_t /*shift*/ = 0) {}

    /// Optional finalisation hook called after diag() completes.
    virtual int  finalize(int info, bool /*bcast*/ = true, bool /*empty*/ = false) { return info; }

    int diag() {
      if (_n == 0) return finalize(0, true, true);
      if (_ntot == 1) {
        zero_eigenapair();
        return finalize(0);
      }

      const int ncv = std::min(_ncv, _ntot);
      const int nev = std::min(_nev, ncv - 1);

#ifdef USE_MPI
      arnoldi::MPIComm arnoldi_comm(comm());
      arnoldi::Arnoldi<arnoldi::Kind::Sym, Prec, arnoldi::MPIComm>
          solver("I", _n, "SA", nev, ncv, arnoldi_comm);
#else
      arnoldi::Arnoldi<arnoldi::Kind::Sym, Prec>
          solver("I", _n, "SA", nev, ncv);
#endif
      solver.tol(Prec(1e-14)).maxiter(1000).mode(1).ishift(1);

      prepare_work_arrays(solver.workd(), std::size_t(2 * _n));

      solver.solve([this](const Prec* x, Prec* y) {
        this->av(const_cast<Prec*>(x), y, _n, /*clear=*/true);
      });

      const int info = solver.info();
      if (info < 0) {
        std::cout << "' '\n' Error with saupd, info = '  " << info << "\n' '" << std::endl;
        return finalize(info);
      }

      const int nconv = solver.num_converged();
      auto r = solver.eigenpairs(/*compute_vectors=*/_eval_only == 0, Prec(0));

      evals.assign(r.values.begin(), r.values.begin() + nconv);
      if (_eval_only == 0) {
        evecs.assign(nconv, std::vector<Prec>(_n, Prec(0)));
        for (int i = 0; i < nconv; ++i) {
          std::memcpy(evecs[i].data(),
                      r.vectors.data() + std::size_t(i) * _n,
                      _n * sizeof(Prec));
        }
      } else {
        evecs.assign(nconv, std::vector<Prec>(1, Prec(0)));
      }

#ifdef USE_MPI
      int myid;
      MPI_Comm_rank(comm(), &myid);
      if (myid == 0) {
#endif
        std::cout << "Here is eigenvalues" << std::endl;
        for (auto e : evals) std::cout << e << "\n";
        if (info == 1) std::cout << "Maximum number of iterations reached.\n";
        else if (info == 3) std::cout << "No shifts could be applied during implicit Arnoldi update, try increasing NCV.\n";
        std::cout << " ========================= \n"
                  << " Size of the matrix is " << _ntot << "\n"
                  << " The number of Ritz values requested is: " << nev << "\n"
                  << " The number of Arnoldi vectors generated: " << ncv << "\n"
                  << " What portion of the spectrum: SA\n"
                  << " The number of converged Ritz values is:  " << nconv << "\n"
                  << " The number of Implicit Arnoldi update iterations taken is: " << solver.num_iterations() << "\n"
                  << " The number of OP*x is: " << solver.num_op_applies() << "\n"
                  << " The convergence criterion is:  " << Prec(1e-14) << "\n"
                  << " ========================= " << std::endl;
#ifdef USE_MPI
      }
#endif
      return finalize(info);
    }

    const std::vector<Prec>&              eigenvalues()  const { return evals; }
    std::vector<Prec>&                    eigenvalues()        { return evals; }
    const std::vector<std::vector<Prec>>& eigenvectors() const { return evecs; }
    std::vector<std::vector<Prec>>&       eigenvectors()       { return evecs; }

    // Generic eigenpair accessors used by Hamiltonian to build the
    // EigenPair set. Default (host) eigenvector type is std::vector<Prec>;
    // device-resident storages shadow `eigenvector_type` and
    // `eigenpair_vector(i)` to hand back a device buffer instead, so the
    // eigenvectors never round-trip through host memory.
    using eigenvector_type = std::vector<Prec>;
    int  num_eigenpairs() const { return static_cast<int>(evals.size()); }
    const Prec& eigenpair_value(int i) const { return evals[i]; }
    const eigenvector_type& eigenpair_vector(int i) const { return evecs[i]; }

#ifdef USE_MPI
    virtual MPI_Comm comm() { return _comm; }
#endif

  protected:
    int& n()    { return _n; }
    int& ntot() { return _ntot; }

#ifdef USE_MPI
    void broadcast_evals(bool empty = false) {
      MPI_Barrier(_comm);
      int nconv = static_cast<int>(evals.size());
      MPI_Bcast(&nconv, 1, MPI_INT, 0, _comm);
      int rank;
      MPI_Comm_rank(_comm, &rank);
      if (rank != 0) {
        evals.resize(nconv);
        if (empty) evecs.assign(nconv, std::vector<Prec>(0, Prec(0)));
      }
      MPI_Bcast(evals.data(), nconv, mpi_type<Prec>(), 0, _comm);
    }
#endif

  private:
    int _ntot = 0;
    int _n    = 0;
    int _nev;
    int _ncv;
    int _eval_only;

    std::vector<Prec>              evals;
    std::vector<std::vector<Prec>> evecs;
#ifdef USE_MPI
    MPI_Comm _comm;
#endif
  };

}

#endif
