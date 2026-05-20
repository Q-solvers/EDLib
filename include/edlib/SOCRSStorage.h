#ifndef EDLIB_SOCRSSTORAGE_H
#define EDLIB_SOCRSSTORAGE_H

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "edlib/HostKernel.h"
#include "edlib/Parameters.h"
#include "edlib/Storage.h"

namespace edlib {

  /**
   * Spin-orbital-coupled CRS storage: holds diagonal explicitly and stores
   * off-diagonal column indices + signs in a compact bitmap. Matrix-vector
   * product re-evaluates state values on the fly.
   */
  template <class ModelType>
  class SOCRSStorage : public Storage<typename ModelType::precision> {
  public:
    using Model = ModelType;
    using prec  = typename ModelType::precision;
    using kernel_type = HostKernel<SOCRSStorage>;
    using Storage<prec>::n;
    using Storage<prec>::ntot;

#ifdef USE_MPI
    SOCRSStorage(const Parameters& p, Model& m, MPI_Comm comm)
        : Storage<prec>(p, comm),
#else
    SOCRSStorage(const Parameters& p, Model& m)
        : Storage<prec>(p),
#endif
          _max_size(p.storage_max_size),
          _max_dim (p.storage_max_dim),
#ifdef _OPENMP
          _nthreads(omp_get_max_threads()),
#else
          _nthreads(1),
#endif
          _row_offset (_nthreads + 1),
          _vind_offset(_nthreads + 1),
          _vind       (_nthreads),
          _vind_byte  (_nthreads),
          _vind_bit   (_nthreads),
          _vind_start (_nthreads),
          _model(m) {
      col_ind.assign(_max_size, 0);
      signs  .assign(static_cast<std::size_t>(std::ceil(double(_max_size) / sizeof(char))), 1);
      dvalues.assign(_max_dim, prec(0));
    }

    void av(prec* v, prec* w, int n_local, bool clear = true) override {
      _model.symmetry().init();
#ifdef _OPENMP
#pragma omp parallel
      {
        int myid = omp_get_thread_num();
#else
        int myid = 0;
#endif
        std::size_t vind      = _vind_offset[myid];
        std::size_t vind_byte = vind / sizeof(char);
        std::size_t vind_bit  = vind % sizeof(char);
        for (int i = _row_offset[myid]; (i < (int)_row_offset[myid + 1]) && (i < n_local); ++i) {
          long long nst = _model.symmetry().state_by_index(i);
          w[i] = dvalues[i] * v[i] + (clear ? prec(0) : w[i]);
          // T_states
          for (std::size_t kkk = 0; kkk < _model.T_states().size(); ++kkk) {
            int test = _model.valid(_model.T_states()[kkk], nst);
            w[i] += test * _model.T_states()[kkk].value()
                  * (1 - 2 * ((signs[vind_byte] >> vind_bit) & 1)) * v[col_ind[vind]];
            vind_bit  += test;
            vind_byte += vind_bit / sizeof(char);
            vind_bit  %= sizeof(char);
            vind      += test;
          }
          // V_states
          for (std::size_t kkk = 0; kkk < _model.V_states().size(); ++kkk) {
            int test = _model.valid(_model.V_states()[kkk], nst);
            w[i] += test * _model.V_states()[kkk].value()
                  * (1 - 2 * ((signs[vind_byte] >> vind_bit) & 1)) * v[col_ind[vind]];
            vind_bit  += test;
            vind_byte += vind_bit / sizeof(char);
            vind_bit  %= sizeof(char);
            vind      += test;
          }
        }
#ifdef _OPENMP
      }
#endif
    }

    void init() { _model.symmetry().init(); }

    void reset() {
      _model.symmetry().init();
      const std::size_t sector_size = _model.symmetry().sector().size();
      if (sector_size > _max_dim) {
        std::stringstream s;
        s << "SOCRSStorage: sector requests more memory than allocated. "
             "Increase storage.MAX_DIM. Requested " << sector_size
          << ", allocated " << _max_dim << ".";
        throw std::runtime_error(s.str());
      }
      const std::size_t nnz_estimate = sector_size
            * (_model.T_states().size() + _model.V_states().size());
      if (nnz_estimate > _max_size) {
        std::stringstream s;
        s << "SOCRSStorage: sector requests more memory than allocated. "
             "Increase storage.MAX_SIZE. Requested " << nnz_estimate
          << ", allocated " << _max_size << ".";
        throw std::runtime_error(s.str());
      }
      for (int t = 0; t < _nthreads; ++t) {
        _vind[t]      = 0;
        _vind_byte[t] = 0;
        _vind_bit[t]  = 0;
      }
      ntot() = static_cast<int>(sector_size);
      n()    = ntot();
    }

    void fill() {
      reset();
      const std::size_t sector_size = _model.symmetry().sector().size();
      const int step = static_cast<int>(std::floor(double(sector_size) / _nthreads));
      for (int i = 0; i <= _nthreads; ++i) _row_offset[i] = step * i;
      int more = static_cast<int>(sector_size) - static_cast<int>(_row_offset[_nthreads]);
      for (int i = 0; i < more; ++i)         _row_offset[i] += i;
      for (int i = more; i <= _nthreads; ++i) _row_offset[i] += more;
      for (int t = 0; t <= _nthreads; ++t)
        _vind_offset[t] = (_model.T_states().size() + _model.V_states().size()) * _row_offset[t];

#ifdef _OPENMP
#pragma omp parallel
      {
        int myid = omp_get_thread_num();
#else
        int myid = 0;
#endif
        _vind[myid]      = _vind_offset[myid];
        _vind_byte[myid] = _vind[myid] / sizeof(char);
        _vind_bit[myid]  = _vind[myid] % sizeof(char);
        for (int i = _row_offset[myid]; i < (int)_row_offset[myid + 1]; ++i) {
          long long nst = _model.symmetry().state_by_index(i);
          addDiagonal(i, _model.diagonal(nst), myid);
          off_diagonal(nst, i, _model.T_states(), myid);
          off_diagonal(nst, i, _model.V_states(), myid);
        }
#ifdef _OPENMP
      }
#endif
    }

    void zero_eigenapair() override {
      this->eigenvalues().resize(1);
      this->eigenvalues()[0] = dvalues[0];
      this->eigenvectors().assign(1, std::vector<prec>(1, prec(1)));
    }

    std::size_t vector_size(typename Model::Sector sector) const { return sector.size(); }

#ifdef USE_MPI
    prec vv(const std::vector<prec>& v, const std::vector<prec>& w, MPI_Comm /*com*/) const {
      return vv(v, w);
    }
#endif
    prec vv(const std::vector<prec>& v, const std::vector<prec>& w) const {
      prec alf = prec(0);
      for (std::size_t k = 0; k < v.size(); ++k) alf += w[k] * v[k];
      return alf;
    }

    void a_adag(int iii, const std::vector<prec>& invec, std::vector<prec>& outvec,
                const typename Model::Sector& next_sec, bool a) {
      long long k;
      int sign;
      int i = 0;
      while (_model.symmetry().next_state()) {
        long long nst = _model.symmetry().state();
        if (_model.checkState(nst, iii, _model.max_total_electrons()) == (a ? 1 : 0)) {
          if (a) _model.a   (iii, nst, k, sign);
          else   _model.adag(iii, nst, k, sign);
          int i1 = _model.symmetry().index(k, next_sec);
          outvec[i1] = sign * invec[i];
        }
        ++i;
      }
    }

#ifdef _OPENMP
    int& nprocs() { return _nthreads; }
#endif

    void constant_shift(prec shift) {
      std::transform(dvalues.begin(), dvalues.end(), dvalues.begin(),
                     [shift](prec x) { return x + shift; });
    }

  private:
    void addDiagonal(int i, prec v, int chunk) {
      dvalues[i] = v;
      _vind_start[chunk] = _vind[chunk];
    }

    void addElement(int i, int j, prec /*t*/, int sign, int chunk) {
      if (i == j) {
        throw std::logic_error("SOCRSStorage::addElement: cannot add diagonal; use addDiagonal");
      }
      for (std::size_t k = _vind_start[chunk]; k < _vind[chunk]; ++k) {
        if (col_ind[k] == j) {
          throw std::logic_error("SOCRSStorage::addElement: column collision");
        }
      }
      if (_vind[chunk] >= _vind_offset[chunk + 1]) {
        throw std::runtime_error("SOCRSStorage: sector requests more memory than allocated. Increase storage.MAX_SIZE.");
      }
      col_ind[_vind[chunk]] = j;
      signs[_vind_byte[chunk]] &= ~(1ll << _vind_bit[chunk]);
      signs[_vind_byte[chunk]] |= (sign < 0) ? (1ll << _vind_bit[chunk]) : 0;
      ++_vind_bit[chunk];
      ++_vind[chunk];
      _vind_byte[chunk] += _vind_bit[chunk] / sizeof(char);
      _vind_bit [chunk] %= sizeof(char);
    }

    template <class TStates>
    void off_diagonal(long long nst, int i, const TStates& states, int chunk) {
      long long k = 0;
      int isign  = 0;
      for (std::size_t kkk = 0; kkk < states.size(); ++kkk) {
        if (_model.valid(states[kkk], nst)) {
          _model.set(states[kkk], nst, k, isign);
          int k_index = _model.symmetry().index(k);
          addElement(i, k_index, states[kkk].value(), isign, chunk);
        }
      }
    }

    std::vector<prec> dvalues;
    std::vector<int>  col_ind;
    std::vector<char> signs;

    std::size_t _max_size;
    std::size_t _max_dim;

    int _nthreads;
    std::vector<std::size_t> _row_offset;
    std::vector<std::size_t> _vind_offset;
    std::vector<std::size_t> _vind;
    std::vector<std::size_t> _vind_byte;
    std::vector<std::size_t> _vind_bit;
    std::vector<std::size_t> _vind_start;

    Model& _model;
  };

}

#endif
