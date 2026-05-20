#ifndef EDLIB_CRSSTORAGE_H
#define EDLIB_CRSSTORAGE_H

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "edlib/HostKernel.h"
#include "edlib/Parameters.h"
#include "edlib/Storage.h"

namespace edlib {

  template <class ModelType>
  class CRSStorage : public Storage<typename ModelType::precision> {
  public:
    using Model = ModelType;
    using prec  = typename ModelType::precision;
    using kernel_type = HostKernel<CRSStorage>;
    using Storage<prec>::n;
    using Storage<prec>::ntot;

#ifdef USE_MPI
    CRSStorage(const Parameters& p, Model& m, MPI_Comm comm)
        : Storage<prec>(p, comm),
          _max_size(p.storage_max_size),
          _max_dim(p.storage_max_dim),
          _model(m) {}
#else
    CRSStorage(const Parameters& p, Model& m)
        : Storage<prec>(p),
          _max_size(p.storage_max_size),
          _max_dim(p.storage_max_dim),
          _model(m) {}
#endif

    void init() { _model.symmetry().init(); }

    void reset() {
      _model.symmetry().init();
      const std::size_t sector_size = _model.symmetry().sector().size();
      if (sector_size > _max_dim) {
        std::stringstream s;
        s << "CRSStorage: sector requests more memory than allocated. "
             "Increase storage.MAX_DIM. Requested " << sector_size
          << ", allocated " << _max_dim << ".";
        throw std::runtime_error(s.str());
      }
      _vind = 0;
      row_ptr.assign(_max_dim + 1, 0);
      col_ind.assign(_max_size, 0);
      values .assign(_max_size, prec(0));
      n()    = 0;
      ntot() = 0;
    }

    void av(prec* v, prec* w, int n_local, bool clear = true) override {
      for (int i = 0; i < n_local; ++i) {
        w[i] = clear ? prec(0) : w[i];
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
          w[i] += values[j] * v[col_ind[j]];
        }
      }
    }

    void fill() {
      reset();
      int i = 0;
      while (_model.symmetry().next_state()) {
        long long nst = _model.symmetry().state();
        addDiagonal(i, _model.diagonal(nst));
        off_diagonal(nst, i, _model.T_states());
        off_diagonal(nst, i, _model.V_states());
        ++i;
      }
      endMatrix();
    }

    void print() const {
      std::cout << std::setprecision(2) << std::fixed << "{";
      for (int i = 0; i < n(); ++i) {
        std::cout << "{";
        for (int j = 0; j < n(); ++j) {
          bool f = true;
          for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
            if (col_ind[k] == j) {
              std::cout << std::setw(6) << values[k] << (j == n() - 1 ? "" : ", ");
              f = false;
            }
          }
          if (f) std::cout << std::setw(6) << 0.0 << (j == n() - 1 ? "" : ", ");
        }
        std::cout << "}" << (i == n() - 1 ? "" : ", \n");
      }
      std::cout << "}" << std::endl;
    }

    void zero_eigenapair() override {
      this->eigenvalues().resize(1);
      this->eigenvalues()[0] = values[0];
      this->eigenvectors().assign(1, std::vector<prec>(1, prec(1)));
    }

    std::size_t vector_size(typename Model::Sector sector) const {
      return sector.size();
    }

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

    void constant_shift(prec shift) {
      for (int i = 0; i < n(); ++i) {
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
          if (col_ind[j] == i) values[j] += shift;
        }
      }
    }

  private:
    void addDiagonal(int i, prec v) {
      row_ptr[i]    = static_cast<int>(_vind);
      col_ind[_vind] = i;
      values [_vind] = v;
      ++_vind;
      ++n();
      ++ntot();
    }

    void addElement(int i, int j, prec t, int sign) {
      bool        hasstate = false;
      std::size_t foundstate = 0;
      for (std::size_t k = row_ptr[i]; k < _vind; ++k) {
        if (col_ind[k] == j) { hasstate = true; foundstate = k; }
      }
      if (hasstate) {
        values[foundstate] += static_cast<prec>(sign) * t;
      } else {
        col_ind[_vind] = j;
        values [_vind] = static_cast<prec>(sign) * t;
        ++_vind;
      }
      if (_vind > _max_size) {
        std::stringstream s;
        s << "CRSStorage: sector requests more memory than allocated. "
             "Increase storage.MAX_SIZE. Requested " << _vind
          << ", allocated " << _max_size << ".";
        throw std::runtime_error(s.str());
      }
    }

    template <class TStates>
    void off_diagonal(long long nst, int i, const TStates& states) {
      long long k = 0;
      int isign = 0;
      for (std::size_t kkk = 0; kkk < states.size(); ++kkk) {
        if (_model.valid(states[kkk], nst)) {
          prec val = _model.set(states[kkk], nst, k, isign);
          int  k_index = _model.symmetry().index(k);
          addElement(i, k_index, val, isign);
        }
      }
    }

    void endMatrix() { row_ptr[n()] = static_cast<int>(_vind); }

    std::vector<prec> values;
    std::vector<int>  row_ptr;
    std::vector<int>  col_ind;
    std::size_t       _max_size;
    std::size_t       _max_dim;
    std::size_t       _vind = 0;
    Model&            _model;
  };

}

#endif
