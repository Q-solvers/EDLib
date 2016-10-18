//
// Created by iskakoff on 28/07/16.
//

#ifndef HUBBARD_SOCRSSTORAGE_H
#define HUBBARD_SOCRSSTORAGE_H

#include <vector>
#include <iomanip>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Storage.h"

namespace EDLib {
  namespace Storage {

    template<typename prec, class Model>
    class SOCRSStorage : public Storage < prec > {
    public:
      using Storage < prec >::n;
      using Storage < prec >::ntot;
#ifdef USE_MPI
      SOCRSStorage(alps::params &p, Model &m, alps::mpi::communicator &comm) : Storage < prec >(p, comm),
#else
      SOCRSStorage(EDParams &p, Model &m) : Storage < prec >(p),
#endif
                                            _max_size(p["storage.MAX_SIZE"]),
#ifdef _OPENMP
                                            _nthreads(omp_get_max_threads()),
#else
                                            _nthreads(1),
#endif
                                            _row_offset(_nthreads + 1), _vind_offset(_nthreads + 1),
                                            _vind(_nthreads), _vind_byte(_nthreads), _vind_bit(_nthreads), _vind_start(_nthreads),
                                            _max_dim(p["storage.MAX_DIM"]), _model(m) {
        /** init what you need from parameters*/
        col_ind.assign(_max_size, 0);
        // XXX I don't trust myself about this one:
        signs.assign(std::ceil(_max_size / sizeof(char)), 1);
        dvalues.assign(_max_dim, prec(0.0));
      };

      virtual void av(prec *v, prec *w, int n, bool clear = true) {
        _model.symmetry().init();
#pragma omp parallel for schedule(static, 1)
        for(int myid = 0; myid < _nthreads; ++myid){
          _vind[myid] = _vind_offset[myid];
          _vind_byte[myid] = _vind[myid] / sizeof(char);
          _vind_bit[myid] = _vind[myid] % sizeof(char);
          // Iteration over rows.
          for(int i = _row_offset[myid]; (i < _row_offset[myid + 1]) && (i < n); ++i){
            long long nst = _model.symmetry().state_by_index(i);
            // Diagonal contribution.
            w[i] = dvalues[i] * v[i] + (clear ? 0.0 : w[i]);
            // Offdiagonal contribution.
            // Iteration over columns(unordered).
            for (int kkk = 0; kkk < _model.T_states().size(); ++kkk) {
              int test = _model.valid(_model.T_states()[kkk], nst);
              // If transition between states corresponding to row and column is possible, calculate the offdiagonal element.
              w[i] += test * _model.T_states()[kkk].value() * (1 - 2 * ((signs[_vind_byte[myid]] >> _vind_bit[myid]) & 1)) * v[col_ind[_vind[myid]]];
              _vind_bit[myid] += test;
              _vind_byte[myid] += _vind_bit[myid] / sizeof(char);
              _vind_bit[myid] %= sizeof(char);
              _vind[myid] += test;
            }
            for (int kkk = 0; kkk < _model.V_states().size(); ++kkk) {
              int test = _model.valid(_model.V_states()[kkk], nst);
              // If transition between states corresponding to row and column is possible, calculate the offdiagonal element.
              w[i] += test * _model.V_states()[kkk].value() * (1 - 2 * ((signs[_vind_byte[myid]] >> _vind_bit[myid]) & 1)) * v[col_ind[_vind[myid]]];
              _vind_bit[myid] += test;
              _vind_byte[myid] += _vind_bit[myid] / sizeof(char);
              _vind_bit[myid] %= sizeof(char);
              _vind[myid] += test;
            }
          }
        }
      }

      void reset(size_t sector_size) {
        if (sector_size > _max_dim) {
          std::stringstream s;
          s << "New sector request more memory than allocated. Increase MAX_DIM parameter. Requested " << sector_size << ", allocated " << _max_dim << ".";
          throw std::runtime_error(s.str().c_str());
        }
        if (sector_size * (_model.T_states().size() + _model.V_states().size()) > _max_size) {
          std::stringstream s;
          s << "New sector request more memory than allocated. Increase MAX_SIZE parameter. Requested " << sector_size * (_model.T_states().size() + _model.V_states().size()) << ", allocated " << _max_size << ".";
          throw std::runtime_error(s.str().c_str());
        }
        for(int myid = 0; myid < _nthreads; ++myid){
          _vind[myid] = 0;
          _vind_byte[myid] = 0;
          _vind_bit[myid] = 0;
        }
        n() = 0;
        ntot() = sector_size;
      }

      void fill() {
        _model.symmetry().init();
        reset(_model.symmetry().sector().size());
        int i = 0;
        long long k = 0;
        int isign = 0;
        // Size chunks equally.
        int step = (int)std::floor(_model.symmetry().sector().size() / _nthreads);
        for (int i = 0; i <= _nthreads; i++){
          _row_offset[i] = step * i;
        }
        // Put the rest into some of the first threads.
        int more = _model.symmetry().sector().size() - _row_offset[_nthreads];
        for (int i = 0; i < more; i++){
          _row_offset[i] += i;
        }
        for (int i = more; i <= _nthreads; i++){
          _row_offset[i] += more;
        }
        for(int myid = 0; myid <= _nthreads; ++myid){
          _vind_offset[myid] = (_model.T_states().size() + _model.V_states().size()) * _row_offset[myid];
        }
#pragma omp parallel for schedule(static, 1)
        for(int myid = 0; myid < _nthreads; ++myid){
          _vind[myid] = _vind_offset[myid];
          _vind_byte[myid] = _vind[myid] / sizeof(char);
          _vind_bit[myid] = _vind[myid] % sizeof(char);
          for(int i = _row_offset[myid]; i < _row_offset[myid + 1]; ++i){
            long long nst = _model.symmetry().state_by_index(i);
            // Compute diagonal element for current i state
            addDiagonal(i, _model.diagonal(nst), myid);
            // non-diagonal terms calculation
            off_diagonal<decltype(_model.T_states())>(nst, i, _model.T_states(), myid);
            off_diagonal<decltype(_model.V_states())>(nst, i, _model.V_states(), myid);
          }
        }
      }

      template<typename T_states>
      inline void off_diagonal(long long nst, int i, T_states& states, int chunk) {
        long long k = 0;
        int isign = 0;
        for (int kkk = 0; kkk < states.size(); ++kkk) {
          if (_model.valid(states[kkk], nst)) {
            _model.set(states[kkk], nst, k, isign);
            int k_index = _model.symmetry().index(k);
            addElement(i, k_index, states[kkk].value(), isign, chunk);
          }
        }
      };

      /**
       * Add diagonal H(i,i) element with value v.
       */
      void inline addDiagonal(const int &i, prec v, int chunk) {
        dvalues[i] = v;
        ++n();
        _vind_start[chunk] = _vind[chunk];
      }

      /**
       * Add off-diagonal H(i,j) element with value t (discarded here, restored in av) and Fermi sign.
       */
      void inline addElement(const int &i, int j, prec t, int sign, int chunk) {
        if (i == j) {
          throw std::logic_error("Attempt to use addElement() to add diagonal element. Use addDiagonal() instead!");
        }
        // It is an error to add the element (i, j) twice.
        for (size_t iii = _vind_start[chunk]; iii < _vind[chunk]; iii++) {
          if (col_ind[iii] == j) {
            throw std::logic_error("Collision. Check a, adag, numState, ninv_value!");
          }
        }
        // Store sign in CRS-like array, one bit per sign.
        col_ind[_vind[chunk]] = j;
        signs[_vind_byte[chunk]] &= ~(1ll << _vind_bit[chunk]);
        signs[_vind_byte[chunk]] |= sign < 0 ? 1ll << _vind_bit[chunk] : 0;
        ++_vind_bit[chunk];
        ++_vind[chunk];
        _vind_byte[chunk] += _vind_bit[chunk] / sizeof(char);
        _vind_bit[chunk] %= sizeof(char);
        if (_vind[chunk] > _max_size) {
// FIXME For OMP, check every chunk!
          std::stringstream s;
          s << "Current sector request more memory than allocated. Increase MAX_SIZE parameter. Requested " << _vind[chunk] << ", allocated " << _max_size << ".";
          throw std::runtime_error(s.str().c_str());
        }
      }

      void endMatrix() {
        // Nothing to do
      }

      void print() {
        // See: av().
        // Each row of the matrix is first restored from the arrays.
        std::vector < prec > line(n(), prec(0.0));
        _model.symmetry().init();
        std::cout << std::setprecision(2) << std::fixed;
        std::cout << "[";
        for (int myid = 0; myid < _nthreads; ++myid) {
          _vind[myid] = _vind_offset[myid];
          _vind_byte[myid] = _vind[myid] / sizeof(char);
          _vind_bit[myid] =  _vind[myid] % sizeof(char);
          for (int i = 0; i < n(); ++i) {
            _model.symmetry().next_state();
            long long nst = _model.symmetry().state();
            std::fill(line.begin(), line.end(), prec(0.0));
            line[i] = dvalues[i];
            for (int kkk = 0; kkk < _model.states().size(); ++kkk) {
              int test = _model.valid(_model.states()[kkk], nst);
              line[col_ind[_vind[myid]]] += test * _model.states()[kkk].value() * (1 - 2 * ((signs[_vind_byte[myid]] >> _vind_bit[myid]) & 1));
              _vind_bit[myid] += test;
              _vind_byte[myid] += _vind_bit[myid] / sizeof(char);
              _vind_bit[myid] %= sizeof(char);
              _vind[myid] += test;
            }
            std::cout << "[";
            for (int j = 0; j < n(); ++j) {
              std::cout << std::setw(6) << line[j] << (j == n() - 1 ? "" : ", ");
            }
            std::cout << "]" << (i == n() - 1 ? "" : ", \n");
          }
          std::cout << "]" << std::endl;
        }
      }

      virtual void zero_eigenapair() {
        Storage < prec >::eigenvalues().resize(1);
        Storage < prec >::eigenvalues()[0] = dvalues[0];
        Storage < prec >::eigenvectors().assign(1, std::vector < prec >(1, prec(1.0)));
      }

      size_t vector_size(typename Model::Sector sector) {
        return sector.size();
      }

      prec vv(const std::vector<prec> & v, const std::vector<prec> & w) {
        prec alf = prec(0.0);
        for (int k = 0; k < v.size(); ++k) {
          alf += w[k] * v[k];
        }
        return alf;
      }

#ifdef _OPENMP
      int &nprocs() { return _nthreads; }
#endif

    private:
      // Internal storage structure
      std::vector < prec > dvalues;
      std::vector < int > col_ind;
      std::vector < char > signs;

      // the maximum sizes of all the objects
      size_t _max_size;
      size_t _max_dim;

      // number of OMP threads, "1" for serial mode.
      int _nthreads;
      // OMP chunks offset
      std::vector < size_t > _row_offset;
      std::vector < size_t > _vind_offset;
      // internal indicies
      std::vector < size_t > _vind;
      // start of current row, used for checks
      std::vector < size_t > _vind_start;
      // bit and byte of the bitmap corresponding to CRS index
      std::vector < size_t > _vind_bit;
      std::vector < size_t > _vind_byte;

      // Hubbard model parameters
      Model &_model;

    };
  }
}
#endif //HUBBARD_SOCRSSTORAGE_H
