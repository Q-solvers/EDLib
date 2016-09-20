//
// Created by iskakoff on 28/07/16.
//

#ifndef HUBBARD_SOCRSSTORAGE_H
#define HUBBARD_SOCRSSTORAGE_H

#include <vector>
#include <iomanip>
#include "Storage.h"

namespace EDLib {
  namespace Storage {

    template<typename prec, class Model>
    class SOCRSStorage : public Storage < prec > {
    public:
      using Storage < prec >::n;
#ifdef ALPS_HAVE_MPI
      SOCRSStorage(EDParams &p, Model &m, alps::mpi::communicator &comm) : Storage < prec >(p, comm),
#else
      SOCRSStorage(EDParams &p, Model &m) : Storage < prec >(p),
#endif
                                            _vind(0), _max_size(p["storage.MAX_SIZE"]),
                                            _max_dim(p["storage.MAX_DIM"]), _model(m) {
        /** init what you need from parameters*/
        col_ind.assign(_max_size, 0);
        signs.assign(_max_size, 0);
        dvalues.assign(_max_dim, prec(0.0));
      };

      virtual void av(prec *v, prec *w, int n, bool clear = true) {
        _model.symmetry().init();
        _vind = 0;
        _vind_byte = 0;
        _vind_bit = 0;
        // Iteration over rows.
        for (int i = 0; i < n; ++i) {
          _model.symmetry().next_state();
          long long nst = _model.symmetry().state();
          // Diagonal contribution.
          w[i] = dvalues[i] * v[i] + (clear ? 0.0 : w[i]);
          // Offdiagonal contribution.
          // Iteration over columns(unordered).
          for (int kkk = 0; kkk < _model.T_states().size(); ++kkk) {
            int test = _model.valid(_model.T_states()[kkk], nst);
            // If transition between states corresponding to row and column is possible, calculate the offdiagonal element.
            w[i] += test * _model.T_states()[kkk].value() * (1 - 2 * ((signs[_vind_byte] >> _vind_bit) & 1)) * v[col_ind[_vind]];
            _vind_bit += test;
            _vind_byte += _vind_bit / sizeof(char);
            _vind_bit %= sizeof(char);
            _vind += test;
          }
          for (int kkk = 0; kkk < _model.V_states().size(); ++kkk) {
            int test = _model.valid(_model.V_states()[kkk], nst);
            // If transition between states corresponding to row and column is possible, calculate the offdiagonal element.
            w[i] += test * _model.V_states()[kkk].value() * (1 - 2 * ((signs[_vind_byte] >> _vind_bit) & 1)) * v[col_ind[_vind]];
            _vind_bit += test;
            _vind_byte += _vind_bit / sizeof(char);
            _vind_bit %= sizeof(char);
            _vind += test;
          }
        }
      }

      void reset(size_t sector_size) {
        if (sector_size > _max_dim) {
          std::stringstream s;
          s << "Current sector request more memory than allocated. Increase MAX_DIM parameter. Requested " << sector_size << ", allocated " << _max_dim << ".";
          throw std::runtime_error(s.str().c_str());
        }
        _vind = 0;
        _vind_byte = 0;
        _vind_bit = 0;
        n() = 0;
      }

      void fill() {
        _model.symmetry().init();
        reset(_model.symmetry().sector().size());
        int i = 0;
        long long k = 0;
        int isign = 0;
        while (_model.symmetry().next_state()) {
          long long nst = _model.symmetry().state();
          // Compute diagonal element for current i state

          addDiagonal(i, _model.diagonal(nst));
          // non-diagonal terms calculation
          off_diagonal<decltype(_model.T_states())>(nst, i, _model.T_states());
          off_diagonal<decltype(_model.V_states())>(nst, i, _model.V_states());
          i++;
        }
      }
      template<typename T_states>
      inline void off_diagonal(long long nst, int i, T_states& states) {
        long long k = 0;
        int isign = 0;
        for (int kkk = 0; kkk < states.size(); ++kkk) {
          if (_model.valid(states[kkk], nst)) {
            _model.set(states[kkk], nst, k, isign);
            int k_index = _model.symmetry().index(k);
            addElement(i, k_index, states[kkk].value(), isign);
          }
        }
      };

      /**
       * Add diagonal H(i,i) element with value v.
       */
      void inline addDiagonal(const int &i, prec v) {
        dvalues[i] = v;
        ++n();
        _vind_start = _vind;
      }

      /**
       * Add off-diagonal H(i,j) element with value t (discarded here, restored in av) and Fermi sign.
       */
      void inline addElement(const int &i, int j, prec t, int sign) {
        if (i == j) {
          throw std::logic_error("Attempt to use addElement() to add diagonal element. Use addDiagonal() instead!");
        }
        // It is an error to add the element (i, j) twice.
        for (size_t iii = _vind_start; iii < _vind; iii++) {
          if (col_ind[iii] == j) {
            throw std::logic_error("Collision. Check a, adag, numState, ninv_value!");
          }
        }
        // Store sign in CRS-like array, one bit per sign.
        col_ind[_vind] = j;
        signs[_vind_byte] &= ~(1ll << _vind_bit);
        signs[_vind_byte] |= sign < 0 ? 1ll << _vind_bit : 0;
        ++_vind_bit;
        ++_vind;
        _vind_byte += _vind_bit / sizeof(char);
        _vind_bit %= sizeof(char);
        if (_vind > _max_size || _vind_byte > _max_size) {
          std::stringstream s;
          s << "Current sector request more memory than allocated. Increase MAX_SIZE parameter. Requested " << _vind << ", allocated " << _max_size << ".";
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
        _vind = 0;
        _vind_byte = 0;
        _vind_bit = 0;
        std::cout << std::setprecision(2) << std::fixed;
        std::cout << "[";
        for (int i = 0; i < n(); ++i) {
          _model.symmetry().next_state();
          long long nst = _model.symmetry().state();
          std::fill(line.begin(), line.end(), prec(0.0));
          line[i] = dvalues[i];
          for (int kkk = 0; kkk < _model.states().size(); ++kkk) {
            int test = _model.valid(_model.states()[kkk], nst);
            line[col_ind[_vind]] += test * _model.states()[kkk].value() * (1 - 2 * ((signs[_vind_byte] >> _vind_bit) & 1));
            _vind_bit += test;
            _vind_byte += _vind_bit / sizeof(char);
            _vind_bit %= sizeof(char);
            _vind += test;
          }
          std::cout << "[";
          for (int j = 0; j < n(); ++j) {
            std::cout << std::setw(6) << line[j] << (j == n() - 1 ? "" : ", ");
          }
          std::cout << "]" << (i == n() - 1 ? "" : ", \n");
        }
        std::cout << "]" << std::endl;
      }

      virtual void zero_eigenapair() {
        Storage < prec >::eigenvalues().resize(1);
        Storage < prec >::eigenvalues()[0] = dvalues[0];
        Storage < prec >::eigenvectors().assign(1, std::vector < prec >(1, prec(1.0)));
      }

    private:
      // Internal storage structure
      std::vector < prec > dvalues;
      std::vector < int > col_ind;
      std::vector < char > signs;

      // the maximum sizes of all the objects
      size_t _max_size;
      size_t _max_dim;

      // internal indicies
      size_t _vind;
      // start of current row, used for checks
      size_t _vind_start;
      // bit and byte of the bitmap corresponding to CRS index
      size_t _vind_bit;
      size_t _vind_byte;

      // Hubbard model parameters
      Model &_model;

    };
  }
}
#endif //HUBBARD_SOCRSSTORAGE_H
