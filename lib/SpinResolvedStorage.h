//
// Created by iskakoff on 23/08/16.
//

#ifndef HUBBARD_SPINRESOLVEDSTORAGE_H
#define HUBBARD_SPINRESOLVEDSTORAGE_H

#include <bitset>
#include <iomanip>
#include <boost/typeof/typeof.hpp>

#include "Storage.h"
#include <SzSymmetry.h>
#include <NSymmetry.h>

namespace EDLib {
  namespace Storage {

    template<typename prec, class Model>
    class SpinResolvedStorage : public Storage < prec > {
    BOOST_STATIC_ASSERT(boost::is_base_of<Symmetry::SzSymmetry, typename Model::SYMMETRY>::value);
    public:
      using Storage < prec >::n;

      template<typename p>
      class CRSMatrix {
      public:
        CRSMatrix() {
        }

        void init(int N) {
          _values.assign(N * N, p(0));
          _col_ind.assign(N * N, 0);
          _row_ptr.assign(N + 1, 0);
          _vind = 0;
        }

        void inline addElement(int i, int j, prec t, int sign) {
          if (i == j) {
            throw std::logic_error("Attempt to use addElement() to add diagonal element. Use addDiagonal() instead!");
          }
          if (std::abs(t) == 0) {
            return;
          }
          // check that there is no any data on the k state
          for (int iii = _row_ptr[i]; iii < _vind; ++iii) {
            if (_col_ind[iii] == j) {
              throw std::logic_error("Collision. Check a, adag, numState, ninv_value!");
            }
          }
          // create new element in CRS arrays
          _col_ind[_vind] = j;
          _values[_vind] = sign * t;
          ++_vind;
        }

        void inline endLine(int i) {
          _row_ptr[i + 1] = _vind;
        }

        std::vector < int > &row_ptr() {
          return _row_ptr;
        };

        std::vector < int > &col_ind() {
          return _col_ind;
        }

        std::vector < prec > &values() {
          return _values;
        }

      private:
        std::vector < prec > _values;
        std::vector < int > _row_ptr;
        std::vector < int > _col_ind;
        int _vind;
      };

      typedef CRSMatrix < prec > Matrix;

      SpinResolvedStorage(EDParams &p, Model &m) : Storage < prec >(p), _model(m), _interaction_size(m.interacting_orbitals()),
                                                   _Ns(p["NSITES"]), _ms(p["NSPINS"]), _up_symmetry(int(p["NSITES"])), _down_symmetry(int(p["NSITES"])),
                                                   _loc_symmetry(m.interacting_orbitals()) {
        long long max_up = 1ll<<_interaction_size;
        // iteration over down spin
//        while(_loc_symmetry.next_sector()) {
//          _loc_symmetry.init();
//          while(_loc_symmetry.next_state()) {
//            long long nst_dn = _loc_symmetry.state();
//            for(long long up = 0; up< max_up; ++up) {
//              long long nst =
//              for (int kkk = 0; kkk < _model.T_states().size(); ++kkk) {
//                if (_model.valid(_model.T_states()[kkk], nst)) {
//                  _model.set(_model.T_states()[kkk], nst, k, isign);
//                  int j = spin_symmetry.index(k >> shift);
//                  spin_matrix.addElement(i, j, _model.T_states()[kkk].value(), isign);
//                }
//              }
//            }
//          }
//        }
      }

      virtual void zero_eigenapair() {
        Storage < prec >::eigenvalues().resize(1);
        Storage < prec >::eigenvalues()[0] = _diagonal[0];
        Storage < prec >::eigenvectors().assign(1, std::vector < prec >(1, prec(1.0)));
      }

      virtual void av(prec *v, prec *w, int n, bool clear = true) {
        // Iteration over rows.
        for (int i = 0; i < n; ++i) {
          // Diagonal contribution.
          w[i] = _diagonal[i] * v[i] + (clear ? 0.0 : w[i]);
        }

        // Interaction part
        for(int i =0; i<n; i+=_up_symmetry.sector().size()) {

        }

        // Offdiagonal contribution.
        // iterate over up spin blocks
        for (int k = 0; k < _up_symmetry.sector().size(); ++k) {
          // Iteration over rows.
          for (int i = 0; i < _down_symmetry.sector().size(); ++i) {
            // Iteration over columns(unordered).
            for (int j = H_down.row_ptr()[i]; j < H_down.row_ptr()[i + 1]; ++j) {
              w[i + k * _down_symmetry.sector().size()] += H_down.values()[j] * v[H_down.col_ind()[j] + k * _down_symmetry.sector().size()];
            }
          }
        }
        //
        // Iteration over rows.
        for (int i = 0; i < _up_symmetry.sector().size(); ++i) {
          // Iteration over columns(unordered).
          for (int j = H_up.row_ptr()[i]; j < H_up.row_ptr()[i + 1]; ++j) {
            for (int k = 0; k < _down_symmetry.sector().size(); ++k) {
              w[i * _down_symmetry.sector().size() + k] += H_up.values()[j] * v[H_up.col_ind()[j] * _down_symmetry.sector().size() + k];
            }
          }
        }
      }

      void reset() {
        const Symmetry::SzSymmetry &symmetry = static_cast<Symmetry::SzSymmetry>(_model.symmetry());
        const Symmetry::SzSymmetry::Sector &sector = symmetry.sector();
        _up_symmetry.set_sector(Symmetry::NSymmetry::Sector(sector.nup(), symmetry.comb().c_n_k(_Ns, sector.nup())));
        _down_symmetry.set_sector(Symmetry::NSymmetry::Sector(sector.ndown(), symmetry.comb().c_n_k(_Ns, sector.ndown())));
        H_up.init(_up_symmetry.sector().size());
        H_down.init(_down_symmetry.sector().size());
#ifdef ALPS_HAVE_MPI
        alps::mpi::communicator comm;
        _nprocs = comm.size();
        _myid = comm.rank();
        _locsize = sector.size()/_nprocs;
        _diagonal.assign(_locsize, prec(0.0));
#else
        _diagonal.assign(sector.size(), prec(0.0));
#endif
      }

      void fill() {
        _model.symmetry().init();
        reset();
        // fill interaction

        // fill off-diagonal matrix for each spin
        fill_spin(_up_symmetry, _Ns, H_up);
        fill_spin(_down_symmetry, 0, H_down);
        // fill diagonal;
        int i = 0;
        while (_model.symmetry().next_state()) {
          long long nst = _model.symmetry().state();
          _diagonal[i] = _model.diagonal(nst);
          ++i;
        }
        n() = _model.symmetry().sector().size();
      }

      void fill_spin(Symmetry::NSymmetry &spin_symmetry, int shift, Matrix &spin_matrix) {
        long long k = 0;
        int isign = 0;
        int i = 0;
        while (spin_symmetry.next_state()) {
          long long nst = spin_symmetry.state();
          for (int kkk = 0; kkk < _model.T_states().size(); ++kkk) {
            if (_model.valid(_model.T_states()[kkk], nst << shift)) {
              _model.set(_model.T_states()[kkk], nst << shift, k, isign);
              int j = spin_symmetry.index(k >> shift);
              spin_matrix.addElement(i, j, _model.T_states()[kkk].value(), isign);
            }
          }
          spin_matrix.endLine(i);
          ++i;
        }
      }

      void print() {
        std::cout << std::setprecision(2) << std::fixed;
        std::cout << "{";
        for (int i = 0; i < _up_symmetry.sector().size(); ++i) {
          std::cout << "{";
          for (int j = 0; j < _up_symmetry.sector().size(); ++j) {
            bool f = true;
            for (int k = H_up.row_ptr()[i]; k < H_up.row_ptr()[i + 1]; ++k) {
              if ((H_up.col_ind()[k]) == j) {
                std::cout << std::setw(6) << H_up.values()[k] << (j == _up_symmetry.sector().size() - 1 ? "" : ", ");
                f = false;
              } /*else {
            std::cout<<"0.0 ";
          }*/
            }
            if (f) {
              std::cout << std::setw(6) << 0.0 << (j == _up_symmetry.sector().size() - 1 ? "" : ", ");
            }
          }
          std::cout << "}" << (i == _up_symmetry.sector().size() - 1 ? "" : ", \n");
        }
        std::cout << "}" << std::endl;
        std::cout << "\n\n{";
        for (int i = 0; i < _down_symmetry.sector().size(); ++i) {
          std::cout << "{";
          for (int j = 0; j < _down_symmetry.sector().size(); ++j) {
            bool f = true;
            for (int k = H_down.row_ptr()[i]; k < H_down.row_ptr()[i + 1]; ++k) {
              if ((H_down.col_ind()[k]) == j) {
                std::cout << std::setw(6) << H_down.values()[k] << (j == _down_symmetry.sector().size() - 1 ? "" : ", ");
                f = false;
              } /*else {
            std::cout<<"0.0 ";
          }*/
            }
            if (f) {
              std::cout << std::setw(6) << 0.0 << (j == _down_symmetry.sector().size() - 1 ? "" : ", ");
            }
          }
          std::cout << "}" << (i == _down_symmetry.sector().size() - 1 ? "" : ", \n");
        }
        std::cout << "}" << std::endl;
      }

      void endMatrix() {
#ifdef ALPS_HAVE_MPI

#endif
      }


    private:
      std::vector < Matrix > H_loc;
      Model &_model;
      Matrix H_up;
      Matrix H_down;

      std::vector < prec > _diagonal;

      Symmetry::NSymmetry _up_symmetry;
      Symmetry::NSymmetry _down_symmetry;
      Symmetry::NSymmetry _loc_symmetry;


      int _interaction_size;
      int _Ns;
      int _ms;

#ifdef ALPS_HAVE_MPI
      size_t _locsize;
      int _nprocs;
      int _myid;
#endif
    };

  }
}
#endif //HUBBARD_SPINRESOLVEDSTORAGE_H
