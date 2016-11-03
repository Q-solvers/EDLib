//
// Created by iskakoff on 23/08/16.
//

#ifndef HUBBARD_SPINRESOLVEDSTORAGE_H
#define HUBBARD_SPINRESOLVEDSTORAGE_H

#include <bitset>
#include <iomanip>
#include <type_traits>

#include "Storage.h"
#include "SzSymmetry.h"
#include "NSymmetry.h"

namespace EDLib {
  namespace Storage {

    template<typename prec, class Model>
    class SpinResolvedStorage : public Storage < prec > {
      static_assert(std::is_base_of<Symmetry::SzSymmetry, typename Model::SYMMETRY>::value, "Model have wrong symmetry.");
    public:
      using Storage < prec >::n;
      using Storage < prec >::ntot;
#ifdef USE_MPI
      using Storage < prec >::comm;
#endif
      using Storage < prec >::prepare_work_arrays;
      using Storage < prec >::finalize;

      template<typename p>
      class CRSMatrix {
      public:
        CRSMatrix() {
        }

        void init(size_t N) {
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

#ifdef USE_MPI
      SpinResolvedStorage(alps::params &p, Model &m, MPI_Comm comm) : Storage < prec >(p, comm), _comm(comm), _model(m), _interaction_size(m.interacting_orbitals()),
                                                   _loc_symmetry(m.interacting_orbitals()), _Ns(p["NSITES"]), _ms(p["NSPINS"]), _up_symmetry(int(p["NSITES"])), _down_symmetry(int(p["NSITES"])) {
        MPI_Comm_size(_comm, &_nprocs);
        MPI_Comm_rank(_comm, &_myid);
      }
#else
      SpinResolvedStorage(alps::params &p, Model &m) : Storage < prec >(p), _model(m), _interaction_size(m.interacting_orbitals()),
                                                   _Ns(p["NSITES"]), _ms(p["NSPINS"]), _up_symmetry(int(p["NSITES"])), _down_symmetry(int(p["NSITES"])),
                                                   _loc_symmetry(m.interacting_orbitals()) {}
#endif

      virtual void zero_eigenapair() {
        Storage < prec >::eigenvalues().resize(1);
        Storage < prec >::eigenvalues()[0] = _diagonal[0];
        Storage < prec >::eigenvectors().assign(1, std::vector < prec >(1, prec(1.0)));
      }

      virtual void av(prec *v, prec *w, int n, bool clear = true) {
#ifdef USE_MPI
        MPI_Win_fence(MPI_MODE_NOPRECEDE, _win);
        for(int i = 0; i<_procs.size(); ++i) {
          MPI_Get(&_vecval[_proc_offset[i]], _proc_size[i], alps::mpi::detail::mpi_type<prec>(), i, 0, _proc_size[i], alps::mpi::detail::mpi_type<prec>(), _win);
        }
#endif
        // Iteration over rows.
        for (int i = 0; i < n; ++i) {
          // Diagonal contribution.
          w[i] = _diagonal[i] * v[i] + (clear ? 0.0 : w[i]);
        }
        // Offdiagonal contribution.
        // iterate over up spin blocks
        for (int k = 0; k < _up_size; ++k) {
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
#ifdef USE_MPI
        MPI_Win_fence(MPI_MODE_NOSUCCEED | MPI_MODE_NOPUT | MPI_MODE_NOSTORE, _win);
#endif
        for (int i = 0; i < _up_size; ++i) {
          // Iteration over columns(unordered).
          for (int j = H_up.row_ptr()[i+_up_shift]; j < H_up.row_ptr()[i + _up_shift + 1]; ++j) {
            for (int k = 0; k < _down_symmetry.sector().size(); ++k) {
#ifdef USE_MPI
              w[i * _down_symmetry.sector().size() + k] += H_up.values()[j] * _vecval[H_up.col_ind()[j] * _down_symmetry.sector().size() + k];
#else
              w[i * _down_symmetry.sector().size() + k] += H_up.values()[j] * v[H_up.col_ind()[j] * _down_symmetry.sector().size() + k];
#endif
            }
          }
        }
      }

      void fill() {
        reset();
        if(n()==0) {
          // Do nothing if the matrix size is zero;
          return;
        }
        // fill off-diagonal matrix for each spin
        fill_spin(_up_symmetry, _Ns, H_up);
        fill_spin(_down_symmetry, 0, H_down);
        // fill local part;
        int isign;
        long long k;
        for(int i =0; i<_locsize; ++i) {
          _model.symmetry().next_state();
          long long nst = _model.symmetry().state();
          _diagonal[i] = _model.diagonal(nst);
//          for (int kkk = 0; kkk < _model.T_states().size(); ++kkk) {
//            if (_model.valid(_model.T_states()[kkk], nst)) {
//              _model.set(_model.V_states()[kkk], nst, k, isign);
//              int j = _model.symmetry().index(k);
//              H_loc[0].addElement(i, j, _model.T_states()[kkk].value(), isign);
//            }
//          }
        }
#ifdef USE_MPI
        find_neighbours();
#endif
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
              }
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
      }


      void reset() {
        _model.symmetry().init();
        const Symmetry::SzSymmetry &symmetry = static_cast<Symmetry::SzSymmetry>(_model.symmetry());
        const Symmetry::SzSymmetry::Sector &sector = symmetry.sector();
        _up_symmetry.set_sector(Symmetry::NSymmetry::Sector(sector.nup(), symmetry.comb().c_n_k(_Ns, sector.nup())));
        _down_symmetry.set_sector(Symmetry::NSymmetry::Sector(sector.ndown(), symmetry.comb().c_n_k(_Ns, sector.ndown())));
        size_t up_size = _up_symmetry.sector().size();
        size_t down_size = _down_symmetry.sector().size();
        H_up.init(up_size);
        H_down.init(down_size);
#ifdef USE_MPI
        MPI_Comm run_comm;
        int color = _myid < up_size ? 1 : MPI_UNDEFINED;
        MPI_Comm_split(_comm, color, _myid, &run_comm);
        if(color == 1) {
          _run_comm = run_comm;
          int myid;
          MPI_Comm_rank(_run_comm,&myid);
          int size;
          MPI_Comm_size(_run_comm,&size);
          int locsize = up_size / size;
          if ((up_size % size) > myid) {
            locsize += 1;
            _offset =  myid * locsize* _down_symmetry.sector().size();
          } else {
            _offset = (myid* locsize + (up_size % size))* _down_symmetry.sector().size();
          }
          _up_size = locsize;
          _up_shift = _offset / down_size;
          _locsize = locsize * down_size;
          _model.symmetry().set_offset(_offset);
          _procs.assign(size, 0);
          _proc_offset.assign(size, 0);
          _proc_size.assign(size, 0);
        } else {
          _up_size = 0;
          _locsize=0;
        }
#else
        _locsize = up_size*down_size;
        _up_size = up_size;
        _up_shift = 0;
#endif
        _diagonal.assign(_locsize, prec(0.0));
        n() = _locsize;
        ntot() = sector.size();
      }

      size_t vector_size(typename Model::Sector sector) {
        size_t sector_size = sector.size();
        int myid,size;
#ifdef USE_MPI
        MPI_Comm_rank(_comm,&myid);
        MPI_Comm_size(_comm,&size);
        size_t up_size = _model.symmetry().comb().c_n_k(_Ns, sector.nup());
        size_t down_size = sector_size / up_size;
        size = up_size>size ? size : up_size;
        if(myid >= size) {
          return 0;
        }
        size_t locsize = up_size / size;
        if ((up_size % size) > myid) {
          locsize += 1;
        }
        std::cout<<"Vec size:"<<locsize * down_size<<std::endl;
        return locsize * down_size;
#else
        return sector_size;
#endif
      }

      prec norm(const std::vector<prec>& vec) {

      }


      void a_adag(int i, const std::vector < prec > &invec, std::vector < prec > &outvec, const typename Model::Sector& next_sec, bool a) {
        size_t locsize = invec.size();
        size_t locsize_max = locsize;
        size_t next_size = next_sec.size();
        size_t up_size = _model.symmetry().comb().c_n_k(_Ns, next_sec.nup());
        size_t down_size = next_size / up_size;
        long long k;
        int sign;
        int ci;
        int cid;
#ifdef USE_MPI
        int myid;
        MPI_Comm_rank(_comm,&myid);
        int size;
        MPI_Comm_size(_comm,&size);
        int t = 0;
        bool fence;

        if(_up_symmetry.sector().size()%size != 0) {
          locsize_max+= _down_symmetry.sector().size();
        }
        std::vector<prec> buff(1000, 0.0);

        MPI_Win eigwin;
        MPI_Win_create(outvec.data(), sizeof(prec) * vector_size(next_sec), sizeof(prec), MPI_INFO_NULL, MPI_COMM_WORLD, &eigwin);
        MPI_Win_fence(MPI_MODE_NOPRECEDE,eigwin);
#endif
        for (int ind = 0; ind < locsize_max; ++ind) {
#ifdef USE_MPI
          if(fence)
            MPI_Win_fence(MPI_MODE_NOPRECEDE,eigwin);
          fence=false;
#endif
          if(ind<locsize) {
            _model.symmetry().next_state();
            long long nst = _model.symmetry().state();
            if (_model.checkState(nst, i, _model.max_total_electrons()) == (a ? 1 : 0)) {
              if (a) _model.a(i, nst, k, sign);
              else _model.adag(i, nst, k, sign);
              int i1 = _model.symmetry().index(k, next_sec);
#ifdef USE_MPI
              int size;
              MPI_Comm_size(_run_comm, &size);
              calcIndex(ci, cid, i1, up_size, down_size, size);
              if(myid == cid) {
                outvec[ci] = sign * invec[ind];
              } else {
                buff[t]=sign*invec[ind];
                MPI_Put(&buff[t],1,alps::mpi::detail::mpi_type<prec>(),cid,ci,1,alps::mpi::detail::mpi_type<prec>(), eigwin);
              }
#else
              outvec[i1] = sign * invec[ind];
#endif

            }
          }
#ifdef USE_MPI
          if((t+1)==buff.size()){fence=true;t=0;}
          if(fence)
            MPI_Win_fence(MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE,eigwin);
#endif
        }
#ifdef USE_MPI
        MPI_Win_fence(MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE, eigwin);
        MPI_Win_free(&eigwin);
#endif
      }


      prec vv(const std::vector<prec> & v, const std::vector<prec> & w) {
#ifdef USE_MPI
        prec alf = prec(0.0);
        prec temp = prec(0.0);
        for (int k = 0; k < v.size(); ++k) {
          temp += w[k] * v[k];
        }
        MPI_Allreduce(&temp, &alf, 1, alps::mpi::detail::mpi_type<prec>(), MPI_SUM, comm());
        return alf;
#else
        prec alf = prec(0.0);
        for (int k = 0; k < v.size(); ++k) {
          alf += w[k] * v[k];
        }
        return alf;
#endif
      }

#ifdef USE_MPI
      virtual void prepare_work_arrays(prec * data, size_t shift = 0) {
        MPI_Info info;
        MPI_Info_create( &info );
        MPI_Info_set( info, "no_locks", "true");
        MPI_Win_create(&data[shift], n() * sizeof(prec), sizeof(prec), info, _run_comm, &_win);
        MPI_Info_free(&info);
      }

      virtual MPI_Comm comm() {
        return _run_comm;
      }

      virtual void finalize(){
        MPI_Win_free(&_win);
        MPI_Comm run_comm = _run_comm;
        MPI_Comm_free(&run_comm);
        _run_comm = Storage<prec>::comm();
      }
#endif

    protected:

    private:
      std::vector < Matrix > H_loc;
      Model &_model;
      Matrix H_up;
      Matrix H_down;

      std::vector < prec > _diagonal;
      std::vector < prec > _vecval;

      Symmetry::NSymmetry _up_symmetry;
      Symmetry::NSymmetry _down_symmetry;
      Symmetry::NSymmetry _loc_symmetry;


      int _interaction_size;
      int _Ns;
      int _ms;


      size_t _up_size;
      size_t _up_shift;
      size_t _locsize;

#ifdef USE_MPI
      MPI_Comm _comm;
      MPI_Comm _run_comm;
      size_t _offset;
      int _nprocs;
      std::vector<int> _proc_offset;
      std::vector<int> _procs;
      std::vector<int> _proc_size;
      int _myid;
      MPI_Win _win;
#endif


#ifdef USE_MPI
      void find_neighbours() {
        int ci, cid;
        for(int i = 0; i< _up_size; ++ i) {
          for (int j = H_up.row_ptr()[i+_up_shift]; j < H_up.row_ptr()[i + _up_shift + 1]; ++j) {
            calcIndex(ci, cid, H_up.col_ind()[j]);
            if(_procs[cid]==0)  {_procs[cid]=1;}
          }
        }
        int oset = 0;
        int nprocs;
        MPI_Comm_size(_run_comm, &nprocs);
        for(int i=0; i < nprocs; i++) {
          if(_procs[i]!=0) {
            _procs[i]=1;
            _proc_offset[i]=oset * _down_symmetry.sector().size();
            int ls=_up_symmetry.sector().size()/nprocs;
            if((_up_symmetry.sector().size()% nprocs) > i) {
              ls++;
//              loc_offset[i] = i * ls;
            }else{
//              loc_offset[i] = i * ls + (Nstates[nup][ndo] % nprocs);
            }
            _proc_size[i]=ls * _down_symmetry.sector().size();
            oset+=ls;
          }
        }
        _vecval.assign(oset * _down_symmetry.sector().size(), prec(0.0));
//        std::cout<<_run_comm.rank()<<"  neighbours: ";
//        for(int i=0; i < nprocs; i++) {
//          std::cout<<" "<<i<<" "<<bool(_procs[i])<<" "<<_proc_size[i]<<" "<<_proc_offset[i]<<" ";
//        }
      }

      void calcIndex(int &ci, int &cid, int i) {
        int size;
        MPI_Comm_size(_run_comm, &size);
        calcIndex(ci, cid, i*_down_symmetry.sector().size(), _up_symmetry.sector().size(), _down_symmetry.sector().size(), size);
      }
      void calcIndex(int &ci, int &cid, int i, size_t u_s, size_t d_s, int nprocs) {
        //       local variables
        int tmp1, tmp2, tmp3, tmp4;
        int i_rest = i % d_s;
        int i_up = i / d_s;
        tmp1 = u_s / nprocs + 1;
        tmp2 = u_s % nprocs;
        tmp3 = u_s / nprocs;
        tmp4 = (i_up) - (tmp1 * tmp2);
        if (i_up > (tmp1 * tmp2)) {
          ci = ((tmp4%tmp3))*d_s + i_rest;
          cid = (i_up - tmp2) / tmp3;
        } else {
          ci = (i_up%(tmp3 + 1))*d_s + i_rest;
          cid = (i_up)/ (tmp3 + 1);
        }
      }
#endif

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
    };
  }
}
#endif //HUBBARD_SPINRESOLVEDSTORAGE_H
