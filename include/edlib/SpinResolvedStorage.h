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
#include "CRSMatrix.h"

namespace EDLib {
  namespace Storage {

    template<class ModelType>
    class SpinResolvedStorage : public Storage < typename ModelType::precision > {
      static_assert(std::is_base_of<Symmetry::SzSymmetry, typename ModelType::SYMMETRY>::value, "Model have wrong symmetry.");
    public:
      typedef ModelType Model;
      typedef typename Model::precision prec;
    public:
      using Storage < prec >::n;
      using Storage < prec >::ntot;
#ifdef USE_MPI
      using Storage < prec >::comm;
      using Storage < prec >::broadcast_evals;
#endif
      using Storage < prec >::prepare_work_arrays;
      using Storage < prec >::finalize;



      typedef CRSMatrix < prec > Matrix;

#ifdef USE_MPI
      SpinResolvedStorage(alps::params &p, Model &m, MPI_Comm comm) : Storage < prec >(p, comm), _comm(comm), _run_comm(MPI_COMM_NULL), _model(m),_interaction_size(m.interacting_orbitals()),
                                                                      _Ns(p["NSITES"].as<int>()), _ms(p["NSPINS"].as<int>()), _up_symmetry(p["NSITES"].as<int>()),
                                                                      _down_symmetry(p["NSITES"].as<int>()) {
        MPI_Comm_size(_comm, &_nprocs);
        MPI_Comm_rank(_comm, &_myid);
      }
#else
      SpinResolvedStorage(alps::params &p, Model &m) : Storage < prec >(p), _model(m), _interaction_size(m.interacting_orbitals()),
                                                   _Ns(p["NSITES"]), _ms(p["NSPINS"]), _up_symmetry(int(p["NSITES"])), _down_symmetry(int(p["NSITES"])) {}
#endif

      virtual void zero_eigenapair() {
        Storage < prec >::eigenvalues().resize(1);
        Storage < prec >::eigenvalues()[0] = _diagonal[0];
        Storage < prec >::eigenvectors().assign(1, std::vector < prec >(1, prec(1.0)));
      }

      virtual void av(prec *v, prec *w, int n, bool clear = true) {
#ifdef USE_MPI
        // Initialize inter-processor communications
        // we collect all data from the remote processes into _vecval array
        MPI_Win_fence(MPI_MODE_NOPRECEDE, _win);
        for(int i = 0; i<_procs.size(); ++i) {
          if(_procs[i]!=0)
          MPI_Get(&_vecval[_proc_offset[i]], _proc_size[i], alps::mpi::detail::mpi_type<prec>(), i, _loc_min[i], _proc_size[i], alps::mpi::detail::mpi_type<prec>(), _win);
        }
#endif
        // Iteration over diagonal contribution.
        for (int i = 0; i < n; ++i) {
          // Diagonal contribution.
          w[i] = _diagonal[i] * v[i] + (clear ? 0.0 : w[i]);
        }
        // Off-diagonal contribution. Process spin-down contribution for each spin-up block
        // iterate over spin-up blocks
        for (int k = 0; k < _up_size; ++k) {
          // Iteration over rows.
          for (int i = 0; i < _down_symmetry.sector().size(); ++i) {
            // Iteration over columns.
            for (int j = H_down.row_ptr()[i]; j < H_down.row_ptr()[i + 1]; ++j) {
              w[i + k * _down_symmetry.sector().size()] += H_down.values()[j] * v[H_down.col_ind()[j] + k * _down_symmetry.sector().size()];
            }
          }
        }
#ifdef USE_MPI
        // Waiting for the data to be received
        MPI_Win_fence(MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE, _win);
#endif
        // Process spin-up hopping contribution
        // Iteration over rows.
        for (int i = 0; i < _up_size; ++i) {
          // Iteration over columns.
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

        // Off-diagonal interaction contribution
        // Check that we have off-diagonal interaction elements
        if(H_loc.row_ptr().size()!=0)
        for (size_t i = _int_start; i < n; ++i) {
          for (int j = H_loc.row_ptr()[i]; j < H_loc.row_ptr()[i + 1]; ++j) {
#ifdef USE_MPI
            w[i] += H_loc.values()[j] * _vecval[H_loc.col_ind()[j]];
#else
            w[i] += H_loc.values()[j] * v[H_loc.col_ind()[j]];
#endif
          }
        }
      }

      /**
       * Fill Hamiltonain matrix for current symmetry sector
       */
      void fill() {
        reset();
        if(n()==0) {
          // Do nothing if the matrix size is zero;
          return;
        }
        // Hopping term
        // fill off-diagonal matrix for each spin
        fill_spin(_up_symmetry, _Ns, H_up);
        fill_spin(_down_symmetry, 0, H_down);
        // fill local part;
        int isign;
        long long k;
        _int_start = _locsize;
        for(size_t i =0; i<_locsize; ++i) {
          _model.symmetry().next_state();
          long long nst = _model.symmetry().state();
          // add diagonal contribution
          _diagonal[i] = _model.diagonal(nst);
          // Add off-diagonal contribution from interaction term
          if(_model.V_states().size() > 0) {
            for (int kkk = 0; kkk < _model.V_states().size(); ++kkk) {
              if (_model.valid(_model.V_states()[kkk], nst)) {
                _int_start = std::min(i, _int_start);
                _model.set(_model.V_states()[kkk], nst, k, isign);
                int j = _model.symmetry().index(k);
                H_loc.addElement(i, j, _model.V_states()[kkk].value(), isign);
              }
            }
            H_loc.endLine(i);
          }
        }
#ifdef USE_MPI
        find_neighbours();
#endif
      }

      void print() {
        // nothing to do
      }

      void init() {
        _model.symmetry().init();
#ifdef USE_MPI
        _model.symmetry().set_offset(_offset);
#endif
      }

      /**
       * Reset storage and symmetry object for the current symmetry sector.
       * Update MPI communicator if necessary, set local dimensions size, setup working arrays size.
       */
      void reset(int t = 0) {
        _model.symmetry().init();
        // For spin-resolved storage we should guarantee that model is spin-symmetric.
        const Symmetry::SzSymmetry &symmetry = static_cast<Symmetry::SzSymmetry>(_model.symmetry());
        // get current symmetry sector
        const Symmetry::SzSymmetry::Sector &sector = symmetry.sector();
        // get symmetries for each spin
        _up_symmetry.set_sector(Symmetry::NSymmetry::Sector(sector.nup(), symmetry.comb().c_n_k(_Ns, sector.nup())));
        _down_symmetry.set_sector(Symmetry::NSymmetry::Sector(sector.ndown(), symmetry.comb().c_n_k(_Ns, sector.ndown())));
        // get each spin dimensions
        size_t up_size = _up_symmetry.sector().size();
        size_t down_size = _down_symmetry.sector().size();
        // init hopping Hamiltonians for each spin
        H_up.init(up_size, 100);
        H_down.init(down_size, 100);
#ifdef USE_MPI
        // release communicator if needed
        if(_comm != _run_comm && _run_comm != MPI_COMM_NULL) {
          MPI_Comm_free(&_run_comm);
        };
        // communicator to be used during diagonalization
        // check that there is data for the current CPU
        int color = _myid < up_size ? 1 : 0;//MPI_UNDEFINED;
        // Create new MPI communicator for the processors with defined color
        MPI_Comm_split(_comm, color, (color == 1 ? _myid : 0), &_run_comm);
        // update working communicator
        if(color == 1) {
          // there is data for current CPU
          // get CPU rank and size for recently created working communicator
          int myid;
          MPI_Comm_rank(_run_comm,&myid);
          int size;
          MPI_Comm_size(_run_comm,&size);
          // compute the size of local arrays and the offset from the beginning
          int locsize = up_size / size;
          if ((up_size % size) > myid) {
            locsize += 1;
            _offset =  myid * locsize* _down_symmetry.sector().size();
          } else {
            _offset = (myid* locsize + (up_size % size))* _down_symmetry.sector().size();
          }
          // local dimension for the spin-up hopping Hamiltonian matrix
          _up_size = locsize;
          // offet in the spin-up channel
          _up_shift = _offset / down_size;
          // local dimension for the whole Hamiltonian matrix
          _locsize = locsize * down_size;
          // apply offset to the symmetry object to generate proper configuration state
          _model.symmetry().set_offset(_offset);
          // inter processor communactions
          // array with flags for each processor
          _procs.assign(size, 0);
          // offset of each processor in the vecval array
          _proc_offset.assign(size, 0);
          // data amount to be received
          _proc_size.assign(size, 0);
          // index of the first element to be received
          _loc_min.assign(size, 0);
        } else {
          _up_size = 0;
          _locsize=0;
        }
#else
        _locsize = up_size*down_size;
        _up_size = up_size;
        _up_shift = 0;
#endif
        // allocate memory for local Hamiltonian
        // density-density contribution
        _diagonal.assign(_locsize, prec(0.0));
        // off-diagonal contribution
        if(_model.V_states().size()>0) {
          H_loc.init(_locsize, 3);
        }
        // local dimension of the Hamiltonian matrix
        n() = _locsize;
        // total dimension of the Hamiltonian matrix
        ntot() = sector.size();
      }

      /**
       * Compute local dimension for the specific sector
       * @param sector -- symmetry sector to compute dimension
       * @return size of local vector for the specific symmetry sector
       */
      size_t vector_size(typename Model::Sector sector) {
        // get the total dimension for the current symmetry sector
        size_t sector_size = sector.size();
#ifdef USE_MPI
        int myid,size;
        // get rank and size for current CPU in the global communicator
        MPI_Comm_rank(_comm,&myid);
        MPI_Comm_size(_comm,&size);
        size_t up_size = _model.symmetry().comb().c_n_k(_Ns, sector.nup());
        size_t down_size = sector_size / up_size;
        size = up_size>size ? size : up_size;
        // there is no data for current CPU
        if(myid >= size) {
          return 0;
        }
        // compute local dimension for the spin-up channel
        size_t locsize = up_size / size;
        if ((up_size % size) > myid) {
          locsize += 1;
        }
        // compute and return the total dimension
        return locsize * down_size;
#else
        return sector_size;
#endif
      }

      /**
       * Perform a_i or a_i^* operation on |invec>
       *
       * @param i -- position to destroy(create) electron
       * @param invec -- input vector
       * @param outvec -- output vector
       * @param next_sec -- the resulting Symmetry sector
       * @param a -- destroy particle if true, create otherwise
       */
      void a_adag(int i, const std::vector < prec > &invec, std::vector < prec > &outvec, const typename Model::Sector& next_sec, bool a) {
        // local dimension of current vector
        size_t locsize = invec.size();
        // maximal local dimension of current vector
        size_t locsize_max = locsize;
        // dimension of the resulting symmetry sector
        size_t next_size = next_sec.size();
        // dimension of spin-up channel Hamiltonian matrix
        size_t up_size = _model.symmetry().comb().c_n_k(_Ns, next_sec.nup());
        // dimension of spin-down channel Hamiltonian matrix
        size_t down_size = next_size / up_size;
        long long k;
        int sign;
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &locsize_max, 1, alps::mpi::detail::mpi_type<size_t>(), MPI_MAX, _comm);
        int ci;
        int cid;
        int myid;
        MPI_Comm_rank(_comm,&myid);
        int size;
        MPI_Comm_size(_comm,&size);
        int t = 0;
        // synchronization flag
        bool fence = false;
        // adjust maximal local dimension
        if(_up_symmetry.sector().size()%size != 0) {
          locsize_max+= _down_symmetry.sector().size();
        }
        // communication buffer
        std::vector<prec> buff(1000, 0.0);
        // communication window
        int rank;
 	MPI_Comm_rank(_run_comm, &rank);
        MPI_Comm_size(_run_comm, &size);
        size = outvec.size() > 0 ? 1 : 0;
        MPI_Allreduce(MPI_IN_PLACE, &size, 1, alps::mpi::detail::mpi_type<int>(), MPI_SUM, _comm);
        MPI_Win eigwin;
        MPI_Win_create(outvec.data(), sizeof(prec) * outvec.size(), sizeof(prec), MPI_INFO_NULL, _comm, &eigwin);
        MPI_Win_fence(MPI_MODE_NOPRECEDE,eigwin);
#endif
        // iterate over local part of vector
        for (int ind = 0; ind < locsize_max; ++ind) {
#ifdef USE_MPI
          // perfrom one-sided communication
          if(fence)
            MPI_Win_fence(MPI_MODE_NOPRECEDE,eigwin);
          fence=false;
#endif
          // try to destroy (create) particle if index is within boundary
          if(ind<locsize) {
            _model.symmetry().next_state();
            long long nst = _model.symmetry().state();
            // check that particle can be destroyed (created)
            if (_model.checkState(nst, i, _model.max_total_electrons()) == (a ? 1 : 0)) {
              if (a) _model.a(i, nst, k, sign);
              else _model.adag(i, nst, k, sign);
              // compute index of the new state
              int i1 = _model.symmetry().index(k, next_sec);
#ifdef USE_MPI
              // compute CPU id and local index of new state
              calcIndex(ci, cid, i1, up_size, down_size, size);
              // avoid intra-CPU MPI communications
              if(myid == cid) {
                outvec[ci] = sign * invec[ind];
              } else {
                // initiate one-sided communication
                buff[t]=sign*invec[ind];
                MPI_Put(&buff[t],1,alps::mpi::detail::mpi_type<prec>(),cid,ci,1,alps::mpi::detail::mpi_type<prec>(), eigwin);
              }
#else
              /// update array
              outvec[i1] = sign * invec[ind];
#endif

            }
          }
#ifdef USE_MPI
          // check buffer boundary
          if((++t)==buff.size()){fence=true;t=0;}
          // synchronize if necessary
          if(fence)
            MPI_Win_fence(MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE,eigwin);
#endif
        }
#ifdef USE_MPI
        if(!fence) MPI_Win_fence(MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE, eigwin);
        MPI_Win_free(&eigwin);
#endif
      }


      /**
       * Compute <v|w> product
       * @param v - bra-state
       * @param w - ket-state
       * @return <v|w> product
       */
      prec vv(const std::vector<prec> & v, const std::vector<prec> & w) {
#ifdef USE_MPI
        return vv(v, w, comm());
#else
        prec alf = prec(0.0);
        for (int k = 0; k < v.size(); ++k) {
          alf += w[k] * v[k];
        }
        return alf;
#endif
      }

#ifdef USE_MPI
      prec vv(const std::vector<prec> & v, const std::vector<prec> & w, MPI_Comm com) {
        prec alf = prec(0.0);
        prec temp = prec(0.0);
        for (int k = 0; k < v.size(); ++k) {
          temp += w[k] * v[k];
        }
        MPI_Allreduce(&temp, &alf, 1, alps::mpi::detail::mpi_type<prec>(), MPI_SUM, com);
        return alf;
      }

      /**
       * Initialize the communication window for the *data object from the specific offset
       * @param data -- input array
       * @param shift -- offset in the input array
       */
      virtual void prepare_work_arrays(prec * data, size_t shift = 0) {
//        if(MPI_WIN_NULL != _win){
//          // TODO: handle already allocated window
//        }
        MPI_Info info;
        MPI_Info_create( &info );
        MPI_Info_set( info, (char *) "no_locks", (char *) "true");
        MPI_Win_create(&data[shift], n() * sizeof(prec), sizeof(prec), MPI_INFO_NULL, _run_comm, &_win);
        MPI_Info_free(&info);
      }

      /**
       *
       * @return current working communicator
       */
      virtual MPI_Comm comm() {
        return _run_comm;
      }

      /**
       * Finalize current MPI execution. Broadcast eigenvalues if necessary. Release window and communicator.
       * Set the global MPI communicator as the current MPI communicator.
       */
      virtual int finalize(int info, bool bcast = true, bool empty = true) {
        MPI_Bcast(&info,1, MPI_INT, 0, Storage<prec>::comm());
        if(info>=0) {
          if(bcast) broadcast_evals(empty);
        }
        if(ntot() > 1 && n() > 0) {
          MPI_Win_free(&_win);
        }
        if(_run_comm != MPI_COMM_WORLD && _run_comm != MPI_COMM_NULL) {
          MPI_Comm_free(&_run_comm);
        }
      	_run_comm = Storage < prec >::comm();
        return info;
      }

      /**
       * @return offset from 0 for current CPU
       */
      size_t offset(){
        return _offset;
      }
#endif

    private:
      /// Current model
      Model &_model;
      /// Off-diagonal part of local Hamiltonian
      Matrix H_loc;
      /// Spin-up hopping
      Matrix H_up;
      /// Spin-down hopping
      Matrix H_down;

      /// diagonal part
      std::vector < prec > _diagonal;
      /// array to store remote processes communication data
      std::vector < prec > _vecval;

      /// symmetries for each spin
      Symmetry::NSymmetry _up_symmetry;
      Symmetry::NSymmetry _down_symmetry;

      ///
      int _interaction_size;
      /// The total number of electon
      int _Ns;
      /// Total number of spins
      int _ms;

      /// size of H_up
      size_t _up_size;
      /// offset in spin_up channel
      size_t _up_shift;
      /// local size on each CPU
      size_t _locsize;
      /// staring index in interaction Hamiltonian
      size_t _int_start;

#ifdef USE_MPI
      /// global communicator
      MPI_Comm _comm;
      /// current working communicator
      MPI_Comm _run_comm;
      /// offset for the current CPU
      size_t _offset;
      int _myid;
      int _nprocs;
      /// Inter-process communication auxiliary arrays
      std::vector<int> _proc_offset;
      std::vector<int> _procs;
      std::vector<int> _loc_min;
      std::vector<int> _proc_size;
      /// MPI communication window
      MPI_Win _win;

      /**
       * Find neighbour CPUs for the current Hamiltonian matrix
       */
      void find_neighbours() {
        int ci, cid;
        // size of the working communicator
        int nprocs;
        MPI_Comm_size(_run_comm, &nprocs);
        std::vector<int> loc_offset(nprocs, 0);
        // Find smallest and largest index in the current Hamiltonian
        std::vector<int> l_loc_max(_loc_min.size(), INT_MIN);
        std::vector<int> l_loc_min(_loc_min.size(), INT_MAX);
        // For the spin-up channel
        for(int i = 0; i< _up_size; ++ i) {
          for (int j = H_up.row_ptr()[i+_up_shift]; j < H_up.row_ptr()[i + _up_shift + 1]; ++j) {
            calcIndex(ci, cid, H_up.col_ind()[j]*_down_symmetry.sector().size(), _up_symmetry.sector().size(), _down_symmetry.sector().size(), nprocs);
            l_loc_max[cid] = std::max(ci, l_loc_max[cid]);
            l_loc_min[cid] = std::min(ci, l_loc_min[cid]);
            if(_procs[cid]==0)  {_procs[cid]=1;}
          }
        }
        // For the off-diagonal interaction term
        if(H_loc.row_ptr().size()!=0) {
          for (size_t i = _int_start; i < _locsize; ++i) {
            for (int j = H_loc.row_ptr()[i]; j < H_loc.row_ptr()[i + 1]; ++j) {
              calcIndex(ci, cid, _down_symmetry.sector().size()*(H_loc.col_ind()[j]/_down_symmetry.sector().size()), _up_symmetry.sector().size(), _down_symmetry.sector().size(), nprocs);
              l_loc_max[cid] = std::max(ci, l_loc_max[cid]);
              l_loc_min[cid] = std::min(ci, l_loc_min[cid]);
              if(_procs[cid]==0)  {_procs[cid]=1;}
            }
          }
        }
        int oset = 0;
        for(int i=0; i < nprocs; i++) {
          if(_procs[i]) {
            _procs[i]=1;
            // calculate offset for i-th CPU
            _proc_offset[i]=oset * _down_symmetry.sector().size() + l_loc_min[i];
            // The index of the first element of the vector to be received from i-th CPU
            _loc_min[i] = l_loc_min[i];
            int ls=_up_symmetry.sector().size()/nprocs;
            if((_up_symmetry.sector().size()% nprocs) > i) {
              ls++;
              loc_offset[i] = (i * ls) - oset;
            }else{
              loc_offset[i] = i * ls + (_up_symmetry.sector().size() % nprocs) - oset;
            }
            // number of elements to be received from the i-th CPU
            _proc_size[i]= l_loc_max[i] - l_loc_min[i] + _down_symmetry.sector().size();
            oset+=ls;
          }
        }
        // alloacte memory for the working array
        _vecval.assign(oset * _down_symmetry.sector().size(), prec(0.0));
        // adjust indexes
        for (int i = 0; i < _up_size; ++i) {
          for (int j = H_up.row_ptr()[i + _up_shift]; j < H_up.row_ptr()[i + _up_shift + 1]; ++j) {
            calcIndex(ci, cid, H_up.col_ind()[j]*_down_symmetry.sector().size(), _up_symmetry.sector().size(), _down_symmetry.sector().size(), nprocs);
            H_up.col_ind()[j] -= loc_offset[cid];
          }
        }
        if(H_loc.row_ptr().size()!=0) {
          for (size_t i = _int_start; i < _locsize; ++i) {
            for (int j = H_loc.row_ptr()[i]; j < H_loc.row_ptr()[i + 1]; ++j) {
              calcIndex(ci, cid, H_loc.col_ind()[j], _up_symmetry.sector().size(), _down_symmetry.sector().size(), nprocs);
              H_loc.col_ind()[j] -= loc_offset[cid]*_down_symmetry.sector().size();
            }
          }
        }
      }

      /// Calculate local index, ci, and CPU id, cid, for the global index i
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

      /**
       * Fill the Hamiltonian matrix for the specific spin
       * @param spin_symmetry -- current
       * @param shift -- spin shift in the configuration state (0 for spin-up, and Ns for spin-down)
       * @param spin_matrix -- matrix to be filled
       */
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
