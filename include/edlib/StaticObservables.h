#ifndef HUBBARD_STATEDESCRIPTION_H
#define HUBBARD_STATEDESCRIPTION_H

#include "EigenPair.h"
#include <vector>
#include <bitset>
#include <string>

#include <alps/params.hpp>

namespace EDLib {
  template<class Hamiltonian>
  class StaticObservables {
  protected:
    typedef typename Hamiltonian::ModelType::precision precision;
    typedef typename Hamiltonian::ModelType::Sector sector;

  public:

    const static std::string _N_;
    const static std::string _N_UP_;
    const static std::string _N_DN_;
    const static std::string _M_;
    const static std::string _D_OCC_;

    StaticObservables(alps::params &p) :
      _beta(p["lanc.BETA"].as<precision>()),
      _cutoff(p["lanc.BOLTZMANN_CUTOFF"])
    {
      if(p["storage.EIGENVALUES_ONLY"] == 1) {
        throw std::logic_error("Eigenvectors have not been computed. StateDescription can not continue.");
      }
#ifdef USE_MPI
      const int nitems=2;
      int blocklengths[nitems] = {1, 1};
      MPI_Datatype types[nitems] = {alps::mpi::detail::mpi_type<size_t>(), alps::mpi::detail::mpi_type<precision>()};
      MPI_Aint offsets[nitems];
      offsets[0] = offsetof(Element, ind);
      offsets[1] = offsetof(Element, val);
      MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_Element);
      MPI_Type_commit(&mpi_Element);
#endif
    };

    /**
     * @brief Print static observables.
     *
     * @param ham - the Hamiltonian
     *
     * Prints all the parameters returned by calculate_static_observables().
     */
    void print_static_observables(Hamiltonian& ham){
      std::map<std::string, std::vector<double>> obs = calculate_static_observables(ham);
#ifdef USE_MPI
      int myid;
      MPI_Comm_rank(ham.comm(), &myid);
      if(!myid)
#endif
      {
        for(auto ivar = obs.begin(); ivar != obs.end(); ++ivar){
          std::cout << "<" << (*ivar).first << "> = {";
          for(int i = 0; i < (*ivar).second.size(); ++i){
            if(i){
              std::cout << ", ";
            }
            std::cout << (*ivar).second[i];
          }
          std::cout << "}" << std::endl;
        }
      }
    }

    /**
     * @brief Compute static observables.
     *
     * @param ham - the Hamiltonian
     *
     * Returns the following parameters for all sites:
     *  n: average number of electrons;
     *  n_up: average number of electrons with the spin up;
     *  n_down: average number of electrons with the spin down;
     *  m: average magnetic moment;
     *  d_occ: average double occupancy;
     */
    std::map<std::string, std::vector<precision>> calculate_static_observables(Hamiltonian& ham){
      std::map<std::string, std::vector<precision>> avg = {
        {_N_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_N_UP_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_N_DN_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_M_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_D_OCC_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)}
      };
      precision sum = 0.0;
      const EigenPair<precision, sector> &groundstate =  *ham.eigenpairs().begin();
      // Loop over all eigenpairs.
      for(auto ipair = ham.eigenpairs().begin(); ipair != ham.eigenpairs().end(); ++ipair){
        const EigenPair<precision, sector>& pair = *ipair;
        // Calculate Boltzmann factor, skip the states with trivial contribution.
        precision boltzmann_f = std::exp(
         -(pair.eigenvalue() - groundstate.eigenvalue()) * _beta
        );
        if(boltzmann_f < _cutoff){
//          std::cout<<"Skipped by Boltzmann factor."<<std::endl;
          continue;
        }
        // Sum the contributions.
        std::map<std::string, std::vector<precision>> contrib = calculate_static_observables_eigenvector(ham, pair);
        for(auto ivar = contrib.begin(); ivar != contrib.end(); ++ivar){
          for(int i = 0; i < (*ivar).second.size(); ++i){
            avg[(*ivar).first][i] += (*ivar).second[i] * boltzmann_f;
          }
        }
        sum += boltzmann_f;
      }
      for(auto ob = avg.begin(); ob != avg.end(); ++ob){
        for(int i = 0; i < (*ob).second.size(); ++i){
          (*ob).second[i] /= sum;
        }
      }
      return std::move(avg);
    }

    std::vector<std::pair<long long, precision>> find_major_electronic_configuration(Hamiltonian &_ham, const EigenPair <precision, sector> &pair, size_t nmax, precision trivial){
      _ham.model().symmetry().set_sector(pair.sector());
      _ham.storage().reset();
      int count = std::min(nmax, pair.eigenvector().size());
      std::vector<size_t> largest = std::vector<size_t>(pair.eigenvector().size());
#ifdef USE_MPI
      int myid;
      int nprocs;
      MPI_Comm_rank(_ham.comm(), &myid);
      MPI_Comm_size(_ham.comm(), &nprocs);
      std::vector<int> counts(nprocs);
      std::vector<int> displs(nprocs + 1);
      MPI_Gather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, _ham.comm());
      if(!myid){
        displs[0] = 0;
        for(size_t i = 0; i < nprocs; i++){
         displs[i + 1] = displs[i] + counts[i];
        }
      }
#endif
      for(size_t i = 0; i < largest.size(); ++i){
        largest[i] = i;
      }
      std::partial_sort(largest.begin(), largest.begin()+count, largest.end(), [&pair] (int a, int b) -> bool {
        return ((std::abs(pair.eigenvector()[a]) > std::abs(pair.eigenvector()[b]))
                || ((std::abs(pair.eigenvector()[a]) == std::abs(pair.eigenvector()[b])) && (a < b))
               );
      });
#ifdef USE_MPI
      std::vector<Element> send(count);
      for(size_t i = 0; i < count; ++i){
       send[i] = Element(largest[i] + _ham.storage().offset(), pair.eigenvector()[largest[i]]);
      }
      std::vector<Element> all(displs[nprocs]);
      MPI_Gatherv(send.data(), count, mpi_Element, all.data(), counts.data(), displs.data(), mpi_Element, 0, _ham.comm());
      if (!myid) {
        nmax = std::min(nmax, all.size());
        std::partial_sort(all.begin(), all.begin()+nmax, all.end(), [] (Element a, Element b) -> bool {return (a > b);});
        std::vector<std::pair<long long, precision>> ret(nmax);
        for(size_t i = 0; i < nmax; ++i){
          long long nst = _ham.model().symmetry().state_by_index(all[i].ind);
          ret[i] = std::pair<long long, precision>(nst, all[i].val);
        }
        return ret;
      } else {
        std::vector<std::pair<long long, precision>> ret(0);
        return ret;
      }
#else
      std::vector<std::pair<long long, precision>> ret(count);
      for(size_t i = 0; i < count; ++i){
        long long nst = _ham.model().symmetry().state_by_index(largest[i]);
        ret[i] = std::pair<long long, precision>(nst, pair.eigenvector()[largest[i]]);
      }
      return ret;
#endif
    }

    void print_major_electronic_configuration(Hamiltonian &_ham, const EigenPair <precision, sector> &pair, size_t nmax, precision trivial){
      std::vector<std::pair<long long, precision>> contribs = find_major_electronic_configuration(_ham, pair, nmax, trivial);
#ifdef USE_MPI
      int myid;
      MPI_Comm_rank(_ham.comm(), &myid);
      if(!myid)
#endif
      {
        std::cout << "Eigenvector components for eigenvalue " << pair.eigenvalue() << " ";
        pair.sector().print();
        std::cout << std::endl;
        for(size_t i = 0; i < contribs.size(); ++i){
          std::cout << contribs[i].second << " * |";
          std::string spin_down = std::bitset< 64 >( contribs[i].first ).to_string().substr(64-  _ham.model().orbitals(), _ham.model().orbitals());
          std::string spin_up   = std::bitset< 64 >( contribs[i].first ).to_string().substr(64-2*_ham.model().orbitals(), _ham.model().orbitals());
          std::cout<<spin_up<< "|"<<spin_down;
          std::cout << ">" << std::endl;
        }
      }
    }

#ifdef USE_MPI
    struct Element{
      Element() {};
      Element(size_t _ind, double _val) : ind(_ind), val(_val) {}
      size_t ind;
      double val;

      bool operator>(const Element &el) const {
        return ((std::abs(val) > std::abs(el.val))
                || ((std::abs(val) == std::abs(el.val)) && (ind < el.ind))
               );
      };

      bool operator<(const Element &el) const {
        return ((std::abs(val) < std::abs(el.val))
                || ((std::abs(val) == std::abs(el.val)) && (ind > el.ind))
               );
      };
    };

     MPI_Datatype mpi_Element;
#endif

  private:

    /**
     * @brief Compute static observables for one eigenvector.
     *
     * @param ham - the Hamiltonian
     * @param pair - the eigenpair
     */
    std::map<std::string, std::vector<precision>> calculate_static_observables_eigenvector(Hamiltonian& ham, const EigenPair<precision, sector>& pair){
      std::vector<precision> n(ham.model().interacting_orbitals(), 0.0);
      std::vector<precision> n_up(ham.model().interacting_orbitals(), 0.0);
      std::vector<precision> n_down(ham.model().interacting_orbitals(), 0.0);
      std::vector<precision> m(ham.model().interacting_orbitals(), 0.0);
      std::vector<precision> d_occ(ham.model().interacting_orbitals(), 0.0);
      ham.model().symmetry().set_sector(pair.sector());
      ham.storage().reset();
      // Loop over basis vectors.
      for(int i = 0; i < pair.eigenvector().size(); ++i){
        precision weight = pair.eigenvector()[i] * pair.eigenvector()[i];
        // Calculate static variables for each orbital.
        ham.model().symmetry().next_state();
        long long nst = ham.model().symmetry().state();
        for(int orb = 0; orb < ham.model().interacting_orbitals(); ++orb){
          int el_up = ham.model().checkState(nst, orb, ham.model().max_total_electrons());
          int el_down = ham.model().checkState(nst, orb + ham.model().orbitals(), ham.model().max_total_electrons());
          n[orb] += (el_up + el_down) * weight;
          n_up[orb] += el_up * weight;
          n_down[orb] += el_down * weight;
          m[orb] += (el_up - el_down) * weight;
          d_occ[orb] += el_up * el_down * weight;
        }
      }
      std::map<std::string, std::vector<precision>> result = {
        {_N_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_N_UP_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_N_DN_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_M_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_D_OCC_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)}
      };
#ifdef USE_MPI
      // Add the sums from all processes.
      MPI_Reduce(n.data(), result[_N_].data(), n.size(), alps::mpi::detail::mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(n_up.data(), result[_N_UP_].data(), n_up.size(), alps::mpi::detail::mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(n_down.data(), result[_N_DN_].data(), n_down.size(), alps::mpi::detail::mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(m.data(), result[_M_].data(), m.size(), alps::mpi::detail::mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(d_occ.data(), result[_D_OCC_].data(), d_occ.size(), alps::mpi::detail::mpi_type<precision>(), MPI_SUM, 0, ham.comm());
#else
      result[_N_] = n;
      result[_N_UP_] = n_up;
      result[_N_DN_] = n_down;
      result[_M_] = m;
      result[_D_OCC_] = d_occ;
#endif
     return std::move(result);
    }

    /// Inverse temperature
    precision _beta;
    /// Boltzmann-factor cutoff
    precision _cutoff;

  };

  template<class Hamiltonian>
  const std::string StaticObservables<Hamiltonian>::_N_ = "N";
  template<class Hamiltonian>
  const std::string StaticObservables<Hamiltonian>::_N_UP_ = "N_up";
  template<class Hamiltonian>
  const std::string StaticObservables<Hamiltonian>::_N_DN_ = "N_dn";
  template<class Hamiltonian>
  const std::string StaticObservables<Hamiltonian>::_M_ = "M";
  template<class Hamiltonian>
  const std::string StaticObservables<Hamiltonian>::_D_OCC_ = "D_occ";

}

#endif
