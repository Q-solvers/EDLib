#ifndef HUBBARD_STATICOBSERVABLES_H
#define HUBBARD_STATICOBSERVABLES_H

#include "EigenPair.h"

#include <alps/params.hpp>

#include <vector>
#include <bitset>
#include <string>
#include <cmath>

namespace EDLib {
  /**
   * Class for calculation of the static observables.
   *
   *
   *
   * @tparam Hamiltonian - type of Hamiltonian object
   */
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
    const static std::string _N_EFF_;
    const static std::string _MI_MJ_;
    const static std::string _E_;

    /**
     * Construct an object of the static observables class
     *
     * @param p - AlpsCore parameter object
     */
    StaticObservables(alps::params &p) :
      _beta(p["lanc.BETA"].as<precision>()),
      _cutoff(p["lanc.BOLTZMANN_CUTOFF"])
    {
      if(p["storage.EIGENVALUES_ONLY"] == 1) {
        throw std::logic_error("Eigenvectors have not been computed. StaticObservables can not continue.");
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
     *  m_i m_j: product of magnetic moments on the i-th and j-th sites;
     *  d_occ: average double occupancy;
     *  N_eff: average effective dimension of Hilbert space;
     */
    std::map<std::string, std::vector<precision>> calculate_static_observables(Hamiltonian& ham){
#ifdef USE_MPI
      int myid;
      MPI_Comm_rank(ham.storage().comm(), &myid);
#endif
      std::map<std::string, std::vector<precision>> avg = {
        {_N_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_N_UP_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_N_DN_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_M_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_D_OCC_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_MI_MJ_, std::vector<precision>(ham.model().interacting_orbitals() * ham.model().interacting_orbitals(), 0.0)},
        {_N_EFF_, std::vector<precision>(1, 0.0)},
        {_E_, std::vector<precision>(1, 0.0)}
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
#ifdef USE_MPI
        if(!myid)
#endif
        std::cout << "Compute static observables contribution for eigenvalue E=" << pair.eigenvalue() << " with Boltzmann factor = " << boltzmann_f << "; for sector" << pair.sector() << std::endl;

        std::map<std::string, std::vector<precision>> contrib = calculate_static_observables_eigenvector(ham, pair);
        for(auto ivar = contrib.begin(); ivar != contrib.end(); ++ivar){
          for(int i = 0; i < (*ivar).second.size(); ++i){
            avg[(*ivar).first][i] += (*ivar).second[i] * boltzmann_f;
          }
        }
        sum += boltzmann_f;
      }
#ifdef USE_MPI
      if(!myid)
#endif
      std::cout << "Statsum: " << sum << std::endl;
      for(auto ob = avg.begin(); ob != avg.end(); ++ob){
        for(int i = 0; i < (*ob).second.size(); ++i){
          (*ob).second[i] /= sum;
        }
      }
      return std::move(avg);
    }

    /**
     * @brief Find largest (by magnitude) coefficients in an eigenvector.
     *
     * @param ham - the Hamiltonian
     * @param pair - the eigenpair
     * @param nmax - maximum number of coefficients to be returned;
     * @param trivial - skip the coefficients smaller than this number.
     *
     * Returns a vector of pairs sorted starting from largest magnitude:
     *  first: symmetry state;
     *  second: coefficient;
     */
    std::vector<std::pair<long long, precision>> find_largest_coefficients(Hamiltonian& ham, const EigenPair<precision, sector>& pair, size_t nmax, precision trivial){
      ham.model().symmetry().set_sector(pair.sector());
      ham.storage().reset();
      int count = std::min(nmax, pair.eigenvector().size());
      std::vector<size_t> largest = std::vector<size_t>(pair.eigenvector().size());
#ifdef USE_MPI
      int myid;
      int nprocs;
      MPI_Comm_rank(ham.comm(), &myid);
      MPI_Comm_size(ham.comm(), &nprocs);
      std::vector<int> counts(nprocs);
      std::vector<int> displs(nprocs + 1);
      MPI_Gather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, ham.comm());
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
       send[i] = Element(largest[i] + ham.storage().offset(), pair.eigenvector()[largest[i]]);
      }
      std::vector<Element> all(displs[nprocs]);
      MPI_Gatherv(send.data(), count, mpi_Element, all.data(), counts.data(), displs.data(), mpi_Element, 0, ham.comm());
      if (!myid) {
        nmax = std::min(nmax, all.size());
        std::partial_sort(all.begin(), all.begin()+nmax, all.end(), [] (Element a, Element b) -> bool {return (a > b);});
        std::vector<std::pair<long long, precision>> ret(nmax);
        for(size_t i = 0; i < nmax; ++i){
          long long nst = ham.model().symmetry().state_by_index(all[i].ind);
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
        long long nst = ham.model().symmetry().state_by_index(largest[i]);
        ret[i] = std::pair<long long, precision>(nst, pair.eigenvector()[largest[i]]);
      }
      return ret;
#endif
    }

    /**
     * @brief Calculate contribution of the states with the same coefficient (classes) to the eigenvector.
     *
     * @param ham - the Hamiltonian
     * @param pair - the eigenpair
     * @param nmax - maximum number of coefficients to be processed;
     * @param trivial - skip the coefficients smaller than this number;
     * @param cumulative - calculate cumulative contribution: this class and all previous classes.
     *
     * Returns a vector of pairs:
     *  first: index of the first state in the class, as returned by find_largest_coefficients();
     *  second: contribution of the class;
     * The vector is sorted the largest coefficient of the class.
     */
    std::vector<std::pair<size_t, precision>> calculate_class_contrib(Hamiltonian& ham, const EigenPair<precision, sector>& pair, size_t nmax, precision trivial, bool cumulative){
      std::vector<std::pair<long long, precision>> coeffs = find_largest_coefficients(ham, pair, nmax, trivial);
      std::vector<std::pair<size_t, precision>> contribs(0);
#ifdef USE_MPI
      int myid;
      MPI_Comm_rank(ham.comm(), &myid);
      if(!myid)
#endif
      {
        for(size_t i = 0; i < coeffs.size(); ++i){
          // Add up squares of coefficients within each class.
          if(!i || std::abs(coeffs[i - 1].second - coeffs[i].second) > trivial){
            contribs.push_back(std::pair<size_t, precision>(i, coeffs[i].second * coeffs[i].second));
            if(i){
             contribs.back().second += (cumulative ? contribs[contribs.size() - 2].second : 0.0);
            }
          }else{
            contribs.back().second += coeffs[i].second * coeffs[i].second;
          }
        }
        // The last contribs may be truncated, remove.
        contribs.pop_back();
      }
      return contribs;
    }

    /**
     * @brief Print largest coefficients of an eigenvector in decreasing order of magnitude.
     *
     * @param ham - the Hamiltonian
     * @param pair - the eigenpair
     * @param nmax - maximum number of coefficients to be processed;
     * @param trivial - skip the coefficients smaller than this number.
     */
    void print_major_electronic_configuration(Hamiltonian& ham, const EigenPair<precision, sector>& pair, size_t nmax, precision trivial){
      std::vector<std::pair<long long, precision>> coeffs = find_largest_coefficients(ham, pair, nmax, trivial);
#ifdef USE_MPI
      int myid;
      MPI_Comm_rank(ham.comm(), &myid);
      if(!myid)
#endif
      {
        std::cout << "Eigenvector components for eigenvalue " << pair.eigenvalue() << " ";
        pair.sector().print();
        std::cout << std::endl;
        for(size_t i = 0; i < coeffs.size(); ++i){
          std::cout << coeffs[i].second << " * |";
          std::string spin_down = std::bitset< 64 >( coeffs[i].first ).to_string().substr(64-  ham.model().orbitals(), ham.model().orbitals());
          std::string spin_up   = std::bitset< 64 >( coeffs[i].first ).to_string().substr(64-2*ham.model().orbitals(), ham.model().orbitals());
          std::cout<<spin_up<< "|"<<spin_down;
          std::cout << ">" << std::endl;
        }
      }
    }

    /**
     * @brief Print contribution of the states with the same coefficient (classes) to the eigenvector.
     *
     * @param ham - the Hamiltonian
     * @param pair - the eigenpair
     * @param nmax - maximum number of coefficients to be processed;
     * @param trivial - skip the coefficients smaller than this number.
     * @param cumulative - calculate cumulative contribution: this class and all previous classes.
     */
    void print_class_contrib(Hamiltonian& ham, const EigenPair<precision, sector>& pair, size_t nmax, precision trivial, bool cumulative){
      std::vector<std::pair<size_t, precision>> contribs = calculate_class_contrib(ham, pair, nmax, trivial, cumulative);
#ifdef USE_MPI
      int myid;
      MPI_Comm_rank(ham.comm(), &myid);
      if(!myid)
#endif
      {
        std::cout << "Contributions of eigenvector component classes for eigenvalue " << pair.eigenvalue() << " ";
        pair.sector().print();
        std::cout << std::endl;
        for(size_t i = 0; i < contribs.size(); ++i){
          std::cout << contribs[i].first << "\t" << contribs[i].second << std::endl;
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
      std::vector<precision> mimj(ham.model().interacting_orbitals() * ham.model().interacting_orbitals(), 0.0);
      std::vector<precision> d_occ(ham.model().interacting_orbitals(), 0.0);
      precision inverse_N_eff = 0.0;
      precision energy = pair.eigenvalue();

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
          for(int orb2 = 0; orb2 < ham.model().interacting_orbitals(); ++orb2){
            int el_up2 = ham.model().checkState(nst, orb2, ham.model().max_total_electrons());
            int el_down2 = ham.model().checkState(nst, orb2 + ham.model().interacting_orbitals(), ham.model().max_total_electrons());
            mimj[ham.model().interacting_orbitals() * orb + orb2] += (el_up - el_down) * (el_up2 - el_down2) * weight;
          }
          d_occ[orb] += el_up * el_down * weight;
        }
        inverse_N_eff += weight * weight;
      }

      std::map<std::string, std::vector<precision>> result = {
        {_N_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_N_UP_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_N_DN_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_M_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_D_OCC_, std::vector<precision>(ham.model().interacting_orbitals(), 0.0)},
        {_MI_MJ_, std::vector<precision>(ham.model().interacting_orbitals() * ham.model().interacting_orbitals(), 0.0)},
        {_N_EFF_, std::vector<precision>(1, 0.0)},
        {_E_, std::vector<precision>(1, 0.0)}
      };
#ifdef USE_MPI
      // Add the sums from all processes.
      MPI_Reduce(n.data(), result[_N_].data(), n.size(), alps::mpi::detail::mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(n_up.data(), result[_N_UP_].data(), n_up.size(), alps::mpi::detail::mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(n_down.data(), result[_N_DN_].data(), n_down.size(), alps::mpi::detail::mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(m.data(), result[_M_].data(), m.size(), alps::mpi::detail::mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(d_occ.data(), result[_D_OCC_].data(), d_occ.size(), alps::mpi::detail::mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(mimj.data(), result[_MI_MJ_].data(), mimj.size(), alps::mpi::detail::mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(&inverse_N_eff, &result[_N_EFF_][0], 1, alps::mpi::detail::mpi_type<precision>(), MPI_SUM, 0, ham.comm());
      MPI_Reduce(&energy, &result[_E_][0], 1, alps::mpi::detail::mpi_type<precision>(), MPI_SUM, 0, ham.comm());
#else
      result[_N_] = n;
      result[_N_UP_] = n_up;
      result[_N_DN_] = n_down;
      result[_M_] = m;
      result[_D_OCC_] = d_occ;
      result[_MI_MJ_] = mimj;
      result[_N_EFF_][0] = inverse_N_eff;
      result[_E_][0] = energy;
#endif
      result[_N_EFF_][0] = 1 / result[_N_EFF_][0];
      return result;
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
  template<class Hamiltonian>
  const std::string StaticObservables<Hamiltonian>::_N_EFF_ = "N_eff";
  template<class Hamiltonian>
  const std::string StaticObservables<Hamiltonian>::_MI_MJ_ = "M_i M_j";
  template<class Hamiltonian>
  const std::string StaticObservables<Hamiltonian>::_E_ = "E";

}

#endif
