//
// Created by iskakoff on 11/01/17.
//

#ifndef HUBBARD_HDF5UTILS_HPP
#define HUBBARD_HDF5UTILS_HPP


/**
 * @brief HDF5Utils class
 *
 * @author iskakoff
 */
#include "SzSymmetry.h"
#include "EigenPair.h"

namespace EDLib {
  namespace hdf5 {
    template<typename T>
    struct HDF5Utils {

      /**
       * Generic function to save objects into hdf5 archive
       * @param t - obeject to be saved
       * @param ar - hdf5 archive
       * @param path - root path in the hdf5 archive
       */
      void save(const T& t,alps::hdf5::archive & ar, const std::string& path);
    };

    /**
     * Implementation for Sz-symmetry sector
     */
    template<>
    void HDF5Utils<typename Symmetry::SzSymmetry::Sector>::save(const typename Symmetry::SzSymmetry::Sector& s, alps::hdf5::archive & ar, const std::string& path) {
      ar[path + "/nup"]<<s.nup();
      ar[path + "/ndown"]<<s.ndown();
      ar[path + "/size"]<<s.size();
    }

    /**
     * Implementation for N-symmetry sector
     */
    template<>
    void HDF5Utils<typename Symmetry::NSymmetry::Sector>::save(const typename Symmetry::NSymmetry::Sector& s, alps::hdf5::archive & ar, const std::string& path) {
      ar[path + "/n"]<<s.n();
      ar[path + "/size"]<<s.size();
    }

    /**
     * Store eigenvalues into hdf5 archive
     *
     * @tparam Ham - type of Hamiltonian instance
     * @param h - Hamiltonian instance
     * @param ar - hdf5 archive
     * @param path - root path in hdf5 archive
     */
    template<typename Ham>
    void save_eigen_pairs(const Ham &h, alps::hdf5::archive & ar, const std::string& path) {
#ifdef USE_MPI
      int rank;
      MPI_Comm_rank(h.comm(), &rank);
      if(!rank){
#endif
        ar[path + "/eigenvalues/N"] << h.eigenpairs().size();
        int i = 0;
        std::vector<typename Ham::ModelType::precision> values;
        for(const EigenPair<typename Ham::ModelType::precision, typename Ham::ModelType::Sector>& e : h.eigenpairs()) {
          values.push_back(e.eigenvalue());
          HDF5Utils<typename Ham::ModelType::Sector>().save(e.sector(), ar, path + "/eigenvalues/sectors/" + boost::lexical_cast<std::string>(i));
          ++i;
        }
        ar[path + "/eigenvalues/data/"]<<values;
#ifdef USE_MPI
      }
#endif
    }

    /**
     * Store static observables into hdf5 archive
     * @param observables - map with static observables
     * @param ar - hdf5 archive
     * @param root_path - root path in hdf5 archive
     */
    void save_static_observables(const std::map < std::string, std::vector < double>> &observables, alps::hdf5::archive& ar, const std::string &root_path) {
      for(auto ob = observables.begin(); ob != observables.end(); ++ob){
        ar[root_path + "/static_observables/" + ob->first]<<ob->second;
      }
    }
  }
}


#endif //HUBBARD_HDF5UTILS_HPP
