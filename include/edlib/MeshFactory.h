//
// Created by iskakoff on 17/01/17.
//

#ifndef HUBBARD_MESHFACTORY_HPP
#define HUBBARD_MESHFACTORY_HPP


#include <alps/params.hpp>
#include <alps/gf/grid.hpp>

/**
 * @brief Factory class for frequency mesh initialization. Creates proper Mesh-object based on the Mesh type .
 *
 * @author iskakoff
 */
namespace EDLib {

  template<typename Mesh, typename ... Args>
  class MeshFactory {
  public:
    static Mesh createMesh(alps::params &p, Args... args);
  };

  template<>
  alps::gf::matsubara_positive_mesh MeshFactory < alps::gf::matsubara_positive_mesh, alps::gf::statistics::statistics_type >::
  createMesh(alps::params &p, alps::gf::statistics::statistics_type type) {
    return std::move(alps::gf::matsubara_positive_mesh(p["lanc.BETA"], p["lanc.NOMEGA"], type));
  }

  template<>
  alps::gf::real_frequency_mesh MeshFactory < alps::gf::real_frequency_mesh >::
  createMesh(alps::params &p) {
    alps::gf::grid::linear_real_frequency_grid g(p["lanc.EMIN"], p["lanc.EMAX"], p["lanc.NOMEGA"]);
    return std::move(alps::gf::real_frequency_mesh(g));
  }

}

#endif //HUBBARD_MESHFACTORY_HPP
