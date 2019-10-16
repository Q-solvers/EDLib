//
// Created by iskakoff on 17/01/17.
//

#ifndef EDLIB_MESHFACTORY_HPP
#define EDLIB_MESHFACTORY_HPP


#include <alps/params.hpp>
#include <alps/gf/grid.hpp>

/**
 * @brief Factory class for frequency mesh initialization. Creates proper Mesh-object based on the Mesh type .
 *
 * @author iskakoff
 */
namespace EDLib {

  class MatsubaraMeshFactory {
  public:
    using MeshType = alps::gf::matsubara_positive_mesh;
    static MeshType createMesh(alps::params &p, alps::gf::statistics::statistics_type type) {
      return std::move(alps::gf::matsubara_positive_mesh(p["lanc.BETA"], p["lanc.NOMEGA"], type));
    }
  };

  class RealFreqMeshFactory {
  public:
    using MeshType = alps::gf::real_frequency_mesh;
    static MeshType createMesh(alps::params &p) {
      alps::gf::grid::linear_real_frequency_grid g(p["lanc.EMIN"], p["lanc.EMAX"], p["lanc.NOMEGA"]);
      return std::move(alps::gf::real_frequency_mesh(g));
    }
  };
}

#endif //EDLIB_MESHFACTORY_HPP
