#ifndef EDLIB_MESHFACTORY_H
#define EDLIB_MESHFACTORY_H

#include "edlib/Mesh.h"
#include "edlib/Parameters.h"

namespace edlib {

  /**
   * Convenience factories. Most call sites can construct meshes directly;
   * these mirror the legacy API for symmetry.
   */
  struct MatsubaraMeshFactory {
    using MeshType = MatsubaraMesh;
    static MatsubaraMesh createMesh(const Parameters& p, Statistics stat) {
      return MatsubaraMesh(p.lanc_beta, p.lanc_nomega, stat);
    }
  };

  struct RealFreqMeshFactory {
    using MeshType = RealFreqMesh;
    static RealFreqMesh createMesh(const Parameters& p) {
      return RealFreqMesh(p.lanc_emin, p.lanc_emax, p.lanc_nomega);
    }
  };

}

#endif
