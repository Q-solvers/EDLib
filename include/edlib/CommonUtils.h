#ifndef EDLIB_COMMONUTILS_H
#define EDLIB_COMMONUTILS_H

#include <cmath>
#include <complex>

#include "edlib/Mesh.h"

namespace edlib {

  inline std::complex<double>
  freq_point(int index, const MatsubaraMesh& mesh, double /*beta*/) {
    return std::complex<double>(0.0, mesh.points()[index]);
  }

  inline std::complex<double>
  freq_point(int index, const RealFreqMesh& mesh, double beta) {
    return std::complex<double>(mesh.points()[index], M_PI / beta);
  }

}

#endif
