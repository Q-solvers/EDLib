//
// Created by iskakoff on 01/02/17.
//

#ifndef EDLIB_COMMONUTILS_H
#define EDLIB_COMMONUTILS_H

#include <complex>
#include <alps/gf/mesh.hpp>


namespace EDLib {
  namespace common {
    std::complex<double> freq_point(int index, const alps::gf::matsubara_positive_mesh & mesh, double beta) {
      return std::complex<double>(0.0, mesh.points()[index]);
    };

    std::complex<double> freq_point(int index, const alps::gf::real_frequency_mesh& mesh, double beta) {
      return std::complex<double>(mesh.points()[index], M_PI/beta);
    };
  }
}

#endif //EDLIB_COMMONUTILS_H
