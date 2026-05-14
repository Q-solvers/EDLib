#ifndef EDLIB_DYSON_H
#define EDLIB_DYSON_H

#include <complex>
#include <stdexcept>

#include <Eigen/Core>
#include <Eigen/LU>

#include "edlib/Gf.h"

namespace edlib {

  /**
   * Solve Dyson's equation Sigma(iw) = G0^{-1}(iw) - G^{-1}(iw) for each
   * frequency and spin slice.
   *
   * All three GFs share shape [n_omega, nsites*nsites, n_spin]. The second
   * axis indexes orbital pairs as I*nsites + J, matching what
   * HubbardModel::bare_greens_function emits.
   */
  inline void solve_dyson(const Gf<std::complex<double>, 3>& bare,
                          const Gf<std::complex<double>, 3>& G,
                          Gf<std::complex<double>, 3>&       sigma,
                          int nsites) {
    if (bare.shape() != G.shape() || G.shape() != sigma.shape()) {
      throw std::invalid_argument("solve_dyson: shape mismatch among bare, G, sigma");
    }
    const int n_omega = bare.shape(0);
    const int n_orb   = bare.shape(1);
    const int n_spin  = bare.shape(2);
    if (n_orb != nsites * nsites) {
      throw std::invalid_argument("solve_dyson: shape(1) must equal nsites*nsites");
    }
    for (int iw = 0; iw < n_omega; ++iw) {
      for (int is = 0; is < n_spin; ++is) {
        Eigen::MatrixXcd b(nsites, nsites);
        Eigen::MatrixXcd g(nsites, nsites);
        for (int I = 0; I < nsites; ++I) {
          for (int J = 0; J < nsites; ++J) {
            b(I, J) = bare(iw, I * nsites + J, is);
            g(I, J) = G   (iw, I * nsites + J, is);
          }
        }
        Eigen::MatrixXcd s = b.inverse() - g.inverse();
        for (int I = 0; I < nsites; ++I) {
          for (int J = 0; J < nsites; ++J) {
            sigma(iw, I * nsites + J, is) = s(I, J);
          }
        }
      }
    }
  }

}

#endif
