//
// Created by iskakoff on 28/08/16.
//

#ifndef HUBBARD_FERMIONICMODEL_H
#define HUBBARD_FERMIONICMODEL_H

#include <alps/params.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

namespace EDLib {
  namespace Model {
/**
 * @brief FermionicModel base class
 *
 * Define common fermionic routines for binary represented state
 *
 * @author iskakoff
 */
    class FermionicModel {
    public:
      FermionicModel(alps::params &p) : _Ns(p["NSITES"]), _ms(p["NSPINS"]), _Ip(int(p["NSPINS"]) * int(p["NSITES"])) {
      }

      /**
       * @brief Check that im state is occupated
       *
       * @param nst - current state
       * @param im - state to check
       * @param Ip - total number of fermionic spins for all sites
       *
       * @return 0 if state is empty, 1 - otherwise
       */
      int inline checkState(long long nst, const int im, int Ip) const {
        return (int) ((nst & (1ll << (Ip - 1 - im))) >> (Ip - 1 - im));
      }
      /**
       * @brief Anihilate particle
       * @param i [in] - site to anihilate particle
       * @param jold [in] - current state
       * @param k [out] - resulting state
       * @param isign [out] - fermionic sign
       */
      void inline a(int i, long long jold, long long &k, int &isign) {
        long long sign = 0;
        for (int ll = 0; ll < i; ++ll) {
          sign += ((jold & (1ll << (_Ip - ll - 1))) != 0) ? 1 : 0;
        }
        isign = (sign % 2) == 0 ? 1 : -1;
        k = jold - (1ll << (_Ip - i - 1));
      }

      /**
       * @brief Create particle
       * \param i [in] - site to create particle
       * \param jold [in] - current state
       * \param k [out] - resulting state
       * \param isign [out] - fermionic sign
       */
      void inline adag(int i, long long jold, long long &k, int &isign) {
        long long sign = 0;
        for (int ll = 0; ll < i; ++ll) {
          sign += ((jold & (1ll << (_Ip - ll - 1))) != 0) ? 1 : 0;
        }
        isign = (sign % 2) == 0 ? 1 : -1;
        k = jold + (1ll << (_Ip - i - 1));
      }


      int orbitals() const {
        return _Ns;
      }

      int max_total_electrons() const {
        return _Ip;
      }

      int spins() const {
        return _ms;
      }
      
      template<typename Mesh>
      void solve_dyson(const alps::gf::three_index_gf<std::complex<double>, Mesh, alps::gf::index_mesh, alps::gf::index_mesh >& bare_gf,
                       const alps::gf::three_index_gf<std::complex<double>, Mesh, alps::gf::index_mesh, alps::gf::index_mesh >& G_ij,
                       alps::gf::three_index_gf<std::complex<double>, Mesh, alps::gf::index_mesh, alps::gf::index_mesh >& sigma) {
        // solve Dyson equation
        for(int iw = 0; iw< bare_gf.mesh1().points().size(); ++iw) {
          typename Mesh::index_type w(iw);
          for (int is : bare_gf.mesh3().points()) {
            Eigen::MatrixXcd bare(_Ns, _Ns);
            Eigen::MatrixXcd bold(_Ns, _Ns);
            Eigen::MatrixXcd sigm(_Ns, _Ns);
            for (int im: bare_gf.mesh2().points()) {
              int I = im / _Ns;
              int J = im % _Ns;
              bare(I, J) = bare_gf(w, alps::gf::index_mesh::index_type(im), alps::gf::index_mesh::index_type(is));
              bold(I, J) = G_ij(w, alps::gf::index_mesh::index_type(im), alps::gf::index_mesh::index_type(is));
            }
            sigm = bare.inverse() - bold.inverse();
            for (int im: bare_gf.mesh2().points()) {
              int I = im / _Ns;
              int J = im % _Ns;
              sigma(w, alps::gf::index_mesh::index_type(im), alps::gf::index_mesh::index_type(is)) = sigm(I, J);
            }
          }
        }
      }

    protected:
      /**
       * _Ns - number of lattice sites
       * _ms - number of electron spins
       * _Ip - maximum number of electrons
       */
      int _Ns;
      int _ms;
      int _Ip;
    };
  }
}

#endif //HUBBARD_FERMIONICMODEL_H
