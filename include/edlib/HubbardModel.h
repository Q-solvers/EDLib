#ifndef EDLIB_HUBBARDMODEL_H
#define EDLIB_HUBBARDMODEL_H

#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>

#include "edlib/CommonUtils.h"
#include "edlib/FermionicModel.h"
#include "edlib/Gf.h"
#include "edlib/Mesh.h"
#include "edlib/Parameters.h"
#include "edlib/SzSymmetry.h"

namespace edlib {

  namespace hubbard {

    template <class Prec>
    class InnerState {
    public:
      InnerState(int ii, int jj, int spin, Prec val)
          : _indicies(ii, jj), _spin(spin), _value(val) {}

      const std::pair<int, int>& indicies() const { return _indicies; }
      const Prec&                value()    const { return _value; }
      int                        spin()     const { return _spin; }

    private:
      std::pair<int, int> _indicies;
      int                 _spin;
      Prec                _value;
    };

  }

  template <class Prec>
  class HubbardModel : public FermionicModel {
  public:
    using precision = Prec;
    using SYMMETRY  = SzSymmetry;
    using St        = hubbard::InnerState<Prec>;
    using Sector    = typename SzSymmetry::Sector;

    /**
     * Per-instance Hubbard parameters. All arrays are caller-supplied; the
     * model performs no file I/O. Sizes are validated against the Parameters
     * passed to the constructor.
     */
    struct ModelData {
      std::vector<std::vector<Prec>> hopping;         ///< [nsites][nsites] -- mandatory
      std::vector<Prec>              U;               ///< [nsites]         -- mandatory
      std::vector<Prec>              mu;              ///< [nsites]         -- mandatory
      std::vector<Prec>              magnetic_field;  ///< [nsites]         -- optional (zeros)
      std::vector<std::vector<Prec>> exchange;        ///< [nsites][nsites] -- optional (zeros)
      std::vector<std::vector<Prec>> site_energy;     ///< [nsites][nspins] -- optional (zeros)
      std::vector<std::array<int,2>> sectors;         ///< optional sector restriction
    };

    HubbardModel(const Parameters& p, const ModelData& bath)
        : FermionicModel(p),
          _symmetry(p, bath.sectors),
          _t(bath.hopping),
          _U(bath.U),
          _xmu(bath.mu),
          _Hmag(bath.magnetic_field.empty()
                ? std::vector<Prec>(p.nsites, Prec(0))
                : bath.magnetic_field),
          _J(bath.exchange.empty()
             ? std::vector<std::vector<Prec>>(p.nsites,
                   std::vector<Prec>(p.nsites, Prec(0)))
             : bath.exchange),
          _Eps(bath.site_energy.empty()
               ? std::vector<std::vector<Prec>>(p.nsites,
                     std::vector<Prec>(p.nspins, Prec(0)))
               : bath.site_energy) {
      validate_model_data(p);
      for (int ii = 0; ii < _Ns; ++ii) {
        for (int jj = 0; jj < _Ns; ++jj) {
          if (std::abs(_t[ii][jj]) > 1e-10) {
            for (int is = 0; is < _ms; ++is) {
              _states.emplace_back(ii, jj, is, _t[ii][jj]);
            }
          }
        }
      }
    }

    inline int valid(const St& state, long long nst) const {
      return checkState(nst, state.indicies().first  + state.spin() * _Ns, _Ip)
           * (1 - checkState(nst, state.indicies().second + state.spin() * _Ns, _Ip));
    }

    inline Prec set(const St& state, long long nst, long long& k, int& sign) const {
      long long k1, k2;
      int isign1, isign2;
      a   (state.indicies().first  + state.spin() * _Ns, nst, k1, isign1);
      adag(state.indicies().second + state.spin() * _Ns, k1,  k2, isign2);
      k    = k2;
      sign = -isign1 * isign2;  // -t c^+ c
      return state.value();
    }

    inline Prec diagonal(long long state) const {
      Prec xtemp = Prec(0);
      for (int im = 0; im < _Ns; ++im) {
        for (int is = 0; is < _ms; ++is) {
          xtemp += (_Eps[im][is] - _xmu[is]) * checkState(state, im + is * _Ns, _Ip);
        }
        xtemp += _U[im]    * checkState(state, im,       _Ip) * checkState(state, im + _Ns, _Ip);
        xtemp += _Hmag[im] * (checkState(state, im + _Ns, _Ip) - checkState(state, im, _Ip));
        for (int im2 = 0; im2 < _Ns; ++im2) {
          xtemp += _J[im][im2]
                 * (checkState(state, im,        _Ip) - checkState(state, im + _Ns, _Ip))
                 * (checkState(state, im2,       _Ip) - checkState(state, im2 + _Ns, _Ip));
        }
      }
      return xtemp;
    }

    /// @deprecated kept for API parity with legacy callers
    inline long long interacting_states(long long nst) const { return nst; }

    const std::vector<St>& T_states() const { return _states; }
    const std::vector<St>& V_states() const { return _V_states; }

    int  interacting_orbitals() const { return _Ns; }

    const SzSymmetry& symmetry() const { return _symmetry; }
    SzSymmetry&       symmetry()       { return _symmetry; }

    /// Read-only accessors used by legacy shims. Stable but not part of the
    /// recommended public API; use bare_greens_function for the standard path.
    const std::vector<std::vector<Prec>>& hopping_matrix() const { return _t; }
    const std::vector<Prec>&              interaction()    const { return _U; }
    const std::vector<Prec>&              chem_potential() const { return _xmu; }
    const std::vector<std::vector<Prec>>& site_energy()    const { return _Eps; }

    /**
     * Compute the non-interacting Green's function on the supplied frequency
     * mesh and fill bare_gf with shape [n_omega, _Ns * _Ns, _ms].
     */
    template <class Mesh>
    void bare_greens_function(Gf<std::complex<double>, 3>& bare_gf,
                              const Mesh& mesh,
                              double beta) const {
      const int n_omega = bare_gf.shape(0);
      const int n_orb   = bare_gf.shape(1);
      const int n_spin  = bare_gf.shape(2);
      if (n_orb != _Ns * _Ns) {
        throw std::invalid_argument("bare_greens_function: shape(1) must equal nsites*nsites");
      }
      for (int iw = 0; iw < n_omega; ++iw) {
        std::complex<double> z = freq_point(iw, mesh, beta);
        for (int is = 0; is < n_spin; ++is) {
          Eigen::MatrixXcd G_inv = Eigen::MatrixXcd::Zero(_Ns, _Ns);
          for (int I = 0; I < _Ns; ++I) {
            G_inv(I, I) = z + _xmu[I] - _Eps[I][is];
            for (int J = 0; J < _Ns; ++J) {
              G_inv(I, J) += _t[I][J];
            }
          }
          Eigen::MatrixXcd G = G_inv.inverse();
          for (int I = 0; I < _Ns; ++I) {
            for (int J = 0; J < _Ns; ++J) {
              bare_gf(iw, I * _Ns + J, is) = G(I, J);
            }
          }
        }
      }
    }

  private:
    void validate_model_data(const Parameters& p) const {
      const int N = p.nsites;
      auto bad = [](const char* what) { throw std::invalid_argument(what); };

      if (static_cast<int>(_t.size())   != N) bad("HubbardModel: hopping must be [nsites][nsites]");
      for (const auto& row : _t)
        if (static_cast<int>(row.size()) != N) bad("HubbardModel: hopping rows must have size nsites");
      if (static_cast<int>(_U.size())   != N) bad("HubbardModel: U must have size nsites");
      if (static_cast<int>(_xmu.size()) != N) bad("HubbardModel: mu must have size nsites");
      if (static_cast<int>(_Hmag.size())!= N) bad("HubbardModel: magnetic_field size mismatch");
      if (static_cast<int>(_J.size())   != N) bad("HubbardModel: exchange must be [nsites][nsites]");
      for (const auto& row : _J)
        if (static_cast<int>(row.size()) != N) bad("HubbardModel: exchange rows must have size nsites");
      if (static_cast<int>(_Eps.size()) != N) bad("HubbardModel: site_energy must be [nsites][nspins]");
      for (const auto& row : _Eps)
        if (static_cast<int>(row.size()) != p.nspins) bad("HubbardModel: site_energy rows must have size nspins");
    }

    SzSymmetry                       _symmetry;
    std::vector<std::vector<Prec>>   _t;
    std::vector<Prec>                _U;
    std::vector<Prec>                _xmu;
    std::vector<Prec>                _Hmag;
    std::vector<std::vector<Prec>>   _J;
    std::vector<std::vector<Prec>>   _Eps;

    std::vector<St>                  _states;
    std::vector<St>                  _V_states;
  };

}

#endif
