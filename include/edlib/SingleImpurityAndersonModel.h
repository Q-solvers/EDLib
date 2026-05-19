#ifndef EDLIB_SINGLEIMPURITYANDERSONMODEL_H
#define EDLIB_SINGLEIMPURITYANDERSONMODEL_H

#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <utility>
#include <vector>

#include "edlib/CommonUtils.h"
#include "edlib/FermionicModel.h"
#include "edlib/Gf.h"
#include "edlib/Mesh.h"
#include "edlib/Parameters.h"
#include "edlib/SzSymmetry.h"

namespace edlib {

  namespace siam {

    /**
     * Base inner state. Concrete derived states encode either a hopping /
     * hybridisation transition (HybridisationInnerState) or a 4-operator
     * interaction term (InteractionInnerState).
     */
    template <class Prec>
    class InnerState {
    public:
      virtual ~InnerState() = default;
      // Non-pure defaults to mirror legacy alpscore-based hierarchy. Concrete
      // derived states (HybridisationInnerState / InteractionInnerState) override
      // them; extensions that add their own dispatch (HolsteinAnderson etc.) can
      // ignore the 4-arg form and provide additional overloads of their own.
      virtual int  valid(long long, int)                       const { return 0; }
      virtual void set  (long long, long long&, int&, int)     const {}
      virtual Prec value()                                     const { return Prec(0); }

    protected:
      static int checkState(long long nst, int im, int Ns) {
        return static_cast<int>((nst & (1ll << (2 * Ns - 1 - im))) >> (2 * Ns - 1 - im));
      }
      static void a(int i, long long jold, long long& k, int& isign, int Ip) {
        long long sign = 0;
        for (int ll = 0; ll < i; ++ll) sign += ((jold & (1ll << (Ip - ll - 1))) != 0) ? 1 : 0;
        isign = (sign % 2) == 0 ? 1 : -1;
        k = jold - (1ll << (Ip - i - 1));
      }
      static void adag(int i, long long jold, long long& k, int& isign, int Ip) {
        long long sign = 0;
        for (int ll = 0; ll < i; ++ll) sign += ((jold & (1ll << (Ip - ll - 1))) != 0) ? 1 : 0;
        isign = (sign % 2) == 0 ? 1 : -1;
        k = jold + (1ll << (Ip - i - 1));
      }
    };

    template <class Prec>
    class HybridisationInnerState : public InnerState<Prec> {
    public:
      HybridisationInnerState(int ii, int jj, int spin, Prec val)
          : _indicies(ii, jj), _spin(spin), _value(val) {}

      const std::pair<int, int>& indicies() const { return _indicies; }
      int                        spin()     const { return _spin; }
      Prec                       value()    const override { return _value; }

      int valid(long long nst, int Ns) const override {
        return InnerState<Prec>::checkState(nst, _indicies.first  + _spin * Ns, Ns)
             * (1 - InnerState<Prec>::checkState(nst, _indicies.second + _spin * Ns, Ns));
      }
      using InnerState<Prec>::set;
      void set(long long nst, long long& k, int& sign, int Ns) const override {
        long long k1, k2;
        int isign1, isign2;
        InnerState<Prec>::a   (_indicies.first  + _spin * Ns, nst, k1, isign1, 2 * Ns);
        InnerState<Prec>::adag(_indicies.second + _spin * Ns, k1,  k2, isign2, 2 * Ns);
        k    = k2;
        sign = isign1 * isign2;
      }

    private:
      std::pair<int, int> _indicies;
      int                 _spin;
      Prec                _value;
    };

    template <class Prec>
    class InteractionInnerState : public InnerState<Prec> {
    public:
      InteractionInnerState(int i, int j, int k, int l, int sigma, int sigma_prime, Prec U)
          : _i(i), _j(j), _k(k), _l(l),
            _sigma(sigma), _sigmaprime(sigma_prime), _U(U) {}

      int  i() const { return _i; }
      int  j() const { return _j; }
      int  k() const { return _k; }
      int  l() const { return _l; }
      Prec U() const { return _U; }

      Prec value() const override { return Prec(0.5) * _U; }

      using InnerState<Prec>::set;
      int valid(long long nst, int Ns) const override {
        int Ip = 2 * Ns;
        if (InnerState<Prec>::checkState(nst, _k + _sigma * Ns, Ns) != 0) {
          long long k3 = nst - (1ll << (Ip - 1 - _k - _sigma * Ns));
          if (InnerState<Prec>::checkState(k3, _l + _sigmaprime * Ns, Ns) != 0) {
            long long k4 = k3 - (1ll << (Ip - 1 - _l - _sigmaprime * Ns));
            if (InnerState<Prec>::checkState(k4, _j + _sigmaprime * Ns, Ns) == 0) {
              long long k2 = k4 | (1ll << (Ip - 1 - _j - _sigmaprime * Ns));
              return (1 - InnerState<Prec>::checkState(k2, _i + _sigma * Ns, Ns));
            }
          }
        }
        return 0;
      }
      void set(long long nst, long long& k, int& sign, int Ns) const override {
        long long k1, k2, k3, k4;
        int isign1, isign2, isign3, isign4;
        InnerState<Prec>::a   (_k + _sigma      * Ns, nst, k3, isign1, 2 * Ns);
        InnerState<Prec>::a   (_l + _sigmaprime * Ns, k3,  k4, isign2, 2 * Ns);
        InnerState<Prec>::adag(_j + _sigmaprime * Ns, k4,  k2, isign3, 2 * Ns);
        InnerState<Prec>::adag(_i + _sigma      * Ns, k2,  k1, isign4, 2 * Ns);
        k    = k1;
        sign = isign1 * isign2 * isign3 * isign4;
      }

    private:
      int  _i, _j, _k, _l;
      int  _sigma, _sigmaprime;
      Prec _U;
    };

  }

  template <class Prec>
  class SingleImpurityAndersonModel : public FermionicModel {
  public:
    using precision = Prec;
    using SYMMETRY  = SzSymmetry;
    using St        = siam::InnerState<Prec>;
    using HSt       = siam::HybridisationInnerState<Prec>;
    using USt       = siam::InteractionInnerState<Prec>;
    using Sector    = typename SzSymmetry::Sector;

    /**
     * Caller-supplied bath / model data. Dimensions checked against
     * Parameters::nsites, ::nspins and ::siam_norbitals at construction.
     *
     *   ml = p.siam_norbitals          (number of impurity orbitals)
     *   Nk = p.nsites - ml             (number of bath levels)
     *
     *   Vk   [ml][Nk][nspins]          impurity-bath hybridisation
     *   H0   [ml][ml][nspins]          non-interacting impurity Hamiltonian
     *   Epsk [Nk][nspins]              bath energies
     *   U    Gf<Prec,6> shape {nspins, nspins, ml, ml, ml, ml}
     */
    struct ModelData {
      std::vector<std::vector<std::vector<Prec>>> Vk;
      std::vector<std::vector<std::vector<Prec>>> H0;
      std::vector<std::vector<Prec>>              Epsk;
      Prec                                         mu = Prec(0);
      Gf<Prec, 6>                                  U;
      std::vector<std::array<int,2>>              sectors;
    };

    SingleImpurityAndersonModel(const Parameters& p, const ModelData& bath)
        : FermionicModel(p),
          _symmetry(p, bath.sectors),
          _ml(p.siam_norbitals),
          _Vk(bath.Vk),
          _H0(bath.H0),
          _Epsk(bath.Epsk),
          _xmu(bath.mu),
          _U(bath.U) {
      if (p.nspins != 2) {
        throw std::invalid_argument("SingleImpurityAndersonModel: NSPINS must be 2");
      }
      if (_ml > _Ns) {
        throw std::invalid_argument("SingleImpurityAndersonModel: siam.NORBITALS exceeds NSITES");
      }
      const int Nk = _Ns - _ml;
      if (static_cast<int>(_Epsk.size()) != Nk) {
        throw std::invalid_argument("SingleImpurityAndersonModel: Epsk size must equal nsites - siam.NORBITALS");
      }
      if (_U.shape(2) != _ml || _U.shape(3) != _ml || _U.shape(4) != _ml || _U.shape(5) != _ml
          || _U.shape(0) != p.nspins || _U.shape(1) != p.nspins) {
        throw std::invalid_argument("SingleImpurityAndersonModel: U must have shape [nspins,nspins,ml,ml,ml,ml]");
      }

      // Inter-orbital hoppings within the impurity cluster
      for (int im = 0; im < _ml; ++im) {
        for (int jm = 0; jm < im; ++jm) {
          for (int is = 0; is < _ms; ++is) {
            if (std::abs(_H0[im][jm][is]) > 1e-10) {
              _T_states.emplace_back(im, jm, is, _H0[im][jm][is]);
              _T_states.emplace_back(jm, im, is, _H0[im][jm][is]);
            }
          }
        }
      }

      // Impurity-bath hybridisation
      for (int im = 0; im < _ml; ++im) {
        for (int ik = 0; ik < static_cast<int>(_Vk[im].size()); ++ik) {
          for (int is = 0; is < _ms; ++is) {
            if (std::abs(_Vk[im][ik][is]) > 1e-10) {
              int imk = ik + _ml;
              _T_states.emplace_back(im,  imk, is, _Vk[im][ik][is]);
              _T_states.emplace_back(imk, im,  is, _Vk[im][ik][is]);
            }
          }
        }
      }

      // Off-diagonal interaction terms (skip density-density which is in diagonal)
      for (int is1 = 0; is1 < _ms; ++is1) {
        for (int is2 = 0; is2 < _ms; ++is2) {
          for (int i = 0; i < _ml; ++i) {
            for (int j = 0; j < _ml; ++j) {
              for (int k = 0; k < _ml; ++k) {
                for (int l = 0; l < _ml; ++l) {
                  if (((i == l) && (j == k) && (is1 == is2)) || ((i == k) && (j == l))) continue;
                  if (std::abs(_U(is1, is2, i, j, k, l)) != Prec(0)) {
                    _V_states.emplace_back(i, j, k, l, is1, is2, _U(is1, is2, i, j, k, l));
                  }
                }
              }
            }
          }
        }
      }
    }

    inline Prec diagonal(long long state) const {
      Prec xtemp = Prec(0);
      for (int is = 0; is < _ms; ++is) {
        for (int ik = 0; ik < static_cast<int>(_Epsk.size()); ++ik) {
          int ikm = ik + _ml;
          xtemp += _Epsk[ik][is] * checkState(state, ikm + is * _Ns, _Ip);
        }
      }
      for (int im = 0; im < _ml; ++im) {
        for (int is = 0; is < _ms; ++is) {
          xtemp += (_H0[im][im][is] - _xmu) * checkState(state, im + is * _Ns, _Ip);
        }
        for (int is = 0; is < _ms; ++is) {
          xtemp += Prec(0.5) * _U(is, is, im, im, im, im)
                 * checkState(state, im, _Ip) * checkState(state, im + _Ns, _Ip);
        }
        for (int jm = 0; jm < _ml; ++jm) {
          for (int is = 0; is < _ms; ++is) {
            if (im != jm) {
              xtemp += Prec(0.5) * (_U(is, is, im, jm, im, jm) - _U(is, is, im, jm, jm, im))
                     * checkState(state, im + is * _Ns, _Ip) * checkState(state, jm + is * _Ns, _Ip);
              xtemp += Prec(0.5) * _U(is, 1 - is, im, jm, im, jm)
                     * checkState(state, im + is * _Ns, _Ip) * checkState(state, jm + (1 - is) * _Ns, _Ip);
            }
          }
        }
      }
      return xtemp;
    }

    inline int  valid(const St& state, long long nst) const { return state.valid(nst, _Ns); }
    inline Prec set  (const St& state, long long nst, long long& k, int& sign) const {
      state.set(nst, k, sign, _Ns);
      return state.value();
    }

    /// @deprecated kept for API parity
    inline long long interacting_states(long long nst) const {
      long long up   = nst >> (_Ip - _ml);
      long long down = (nst & ((1ll << _Ns) - 1)) >> (_Ns - _ml);
      return (up << _ml) + down;
    }

    SzSymmetry&       symmetry()       { return _symmetry; }
    const SzSymmetry& symmetry() const { return _symmetry; }

    int  interacting_orbitals() const { return _ml; }

    const std::vector<HSt>& T_states() const { return _T_states; }
    const std::vector<USt>& V_states() const { return _V_states; }

    /// Read-only accessors used by legacy shims.
    const std::vector<std::vector<std::vector<Prec>>>& hybridisation() const { return _Vk; }
    const std::vector<std::vector<std::vector<Prec>>>& H0()            const { return _H0; }
    const std::vector<std::vector<Prec>>&              bath_energies() const { return _Epsk; }
    Prec                                               chem_potential() const { return _xmu; }

    template <class Mesh>
    void bare_greens_function(Gf<std::complex<double>, 3>& bare_gf,
                              const Mesh& mesh,
                              double beta) const {
      const int n_omega = bare_gf.shape(0);
      const int n_orb   = bare_gf.shape(1);
      const int n_spin  = bare_gf.shape(2);
      for (int iw = 0; iw < n_omega; ++iw) {
        std::complex<double> z = freq_point(iw, mesh, beta);
        for (int im = 0; im < n_orb; ++im) {
          for (int is = 0; is < n_spin; ++is) {
            std::complex<double> delta(0.0, 0.0);
            for (int ik = 0; ik < static_cast<int>(_Epsk.size()); ++ik) {
              delta += static_cast<double>(_Vk[im][ik][is]) * static_cast<double>(_Vk[im][ik][is])
                     / (z - static_cast<double>(_Epsk[ik][is]));
            }
            bare_gf(iw, im, is) = 1.0
                / (z - static_cast<double>(_H0[im][im][is]) - delta);
          }
        }
      }
    }

  private:
    SzSymmetry                                    _symmetry;
    int                                           _ml;
    std::vector<std::vector<std::vector<Prec>>>   _Vk;
    std::vector<std::vector<std::vector<Prec>>>   _H0;
    std::vector<std::vector<Prec>>                _Epsk;
    Prec                                          _xmu;
    Gf<Prec, 6>                                   _U;

    std::vector<HSt>                              _T_states;
    std::vector<USt>                              _V_states;
  };

}

#endif
