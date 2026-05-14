#ifndef EDLIB_EXT_HOLSTEINANDERSONMODEL_H
#define EDLIB_EXT_HOLSTEINANDERSONMODEL_H

#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <utility>
#include <vector>

#include "edlib/FermionicModel.h"
#include "edlib/Gf.h"
#include "edlib/Parameters.h"
#include "edlib/SingleImpurityAndersonModel.h"
#include "ext/SzSymmetryWithBoson.h"

namespace edlib { namespace ext {

  /**
   * Extra parameters for Holstein-Anderson model. Mirrors the legacy
   * EDLib::Ext::define_parameters() NBBITS/NBLEVEL keys plus the impurity
   * orbital count (legacy: NORBITALS). Embedded in the model's ModelData struct
   * so the model fits the standard (Parameters, ModelData) ctor signature.
   */
  struct HolsteinAndersonParameters {
    int nbbits   = 3;  ///< bits per bosonic mode
    int nblevel  = 1;  ///< number of bosonic modes
    int ml       = 1;  ///< number of impurity orbitals
    double avg   = 0;  ///< average orbital occupation reference
  };

  namespace holstein {

    template <class Prec>
    class InnerState : public edlib::siam::InnerState<Prec> {
    public:
      using edlib::siam::InnerState<Prec>::a;
      using edlib::siam::InnerState<Prec>::adag;
      using edlib::siam::InnerState<Prec>::checkState;

      int valid(long long, int) const override { return 0; }
      void set(long long, long long&, int&, int) const override {}
      Prec value() const override { return Prec(0); }

      virtual int  valid(long long, int, int) const { return 0; }
      virtual Prec set  (long long, long long&, int&, int, int) const { return Prec(0); }
    };

    template <class Prec>
    class HybridisationInnerState : public InnerState<Prec> {
    public:
      HybridisationInnerState(int i1, int i2, int is, Prec val)
          : _indicies(i1, i2), _spin(is), _value(val) {}

      const std::pair<int, int>& indicies() const { return _indicies; }
      Prec                       value()    const override { return _value; }
      int                        spin()     const { return _spin; }

      int valid(long long nst, int Ns, int Nb) const override {
        return this->checkState(nst >> Nb, _indicies.first  + _spin * Ns, Ns)
             * (1 - this->checkState(nst >> Nb, _indicies.second + _spin * Ns, Ns));
      }
      Prec set(long long nst, long long& k, int& sign, int Ns, int Nb) const override {
        long long k1, k2;
        int isign1, isign2;
        long long fnst = nst >> Nb;
        long long bnst = nst & ((1ll << Nb) - 1);
        this->a   (_indicies.first  + _spin * Ns, fnst, k1, isign1, 2 * Ns);
        this->adag(_indicies.second + _spin * Ns, k1,   k2, isign2, 2 * Ns);
        k    = (k2 << Nb) + bnst;
        sign = isign1 * isign2;
        return _value;
      }

    private:
      std::pair<int, int> _indicies;
      int                 _spin;
      Prec                _value;
    };

    template <class Prec>
    class BosonInnerState : public InnerState<Prec> {
    public:
      BosonInnerState(int ib, int i, Prec value, int bit_cutoff, Prec avg, bool dag)
          : _b(ib), _i(i),
            _bit_cutoff(bit_cutoff),
            _cutoff((1ll << bit_cutoff) - 1),
            _dag(dag), _avg(avg), _value(value) {}

      int valid(long long nst, int Ns, int Nb) const override {
        long long bnst = nst & ((1ll << Nb) - 1);
        long long cbos = (bnst >> (_bit_cutoff * _b)) & _cutoff;
        return (std::abs(this->checkState(nst >> Nb, _i,       Ns)
                       + this->checkState(nst >> Nb, Ns + _i, Ns)
                       - _avg) > Prec(1e-9))
             * (_dag ? cbos < _cutoff : cbos > 0);
      }
      Prec set(long long nst, long long& k, int& sign, int Ns, int Nb) const override {
        Prec N = (this->checkState(nst >> Nb, _i,       Ns)
                + this->checkState(nst >> Nb, Ns + _i, Ns)
                - _avg);
        long long bnst = nst & ((1ll << Nb) - 1);
        long long cbos = (bnst >> (_bit_cutoff * _b)) & _cutoff;
        if (_dag) {
          k    = nst + (1ll << (_bit_cutoff * _b));
          sign = 1;
          return N * _value * static_cast<Prec>(std::sqrt(static_cast<double>(cbos + 1)));
        } else {
          k    = nst - (1ll << (_bit_cutoff * _b));
          sign = 1;
          return N * _value * static_cast<Prec>(std::sqrt(static_cast<double>(cbos)));
        }
      }

    private:
      int       _b;
      int       _i;
      int       _bit_cutoff;
      long long _cutoff;
      bool      _dag;
      Prec      _avg;
      Prec      _value;
    };

  }

  template <class Prec>
  class HolsteinAndersonModel : public edlib::FermionicModel {
  public:
    using precision = Prec;
    using SYMMETRY  = SzSymmetryWithBoson;
    using St        = holstein::InnerState<Prec>;
    using HSt       = holstein::HybridisationInnerState<Prec>;
    using BSt       = holstein::BosonInnerState<Prec>;
    using Sector    = typename SzSymmetryWithBoson::Sector;

    /**
     * Caller-supplied Holstein-Anderson bath / model data.
     *
     *   ml         = ep.ml         (impurity orbitals)
     *   Nk         = p.nsites - ml (fermionic bath levels)
     *   nblevel    = ep.nblevel    (bosonic modes per orbital)
     *
     *   tk     [ml][ml]                  intracluster hopping
     *   U      [ml][ml]                  on-site/inter-orbital interaction
     *   Eps    [ml][nspins]              impurity site energies
     *   Vk     [ml][Nk][nspins]          fermionic hybridisation
     *   Epsk   [Nk][nspins]              fermionic bath levels
     *   w0     [nblevel]                 bosonic mode frequencies
     *   W      [ml][nblevel]             bosonic couplings
     */
    struct ModelData {
      HolsteinAndersonParameters                  ep;
      std::vector<std::vector<Prec>>              tk;
      std::vector<std::vector<Prec>>              U;
      std::vector<std::vector<Prec>>              Eps;
      std::vector<std::vector<std::vector<Prec>>> Vk;
      std::vector<std::vector<Prec>>              Epsk;
      std::vector<Prec>                           w0;
      std::vector<std::vector<Prec>>              W;
      Prec                                         mu = Prec(0);
      std::vector<std::array<int,2>>              sectors;
    };

    HolsteinAndersonModel(const Parameters& p, const ModelData& bath)
        : HolsteinAndersonModel(p, bath.ep, bath) {}

    HolsteinAndersonModel(const Parameters& p,
                          const HolsteinAndersonParameters& ep,
                          const ModelData& bath)
        : edlib::FermionicModel(p),
          _symmetry(p, ep.nbbits * ep.nblevel, bath.sectors),
          _Nb(ep.nbbits * ep.nblevel),
          _ml(ep.ml),
          _xmu(bath.mu),
          _U(bath.U),
          _Eps(bath.Eps),
          _tk(bath.tk),
          _Vk(bath.Vk),
          _Epsk(bath.Epsk),
          _w0(bath.w0),
          _W(bath.W) {
      if (p.nspins != 2) {
        throw std::invalid_argument("HolsteinAndersonModel: NSPINS must be 2");
      }
      if (static_cast<int>(_w0.size()) != ep.nblevel) {
        throw std::invalid_argument("HolsteinAndersonModel: w0.size() must equal nblevel");
      }
      for (int i = 0; i < _ml; ++i) {
        if (_Vk[i].size() != _Epsk.size()) {
          throw std::invalid_argument("HolsteinAndersonModel: Vk[i].size() must equal Nk = Epsk.size()");
        }
        if (_W[i].size() != _w0.size()) {
          throw std::invalid_argument("HolsteinAndersonModel: W[i].size() must equal nblevel = w0.size()");
        }
      }

      // Fermionic bath - intracluster hopping
      for (int im = 0; im < _ml; ++im) {
        for (int jm = 0; jm < im; ++jm) {
          for (int is = 0; is < _ms; ++is) {
            if (std::abs(_tk[im][jm]) > 1e-10) {
              _F_states.emplace_back(im, jm, is, _tk[im][jm]);
              _F_states.emplace_back(jm, im, is, _tk[im][jm]);
            }
          }
        }
      }
      // Impurity-bath hybridisation
      for (int i = 0; i < _ml; ++i) {
        for (int ik = 0; ik < static_cast<int>(_Vk[i].size()); ++ik) {
          for (int is = 0; is < _ms; ++is) {
            if (std::abs(_Vk[i][ik][is]) > 1e-10) {
              int imk = _ml + ik;
              _F_states.emplace_back(i,   imk, is, _Vk[i][ik][is]);
              _F_states.emplace_back(imk, i,   is, _Vk[i][ik][is]);
            }
          }
        }
        // Bosonic bath
        for (int ib = 0; ib < static_cast<int>(_W[i].size()); ++ib) {
          if (std::abs(_W[i][ib]) > 1e-10) {
            _B_states.emplace_back(ib, i, _W[i][ib], ep.nbbits, ep.avg, false);
            _B_states.emplace_back(ib, i, _W[i][ib], ep.nbbits, ep.avg, true);
          }
        }
      }
    }

    inline Prec diagonal(long long full_state) const {
      Prec xtemp = Prec(0);
      long long bosons = full_state & ((1ll << _Nb) - 1);

      for (int is = 0; is < _ms; ++is) {
        for (int ik = 0; ik < static_cast<int>(_Epsk.size()); ++ik) {
          int ikm = _ml + ik;
          xtemp += _Epsk[ik][is] * checkState(full_state, ikm + is * _Ns, _Ip);
        }
        for (int i = 0; i < _ml; ++i) {
          xtemp += (_Eps[i][is] - _xmu) * checkState(full_state, i + is * _Ns, _Ip);
        }
      }

      for (int i = 0; i < _ml; ++i) {
        xtemp += _U[i][i] * checkState(full_state, i, _Ip) * checkState(full_state, _Ns + i, _Ip);
        for (int j = 0; j < _ml; ++j) {
          if (i != j) {
            xtemp += Prec(0.5)
                * (_U[i][j] * checkState(full_state, i,       _Ip) * checkState(full_state, _Ns + j, _Ip)
                 + _U[i][j] * checkState(full_state, _Ns + i, _Ip) * checkState(full_state, j,       _Ip)
                 + _U[i][j] * checkState(full_state, i,       _Ip) * checkState(full_state, j,       _Ip)
                 + _U[i][j] * checkState(full_state, _Ns + i, _Ip) * checkState(full_state, _Ns + j, _Ip));
          }
        }
      }

      int bit_cutoff = _Nb > 0 ? (_Nb / static_cast<int>(_w0.size())) : 0;
      long long cutoff = (1ll << bit_cutoff) - 1;
      for (int i = 0; i < static_cast<int>(_w0.size()); ++i) {
        long long cbos = (bosons >> (bit_cutoff * i)) & cutoff;
        xtemp += static_cast<Prec>(cbos) * _w0[i];
      }
      return xtemp;
    }

    inline int  valid(const St& state, long long nst) const { return state.valid(nst, _Ns, _Nb); }
    inline Prec set  (const St& state, long long nst, long long& k, int& sign) const {
      return state.set(nst, k, sign, _Ns, _Nb);
    }

    inline int checkState(long long nst, int im, int Ip) const {
      return static_cast<int>(((nst >> _Nb) & (1ll << (Ip - 1 - im))) >> (Ip - 1 - im));
    }

    inline void a(int i, long long jold, long long& k, int& isign) const {
      long long sign = 0;
      long long bos  = jold & ((1ll << _Nb) - 1);
      jold >>= _Nb;
      for (int ll = 0; ll < i; ++ll) sign += ((jold & (1ll << (_Ip - ll - 1))) != 0) ? 1 : 0;
      isign = (sign % 2) == 0 ? 1 : -1;
      k = ((jold - (1ll << (_Ip - i - 1))) << _Nb) + bos;
    }
    inline void adag(int i, long long jold, long long& k, int& isign) const {
      long long sign = 0;
      long long bos  = jold & ((1ll << _Nb) - 1);
      jold >>= _Nb;
      for (int ll = 0; ll < i; ++ll) sign += ((jold & (1ll << (_Ip - ll - 1))) != 0) ? 1 : 0;
      isign = (sign % 2) == 0 ? 1 : -1;
      k = ((jold + (1ll << (_Ip - i - 1))) << _Nb) + bos;
    }

    std::size_t bos_dim() const { return _w0.size(); }

    int number_of_bosons(long long nst, int orb) const {
      long long bosons = nst & ((1ll << _Nb) - 1);
      int bit_cutoff = _Nb > 0 ? (_Nb / static_cast<int>(_w0.size())) : 0;
      long long cutoff = (1ll << bit_cutoff) - 1;
      long long cbos = (bosons >> (bit_cutoff * orb)) & cutoff;
      return static_cast<int>(cbos);
    }

    const std::vector<HSt>& T_states() const { return _F_states; }
    const std::vector<BSt>& V_states() const { return _B_states; }

    int  interacting_orbitals() const { return _ml; }

    SzSymmetryWithBoson&       symmetry()       { return _symmetry; }
    const SzSymmetryWithBoson& symmetry() const { return _symmetry; }

    template <class Mesh>
    void bare_greens_function(edlib::Gf<std::complex<double>, 3>& /*bare_gf*/,
                              const Mesh& /*mesh*/,
                              double /*beta*/) const {
      // Mirror of legacy: stub. To be implemented when needed.
    }

  private:
    SzSymmetryWithBoson                           _symmetry;
    int                                           _Nb;
    int                                           _ml;
    Prec                                          _xmu;
    std::vector<std::vector<Prec>>                _U;
    std::vector<std::vector<Prec>>                _Eps;
    std::vector<std::vector<Prec>>                _tk;
    std::vector<std::vector<std::vector<Prec>>>   _Vk;
    std::vector<std::vector<Prec>>                _Epsk;
    std::vector<Prec>                             _w0;
    std::vector<std::vector<Prec>>                _W;

    std::vector<HSt>                              _F_states;
    std::vector<BSt>                              _B_states;
  };

}}  // namespace edlib::ext

#endif
