//
// Created by iskakoff on 20/04/17.
//

#ifndef HUBBARD_HOLSTEINANDERSONMODEL_H
#define HUBBARD_HOLSTEINANDERSONMODEL_H


#include <edlib/FermionicModel.h>
#include <edlib/CommonUtils.h>
#include <alps/gf/gf.hpp>
#include <edlib/SingleImpurityAndersonModel.h>
#include "SzSymmetryWithBoson.h"

namespace EDLib {
  namespace Ext {
    namespace Model {
      namespace HolsteinAnderson {
        template<typename prec>
        class InnerState : public EDLib::Model::SingleImpurityAnderson::InnerState < prec > {
        public:
          using EDLib::Model::SingleImpurityAnderson::InnerState<prec>::a;
          using EDLib::Model::SingleImpurityAnderson::InnerState<prec>::adag;
          using EDLib::Model::SingleImpurityAnderson::InnerState<prec>::checkState;
          virtual int valid(long long, int, int) const { return 0; };

          virtual prec set(long long, long long &, int &, int, int) const {return 0.0;};
        };

        template<typename prec>
        class HybridisationInnerState : public InnerState < prec > {

        public:
          HybridisationInnerState(int i1, int i2, int is, prec value) : _indicies(i1, i2), _spin(is), _value(value) {}

          const inline std::pair < int, int > &indicies() const { return _indicies; }

          virtual inline prec value() const { return _value; }

          inline int spin() const { return _spin; }

          virtual int valid(long long nst, int Ns, int Nb) const {
            return (this->checkState(nst >> Nb, _indicies.first + _spin * Ns, Ns) * (1 - this->checkState(nst >> Nb, _indicies.second + _spin * Ns, Ns)));
          }

          virtual prec set(long long nst, long long &k, int &sign, int Ns, int Nb) const {
            long long k1, k2;
            int isign1, isign2;
            long long fnst = nst >> Nb;
            long long bnst = nst & ((1 << Nb) - 1);
            this->a(_indicies.first + _spin * Ns, fnst, k1, isign1, 2 * Ns);
            this->adag(_indicies.second + _spin * Ns, k1, k2, isign2, 2 * Ns);
            k = (k2 << Nb) + bnst;
            sign = isign1 * isign2;
            return _value;
          }

        private:
          std::pair < int, int > _indicies;
          int _spin;
          prec _value;
        };

        template<typename prec>
        class BosonInnerState : public InnerState < prec > {

        public:
          BosonInnerState(int ib, int i, prec value, int cutoff, prec avg, bool dag) : _b(ib), _i(i), _value(value), _bit_cutoff(cutoff), _cutoff((1 << cutoff) - 1), _avg(avg), _dag(dag) {}

          virtual int valid(long long nst, int Ns, int Nb) const {
            // extract bosonic part of state
            long long int bnst = nst & ((1 << Nb) - 1);
            // extract current boson from N-dimensional representation
            long long cbos = ((bnst >> (_bit_cutoff * _b)) & _cutoff);
            return (std::abs(this->checkState(nst >> Nb, _i, Ns) + this->checkState(nst >> Nb, Ns+_i, Ns) - _avg) >1e-9) * (_dag ? (cbos) < _cutoff : cbos > 0);
          }

          virtual prec set(long long nst, long long &k, int &sign, int Ns, int Nb) const {
            double N = (this->checkState(nst>>Nb, _i, Ns) + this->checkState(nst>>Nb, Ns+_i, Ns) - _avg);
            long long int bnst = nst & ((1 << Nb) - 1);
            long long cbos = ((bnst >> (_bit_cutoff * _b)) & _cutoff);
            if (_dag) {
              // create boson
              k = nst + (1 << (_bit_cutoff * _b));
              sign = 1;
              return N*_value * std::sqrt(cbos + 1);
            } else {
              // destroy boson
              k = nst - (1 << (_bit_cutoff * _b));
              sign = 1;
              return N * _value * std::sqrt(cbos);
            }
          }
        private:
          int _b;
          int _i;
          // number of bits per boson
          int _bit_cutoff;
          long long _cutoff;
          bool _dag;
          prec _avg;
          prec _value;
        };
      }
      /**
       * @brief HolsteinAndersonModel class
       *
       * @author iskakoff
       */
      template<typename prec>
      class HolsteinAndersonModel : public EDLib::Model::FermionicModel {
      public:
        typedef prec precision;
        typedef typename Symmetry::SzSymmetryWithBoson SYMMETRY;
        typedef typename HolsteinAnderson::InnerState < precision > St;
        typedef typename HolsteinAnderson::HybridisationInnerState < precision > HSt;
        typedef typename HolsteinAnderson::BosonInnerState < precision > BSt;
        typedef typename Symmetry::SzSymmetryWithBoson::Sector Sector;


        HolsteinAndersonModel(alps::params &p) : FermionicModel(p), _symmetry(p), _Nb(p["NBBITS"].as<int>() * p["NBLEVEL"].as<int>()), _ml(p["NORBITALS"].as<int>()),
               _U(p["NORBITALS"].as<int>(), std::vector<precision>(p["NORBITALS"].as<int>(), 0.0) ) {
          std::string input = p["INPUT_FILE"];
          alps::hdf5::archive input_data(input.c_str(), "r");
          if (_ms != 2) {
            throw std::invalid_argument("Incorrect values for the number of spins. Please check input file.");
          }
          _Ip = _ms * _Ns;
          std::stringstream s;
          s << "Bath/Vk" << "/values";
          input_data >> alps::make_pvp(s.str().c_str(), _Vk);
          s.str("");
          s << "Bath/Epsk" << "/values";
          input_data >> alps::make_pvp(s.str().c_str(), _Epsk);
          s.str("");
          s << "Bath/w0" << "/values";
          input_data >> alps::make_pvp(s.str().c_str(), _w0);
          s.str("");
          s << "Bath/W" << "/values";
          input_data >> alps::make_pvp(s.str().c_str(), _W);
          s.str("");
          s<<"tk/values";
          input_data >> alps::make_pvp(s.str().c_str(), _tk);
          input_data >> alps::make_pvp("Eps0/values", _Eps);
          input_data >> alps::make_pvp("mu", _xmu);
          precision U;
          input_data >> alps::make_pvp("U", _U);

          std::cout<<"ml:"<<_ml<<"\n";
          precision avg;
          input_data >> alps::make_pvp("AVG", avg);

          int cutoff = p["NBBITS"].as<int>();
          input_data.close();
          for (int i = 0; i< _ml; ++i) {
            if (_Vk[i].size() != _Epsk.size()) {
              throw std::invalid_argument("Incorrect bosonic bath. Please check input file.");
            }
            if (_w0.size() != _W[i].size()) {
              throw std::invalid_argument("Incorrect bosonic bath. Please check input file.");
            }
          }
          if(_w0.size() != p["NBLEVEL"].as<int>()) {
            throw std::invalid_argument("Incorrect bosonic bath. Number of bosonic levels in input file inconsistent with parameters. Please check input file.");
          }
          // Fermionic bath
          // interorbital hoppings
          for (int im = 0; im < _ml; ++im) {
            for (int jm = 0; jm < im; ++jm) {
              for(int is = 0; is< _ms; ++is) {
                if (std::abs(_tk[im][jm]) > 1e-10) {
                  _F_states.push_back(HSt(im, jm, is, _tk[im][jm]));
                  _F_states.push_back(HSt(jm, im, is, _tk[im][jm]));
                }
                std::cout<<"hopping\n";
              }
            }
          }
          for (int i = 0; i< _ml; ++i) {
            for (int ik = 0; ik < _Vk[i].size(); ++ik) {
              for (int is = 0; is < _ms; ++is) {
                if (std::abs(_Vk[i][ik][is]) > 1e-10) {
                  int imk = _ml + ik;
                  // hoppings
                  // from impurity to bath
                  _F_states.push_back(HSt(i, imk, is, _Vk[i][ik][is]));
                  // from bath to impurity
                  _F_states.push_back(HSt(imk, i, is, _Vk[i][ik][is]));
                }
              }
            }
            // Bosonic bath
            for (int ib = 0; ib < _W[i].size(); ++ib) {
              std::cout<<"Wb["<<i<<"]["<<ib<<"]="<<_W[i][ib]<<"\n";
              if (std::abs(_W[i][ib]) > 1e-10) {
                // destroy boson
                _B_states.push_back(BSt(ib, i, _W[i][ib], cutoff, avg, false));
                // create boson
                _B_states.push_back(BSt(ib, i, _W[i][ib], cutoff, avg, true));
              }
            }
          }
        }

        inline const precision diagonal(long long full_state) const {
          precision xtemp = 0.0;
          long long bosons = full_state & ((1 << _Nb) - 1);
          // density term
          for (int is = 0; is < _ms; ++is) {
            for (int ik = 0; ik < _Epsk.size(); ++ik) {
              int ikm = _ml + ik;
              xtemp += (_Epsk[ik][is] * checkState(full_state, ikm + is * _Ns, _Ip));
            }
            for( int i = 0; i < _ml; ++i) {
              xtemp += (_Eps[i][is] - _xmu) * checkState(full_state, i + is * _Ns, _Ip);
            }
          }
          // interaction term
          for( int i = 0; i < _ml; ++i) {
            xtemp += _U[i][i] * checkState(full_state, i, _Ip) * checkState(full_state, _Ns + i, _Ip);
            for( int j = 0; j < _ml; ++j) {
              if(i!=j) {
                xtemp += 0.5*(_U[i][j] * checkState(full_state, i, _Ip) * checkState(full_state, _Ns + j, _Ip) +
                     _U[i][j] * checkState(full_state, _Ns + i, _Ip) * checkState(full_state, j, _Ip) +
                     _U[i][j] * checkState(full_state, i, _Ip) * checkState(full_state, j, _Ip) +
                     _U[i][j] * checkState(full_state, _Ns + i, _Ip) * checkState(full_state, _Ns + j, _Ip) );
              }
            }
          }

          int bit_cutoff = _Nb > 0 ? (_Nb / _w0.size()) :0;
          long long cutoff = (1 << bit_cutoff) - 1;
          for (int i = 0; i < _w0.size(); ++i) {
            long long cbos = ((bosons >> (bit_cutoff * i)) & cutoff);
            xtemp += cbos * (_w0[i]);
          }
          return xtemp;
        }

        inline int valid(const St &state, long long nst) {
          return state.valid(nst, _Ns, _Nb);
        }

        inline prec set(const St &state, long long nst, long long &k, int &sign) {
          return state.set(nst, k, sign, _Ns, _Nb);
        }

        int inline checkState(long long nst, const int im, int Ip) const {
          return (int) (((nst>>_Nb) & (1ll << (Ip - 1 - im))) >> (Ip - 1 - im));
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
          long long bos = jold & ((1<<_Nb) - 1);
          jold >>= _Nb;
          for (int ll = 0; ll < i; ++ll) {
            sign += ((jold & (1ll << (_Ip - ll - 1))) != 0) ? 1 : 0;
          }
          isign = (sign % 2) == 0 ? 1 : -1;
          k = ((jold - (1ll << (_Ip - i - 1))) << _Nb) + bos;
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
          long long bos = jold & ((1<<_Nb) - 1);
          jold >>= _Nb;
          for (int ll = 0; ll < i; ++ll) {
            sign += ((jold & (1ll << (_Ip - ll - 1))) != 0) ? 1 : 0;
          }
          isign = (sign % 2) == 0 ? 1 : -1;
          k = ((jold + (1ll << (_Ip - i - 1))) << _Nb) + bos;
        }

        size_t bos_dim() const {
          return _w0.size();
        }

        int number_of_bosons(long long nst, int orb) const {
          long long bosons = nst & ((1 << _Nb) - 1);
          int bit_cutoff = _Nb > 0 ? (_Nb / _w0.size()) :0;
          long long cutoff = (1 << bit_cutoff) - 1;
          long long cbos = ((bosons >> (bit_cutoff * orb)) & cutoff);
          return int(cbos);
        }

        inline const std::vector < HSt > &T_states() const {
          return _F_states;
        }

        inline const std::vector < BSt > &V_states() const {
          return _B_states;
        }

        int interacting_orbitals() const {
          return _ml;
        }

        template<typename Mesh>
        void bare_greens_function(alps::gf::three_index_gf < std::complex < double >, Mesh, alps::gf::index_mesh, alps::gf::index_mesh > &bare_gf, double beta) {
        }

        inline const Symmetry::SzSymmetryWithBoson &symmetry() const {
          return _symmetry;
        }

        inline Symmetry::SzSymmetryWithBoson &symmetry() {
          return _symmetry;
        }

      private:
        Symmetry::SzSymmetryWithBoson _symmetry;
        int _Nb;
        // impurity parameters
        // number of sites in cluster impurity
        int _ml;
        // chemical potential
        precision _xmu;
        // local interaction
        std::vector< std::vector< precision > > _U;
        std::vector< std::vector< precision > > _Eps;
        // intracluster hoppings
        std::vector< std::vector < precision > > _tk;
        // fermionic bath
        std::vector< std::vector< std::vector < precision > > > _Vk;
        std::vector< std::vector < precision > > _Epsk;
        // bosonic bath
        std::vector < precision > _w0;
        std::vector< std::vector < precision > > _W;


        // Fermionic non-digaonal part of Hamiltonian
        std::vector < HSt > _F_states;
        // Bosonic non-digaonal part of Hamiltonian
        std::vector < BSt > _B_states;
      };

    }
  }
}

#endif //HUBBARD_HOLSTEINANDERSONMODEL_H
