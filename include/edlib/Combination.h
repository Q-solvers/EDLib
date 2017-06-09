//
// Created by iskakoff on 21/08/16.
//

#ifndef HUBBARD_COMBINATION_H
#define HUBBARD_COMBINATION_H


#include <vector>

namespace EDLib {

  class Combination {
  public:
    Combination(int N) : _c_n_k(N + 1, std::vector < int >(N + 1, 0)) {
      for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
          _c_n_k[i][j] = C_n_k_i(i, j);
        }
      }
    }

    inline int c_n_k(int n, int k) const {
      return _c_n_k[n][k];
    }

    /**
     * reset to initial state
     */
    inline void init_state(int ik, std::vector < int > &vec) {
      for (int i = 0; i < ik; i++) {
        vec[i] = i;
      }
    }

    /**
     * compute next combination in lexicographicaly ordered basis
     */
    inline bool next_combination(int n, int k, std::vector < int > &old) {
      for (int i = k - 1; i >= 0; i--) {
        if (old[i] < (n - 1 - k + (i + 1))) {
          old[i] += 1;
          for (int j = i + 1; j < k; j++) {
            old[j] = old[j - 1] + 1;
          }
          return true;
        }
      }
      return false;
    }

  private:
    /// chached values for combination k of n
    std::vector < std::vector < int > > _c_n_k;

    /**
     * Calculate number of combinations:
     * C_k^n = n!/(k!*(n-k)!)
     */
    int C_n_k_i(int n, int k) {
      if ((n - k) > k) {
        return variation(n - k + 1, n) / variation(1, k);
      }
      return variation(k + 1, n) / variation(1, n - k);
    }

    /**
     * Calculate n2!/n1!
     */
    int variation(int n1, int n2) {
      int result = 1;
      for (int i = n1; i <= n2; i++) {
        result *= i;
      }
      return result;
    }
  };
}
#endif //HUBBARD_COMBINATION_H
