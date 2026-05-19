#ifndef EDLIB_COMBINATION_H
#define EDLIB_COMBINATION_H

#include <cstdint>
#include <vector>

namespace edlib {

  class Combination {
  public:
    explicit Combination(int N) : _c_n_k(N + 1, std::vector<int>(N + 1, 0)) {
      for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
          _c_n_k[i][j] = static_cast<int>(C_n_k_i(i, j));
        }
      }
    }

    inline int c_n_k(int n, int k) const { return _c_n_k[n][k]; }

    inline void init_state(int ik, std::vector<int>& vec) {
      for (int i = 0; i < ik; ++i) vec[i] = i;
    }

    inline bool next_combination(int n, int k, std::vector<int>& old) {
      for (int i = k - 1; i >= 0; --i) {
        if (old[i] < (n - 1 - k + (i + 1))) {
          old[i] += 1;
          for (int j = i + 1; j < k; ++j) old[j] = old[j - 1] + 1;
          return true;
        }
      }
      return false;
    }

  private:
    std::vector<std::vector<int>> _c_n_k;

    static std::int64_t variation(std::int64_t n1, std::int64_t n2) {
      std::int64_t r = 1;
      for (std::int64_t i = n1; i <= n2; ++i) r *= i;
      return r;
    }

    static std::int64_t C_n_k_i(int n, int k) {
      if ((n - k) > k) return variation(n - k + 1, n) / variation(1, k);
      return variation(k + 1, n) / variation(1, n - k);
    }
  };

}

#endif
