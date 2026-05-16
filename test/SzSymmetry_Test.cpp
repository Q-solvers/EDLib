// Core-only SzSymmetry parity tests.

#include <edlib/Parameters.h>
#include <edlib/SzSymmetry.h>

#include <gtest/gtest.h>

TEST(SzSymmetryCore, Combinatorics) {
  edlib::Parameters p;
  edlib::SzSymmetry sym(p);
  sym.init();
  ASSERT_EQ(sym.comb().c_n_k(3, 2), 3);
}

TEST(SzSymmetryCore, Construct) {
  edlib::Parameters p;
  edlib::SzSymmetry sym(p);
  sym.init();
}

TEST(SzSymmetryCore, IndexMatchesIteration) {
  edlib::Parameters p;
  edlib::SzSymmetry sym(p);
  while (sym.next_sector()) {
    sym.init();
    int i = 0;
    while (sym.next_state()) {
      ASSERT_EQ(i, sym.index(sym.state()));
      ++i;
    }
  }
}
