// Core-only NSymmetry parity test.

#include <edlib/NSymmetry.h>
#include <edlib/Parameters.h>

#include <gtest/gtest.h>

TEST(NSymmetryCore, IndexMatchesIteration) {
  edlib::Parameters p;
  edlib::NSymmetry sym(p);
  while (sym.next_sector()) {
    sym.init();
    int i = 0;
    while (sym.next_state()) {
      ASSERT_EQ(i, sym.index(sym.state()));
      ++i;
    }
  }
}
