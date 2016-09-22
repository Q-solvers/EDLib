//
// Created by iskakoff on 22/08/16.
//

#include "gtest/gtest.h"

#include "edlib/NSymmetryWithBoson.h"

TEST(NSymmetryWithBosonTest, Initialization) {
  EDLib::EDParams p;
  EDLib::Symmetry::NSymmetryWithBoson sym(p);
  while(sym.next_sector()) {
    sym.init();
    int i = 0;
    while (sym.next_state()) {
      ASSERT_EQ(i, sym.index(sym.state()));
      ++i;
    }
  }
}