//
// Created by iskakoff on 21/08/16.
//

#include "gtest/gtest.h"

#include "edlib/NSymmetry.h"
#include "edlib/EDParams.h"

TEST(NSymmetryTest, Initialization) {
  alps::params p;
  EDLib::define_parameters(p);
  EDLib::Symmetry::NSymmetry sym(p);
  while(sym.next_sector()) {
    sym.init();
    int i = 0;
    while (sym.next_state()) {
      ASSERT_EQ(i, sym.index(sym.state()));
      ++i;
    }
  }
}
