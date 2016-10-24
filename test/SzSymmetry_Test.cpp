//
// Created by iskakoff on 22/07/16.
//

#include "gtest/gtest.h"

#include "edlib/SzSymmetry.h"
#include "edlib/EDParams.h"


TEST(SzSymmetryTest, Combinatorics) {
  alps::params p;
  EDLib::Symmetry::SzSymmetry sym(p);
  sym.init();

  ASSERT_EQ(sym.comb().c_n_k(3, 2), 3);
}


TEST(SzSymmetryTest, States) {
  alps::params p;
  EDLib::Symmetry::SzSymmetry sym(p);
  sym.init();
}

TEST(SzSymmetryTest, Initialization) {
  alps::params p;
  EDLib::define_parameters(p);
  EDLib::Symmetry::SzSymmetry sym(p);
  while(sym.next_sector()) {
    sym.init();
    int i = 0;
    while (sym.next_state()) {
      ASSERT_EQ(i, sym.index(sym.state()));
      ++i;
    }
  }
}
