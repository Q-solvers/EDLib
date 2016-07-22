//
// Created by iskakoff on 22/07/16.
//

#include <alps/params.hpp>
#include "gtest/gtest.h"

#include "SzSymmetry.h"


TEST(SzSymmetryTest, Combinatorics) {
  alps::params p;
  p.define("NSITES", 4, "");
  p.define("NSPINS", 2, "");
  SzSymmetry sym(p);
  sym.init();

  ASSERT_EQ(sym.C_n_k_i(3, 2), 3);
}


TEST(SzSymmetryTest, States) {
  alps::params p;
  p.define("NSITES", 2, "");
  p.define("NSPINS", 2, "");
  SzSymmetry sym(p);
  sym.init();
}