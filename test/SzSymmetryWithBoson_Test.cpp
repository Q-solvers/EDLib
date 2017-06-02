//
// Created by iskakoff on 02/06/17.
//

#include <gtest/gtest.h>

#include <alps/params.hpp>
#include <ext/HolsteinAndersonParameter.h>
#include "edlib/EDParams.h"
#include "ext/SzSymmetryWithBoson.h"

TEST(SzSymmetryWithBosonTest, Indexing) {
  alps::params p;
  EDLib::define_parameters(p);
  EDLib::Ext::define_parameters(p);
  p["NBLEVEL"] = 2;
  EDLib::Symmetry::SzSymmetryWithBoson sym(p);
  while(sym.next_sector()) {
    sym.init();
    int i = 0;
    while (sym.next_state()) {
      ASSERT_EQ(i, sym.index(sym.state()));
      ++i;
    }
  }
}
