//
// Created by iskakoff on 23/08/16.
//

#include <gtest/gtest.h>
#include <HubbardModel.h>

#include "SpinResolvedStorage.h"
#include "SingleImpurityAndersonModel.h"


TEST(SpinResolvedStorageTest, Init) {
  EDParams p;
  HubbardModel<double> m(p);
  SpinResolvedStorage<double, HubbardModel<double> > storage(p, m);
  while (m.symmetry().next_sector()) {
    m.symmetry().init();
    storage.fill();
    storage.diag();
  }
}