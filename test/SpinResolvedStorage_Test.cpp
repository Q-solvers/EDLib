//
// Created by iskakoff on 23/08/16.
//

#include <gtest/gtest.h>
#include "edlib/HubbardModel.h"

#include "edlib/SpinResolvedStorage.h"
#include "edlib/SingleImpurityAndersonModel.h"
#include "edlib/CRSStorage.h"

TEST(SpinResolvedStorageTest, Init) {
  alps::params p;
  p["NSITES"]=8;
  p["storage.MAX_SIZE"] = 80000;
  p["storage.MAX_DIM"] = 4900;
  EDLib::Model::HubbardModel<double> m(p);
  EDLib::Model::HubbardModel<double> m2(p);
  EDLib::Storage::SpinResolvedStorage<double, EDLib::Model::HubbardModel<double> > storage(p, m);
  EDLib::Storage::CRSStorage<double, EDLib::Model::HubbardModel<double> > storage2(p, m2);
  while (m.symmetry().next_sector() && m2.symmetry().next_sector()) {
    std::vector<double> v(m.symmetry().sector().size(), 1.0);
    std::vector<double> w(m.symmetry().sector().size(), 0.0);
    std::vector<double> w2(m.symmetry().sector().size(), 0.0);
    m.symmetry().init();
    m2.symmetry().init();
    storage.fill();
    storage2.fill();
    storage.av(v.data(), w.data(), m.symmetry().sector().size());
    storage2.av(v.data(), w2.data(), m.symmetry().sector().size());
    for(int i = 0; i< w.size(); ++i) {
      ASSERT_EQ(w[i], w2[i]);
    }
  }
}
