//
// Created by iskakoff on 23/08/16.
//

#include <gtest/gtest.h>
#include <HubbardModel.h>

#include "SpinResolvedStorage.h"
#include "SingleImpurityAndersonModel.h"


#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <CRSStorage.h>

TEST(SpinResolvedStorageTest, Init) {
  EDParams p;
  p["NSITES"]=8;
  HubbardModel<double> m(p);
  SpinResolvedStorage<double, HubbardModel<double> > storage(p, m);
//  while (m.symmetry().next_sector()) {
//    m.symmetry().init();
//    storage.fill();
//    storage.diag();
//  }
  m.symmetry().set_sector(SzSymmetry::Sector(1, 6, 224));
  m.symmetry().init();
  storage.fill();
  std::vector<double> v(224, 1.0);
  std::vector<double> w(224, 0.0);
  std::vector<double> w2(224, 0.0);
  storage.av(v.data(), w.data(), 224);
//  std::cout<<"========"<<std::endl;
  CRSStorage<double, HubbardModel<double> > storage2(p, m);
  m.symmetry().init();
  storage2.fill();
  storage2.av(v.data(), w2.data(), 224);
//  std::cout<< boost::algorithm::join( v | boost::adaptors::transformed( static_cast<std::string(*)(double)>(std::to_string) ), ", " )<<std::endl;
//  std::cout<< boost::algorithm::join( w | boost::adaptors::transformed( static_cast<std::string(*)(double)>(std::to_string) ), ", " )<<std::endl;
//  std::cout<< boost::algorithm::join( w2 | boost::adaptors::transformed( static_cast<std::string(*)(double)>(std::to_string) ), ", " )<<std::endl;
//  storage.print();
//  storage2.print();
  for(int i = 0; i< w.size(); ++i) {
    ASSERT_EQ(w[i], w2[i]);
  }
}