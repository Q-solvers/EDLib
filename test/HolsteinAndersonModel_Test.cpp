//
// Created by iskakoff on 22/08/16.
//

#include <gtest/gtest.h>
#include "edlib/Hamiltonian.h"
#include "edlib/HolsteinAndersonModel.h"
#include "edlib/Storage.h"

TEST(HubbardModelTest, Filling) {
//  const char *string = "--NSITES=16 --INPUT_FILE=../input/input.h5";
  EDLib::EDParams p;
  p["INPUT_FILE"] = "./input/input.h5";
  EDLib::Model::HolsteinAndersonModel<double> model(p);
  model.diagonal(0);
}