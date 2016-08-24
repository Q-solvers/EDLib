#include <iostream>


#include <EDParams.h>
#include <Hamiltonian.h>
#include <SzSymmetry.h>
#include <SOCRSStorage.h>
#include <CRSStorage.h>
#include <HubbardModel.h>
#include <GreensFunction.h>
#include <SpinResolvedStorage.h>
#include "HolsteinAndersonModel.h"


int main(int argc, const char ** argv) {
  EDParams params(argc, argv);
  if(params.help_requested(std::cout)) {
    exit(0);
  }
//  CSRHubbardHamiltonian_float ham(params);
  Hamiltonian<float, SpinResolvedStorage<float, HubbardModel<float> > , HubbardModel<float> > ham(params);
//  SOCSRHubbardHamiltonian_float ham(params);
  ham.diag();
//  GreensFunction<float, CSRHubbardHamiltonian_float > greensFunction(params, ham);
  GreensFunction<float, Hamiltonian<float, SpinResolvedStorage<float, HubbardModel<float> > , HubbardModel<float> > > greensFunction(params, ham);
  greensFunction.compute();
  return 0;
}
