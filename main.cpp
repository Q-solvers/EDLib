#include <iostream>


#include <EDParams.h>
#include <Hamiltonian.h>
#include <SzSymmetry.h>
#include <SOCRSStorage.h>
#include <CRSStorage.h>
#include <HubbardModel.h>
#include <GreensFunction.h>


int main(int argc, const char ** argv) {
  EDParams params(argc, argv);
  if(params.help_requested(std::cout)) {
    exit(0);
  }
  Hamiltonian<double, SOCRSStorage<double, HubbardModel<double> > , HubbardModel<double> > ham(params);
//  Hamiltonian<double, CRSStorage<double,HubbardModel<double> > , HubbardModel<double> > ham(params);
  ham.diag();
  GreensFunction<double, Hamiltonian<double, SOCRSStorage<double, HubbardModel<double> > , HubbardModel<double> > > greensFunction(params, ham);
  greensFunction.compute();
  return 0;
}
