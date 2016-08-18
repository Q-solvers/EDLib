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
//  CSRHubbardHamiltonian ham(params);
  SOCSRHubbardHamiltonian ham(params);
  ham.diag();
//  GreensFunction<double, CSRHubbardHamiltonian > greensFunction(params, ham);
  GreensFunction<double, SOCSRHubbardHamiltonian > greensFunction(params, ham);
  greensFunction.compute();
  return 0;
}
