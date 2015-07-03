/*! 
  @file Integrator_base.cc
  @brief utilities for MD including funcs to generate initial HMC momentum
 */


#include <Grid.h>

static const double sq3i = 1.0/sqrt(3.0);

namespace Grid{
  namespace QCD{

    

    void MDutils::generate_momenta_su3(LatticeColourMatrix& P,GridParallelRNG& pRNG){
      SU3::GaussianLieAlgebraMatrix(pRNG, P);
    }


  }
}
