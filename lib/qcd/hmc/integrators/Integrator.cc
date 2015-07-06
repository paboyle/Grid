/*! 
  @file Integrator_base.cc
  @brief utilities for MD including funcs to generate initial HMC momentum
 */

#include <Grid.h>

namespace Grid{
  namespace QCD{

    void MDutils::generate_momenta(LatticeLorentzColourMatrix& P,GridParallelRNG& pRNG){
      // for future support of different groups
      MDutils::generate_momenta_su3(P, pRNG);
    }

    void MDutils::generate_momenta_su3(LatticeLorentzColourMatrix& P,GridParallelRNG& pRNG){
      LatticeColourMatrix Pmu(P._grid);
      Pmu = zero;
      for(int mu=0;mu<Nd;mu++){
	SU3::GaussianLieAlgebraMatrix(pRNG, Pmu);
	pokeLorentz(P, Pmu, mu);
      }
      
    }


  }
}
