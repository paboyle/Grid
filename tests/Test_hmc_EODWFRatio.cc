#include "Grid.h"


using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls = 8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  GridParallelRNG  pRNG(UGrid);
  pRNG.SeedRandomDevice();

  LatticeLorentzColourMatrix     U(UGrid);

  SU3::HotConfiguration(pRNG, U);

  // simplify template declaration? Strip the lorentz from the second template
  WilsonGaugeAction<LatticeLorentzColourMatrix, LatticeColourMatrix> Waction(5.6);

  Real mass=0.04;
  Real pv  =1.0;
  RealD M5=1.5;

  DomainWallFermionR DenOp(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  DomainWallFermionR NumOp(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,pv,M5);
  
  ConjugateGradient<LatticeFermion>  CG(1.0e-8,10000);
  TwoFlavourEvenOddRatioPseudoFermionAction<WilsonImplR> Nf2(NumOp, DenOp,CG,CG);
  
  //Collect actions
  ActionLevel Level1;
  Level1.push_back(&Nf2);
  Level1.push_back(&Waction);
  ActionSet FullSet;
  FullSet.push_back(Level1);

  // Create integrator
  //  typedef LeapFrog  IntegratorAlgorithm;// change here to change the algorithm
  typedef MinimumNorm2  IntegratorAlgorithm;// change here to change the algorithm
  IntegratorParameters MDpar(12,20,1.0);
  Integrator<IntegratorAlgorithm> MDynamics(UGrid,MDpar, FullSet);

  // Create HMC
  HMCparameters HMCpar;
  HybridMonteCarlo<IntegratorAlgorithm>  HMC(HMCpar, MDynamics);

  // Run it
  HMC.evolve(U);

}
