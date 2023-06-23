#include <Grid/Grid.h>

int main(int argc, char **argv) {
  using namespace Grid;
   
  typedef Representations< SpFundamentalRepresentation > TheRepresentations;

  Grid_init(&argc, &argv);
    
  typedef GenericSpHMCRunnerHirep<TheRepresentations, MinimumNorm2> HMCWrapper; // ok
  typedef SpWilsonImplR FermionImplPolicy;                                    // ok
  typedef SpWilsonFermionD FermionAction;                                     // ok
  typedef typename FermionAction::FermionField FermionField;                  // ok?
 
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   
  HMCWrapper TheHMC;
    
  TheHMC.Resources.AddFourDimGrid("gauge");
    
  // Checkpointer definition
  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = 100;
  CPparams.format = "IEEE64BIG";
    
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
    
  RealD beta = 7.2 ;
    
  SpWilsonGaugeActionR Waction(beta);
    
  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
    
    SpFundamentalRepresentation::LatticeField U(GridPtr);
    //LatticeGaugeField U(GridPtr);
    
  RealD mass = -0.76;

  FermionAction FermOp(U, *GridPtr, *GridRBPtr, mass);
    
 
  ConjugateGradient<FermionField> CG(1.0e-8, 2000, false);
  
  TwoFlavourPseudoFermionAction<FermionImplPolicy> Nf2(FermOp, CG, CG);
    
  Nf2.is_smeared = false;
  
    ActionLevel<HMCWrapper::Field, TheRepresentations > Level1(1);
    Level1.push_back(&Nf2);
    
    ActionLevel<HMCWrapper::Field, TheRepresentations > Level2(4);
    Level2.push_back(&Waction);
    
    TheHMC.TheAction.push_back(Level1);
    TheHMC.TheAction.push_back(Level2);
    
    TheHMC.Parameters.MD.MDsteps = 36;
    TheHMC.Parameters.MD.trajL   = 1.0;
    
    TheHMC.ReadCommandLine(argc, argv);
    TheHMC.Run();
    
    
  Grid_finalize();
}

