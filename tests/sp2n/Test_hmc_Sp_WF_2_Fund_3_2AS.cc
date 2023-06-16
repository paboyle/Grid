#include <Grid/Grid.h>

int main(int argc, char **argv) {
  using namespace Grid;
   
  typedef Representations< SpFundamentalRepresentation, SpTwoIndexAntiSymmetricRepresentation > TheRepresentations;

  Grid_init(&argc, &argv);
    
  typedef GenericSp2nHMCRunnerHirep<TheRepresentations, MinimumNorm2> HMCWrapper;
    
  typedef SpWilsonTwoIndexAntiSymmetricImplR TwoIndexFermionImplPolicy;
  typedef SpWilsonTwoIndexAntiSymmetricFermionD TwoIndexFermionAction;
  typedef typename TwoIndexFermionAction::FermionField TwoIndexFermionField;
    
  typedef SpWilsonImplR FundFermionImplPolicy;                                    // ok
  typedef SpWilsonFermionD FundFermionAction;                                     // ok
  typedef typename FundFermionAction::FermionField FundFermionField;
 
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   
  HMCWrapper TheHMC;
    
  TheHMC.Resources.AddFourDimGrid("gauge");
    
  // Checkpointer definition
  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = 5;
  CPparams.format = "IEEE64BIG";
    
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
    
      typedef PolyakovMod<HMCWrapper::ImplPolicy> PolyakovObs;
      TheHMC.Resources.AddObservable<PolyakovObs>();
    

    
  RealD beta = 6 ;
    
  SpWilsonGaugeActionR Waction(beta);
    
  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
    
  SpFundamentalRepresentation::LatticeField fundU(GridPtr);
  SpTwoIndexAntiSymmetricRepresentation::LatticeField asU(GridPtr);
    //LatticeGaugeField U(GridPtr);
    
  RealD Fundmass = -0.71;
  RealD ASmass = -0.71;
  std::vector<Complex> boundary = {-1,-1,-1,-1};

  FundFermionAction::ImplParams bc(boundary);
  TwoIndexFermionAction::ImplParams bbc(boundary);
    
  FundFermionAction FundFermOp(fundU, *GridPtr, *GridRBPtr, Fundmass, bbc);
  TwoIndexFermionAction TwoIndexFermOp(asU, *GridPtr, *GridRBPtr, ASmass, bbc);
  ConjugateGradient<FundFermionField> fCG(1.0e-8, 2000, false);
  ConjugateGradient<TwoIndexFermionField> asCG(1.0e-8, 2000, false);
  OneFlavourRationalParams Params(1.0e-6, 64.0, 2000, 1.0e-6, 16);
    
    
  TwoFlavourPseudoFermionAction<FundFermionImplPolicy> fundNf2(FundFermOp, fCG, fCG);
  TwoFlavourPseudoFermionAction<TwoIndexFermionImplPolicy> asNf2(TwoIndexFermOp, asCG, asCG);
  OneFlavourRationalPseudoFermionAction<TwoIndexFermionImplPolicy> asNf1(TwoIndexFermOp,Params);
    
  fundNf2.is_smeared = false;
  asNf2.is_smeared = false;
  asNf1.is_smeared = false;
  
    ActionLevel<HMCWrapper::Field, TheRepresentations > Level1(1);
    Level1.push_back(&fundNf2);
    Level1.push_back(&asNf2);
    Level1.push_back(&asNf1);
    
    ActionLevel<HMCWrapper::Field, TheRepresentations > Level2(4);
    Level2.push_back(&Waction);
    
    TheHMC.TheAction.push_back(Level1);
    TheHMC.TheAction.push_back(Level2);
    
    TheHMC.Parameters.MD.MDsteps = 28;
    TheHMC.Parameters.MD.trajL   = 1.0;
    
    TheHMC.ReadCommandLine(argc, argv);
    TheHMC.Run();
    
    
  Grid_finalize();
}

