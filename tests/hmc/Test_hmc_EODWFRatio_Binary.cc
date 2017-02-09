/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_EODWFRatio.cc

Copyright (C) 2015-2016

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Guido Cossu <guido.cossu@ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>



int main(int argc, char **argv) {
  using namespace Grid;
  using namespace Grid::QCD;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  typedef WilsonImplR FermionImplPolicy;
  typedef DomainWallFermionR FermionAction;
  typedef typename FermionAction::FermionField FermionField;


  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  HMCWrapper TheHMC;

  // Grid from the command line
  TheHMC.Resources.AddFourDimGrid("gauge");
  // Possibile to create the module by hand 
  // hardcoding parameters or using a Reader


  // Checkpointer definition
  CheckpointerParameters CPparams;  
  CPparams.config_prefix = "ckpoint_EODWF_lat";
  CPparams.rng_prefix = "ckpoint_EODWF_rng";
  CPparams.saveInterval = 5;
  CPparams.format = "IEEE64BIG";
  
  TheHMC.Resources.LoadBinaryCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.SerialSeed = {1,2,3,4,5};
  RNGpar.ParallelSeed = {6,7,8,9,10};
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  // here there is too much indirection 
  PlaquetteObsParameters PlPar;
  PlPar.output_prefix = "Plaquette";
  PlaquetteMod<HMCWrapper::ImplPolicy> PlaqModule(PlPar);
  TheHMC.Resources.AddObservable(&PlaqModule);
  //////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes 
  // that have a complex construction
  // standard
  RealD beta = 5.6 ;
  WilsonGaugeActionR Waction(beta);
    
  auto GridPtr   = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
  auto FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtr);
  auto FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtr);


  // temporarily need a gauge field
  LatticeGaugeField U(GridPtr);

  Real mass = 0.04;
  Real pv   = 1.0;
  RealD M5  = 1.5;

  FermionAction DenOp(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  FermionAction NumOp(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,pv,  M5);

  double StoppingCondition = 1.0e-8;
  double MaxCGIterations = 10000;
  ConjugateGradient<FermionField>  CG(StoppingCondition,MaxCGIterations);
  TwoFlavourEvenOddRatioPseudoFermionAction<ImplPolicy> Nf2(NumOp, DenOp,CG,CG);

    // Set smearing (true/false), default: false
  Nf2.is_smeared = true;

  // Collect actions
  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Nf2);

  ActionLevel<HMCWrapper::Field> Level2(4);
  Level2.push_back(&Waction);

  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  /////////////////////////////////////////////////////////////

  /*
    double rho = 0.1;  // smearing parameter
    int Nsmear = 2;    // number of smearing levels
    Smear_Stout<HMCWrapper::ImplPolicy> Stout(rho);
    SmearedConfiguration<HMCWrapper::ImplPolicy> SmearingPolicy(
        UGrid, Nsmear, Stout);
  */

  // HMC parameters are serialisable 
  TheHMC.Parameters.MD.MDsteps = 20;
  TheHMC.Parameters.MD.trajL   = 1.0;

  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file

  // Reset performance counters 
  NumOp.ZeroCounters();
  DenOp.ZeroCounters();
  TheHMC.Run();  // no smearing
  // TheHMC.Run(SmearingPolicy); // for smearing

  std::cout << GridLogMessage << "Numerator report, Pauli-Villars term         : " << std::endl;
  NumOp.Report();
  std::cout << GridLogMessage << "Denominator report, Dw(m) term (includes CG) : " << std::endl;
  DenOp.Report();




  Grid_finalize();




} // main




// Derive from the BinaryHmcRunner (templated for gauge fields)
class HmcRunner : public BinaryHmcRunner {
public:
  void BuildTheAction(int argc, char **argv)

  {
    typedef WilsonImplR ImplPolicy;
    //typedef ScaledShamirFermion<ImplPolicy> FermionAction;
    typedef DomainWallFermionR FermionAction;
    typedef typename FermionAction::FermionField FermionField;

    const int Ls = 8;

    UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
    UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
    
    FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
    FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

    // temporarily need a gauge field
    LatticeGaugeField  U(UGrid);

    // Gauge action
    double beta = 5.6;
    WilsonGaugeActionR Waction(beta);

    Real mass = 0.04;
    Real pv   = 1.0;
    RealD M5  = 1.5;
    RealD scale = 2.0;
    /*
    FermionAction DenOp(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,scale);
    FermionAction NumOp(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,pv,M5,scale);
    */
    FermionAction DenOp(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
    FermionAction NumOp(U,*FGrid,*FrbGrid,*UGrid,*UrbGrid,pv,M5);

    double StoppingCondition = 1.0e-8;
    double MaxCGIterations = 10000;
    ConjugateGradient<FermionField>  CG(StoppingCondition,MaxCGIterations);
    TwoFlavourEvenOddRatioPseudoFermionAction<ImplPolicy> Nf2(NumOp, DenOp,CG,CG);

    // Set smearing (true/false), default: false
    Nf2.is_smeared = true;

    // Collect actions
    // here an example of 2 level integration
    ActionLevel<Field> Level1(1);
    Level1.push_back(&Nf2);
    Level1.push_back(&Waction);

    // this level will integrate with a
    // step that is 4 times finer
    // than the previous level
    //ActionLevel<Field> Level2(4);
    

    TheAction.push_back(Level1);
    //TheAction.push_back(Level2);

    // Add observables
    int SaveInterval = 1;
    std::string format = std::string("IEEE64BIG");
    std::string conf_prefix = std::string("DWF_ckpoint_lat");
    std::string rng_prefix = std::string("DWF_ckpoint_rng");
    BinaryHmcCheckpointer<BinaryHmcRunner::ImplPolicy> Checkpoint(conf_prefix, rng_prefix, SaveInterval, format);
    // Can implement also a specific function in the hmcrunner
    // AddCheckpoint (...) that takes the same parameters + a string/tag
    // defining the type of the checkpointer
    // with tags can be implemented by overloading and no ifs
    // Then force all checkpoint to have few common functions
    // return an object that is then passed to the Run function

    PlaquetteLogger<BinaryHmcRunner::ImplPolicy> PlaqLog(std::string("Plaquette"));
    ObservablesList.push_back(&PlaqLog);
    ObservablesList.push_back(&Checkpoint);

    // Smearing section, omit if not needed
    double rho = 0.1;  // smearing parameter
    int Nsmear = 2;    // number of smearing levels
    Smear_Stout<BinaryHmcRunner::ImplPolicy> Stout(rho);
    SmearedConfiguration<BinaryHmcRunner::ImplPolicy> SmearingPolicy(
        UGrid, Nsmear, Stout);
    ///////////////////

    NumOp.ZeroCounters();
    DenOp.ZeroCounters();
    //Run(argc, argv, Checkpoint, SmearingPolicy); 
    Run(argc, argv, Checkpoint);  // no smearing



    std::cout << GridLogMessage << "Numerator report, Pauli-Villars term         : " << std::endl;
    NumOp.Report();
    std::cout << GridLogMessage << "Denominator report, Dw(m) term (includes CG) : " << std::endl;
    DenOp.Report();
};
};
}
}
