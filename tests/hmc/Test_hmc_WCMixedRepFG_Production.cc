/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_WilsonAdjointFermionGauge.cc

Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#include "Grid/Grid.h"


namespace Grid{
  struct FermionParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(FermionParameters,
            double, mass,
            double, csw,
				    double, StoppingCondition,
				    int, MaxCGIterations,
				    bool, ApplySmearing);
  };

  struct WilsonCloverHMCParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(WilsonCloverHMCParameters,
				  double, gauge_beta,
          FermionParameters, WilsonCloverFund, 
          FermionParameters, WilsonCloverAS)

  template <class ReaderClass >
  WilsonCloverHMCParameters(Reader<ReaderClass>& Reader){
    read(Reader, "Action", *this);
    }
  };

  struct SmearingParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(SmearingParameters,
				    double, rho,
				    Integer, Nsmear)

    template <class ReaderClass >
    SmearingParameters(Reader<ReaderClass>& Reader){
      read(Reader, "StoutSmearing", *this);
    }

  }; 
}


int main(int argc, char **argv) {
  using namespace Grid;
  using namespace Grid::QCD;

  // Here change the allowed (higher) representations
  typedef Representations< FundamentalRepresentation, TwoIndexAntiSymmetricRepresentation> TheRepresentations;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef GenericHMCRunnerHirep<TheRepresentations, MinimumNorm2> HMCWrapper;

  typedef WilsonImplR FundImplPolicy;
  typedef WilsonCloverFermionR FundFermionAction; 
  typedef typename FundFermionAction::FermionField FundFermionField;

  typedef WilsonTwoIndexAntiSymmetricImplR ASymmImplPolicy; 
  typedef WilsonCloverTwoIndexAntiSymmetricFermionR ASymmFermionAction; 
  typedef typename ASymmFermionAction::FermionField ASymmFermionField;

  typedef Grid::JSONReader Serialiser;
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  HMCWrapper TheHMC;
  
    // Grid from the command line
    TheHMC.ReadCommandLine(argc, argv); 
    if (TheHMC.ParameterFile.empty()){
      std::cout << "Input file not specified."
                << "Use --ParameterFile option in the command line.\nAborting" 
                << std::endl;
      exit(1);
    }
    Serialiser Reader(TheHMC.ParameterFile);
    WilsonCloverHMCParameters MyParams(Reader);  
  
    // Apply smearing to the fermionic action
    bool ApplySmearingFund = MyParams.WilsonCloverFund.ApplySmearing;
    bool ApplySmearingAS = MyParams.WilsonCloverAS.ApplySmearing;
    

    TheHMC.Resources.AddFourDimGrid("gauge");
  
    // Checkpointer definition
    CheckpointerParameters CPparams(Reader);
    
    /*
    CPparams.config_prefix = "ckpoint_lat";
    CPparams.rng_prefix = "ckpoint_rng";
    CPparams.saveInterval = 5;
    CPparams.format = "IEEE64BIG";
    */
    
    TheHMC.Resources.LoadNerscCheckpointer(CPparams);
  
    RNGModuleParameters RNGpar(Reader);
    /*
    RNGpar.serial_seeds = "1 2 3 4 5";
    RNGpar.parallel_seeds = "6 7 8 9 10";
    TheHMC.Resources.SetRNGSeeds(RNGpar);
    */
    TheHMC.Resources.SetRNGSeeds(RNGpar);
  
    // Construct observables
    typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
    TheHMC.Resources.AddObservable<PlaqObs>();
  
    typedef PolyakovMod<HMCWrapper::ImplPolicy> PolyakovObs;
    TheHMC.Resources.AddObservable<PolyakovObs>();
  
    typedef TopologicalChargeMod<HMCWrapper::ImplPolicy> QObs;
    TopologyObsParameters TopParams(Reader);
    TheHMC.Resources.AddObservable<QObs>(TopParams);
    //////////////////////////////////////////////
  
    /////////////////////////////////////////////////////////////
    // Collect actions, here use more encapsulation
    // need wrappers of the fermionic classes
    // that have a complex construction
    // standard
    
    //RealD beta = 5.6;
    WilsonGaugeActionR Waction(MyParams.gauge_beta);
    
      auto GridPtr   = TheHMC.Resources.GetCartesian();
      auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
    
      // temporarily need a gauge field
      FundamentalRepresentation::LatticeField UF(GridPtr);
      TwoIndexAntiSymmetricRepresentation::LatticeField UAS(GridPtr);


      Real Fundmass = MyParams.WilsonCloverFund.mass;
      Real Fundcsw = MyParams.WilsonCloverFund.csw;
      Real ASmass = MyParams.WilsonCloverAS.mass;
      Real AScsw = MyParams.WilsonCloverAS.csw;

      

  std::cout << "Fund: mass and csw" << Fundmass << " and " << Fundcsw << std::endl; 
  std::cout << "AS  : mass and csw" << ASmass << " and " << AScsw << std::endl; 
  
  
  FundFermionAction FundFermOp(UF, *GridPtr, *GridRBPtr, Fundmass, Fundcsw, Fundcsw);
  ConjugateGradient<FundFermionField> CG_Fund(MyParams.WilsonCloverFund.StoppingCondition, MyParams.WilsonCloverFund.MaxCGIterations);
  TwoFlavourPseudoFermionAction<FundImplPolicy> Nf2_Fund(FundFermOp, CG_Fund, CG_Fund);

  ASymmFermionAction ASFermOp(UAS, *GridPtr, *GridRBPtr, ASmass, AScsw, AScsw);
  ConjugateGradient<ASymmFermionField> CG_AS(MyParams.WilsonCloverAS.StoppingCondition, MyParams.WilsonCloverAS.MaxCGIterations);
  TwoFlavourPseudoFermionAction<ASymmImplPolicy> Nf2_AS(ASFermOp, CG_AS, CG_AS);

  Nf2_Fund.is_smeared = ApplySmearingFund;
  Nf2_AS.is_smeared   = ApplySmearingAS;
  

  // Collect actions
  ActionLevel<HMCWrapper::Field, TheRepresentations > Level1(1);
  Level1.push_back(&Nf2_Fund);
  Level1.push_back(&Nf2_AS);


  ActionLevel<HMCWrapper::Field, TheRepresentations > Level2(4);
  Level2.push_back(&Waction);

  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);

  TheHMC.Parameters.initialize(Reader);
  //TheHMC.Parameters.MD.MDsteps = 20;
  //TheHMC.Parameters.MD.trajL = 1.0;
/*
  if (ApplySmearingFund || ApplySmearingAS){
    SmearingParameters SmPar(Reader);
    //double rho = 0.1;  // smearing parameter
    //int Nsmear = 3;    // number of smearing levels
    Smear_Stout<HMCWrapper::ImplPolicy> Stout(SmPar.rho);
    SmearedConfiguration<HMCWrapper::ImplPolicy> SmearingPolicy(GridPtr, SmPar.Nsmear, Stout);
    TheHMC.Run(SmearingPolicy); // for smearing
  } else {
    TheHMC.Run();  // no smearing
  }
*/
  TheHMC.Run(); 


  //TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  //TheHMC.Run();                       // no smearing
  // TheHMC.Run(SmearingPolicy); // for smearing

  Grid_finalize();

} // main
