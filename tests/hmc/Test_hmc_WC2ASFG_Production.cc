/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_WilsonFermionGauge.cc

Copyright (C) 2017

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
				  FermionParameters, WilsonClover)

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

int main(int argc, char **argv)
{
  using namespace Grid;
  using namespace Grid::QCD;

  typedef Representations< FundamentalRepresentation, TwoIndexAntiSymmetricRepresentation > TheRepresentations;  

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  // Typedefs to simplify notation
  typedef GenericHMCRunnerHirep<TheRepresentations, MinimumNorm2> HMCWrapper; // Uses the default minimum norm
  typedef WilsonTwoIndexAntiSymmetricImplR FermionImplPolicy; // gauge field implemetation for the pseudofermions
  typedef WilsonCloverTwoIndexAntiSymmetricFermionR FermionAction; // type of lattice fermions (Wilson, DW, ...)
  typedef typename FermionAction::FermionField FermionField;
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
  bool ApplySmearing = MyParams.WilsonClover.ApplySmearing;

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

  //typedef TopologicalChargeMod<HMCWrapper::ImplPolicy> QObs;
  //TopologyObsParameters TopParams(Reader);
  //TheHMC.Resources.AddObservable<QObs>(TopParams);
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
  TwoIndexAntiSymmetricRepresentation::LatticeField U(GridPtr);

  //Real mass = 0.01;
  //Real csw = 1.0;

  Real mass = MyParams.WilsonClover.mass;
  Real csw = MyParams.WilsonClover.csw;

  std::cout << "mass and csw" << mass << " and " << csw << std::endl; 

  FermionAction FermOp(U, *GridPtr, *GridRBPtr, mass, csw, csw);
  ConjugateGradient<FermionField> CG(MyParams.WilsonClover.StoppingCondition, MyParams.WilsonClover.MaxCGIterations);
  TwoFlavourPseudoFermionAction<FermionImplPolicy> Nf2(FermOp, CG, CG);

  // Set smearing (true/false), default: false
  Nf2.is_smeared = ApplySmearing;

  // Collect actions
  ActionLevel<HMCWrapper::Field, TheRepresentations> Level1(1);
  Level1.push_back(&Nf2);

  ActionLevel<HMCWrapper::Field, TheRepresentations> Level2(4);
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

  TheHMC.Parameters.initialize(Reader);
  //TheHMC.Parameters.MD.MDsteps = 20;
  //TheHMC.Parameters.MD.trajL = 1.0;

  if (ApplySmearing){
    SmearingParameters SmPar(Reader);
    //double rho = 0.1;  // smearing parameter
    //int Nsmear = 3;    // number of smearing levels
    Smear_Stout<HMCWrapper::ImplPolicy> Stout(SmPar.rho);
    SmearedConfiguration<HMCWrapper::ImplPolicy> SmearingPolicy(GridPtr, SmPar.Nsmear, Stout);
    TheHMC.Run(SmearingPolicy); // for smearing
  } else {
    TheHMC.Run();  // no smearing
  }

  //TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  //TheHMC.Run();                       // no smearing
  // TheHMC.Run(SmearingPolicy); // for smearing

  Grid_finalize();

} // main

