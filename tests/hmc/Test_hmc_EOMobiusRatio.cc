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

namespace Grid{
  struct FermionParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(FermionParameters,
				    int, Ls,
				    double, mass,
				    double, M5,
				    double, b,
				    double, c,
				    double, StoppingCondition,
				    int, MaxCGIterations,
				    bool, ApplySmearing);
  };

  
  struct MobiusHMCParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(MobiusHMCParameters,
				  double, gauge_beta,
				  FermionParameters, Mobius)

  template <class ReaderClass >
  MobiusHMCParameters(Reader<ReaderClass>& Reader){
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

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  typedef WilsonImplR FermionImplPolicy;
  typedef MobiusFermionR FermionAction;
  typedef typename FermionAction::FermionField FermionField;
  // Serialiser
  //typedef Grid::XmlReader       Serialiser;
  typedef Grid::JSONReader       Serialiser;
  
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  HMCWrapper TheHMC;
  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
 
  // Reader, file should come from command line
  if (TheHMC.ParameterFile.empty()){
    std::cout << "Input file not specified."
              << "Use --ParameterFile option in the command line.\nAborting" 
              << std::endl;
    exit(1);
  }
  Serialiser Reader(TheHMC.ParameterFile);

  MobiusHMCParameters MyParams(Reader);  
  // Apply smearing to the fermionic action
  bool ApplySmearing = MyParams.Mobius.ApplySmearing;
  
  
  // Use this if you want to tweak the default decomposition
  // commented out as very architecture speficic
  
  //std::vector<int> simd_lanes({2,2,1,1});

  // Grid from the command line arguments --grid and --mpi
  // drop the simd_lanes argument to fall back to the default decomposition for the SIMD lanes
  
  //TheHMC.Resources.AddFourDimGrid("gauge", simd_lanes); // tweak the SIMD lanes
  TheHMC.Resources.AddFourDimGrid("gauge"); // use default simd lanes decomposition
  
  // Possibile to create the module by hand 
  // hardcoding parameters or using a Reader


  // Checkpointer definition (Name: Checkpointer)
  CheckpointerParameters CPparams(Reader);
  // Commenting out since we are using the reader 
  /*
  CPparams.config_prefix = "ckpoint_EODWF_lat";
  CPparams.rng_prefix = "ckpoint_EODWF_rng";
  CPparams.saveInterval = 5;
  CPparams.format = "IEEE64BIG";
  */
  
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  // RNG definition (Name: RandomNumberGenerator)
  RNGModuleParameters RNGpar(Reader);
  // Commenting out since we are using the reader
  /*
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  */
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  // here there is too much indirection 
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  //////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes 
  // that have a complex construction
  // standard

  //RealD beta = 5.6 ;
  WilsonGaugeActionR Waction(MyParams.gauge_beta);
    

  //const int Ls = 8;
  const int Ls = MyParams.Mobius.Ls;
  auto GridPtr   = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
  auto FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtr);
  auto FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtr);


  // temporarily need a gauge field
  LatticeGaugeField U(GridPtr);

  Real mass = MyParams.Mobius.mass; //0.04;
  Real pv   = 1.0;
  RealD M5  = MyParams.Mobius.M5; //1.5;
  // Note: IroIro and Grid notation for b and c differ
  RealD b   = MyParams.Mobius.b; //  3./2.;
  RealD c   = MyParams.Mobius.c; //  1./2.;

  // These lines are unecessary if BC are all periodic
  std::vector<Complex> boundary = {1,1,1,-1};
  FermionAction::ImplParams Params(boundary);
  
  FermionAction DenOp(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,mass,M5,b,c, Params);
  FermionAction NumOp(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,pv,  M5,b,c, Params);

  //double StoppingCondition = 1e-8;
  //double MaxCGIterations = 10000;
  ConjugateGradient<FermionField>  CG(MyParams.Mobius.StoppingCondition,MyParams.Mobius.MaxCGIterations);
  TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy> Nf2(NumOp, DenOp,CG,CG);

  // Set smearing (true/false), default: false
  Nf2.is_smeared = ApplySmearing;
  
  // Collect actions
  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Nf2);

  ActionLevel<HMCWrapper::Field> Level2(4);
  Level2.push_back(&Waction);

  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);

  /////////////////////////////////////////////////////////////
  // HMC parameters are serialisable
  TheHMC.Parameters.initialize(Reader);
  /*
  TheHMC.Parameters.MD.MDsteps = 20;
  TheHMC.Parameters.MD.trajL   = 1.0;
  */

  // Reset performance counters 
  NumOp.ZeroCounters();
  DenOp.ZeroCounters();

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

  std::cout << GridLogMessage << "Numerator report, Pauli-Villars term         : " << std::endl;
  NumOp.Report();
  std::cout << GridLogMessage << "Denominator report, Dw(m) term (includes CG) : " << std::endl;
  DenOp.Report();

  Grid_finalize();
} // main



/* Examples for input files

JSON

{
    "Checkpointer": {
    "config_prefix": "ckpoint_json_lat",
    "rng_prefix": "ckpoint_json_rng",
    "saveInterval": 1,
    "format": "IEEE64BIG"
    },
    "RandomNumberGenerator": {
    "serial_seeds": "1 2 3 4 6",
    "parallel_seeds": "6 7 8 9 11"
    },
    "Action":{
    "gauge_beta": 5.6,
    "Mobius": {
        "Ls"  : 10,
	"mass": 0.01,
	"M5"  : 1.0,
	"b"   : 1.5,
	"c"   : 0.5,
	"StoppingCondition": 1e-8,
	"MaxCGIterations": 10000,
	"ApplySmearing": true
       }
    },
    "HMC":{
    "StartTrajectory": 0,
    "Trajectories": 100,
    "MetropolisTest": true,
    "NoMetropolisUntil": 10,
    "StartingType": "HotStart",
    "MD":{
        "name": "MinimumNorm2",
	"MDsteps": 15,
	"trajL": 2.0
	}
    },
    "StoutSmearing":{
    "rho": 0.1,
    "Nsmear": 3
    }
}


XML example not provided yet

*/
