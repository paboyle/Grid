/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./HMC/DWF2p1fIwasakiGparity.cc

Copyright (C) 2015-2016

Author: Christopher Kelly <ckelly@bnl.gov>
Author: Peter Boyle <pabobyle@ph.ed.ac.uk>


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

using namespace Grid;

//2+1f DWF+I ensemble with G-parity BCs
//designed to reproduce ensembles in https://arxiv.org/pdf/1908.08640.pdf
struct RatQuoParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(RatQuoParameters,
				  double, bnd_lo,
				  double, bnd_hi,
				  Integer, action_degree,
				  double, action_tolerance,
				  Integer, md_degree,
				  double, md_tolerance,
				  Integer, reliable_update_freq,
				  Integer, bnd_check_freq);
  RatQuoParameters() { 
    bnd_lo = 1e-2;
    bnd_hi = 30;
    action_degree = 10;
    action_tolerance = 1e-10;
    md_degree = 10;
    md_tolerance = 1e-8;
    bnd_check_freq = 20;
    reliable_update_freq = 50;
  }

  void Export(RationalActionParams &into) const{
    into.lo = bnd_lo;
    into.hi = bnd_hi;
    into.action_degree = action_degree;
    into.action_tolerance = action_tolerance;
    into.md_degree = md_degree;
    into.md_tolerance = md_tolerance;
    into.BoundsCheckFreq = bnd_check_freq;
  }
};


struct EvolParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(EvolParameters,
                                  Integer, StartTrajectory,
                                  Integer, Trajectories,
				  Integer, SaveInterval,
                                  bool, MetropolisTest,
				  std::string, StartingType,
				  std::vector<Integer>, GparityDirs,
				  RatQuoParameters, rat_quo_l,
				  RatQuoParameters, rat_quo_s);

  EvolParameters() {
    //For initial thermalization; afterwards user should switch Metropolis on and use StartingType=CheckpointStart
    MetropolisTest    = false;
    StartTrajectory   = 0;
    Trajectories      = 50;
    SaveInterval = 5;
    StartingType      = "ColdStart";
    GparityDirs.resize(3, 1); //1 for G-parity, 0 for periodic
  }
};

bool fileExists(const std::string &fn){
  std::ifstream f(fn);
  return f.good();
}

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  std::string param_file = "params.xml";
  for(int i=1;i<argc;i++){
    if(std::string(argv[i]) == "--param_file"){
      assert(i!=argc-1);
      param_file = argv[i+1];
      break;
    }
  }

  //Read the user parameters
  EvolParameters user_params;
  
  if(fileExists(param_file)){
    std::cout << GridLogMessage << " Reading " << param_file << std::endl;
    Grid::XmlReader rd(param_file);
    read(rd, "Params", user_params);
  }else if(!GlobalSharedMemory::WorldRank){
    std::cout << GridLogMessage << " File " << param_file << " does not exist" << std::endl;
    std::cout << GridLogMessage << " Writing xml template to " << param_file << ".templ" << std::endl;
    Grid::XmlWriter wr(param_file + ".templ");
    write(wr, "Params", user_params);

    std::cout << GridLogMessage << " Done" << std::endl;
    Grid_finalize();
    return 0;
  }

  //Check the parameters
  if(user_params.GparityDirs.size() != Nd-1){
    std::cerr << "Error in input parameters: expect GparityDirs to have size = " << Nd-1 << std::endl;
    exit(1);
  }
  for(int i=0;i<Nd-1;i++)
    if(user_params.GparityDirs[i] != 0 && user_params.GparityDirs[i] != 1){
      std::cerr << "Error in input parameters: expect GparityDirs values to be 0 (periodic) or 1 (G-parity)" << std::endl;
      exit(1);
    }

   // Typedefs to simplify notation
  typedef GparityDomainWallFermionD FermionActionD;
  typedef typename FermionActionD::Impl_t FermionImplPolicyD;
  typedef typename FermionActionD::FermionField FermionFieldD;

  typedef GparityDomainWallFermionF FermionActionF;
  typedef typename FermionActionF::Impl_t FermionImplPolicyF;
  typedef typename FermionActionF::FermionField FermionFieldF;

  typedef GeneralEvenOddRatioRationalMixedPrecPseudoFermionAction<FermionImplPolicyD,FermionImplPolicyF> MixedPrecRHMC;

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IntegratorParameters MD;
  typedef ConjugateHMCRunnerD<MinimumNorm2> HMCWrapper; //NB: This is the "Omelyan integrator"
  MD.name    = std::string("MinimumNorm2");
  MD.MDsteps = 5; //5 steps of 0.2 for GP* ensembles
  MD.trajL   = 1.0;

  HMCparameters HMCparams;
  HMCparams.StartTrajectory  = user_params.StartTrajectory;
  HMCparams.Trajectories     = user_params.Trajectories;
  HMCparams.NoMetropolisUntil= 0;
  HMCparams.StartingType     = user_params.StartingType;
  HMCparams.MetropolisTest = user_params.MetropolisTest;
  HMCparams.MD = MD;
  HMCWrapper TheHMC(HMCparams);

  // Grid from the command line arguments --grid and --mpi
  TheHMC.Resources.AddFourDimGrid("gauge"); // use default simd lanes decomposition

  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix    = "ckpoint_rng";
  CPparams.saveInterval  = user_params.SaveInterval;
  CPparams.format        = "IEEE64BIG";
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  //Note that checkpointing saves the RNG state so that this initialization is required only for the very first configuration
  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  //////////////////////////////////////////////

  const int Ls      = 16;
  Real beta         = 2.13;
  Real light_mass   = 0.01;
  Real strange_mass = 0.032;
  Real pv_mass      = 1.0;
  RealD M5  = 1.8;

  //Setup the Grids
  auto GridPtrD   = TheHMC.Resources.GetCartesian();
  auto GridRBPtrD = TheHMC.Resources.GetRBCartesian();
  auto FGridD     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtrD);
  auto FrbGridD   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtrD);

  GridCartesian* GridPtrF = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexF::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* GridRBPtrF = SpaceTimeGrid::makeFourDimRedBlackGrid(GridPtrF);
  auto FGridF     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtrF);
  auto FrbGridF   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtrF);

  ConjugateIwasakiGaugeActionD GaugeAction(beta);

  // temporarily need a gauge field
  LatticeGaugeFieldD Ud(GridPtrD);
  LatticeGaugeFieldF Uf(GridPtrF);

  
  //Setup the BCs
  FermionActionD::ImplParams Params;
  for(int i=0;i<Nd-1;i++) Params.twists = user_params.GparityDirs[i]; //G-parity directions
  Params.twists[Nd-1] = 1; //APBC in time direction

  std::vector<int> dirs4(Nd);
  for(int i=0;i<Nd-1;i++) dirs4[i] = user_params.GparityDirs[i];
  dirs4[Nd-1] = 0; //periodic gauge BC in time

  ConjugateGimplD::setDirections(dirs4); //gauge BC


  ////////////////////////////////////
  // Collect actions
  ////////////////////////////////////
  ActionLevel<HMCWrapper::Field> Level1(1); //light quark 
  ActionLevel<HMCWrapper::Field> Level2(1); //strange quark
  ActionLevel<HMCWrapper::Field> Level3(8); //gauge (8 increments per step)

  ////////////////////////////////////
  // Strange action
  ////////////////////////////////////

  //Use same parameters as used for 16GPX ensembles
  RationalActionParams rat_act_params_s;

  rat_act_params_s.inv_pow  = 4; // (M^dag M)^{1/4}
  rat_act_params_s.precision= 60;
  rat_act_params_s.MaxIter  = 10000;
  user_params.rat_quo_s.Export(rat_act_params_s);

  //For the 16GPX ensembles we used Hasenbusch mass splitting:
  //  det[ (M^dag(0.032) M(0.032)) / (M^dag(1.0) M(1.0)) ]^{1/4}    * det[ (M^dag(0.01) M(0.01)) / (M^dag(1.0) M(1.0)) ]^{1/2}
  //=
  // [ det[ (M^dag(0.032) M(0.032)) / (M^dag(1.0) M(1.0)) ]^{1/4}  ]^3    * det[ (M^dag(0.01) M(0.01)) / (M^dag(0.032) M(0.032)) ]^{1/2} 

  //I don't know if it's actually necessary for the action objects to be independent instances...
  int n_hasenbusch_s = 3;
  std::vector<FermionActionD*> Numerators_sD(n_hasenbusch_s);
  std::vector<FermionActionD*> Denominators_sD(n_hasenbusch_s);
  std::vector<FermionActionF*> Numerators_sF(n_hasenbusch_s);
  std::vector<FermionActionF*> Denominators_sF(n_hasenbusch_s);

  std::vector<MixedPrecRHMC*> Quotients_s(n_hasenbusch_s);

  for(int h=0;h<n_hasenbusch_s;h++){
    Numerators_sD[h] = new FermionActionD(Ud,*FGridD,*FrbGridD,*GridPtrD,*GridRBPtrD,strange_mass,M5,Params);
    Denominators_sD[h] = new FermionActionD(Ud,*FGridD,*FrbGridD,*GridPtrD,*GridRBPtrD,pv_mass,  M5,Params);   

    Numerators_sF[h] = new FermionActionF(Uf,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,strange_mass,M5,Params);
    Denominators_sF[h] = new FermionActionF(Uf,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,pv_mass,  M5,Params);   

    //note I define the numerator operator wrt how they appear in the determinant    
    Quotients_s[h] = new MixedPrecRHMC(*Denominators_sD[h], *Numerators_sD[h], *Denominators_sF[h], *Numerators_sF[h], rat_act_params_s, user_params.rat_quo_s.reliable_update_freq); 
    Level2.push_back(Quotients_s[h]);
  }

  /////////////////////////////////////////////////////////////
  // Light action
  /////////////////////////////////////////////////////////////
 
  //We don't Hasenbusch the light quark directly, instead the denominator mass is set equal to the strange mass; cf above
  FermionActionD Numerator_lD(Ud,*FGridD,*FrbGridD,*GridPtrD,*GridRBPtrD, light_mass,M5,Params);
  FermionActionD Denominator_lD(Ud,*FGridD,*FrbGridD,*GridPtrD,*GridRBPtrD, strange_mass,M5,Params);

  FermionActionF Numerator_lF(Uf,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF, light_mass,M5,Params);
  FermionActionF Denominator_lF(Uf,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF, strange_mass,M5,Params);

  RationalActionParams rat_act_params_l;
  rat_act_params_l.inv_pow  = 2; // (M^dag M)^{1/2}
  rat_act_params_l.precision= 60;
  rat_act_params_l.MaxIter  = 10000;
  user_params.rat_quo_l.Export(rat_act_params_l);
 
  MixedPrecRHMC Quotient_l(Denominator_lD, Numerator_lD, Denominator_lF, Numerator_lF, rat_act_params_l, user_params.rat_quo_l.reliable_update_freq);
  Level1.push_back(&Quotient_l);

  /////////////////////////////////////////////////////////////
  // Gauge action
  /////////////////////////////////////////////////////////////
  Level3.push_back(&GaugeAction);
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  TheHMC.TheAction.push_back(Level3);
  std::cout << GridLogMessage << " Action complete "<< std::endl;

  /////////////////////////////////////////////////////////////
  // HMC parameters are serialisable
  
  if(0){
    TheHMC.Resources.AddRNGs();
    ConjugateGimplR::HotConfiguration(TheHMC.Resources.GetParallelRNG(), Ud);
    Quotient_l.refresh(Ud, TheHMC.Resources.GetParallelRNG());    
    LatticeGaugeFieldD out(Ud);
    std::cout << GridLogMessage << " Running the derivative "<< std::endl;
    Quotient_l.deriv(Ud,out);    
    std::cout << GridLogMessage << " Finished running the derivative "<< std::endl;
    Numerator_lD.Report();
    Denominator_lD.Report();
  }


  if(1){
    std::cout << GridLogMessage << " Running the HMC "<< std::endl;
    TheHMC.Run();  // no smearing
  }


  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
} // main

