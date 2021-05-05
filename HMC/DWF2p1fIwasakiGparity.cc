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
				  Integer, Steps,
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
    Steps = 5;
  }
};

bool fileExists(const std::string &fn){
  std::ifstream f(fn);
  return f.good();
}




struct LanczosParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParameters,
				  double, alpha,
				  double, beta,
				  double, mu,
				  int, ord,
				  int, n_stop,
				  int, n_want,
				  int, n_use,
				  double, tolerance);

  LanczosParameters() {
    alpha = 35;
    beta = 5;
    mu = 0;
    ord = 100;
    n_stop = 10;
    n_want = 10;
    n_use = 15;
    tolerance = 1e-6;
  }
};



template<typename FermionActionD, typename FermionFieldD>
void computeEigenvalues(std::string param_file,
			GridCartesian* Grid, GridRedBlackCartesian* rbGrid, const LatticeGaugeFieldD &latt,  //expect lattice to have been initialized to something
			FermionActionD &action, GridParallelRNG &rng){
  
  LanczosParameters params;
  if(fileExists(param_file)){
    std::cout << GridLogMessage << " Reading " << param_file << std::endl;
    Grid::XmlReader rd(param_file);
    read(rd, "LanczosParameters", params);
  }else if(!GlobalSharedMemory::WorldRank){
    std::cout << GridLogMessage << " File " << param_file << " does not exist" << std::endl;
    std::cout << GridLogMessage << " Writing xml template to " << param_file << ".templ" << std::endl;
    Grid::XmlWriter wr(param_file + ".templ");
    write(wr, "LanczosParameters", params);
  }

  FermionFieldD gauss_o(rbGrid);
  FermionFieldD gauss(Grid);
  gaussian(rng, gauss);
  pickCheckerboard(Odd, gauss_o, gauss);

  action.ImportGauge(latt);

  SchurDiagMooeeOperator<FermionActionD, FermionFieldD> hermop(action);
  PlainHermOp<FermionFieldD> hermop_wrap(hermop);
  //ChebyshevLanczos<FermionFieldD> Cheb(params.alpha, params.beta, params.mu, params.ord);
  assert(params.mu == 0.0);

  Chebyshev<FermionFieldD> Cheb(params.beta*params.beta, params.alpha*params.alpha, params.ord+1);
  FunctionHermOp<FermionFieldD> Cheb_wrap(Cheb, hermop);

  std::cout << "IRL: alpha=" << params.alpha << " beta=" << params.beta << " mu=" << params.mu << " ord=" << params.ord << std::endl;
  ImplicitlyRestartedLanczos<FermionFieldD> IRL(Cheb_wrap, hermop_wrap, params.n_stop, params.n_want, params.n_use, params.tolerance, 10000);

  std::vector<RealD> eval(params.n_use);
  std::vector<FermionFieldD> evec(params.n_use, rbGrid);
  int Nconv;
  IRL.calc(eval, evec, gauss_o, Nconv);

  std::cout << "Eigenvalues:" << std::endl;
  for(int i=0;i<params.n_want;i++){
    std::cout << i << " " << eval[i] << std::endl;
  }
}


//Check the quality of the RHMC approx
template<typename FermionActionD, typename FermionFieldD, typename RHMCtype>
void checkRHMC(GridCartesian* Grid, GridRedBlackCartesian* rbGrid, const LatticeGaugeFieldD &latt,  //expect lattice to have been initialized to something
	       FermionActionD &numOp, FermionActionD &denOp, RHMCtype &rhmc, GridParallelRNG &rng,
	       int inv_pow, const std::string &quark_descr){

  FermionFieldD gauss_o(rbGrid);
  FermionFieldD gauss(Grid);
  gaussian(rng, gauss);
  pickCheckerboard(Odd, gauss_o, gauss);

  numOp.ImportGauge(latt);
  denOp.ImportGauge(latt);

  typedef typename FermionActionD::Impl_t FermionImplPolicyD;
  SchurDifferentiableOperator<FermionImplPolicyD> MdagM(numOp);
  SchurDifferentiableOperator<FermionImplPolicyD> VdagV(denOp);
      
  std::cout << "Starting: Checking quality of RHMC action approx for " << quark_descr << " quark numerator and power -1/" << inv_pow << std::endl;
  InversePowerBoundsCheck(inv_pow, 10000, 1e16, MdagM,gauss_o, rhmc.ApproxNegPowerAction); //use large tolerance to prevent exit on fail; we are trying to tune here!
  std::cout << "Finished: Checking quality of RHMC action approx for " << quark_descr << " quark numerator and power -1/" << inv_pow << std::endl;

  std::cout << "Starting: Checking quality of RHMC action approx for " << quark_descr << " quark numerator and power -1/" << 2*inv_pow << std::endl;
  InversePowerBoundsCheck(2*inv_pow, 10000, 1e16, MdagM,gauss_o, rhmc.ApproxNegHalfPowerAction);
  std::cout << "Finished: Checking quality of RHMC action approx for " << quark_descr << " quark numerator and power -1/" << 2*inv_pow << std::endl;

  std::cout << "Starting: Checking quality of RHMC action approx for " << quark_descr << " quark denominator and power -1/" << inv_pow << std::endl;
  InversePowerBoundsCheck(inv_pow, 10000, 1e16, VdagV,gauss_o, rhmc.ApproxNegPowerAction);
  std::cout << "Finished: Checking quality of RHMC action approx for " << quark_descr << " quark denominator and power -1/" << inv_pow << std::endl;

  std::cout << "Starting: Checking quality of RHMC action approx for " << quark_descr << " quark denominator and power -1/" << 2*inv_pow << std::endl;
  InversePowerBoundsCheck(2*inv_pow, 10000, 1e16, VdagV,gauss_o, rhmc.ApproxNegHalfPowerAction);
  std::cout << "Finished: Checking quality of RHMC action approx for " << quark_descr << " quark denominator and power -1/" << 2*inv_pow << std::endl;

  std::cout << "-------------------------------------------------------------------------------" << std::endl;

  std::cout << "Starting: Checking quality of RHMC MD approx for " << quark_descr << " quark numerator and power -1/" << inv_pow << std::endl;
  InversePowerBoundsCheck(inv_pow, 10000, 1e16, MdagM,gauss_o, rhmc.ApproxNegPowerMD); 
  std::cout << "Finished: Checking quality of RHMC MD approx for " << quark_descr << " quark numerator and power -1/" << inv_pow << std::endl;

  std::cout << "Starting: Checking quality of RHMC MD approx for " << quark_descr << " quark numerator and power -1/" << 2*inv_pow << std::endl;
  InversePowerBoundsCheck(2*inv_pow, 10000, 1e16, MdagM,gauss_o, rhmc.ApproxNegHalfPowerMD);
  std::cout << "Finished: Checking quality of RHMC MD approx for " << quark_descr << " quark numerator and power -1/" << 2*inv_pow << std::endl;

  std::cout << "Starting: Checking quality of RHMC MD approx for " << quark_descr << " quark denominator and power -1/" << inv_pow << std::endl;
  InversePowerBoundsCheck(inv_pow, 10000, 1e16, VdagV,gauss_o, rhmc.ApproxNegPowerMD);
  std::cout << "Finished: Checking quality of RHMC MD approx for " << quark_descr << " quark denominator and power -1/" << inv_pow << std::endl;

  std::cout << "Starting: Checking quality of RHMC MD approx for " << quark_descr << " quark denominator and power -1/" << 2*inv_pow << std::endl;
  InversePowerBoundsCheck(2*inv_pow, 10000, 1e16, VdagV,gauss_o, rhmc.ApproxNegHalfPowerMD);
  std::cout << "Finished: Checking quality of RHMC MD approx for " << quark_descr << " quark denominator and power -1/" << 2*inv_pow << std::endl;
}









int main(int argc, char **argv) {
  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  std::string param_file = "params.xml";
  bool file_load_check = false;
  for(int i=1;i<argc;i++){
    std::string sarg(argv[i]);
    if(sarg == "--param_file"){
      assert(i!=argc-1);
      param_file = argv[i+1];
    }else if(sarg == "--read_check"){ //check the fields load correctly and pass checksum/plaquette repro
      file_load_check = true;
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
  typedef GeneralEvenOddRatioRationalPseudoFermionAction<FermionImplPolicyD> DoublePrecRHMC;

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IntegratorParameters MD;
  typedef ConjugateHMCRunnerD<MinimumNorm2> HMCWrapper; //NB: This is the "Omelyan integrator"
  typedef HMCWrapper::ImplPolicy GaugeImplPolicy;
  MD.name    = std::string("MinimumNorm2");
  MD.MDsteps = user_params.Steps;
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

  typedef PlaquetteMod<GaugeImplPolicy> PlaqObs;
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
  for(int i=0;i<Nd-1;i++) Params.twists[i] = user_params.GparityDirs[i]; //G-parity directions
  Params.twists[Nd-1] = 1; //APBC in time direction

  std::vector<int> dirs4(Nd);
  for(int i=0;i<Nd-1;i++) dirs4[i] = user_params.GparityDirs[i];
  dirs4[Nd-1] = 0; //periodic gauge BC in time

  GaugeImplPolicy::setDirections(dirs4); //gauge BC

  //Run optional gauge field checksum checker and exit
  if(file_load_check){
    TheHMC.initializeGaugeFieldAndRNGs(Ud);
    std::cout << GridLogMessage << " Done" << std::endl;
    Grid_finalize();
    return 0;
  }


  ////////////////////////////////////
  // Collect actions
  ////////////////////////////////////
  ActionLevel<HMCWrapper::Field> Level1(1); //light quark + strange quark
  ActionLevel<HMCWrapper::Field> Level2(8); //gauge (8 increments per step)


  /////////////////////////////////////////////////////////////
  // Light action
  /////////////////////////////////////////////////////////////

  FermionActionD Numerator_lD(Ud,*FGridD,*FrbGridD,*GridPtrD,*GridRBPtrD, light_mass,M5,Params);
  FermionActionD Denominator_lD(Ud,*FGridD,*FrbGridD,*GridPtrD,*GridRBPtrD, pv_mass,M5,Params);

  FermionActionF Numerator_lF(Uf,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF, light_mass,M5,Params);
  FermionActionF Denominator_lF(Uf,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF, pv_mass,M5,Params);

  RationalActionParams rat_act_params_l;
  rat_act_params_l.inv_pow  = 2; // (M^dag M)^{1/2}
  rat_act_params_l.precision= 60;
  rat_act_params_l.MaxIter  = 10000;
  user_params.rat_quo_l.Export(rat_act_params_l);
  std::cout << GridLogMessage << " Light quark bounds check every " << rat_act_params_l.BoundsCheckFreq << " trajectories (avg)" << std::endl;
 
  MixedPrecRHMC Quotient_l(Denominator_lD, Numerator_lD, Denominator_lF, Numerator_lF, rat_act_params_l, user_params.rat_quo_l.reliable_update_freq);
  //DoublePrecRHMC Quotient_l(Denominator_lD, Numerator_lD, rat_act_params_l);
  Level1.push_back(&Quotient_l);


  ////////////////////////////////////
  // Strange action
  ////////////////////////////////////
  FermionActionD Numerator_sD(Ud,*FGridD,*FrbGridD,*GridPtrD,*GridRBPtrD,strange_mass,M5,Params);
  FermionActionD Denominator_sD(Ud,*FGridD,*FrbGridD,*GridPtrD,*GridRBPtrD, pv_mass,M5,Params);

  FermionActionF Numerator_sF(Uf,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF,strange_mass,M5,Params);
  FermionActionF Denominator_sF(Uf,*FGridF,*FrbGridF,*GridPtrF,*GridRBPtrF, pv_mass,M5,Params);

  RationalActionParams rat_act_params_s;
  rat_act_params_s.inv_pow  = 4; // (M^dag M)^{1/4}
  rat_act_params_s.precision= 60;
  rat_act_params_s.MaxIter  = 10000;
  user_params.rat_quo_s.Export(rat_act_params_s);
  std::cout << GridLogMessage << " Heavy quark bounds check every " << rat_act_params_l.BoundsCheckFreq << " trajectories (avg)" << std::endl;

  MixedPrecRHMC Quotient_s(Denominator_sD, Numerator_sD, Denominator_sF, Numerator_sF, rat_act_params_s, user_params.rat_quo_s.reliable_update_freq); 
  //DoublePrecRHMC Quotient_s(Denominator_sD, Numerator_sD, rat_act_params_s); 
  Level1.push_back(&Quotient_s);  


  /////////////////////////////////////////////////////////////
  // Gauge action
  /////////////////////////////////////////////////////////////
  Level2.push_back(&GaugeAction);
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  std::cout << GridLogMessage << " Action complete "<< std::endl;


  //Action tuning
  bool tune_rhmc_l=false, tune_rhmc_s=false, eigenrange_l=false, eigenrange_s=false; 
  std::string lanc_params_l, lanc_params_s;
  for(int i=1;i<argc;i++){
    std::string sarg(argv[i]);
    if(sarg == "--tune_rhmc_l") tune_rhmc_l=true;
    else if(sarg == "--tune_rhmc_s") tune_rhmc_s=true;
    else if(sarg == "--eigenrange_l"){
      assert(i < argc-1);
      eigenrange_l=true;
      lanc_params_l = argv[i+1];
    }
    else if(sarg == "--eigenrange_s"){
      assert(i < argc-1);
      eigenrange_s=true;
      lanc_params_s = argv[i+1];
    }
  }
  if(tune_rhmc_l || tune_rhmc_s || eigenrange_l || eigenrange_s){
    TheHMC.initializeGaugeFieldAndRNGs(Ud);
    if(eigenrange_l) computeEigenvalues<FermionActionD, FermionFieldD>(lanc_params_l, FGridD, FrbGridD, Ud, Numerator_lD, TheHMC.Resources.GetParallelRNG());
    if(eigenrange_s) computeEigenvalues<FermionActionD, FermionFieldD>(lanc_params_s, FGridD, FrbGridD, Ud, Numerator_sD, TheHMC.Resources.GetParallelRNG());
    if(tune_rhmc_l) checkRHMC<FermionActionD, FermionFieldD, decltype(Quotient_l)>(FGridD, FrbGridD, Ud, Numerator_lD, Denominator_lD, Quotient_l, TheHMC.Resources.GetParallelRNG(), 2, "light");
    if(tune_rhmc_s) checkRHMC<FermionActionD, FermionFieldD, decltype(Quotient_s)>(FGridD, FrbGridD, Ud, Numerator_sD, Denominator_sD, Quotient_s, TheHMC.Resources.GetParallelRNG(), 4, "strange");

    std::cout << GridLogMessage << " Done" << std::endl;
    Grid_finalize();
    return 0;
  }


  //Run the HMC
  std::cout << GridLogMessage << " Running the HMC "<< std::endl;
  TheHMC.Run();

  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
} // main

