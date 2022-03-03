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

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();

   // Typedefs to simplify notation
  typedef WilsonImplR FermionImplPolicy;
  typedef MobiusFermionR FermionAction;
  typedef typename FermionAction::FermionField FermionField;

  typedef Grid::XmlReader       Serialiser;

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IntegratorParameters MD;
  //  typedef GenericHMCRunner<LeapFrog> HMCWrapper;
  //  MD.name    = std::string("Leap Frog");
  //  typedef GenericHMCRunner<ForceGradient> HMCWrapper;
  //  MD.name    = std::string("Force Gradient");
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;
  MD.name    = std::string("MinimumNorm2");
  MD.MDsteps =  4;
  MD.trajL   = 1.0;

  HMCparameters HMCparams;
  HMCparams.StartTrajectory  = 8;
  HMCparams.Trajectories     = 200;
  HMCparams.NoMetropolisUntil=  0;
  // "[HotStart, ColdStart, TepidStart, CheckpointStart]\n";
  //  HMCparams.StartingType     =std::string("ColdStart");
  HMCparams.StartingType     =std::string("CheckpointStart");
  HMCparams.MD = MD;
  HMCWrapper TheHMC(HMCparams);

  // Grid from the command line arguments --grid and --mpi
  TheHMC.Resources.AddFourDimGrid("gauge"); // use default simd lanes decomposition

  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_EODWF_lat";
  CPparams.rng_prefix    = "ckpoint_EODWF_rng";
  CPparams.saveInterval  = 1;
  CPparams.format        = "IEEE64BIG";
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  // here there is too much indirection
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  //////////////////////////////////////////////

  const int Ls      = 16;
  Real beta         = 2.13;
  Real light_mass   = 0.01;
  Real strange_mass = 0.04;
  Real pv_mass      = 1.0;
  RealD M5  = 1.8;
  RealD b   = 1.0;
  RealD c   = 0.0;

  // FIXME:
  // Same in MC and MD
  // Need to mix precision too
  OneFlavourRationalParams OFRp;
  OFRp.lo       = 4.0e-3;
  OFRp.hi       = 30.0;
  OFRp.MaxIter  = 10000;
  OFRp.tolerance= 1.0e-10;
  OFRp.degree   = 16;
  OFRp.precision= 50;

  std::vector<Real> hasenbusch({ 0.01, 0.04, 0.2 , pv_mass });
  std::vector<bool> dirichlet ({ true, true, true });

  auto GridPtr   = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();

  ////////////////////////////////////////////////////////////////
  // Domain decomposed
  ////////////////////////////////////////////////////////////////
  Coordinate latt4  = GridPtr->GlobalDimensions();
  Coordinate mpi    = GridPtr->ProcessorGrid();
  Coordinate shm;

  GlobalSharedMemory::GetShmDims(mpi,shm);
  
  Coordinate CommDim(Nd);
  for(int d=0;d<Nd;d++) CommDim[d]= (mpi[d]/shm[d])>1 ? 1 : 0;

  Coordinate Dirichlet(Nd+1,0);
  Dirichlet[1] = CommDim[0]*latt4[0]/mpi[0] * shm[0];
  Dirichlet[2] = CommDim[1]*latt4[1]/mpi[1] * shm[1];
  Dirichlet[3] = CommDim[2]*latt4[2]/mpi[2] * shm[2];
  Dirichlet[4] = CommDim[3]*latt4[3]/mpi[3] * shm[3];

  Coordinate Block4(Nd);
  Block4[0] = Dirichlet[1];
  Block4[1] = Dirichlet[2];
  Block4[2] = Dirichlet[3];
  Block4[3] = Dirichlet[4];
  TheHMC.Resources.SetMomentumFilter(new DDHMCFilter<WilsonImplR::Field>(Block4));

  //////////////////////////
  // Fermion Grid
  //////////////////////////
  auto FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtr);
  auto FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtr);

  IwasakiGaugeActionR GaugeAction(beta);

  // temporarily need a gauge field
  LatticeGaugeField U(GridPtr);

  // These lines are unecessary if BC are all periodic
  std::vector<Complex> boundary = {1,1,1,-1};
  FermionAction::ImplParams Params(boundary);

  double StoppingCondition = 1e-10;
  double MaxCGIterations = 30000;
  ConjugateGradient<FermionField>  CG(StoppingCondition,MaxCGIterations);

  ////////////////////////////////////
  // Collect actions
  ////////////////////////////////////
  ActionLevel<HMCWrapper::Field> Level1(1);
  ActionLevel<HMCWrapper::Field> Level2(2);
  ActionLevel<HMCWrapper::Field> Level3(8);

  ////////////////////////////////////
  // Strange action
  ////////////////////////////////////
  FermionAction StrangeOp (U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,strange_mass,M5,b,c, Params);
  FermionAction StrangePauliVillarsOp(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,pv_mass,  M5,b,c, Params);

  OneFlavourEvenOddRatioRationalPseudoFermionAction<FermionImplPolicy> StrangePseudoFermion(StrangePauliVillarsOp,StrangeOp,OFRp);
  //  Level1.push_back(&StrangePseudoFermion);

  ////////////////////////////////////
  // up down action
  ////////////////////////////////////
  std::vector<Real> light_den;
  std::vector<Real> light_num;
  std::vector<int> dirichlet_den;
  std::vector<int> dirichlet_num;

  int n_hasenbusch = hasenbusch.size();
  light_den.push_back(light_mass);
  dirichlet_den.push_back(0);
  for(int h=0;h<n_hasenbusch;h++){
    light_den.push_back(hasenbusch[h]);
    light_num.push_back(hasenbusch[h]);
    dirichlet_num.push_back(1);
    dirichlet_den.push_back(1);
  }
  light_num.push_back(pv_mass);
  dirichlet_num.push_back(0);

  std::vector<FermionAction *> Numerators;
  std::vector<FermionAction *> Denominators;
  std::vector<TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy> *> Quotients;
  
  for(int h=0;h<n_hasenbusch+1;h++){
    std::cout << GridLogMessage << " 2f quotient Action  "<< light_num[h]<< " (" << dirichlet_num[h]
	      <<") / " << light_den[h]<< " (" << dirichlet_den[h]<<")"<< std::endl;
    Numerators.push_back  (new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_num[h],M5,b,c, Params));
    Denominators.push_back(new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_den[h],M5,b,c, Params));
    Quotients.push_back   (new TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy>(*Numerators[h],*Denominators[h],CG,CG));
    if ( dirichlet_den[h]==1) Denominators[h]->DirichletBlock(Dirichlet);
    if ( dirichlet_num[h]==1) Numerators[h]->DirichletBlock(Dirichlet);
  }

  int nquo=Quotients.size();
  Level1.push_back(Quotients[0]);
  Level1.push_back(Quotients[nquo-1]);
  for(int h=1;h<nquo-1;h++){
    Level2.push_back(Quotients[h]);
  }

  /////////////////////////////////////////////////////////////
  // Gauge action
  /////////////////////////////////////////////////////////////
  Level3.push_back(&GaugeAction);
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  TheHMC.TheAction.push_back(Level3);
  std::cout << GridLogMessage << " Action complete "<< std::endl;

  /////////////////////////////////////////////////////////////

  std::cout << GridLogMessage << " Running the HMC "<< std::endl;
  TheHMC.Run();  // no smearing

  Grid_finalize();
} // main



