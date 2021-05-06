    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_rhmc_EOWilsonRatio_doubleVsMixedPrec.cc

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@bnl.gov>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid/Grid.h>

//This test ensures the mixed precision RHMC gives the same result as the regular double precision
int main(int argc, char **argv) {
  using namespace Grid;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm

  typedef WilsonImplD FermionImplPolicyD;
  typedef WilsonFermionD FermionActionD;
  typedef typename FermionActionD::FermionField FermionFieldD;

  typedef WilsonImplF FermionImplPolicyF;
  typedef WilsonFermionF FermionActionF;
  typedef typename FermionActionF::FermionField FermionFieldF;

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  HMCWrapper TheHMC;
  TheHMC.Resources.AddFourDimGrid("gauge");

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  auto GridPtrD = TheHMC.Resources.GetCartesian();
  auto GridRBPtrD = TheHMC.Resources.GetRBCartesian();

  GridCartesian* GridPtrF = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexF::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* GridRBPtrF = SpaceTimeGrid::makeFourDimRedBlackGrid(GridPtrF);

  // temporarily need a gauge field
  LatticeGaugeFieldD Ud(GridPtrD);
  LatticeGaugeFieldF Uf(GridPtrF);

  Real mass = -0.77;
  Real pv   = 0.0;

  FermionActionD DenOpD(Ud, *GridPtrD, *GridRBPtrD, mass);
  FermionActionD NumOpD(Ud, *GridPtrD, *GridRBPtrD, pv);

  FermionActionF DenOpF(Uf, *GridPtrF, *GridRBPtrF, mass);
  FermionActionF NumOpF(Uf, *GridPtrF, *GridRBPtrF, pv);

  TheHMC.Resources.AddRNGs();
  PeriodicGimplR::HotConfiguration(TheHMC.Resources.GetParallelRNG(), Ud);

  std::string seed_string = "the_seed";

  //Setup the pseudofermion actions
  RationalActionParams GenParams;
  GenParams.inv_pow = 2;
  GenParams.lo = 1e-2;
  GenParams.hi = 64.0;
  GenParams.MaxIter = 1000;
  GenParams.action_tolerance = GenParams.md_tolerance = 1e-6;
  GenParams.action_degree = GenParams.md_degree = 6;
  GenParams.precision = 64;
  GenParams.BoundsCheckFreq = 20;

  GeneralEvenOddRatioRationalPseudoFermionAction<FermionImplPolicyD> GenD(NumOpD,DenOpD,GenParams);
  GeneralEvenOddRatioRationalMixedPrecPseudoFermionAction<FermionImplPolicyD, FermionImplPolicyF> GenFD(NumOpD, DenOpD, 
													NumOpF, DenOpF, 
													GenParams, 50);
  TheHMC.Resources.GetParallelRNG().SeedUniqueString(seed_string);
  GenD.refresh(Ud, TheHMC.Resources.GetSerialRNG(), TheHMC.Resources.GetParallelRNG());    
  RealD Sd = GenD.S(Ud);
  LatticeGaugeField derivD(Ud);
  GenD.deriv(Ud,derivD);   

  TheHMC.Resources.GetParallelRNG().SeedUniqueString(seed_string);
  GenFD.refresh(Ud, TheHMC.Resources.GetSerialRNG(), TheHMC.Resources.GetParallelRNG());    
  RealD Sfd = GenFD.S(Ud);
  LatticeGaugeField derivFD(Ud);
  GenFD.deriv(Ud,derivFD);   

  //Compare
  std::cout << "Action : " << Sd << " " << Sfd << " reldiff " << (Sd - Sfd)/Sd << std::endl;
  
  LatticeGaugeField diff(Ud);
  axpy(diff, -1.0, derivD, derivFD);
  std::cout << "Norm of difference in deriv " << sqrt(norm2(diff)) << std::endl;

  Grid_finalize();
  return 0;
}

