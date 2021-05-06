    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_rhmc_EOWilsonRatio_genericVsOneFlavor.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Christopher Kelly <ckelly@bnl.gov>

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

//This test ensures that the OneFlavourEvenOddRatioRationalPseudoFermionAction and GeneralEvenOddRatioRationalPseudoFermionAction action (with parameters set appropriately0
//give the same results

int main(int argc, char **argv) {
  using namespace Grid;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  typedef WilsonImplR FermionImplPolicy;
  typedef WilsonFermionR FermionAction;
  typedef typename FermionAction::FermionField FermionField;

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  HMCWrapper TheHMC;
  TheHMC.Resources.AddFourDimGrid("gauge");


  // // Checkpointer definition
  // CheckpointerParameters CPparams;  
  // CPparams.config_prefix = "ckpoint_lat";
  // CPparams.rng_prefix = "ckpoint_rng";
  // CPparams.saveInterval = 5;
  // CPparams.format = "IEEE64BIG";
  
  // TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();

  // temporarily need a gauge field
  LatticeGaugeField U(GridPtr);

  Real mass = -0.77;
  Real pv   = 0.0;

  FermionAction DenOp(U, *GridPtr, *GridRBPtr, mass);
  FermionAction NumOp(U, *GridPtr, *GridRBPtr, pv);

  TheHMC.Resources.AddRNGs();
  PeriodicGimplR::HotConfiguration(TheHMC.Resources.GetParallelRNG(), U);

  std::string seed_string = "the_seed";

  //1f action
  OneFlavourRationalParams OneFParams(1.0e-2,64.0,1000,1.0e-6,6); 

  OneFlavourEvenOddRatioRationalPseudoFermionAction<FermionImplPolicy> OneF(NumOp,DenOp,OneFParams);
  TheHMC.Resources.GetParallelRNG().SeedUniqueString(seed_string);
  OneF.refresh(U, TheHMC.Resources.GetParallelRNG());    
  RealD OneFS = OneF.S(U);
  LatticeGaugeField OneFderiv(U);
  OneF.deriv(U,OneFderiv);    
  
  //general action
  RationalActionParams GenParams;
  GenParams.inv_pow = 2;
  GenParams.lo = OneFParams.lo;
  GenParams.hi = OneFParams.hi;
  GenParams.MaxIter = OneFParams.MaxIter;
  GenParams.action_tolerance = GenParams.md_tolerance = OneFParams.tolerance;
  GenParams.action_degree = GenParams.md_degree = OneFParams.degree;
  GenParams.precision = OneFParams.precision;
  GenParams.BoundsCheckFreq = OneFParams.BoundsCheckFreq;

  GeneralEvenOddRatioRationalPseudoFermionAction<FermionImplPolicy> Gen(NumOp,DenOp,GenParams);
  TheHMC.Resources.GetParallelRNG().SeedUniqueString(seed_string);
  Gen.refresh(U, TheHMC.Resources.GetParallelRNG());    
  RealD GenS = Gen.S(U);
  LatticeGaugeField Genderiv(U);
  Gen.deriv(U,Genderiv);   


  //Compare
  std::cout << "Action : " << OneFS << " " << GenS << " reldiff " << (OneFS - GenS)/OneFS << std::endl;
  
  LatticeGaugeField diff(U);
  axpy(diff, -1.0, Genderiv, OneFderiv);
  std::cout << "Norm of difference in deriv " << sqrt(norm2(diff)) << std::endl;

  Grid_finalize();
  return 0;
}

