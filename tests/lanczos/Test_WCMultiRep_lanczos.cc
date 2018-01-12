/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_dwf_lanczos.cc

Copyright (C) 2015

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


//typedef WilsonCloverFermionR FermionOp;
//typedef typename WilsonFermionR::FermionField FermionField;

typedef WilsonImplR FundImplPolicy;
typedef WilsonCloverFermionR FundFermionAction; 
typedef typename FundFermionAction::FermionField FundFermionField;

typedef WilsonTwoIndexAntiSymmetricImplR ASymmImplPolicy; 
typedef WilsonCloverTwoIndexAntiSymmetricFermionR ASymmFermionAction; 
typedef typename ASymmFermionAction::FermionField ASymmFermionField;


RealD AllZero(RealD x) { return 0.; }

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()),
      GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid =
      SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian* FGrid = UGrid;
  GridRedBlackCartesian* FrbGrid = UrbGrid;
  printf("UGrid=%p UrbGrid=%p FGrid=%p FrbGrid=%p\n", UGrid, UrbGrid, FGrid,
         FrbGrid);

  std::vector<int> seeds4({1, 2, 3, 4});
  std::vector<int> seeds5({5, 6, 7, 8});
  GridParallelRNG RNG5(FGrid);
  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG RNG4(UGrid);
  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG RNG5rb(FrbGrid);
  RNG5.SeedFixedIntegers(seeds5);

  GridParallelRNG          pRNG(UGrid); 
  GridSerialRNG            sRNG; 

  FundamentalRepresentation::LatticeField Umu(UGrid);
  
  TwoIndexAntiSymmetricRepresentation HiRep(UGrid);
  TwoIndexAntiSymmetricRepresentation::LatticeField UmuAS(UGrid);

  
  CheckpointerParameters CPparams;
  
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.format = "IEEE64BIG";

//NerscHmcCheckpointer<PeriodicGimplR> Checkpoint(std::string("ckpoint_lat"),
 //                                                 std::string("ckpoint_rng"), 1);
 
NerscHmcCheckpointer<PeriodicGimplR> Checkpoint(CPparams);

  int CNFGSTART=1;
  int CNFGEND=2;
  int CNFGSTEP=1;

  Real Fundmass = -0.1;
  Real Fundcsw  =  1.0;
  Real ASmass   = -0.1;
  Real AScsw    =  1.0;
    
  std::cout << "Fund: mass and csw" << Fundmass << " and " << Fundcsw << std::endl; 
  std::cout << "AS  : mass and csw" << ASmass << " and " << AScsw << std::endl; 

  const int Nstop = 30;
  const int Nk = 40;
  const int Np = 40;
  const int Nm = Nk + Np;
  const int MaxIt = 10000;
  RealD resid = 1.0e-8;

    for (int cnfg=CNFGSTART;cnfg<=CNFGEND;cnfg+=CNFGSTEP){
      Checkpoint.CheckpointRestore(cnfg,Umu, sRNG, pRNG);

  //SU4::HotConfiguration(RNG4, Umu); // temporary, then read.
  
  HiRep.update_representation(Umu);
  UmuAS = HiRep.U;

  FundFermionAction FundFermOp(Umu,*FGrid,*FrbGrid, Fundmass, Fundcsw, Fundcsw);
  MdagMLinearOperator<FundFermionAction,FundFermionField> HermOpFund(FundFermOp); /// <-----
  
  ASymmFermionAction ASFermOp(UmuAS,*FGrid,*FrbGrid, ASmass, AScsw, AScsw);
  MdagMLinearOperator<ASymmFermionAction,ASymmFermionField> HermOpAS(ASFermOp); /// <-----
  
  std::vector<double> Coeffs{0, -1.};
  Polynomial<FundFermionField> FundPolyX(Coeffs);
  //Chebyshev<FundFermionField> FundCheb(0.0, 10., 12);
  
  FunctionHermOp<FundFermionField> FundPolyXOp(FundPolyX,HermOpFund);
  PlainHermOp<FundFermionField>    FundOp     (HermOpFund);

  ImplicitlyRestartedLanczos<FundFermionField> IRL_Fund(FundOp, FundPolyXOp, Nstop, Nk, Nm,
                                               resid, MaxIt);
  
  Polynomial<ASymmFermionField> ASPolyX(Coeffs);
  //Chebyshev<ASymmFermionField> ASCheb(0.0, 10., 12);

  FunctionHermOp<ASymmFermionField> ASPolyXOp(ASPolyX,HermOpAS);
  PlainHermOp<ASymmFermionField>    ASOp     (HermOpAS);

  ImplicitlyRestartedLanczos<ASymmFermionField> IRL_AS(ASOp, ASPolyXOp, Nstop, Nk, Nm,
                                               resid, MaxIt);
                                               
  std::vector<RealD> Fundeval(Nm);
  std::vector<RealD> ASeval(Nm);

  FundFermionField Fundsrc(FGrid);
  ASymmFermionField   ASsrc(FGrid);
  
  gaussian(RNG5, Fundsrc);
  gaussian(RNG5, ASsrc);

  std::vector<FundFermionField> Fundevec(Nm, FGrid);
  std::vector<ASymmFermionField>   ASevec(Nm, FGrid);
  
  for (int i = 0; i < 1; i++) {
    std::cout << i << " / " << Nm << "Fund: grid pointer " << Fundevec[i]._grid
              << std::endl;
  };
  for (int i = 0; i < 1; i++) {
    std::cout << i << " / " << Nm << "AS: grid pointer " << ASevec[i]._grid
              << std::endl;
  };
  
  int FundNconv, ASNconv;
  IRL_Fund.calc(Fundeval, Fundevec, Fundsrc, FundNconv);
  IRL_AS.calc(ASeval, ASevec, ASsrc, ASNconv);

      for (int i=0;i<FundNconv;i++){
      std::cout << "Fund: eval[" << i << "] = " << Fundeval[i] << std::endl;
    }  
    for (int i=0;i<ASNconv;i++){
      std::cout << "2Index: eval[" << i << "] = " << ASeval[i] << std::endl;
    }  
    }

  Grid_finalize();
}
