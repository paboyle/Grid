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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

typedef typename GparityDomainWallFermionR::FermionField FermionField;

RealD AllZero(RealD x){ return 0.;}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);
  printf("UGrid=%p UrbGrid=%p FGrid=%p FrbGrid=%p\n",UGrid,UrbGrid,FGrid,FrbGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5rb(FrbGrid);  RNG5.SeedFixedIntegers(seeds5);

  LatticeGaugeField Umu(UGrid); 
  SU3::TepidConfiguration(RNG4, Umu);

  std::vector<LatticeColourMatrix> U(4,UGrid);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }
  
  RealD mass=0.01;
  RealD M5=1.8;
  RealD mob_b=1.5;
//  DomainWallFermionR Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  GparityMobiusFermionD ::ImplParams params;
  std::vector<int> twists({1,1,1,0});
  params.twists = twists;
  GparityMobiusFermionR  Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,mob_b,mob_b-1.,params);

//  MdagMLinearOperator<DomainWallFermionR,LatticeFermion> HermOp(Ddwf);
//  SchurDiagTwoOperator<DomainWallFermionR,LatticeFermion> HermOp(Ddwf);
  SchurDiagTwoOperator<GparityMobiusFermionR,FermionField> HermOp(Ddwf);
//  SchurDiagMooeeOperator<DomainWallFermionR,LatticeFermion> HermOp(Ddwf);

  const int Nstop = 30;
  const int Nk = 40;
  const int Np = 40;
  const int Nm = Nk+Np;
  const int MaxIt= 10000;
  RealD resid = 1.0e-8;

  std::vector<double> Coeffs { 0.,-1.};
  Polynomial<FermionField> PolyX(Coeffs);
  Chebyshev<FermionField> Cheb(0.2,5.,11);
//  ChebyshevLanczos<LatticeFermion> Cheb(9.,1.,0.,20);
//  Cheb.csv(std::cout);
//  exit(-24);
  ImplicitlyRestartedLanczos<FermionField> IRL(HermOp,Cheb,Nstop,Nk,Nm,resid,MaxIt);

  
  std::vector<RealD>          eval(Nm);
  FermionField    src(FrbGrid); gaussian(RNG5rb,src);
  std::vector<FermionField> evec(Nm,FrbGrid);
  for(int i=0;i<1;i++){
    std::cout << i<<" / "<< Nm<< " grid pointer "<<evec[i]._grid<<std::endl;
  };

  int Nconv;
  IRL.calc(eval,evec,
	   src,
	   Nconv);


  Grid_finalize();
}
