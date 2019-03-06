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

typedef typename ImprovedStaggeredFermionR::FermionField FermionField;
typename ImprovedStaggeredFermionR::ImplParams params;

RealD AllZero(RealD x) { return 0.; }

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  GridCartesian* UGrid =
    SpaceTimeGrid::makeFourDimGrid(
                                   GridDefaultLatt(),
                                   GridDefaultSimd(Nd, vComplex::Nsimd()),
                                   GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian* FGrid = UGrid;
  GridRedBlackCartesian* FrbGrid = UrbGrid;
  //printf("UGrid=%p UrbGrid=%p FGrid=%p FrbGrid=%p\n", UGrid, UrbGrid, FGrid, FrbGrid);

  std::vector<int> seeds4({1, 2, 3, 4});
  GridParallelRNG RNG4rb(UrbGrid);
  RNG4rb.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid);
    
  FieldMetaData header;
  std::string file("../lat.sample.l4444.ildg");
  IldgReader IR;
  IR.open(file);
  IR.readConfiguration(Umu,header);
  IR.close();

  std::vector<LatticeColourMatrix> U(4, UGrid);
  for (int mu = 0; mu < Nd; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
  }

  RealD mass=0.1;
  RealD c1=9.0/8.0;
  RealD c2=-1.0/24.0;
  RealD u0=1.0;
  ImprovedStaggeredFermionR Ds(Umu,Umu,*UGrid,*UrbGrid,mass,c1,c2,u0,params);

  SchurStagOperator < ImprovedStaggeredFermionR, FermionField > HermOp (Ds);
  Chebyshev<FermionField> Cheby(0.12,15.,41);
  FunctionHermOp<FermionField> OpCheby(Cheby,HermOp);
  PlainHermOp<FermionField> Op     (HermOp);
    
  const int Nstop = 20;
  const int Nk = 20;
  const int Np = 10;
  const int Nm = Nk+Np;
  const int MaxIt= 4;
  RealD resid = 1.0e-6;
  ImplicitlyRestartedLanczos<FermionField> IRL(OpCheby,Op,Nstop,Nk,Nm,resid,MaxIt);
  std::vector<RealD> eval(Nm);
  FermionField v0(UrbGrid);
  gaussian(RNG4rb,v0);
  std::vector<ImprovedStaggeredFermionR::FermionField> evec(Nm,UrbGrid);
  for(int i=0;i<1;i++){
    evec[i].checkerboard = Even;
    std::cout << GridLogMessage << i <<" / "<< Nm << " grid pointer " << evec[i]._grid<< std::endl;
  };
    
  int Nconv;
  IRL.calc(eval,evec,v0,Nconv);

  Grid_finalize();
}
