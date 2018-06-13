/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/solver/Test_staggered_fgmres_prec.cc

Copyright (C) 2015-2018

Author: Daniel Richtmann <daniel.richtmann@ur.de>

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
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  typedef typename ImprovedStaggeredFermionR::FermionField FermionField;
  typedef typename ImprovedStaggeredFermionR::ComplexField ComplexField;
  typename ImprovedStaggeredFermionR::ImplParams params;

  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(&Grid);

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&Grid);  pRNG.SeedFixedIntegers(seeds);

  FermionField src(&Grid); random(pRNG,src);
  RealD nrm = norm2(src);
  FermionField result(&Grid); result=zero;
  LatticeGaugeField Umu(&Grid); SU3::HotConfiguration(pRNG,Umu);

  double volume=1;
  for(int mu=0;mu<Nd;mu++){
    volume=volume*latt_size[mu];
  }

  RealD mass=0.5;
  ImprovedStaggeredFermionR Ds(Umu,Umu,Grid,RBGrid,mass);

  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> HermOp(Ds);

  TrivialPrecon<FermionField> simple;

  FlexibleGeneralisedMinimalResidual<FermionField> FGMRES(1.0e-8, 10000, simple, 25);
  FGMRES(HermOp,src,result);

  Grid_finalize();
}
