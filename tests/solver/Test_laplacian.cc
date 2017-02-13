    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_wilson_cg_prec.cc

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
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);


  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(latt_size,simd_layout,mpi_layout);

  std::vector<int> seeds({1,2,3,4,5});
  GridParallelRNG          pRNG(&Grid);  pRNG.SeedFixedIntegers(seeds);

  LatticeGaugeField Umu(&Grid); 
  SU<Nc>::HotConfiguration(pRNG,Umu);

  typedef SU<Nc>::LatticeAlgebraVector AVector;
  // Source and result in the algebra
  AVector src(&Grid); random(pRNG,src);
  AVector result(&Grid); result=zero;

  LaplacianAdjointField<PeriodicGimplR> Laplacian(&Grid);
  Laplacian.ImportGauge(Umu);

  HermitianLinearOperator<LaplacianAdjointField<PeriodicGimplR>,AVector> HermOp(Laplacian);
  ConjugateGradient<AVector> CG(1.0e-8,10000);
  CG(HermOp,src,result);
  
  Grid_finalize();
}
