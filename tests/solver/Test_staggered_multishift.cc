    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_wilson_cg_unprec.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
#include <Grid/algorithms/iterative/BlockConjugateGradient.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

template<class d>
struct scal {
  d internal;
};

  Gamma::Algebra Gmu [] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT
  };

int main (int argc, char ** argv)
{
  typedef typename ImprovedStaggeredFermionR::FermionField FermionField; 
  typename ImprovedStaggeredFermionR::ImplParams params; 

  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();

  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(&Grid);

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&Grid);  pRNG.SeedFixedIntegers(seeds);


  LatticeGaugeField Umu(&Grid); SU3::HotConfiguration(pRNG,Umu);

  double volume=1;
  for(int mu=0;mu<Nd;mu++){
    volume=volume*latt_size[mu];
  }  

  ////////////////////////////////////////
  // sqrt 
  ////////////////////////////////////////
  double     lo=0.001;
  double     hi=1.0;
  int precision=64;
  int    degree=10;
  AlgRemez remez(lo,hi,precision);
  remez.generateApprox(degree,1,2);
  MultiShiftFunction Sqrt(remez,1.0e-6,false);
  std::cout<<GridLogMessage << "Generating degree "<<degree<<" for x^(1/2)"<<std::endl;


  ////////////////////////////////////////////
  // Setup staggered
  ////////////////////////////////////////////
  RealD mass=0.003;
  RealD c1=9.0/8.0;
  RealD c2=-1.0/24.0;
  RealD u0=1.0;

  ImprovedStaggeredFermionR Ds(Umu,Umu,Grid,RBGrid,mass,c1,c2,u0);
  SchurStaggeredOperator<ImprovedStaggeredFermionR,FermionField> HermOpEO(Ds);

  FermionField src(&Grid); random(pRNG,src);
  FermionField src_o(&RBGrid); 
  pickCheckerboard(Odd,src_o,src);


  /////////////////////////////////
  //Multishift CG
  /////////////////////////////////
  std::vector<FermionField> result(degree,&RBGrid);
  ConjugateGradientMultiShift<FermionField> MSCG(10000,Sqrt);

  double deodoe_flops=(1205+15*degree)*volume; // == 66*16 +  == 1146

  double t1=usecond();
  MSCG(HermOpEO,src_o,result);
  double t2=usecond();
  double ncall=MSCG.IterationsToComplete;
  double flops = deodoe_flops * ncall;
  std::cout<<GridLogMessage << "usec    =   "<< (t2-t1)<<std::endl;
  std::cout<<GridLogMessage << "flops   =   "<< flops<<std::endl;
  std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t2-t1)<<std::endl;
  //  HermOpEO.Report();

  Grid_finalize();
}
