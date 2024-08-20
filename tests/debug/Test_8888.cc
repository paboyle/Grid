/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_general_coarse_hdcg.cc

    Copyright (C) 2023

Author: Peter Boyle <pboyle@bnl.gov>

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
#include <Grid/algorithms/iterative/ImplicitlyRestartedBlockLanczos.h>
#include <Grid/algorithms/iterative/ImplicitlyRestartedBlockLanczosCoarse.h>
#include <Grid/algorithms/iterative/AdefMrhs.h>

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=8;
  const int nbasis = 40;
  const int cb = 0 ;
  RealD mass=0.01;
  RealD M5=1.8;
  RealD b=1.0;
  RealD c=0.0;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
								   GridDefaultSimd(Nd,vComplex::Nsimd()),
								   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  ///////////////////////// RNGs /////////////////////////////////
  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  std::vector<int> cseeds({5,6,7,8});

  GridParallelRNG          RNG5(FGrid);   RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);   RNG4.SeedFixedIntegers(seeds4);

  ///////////////////////// Configuration /////////////////////////////////
  LatticeGaugeField Umu(UGrid);

  FieldMetaData header;
  std::string file("ckpoint_EODWF_lat.125");
  NerscIO::readConfiguration(Umu,header,file);

  //////////////////////// Fermion action //////////////////////////////////
  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);

  MdagMLinearOperator<MobiusFermionD, LatticeFermion> HermOp(Ddwf);

  
  std::cout << "**************************************"<<std::endl;
  std::cout << "         Fine Power method            "<<std::endl;
  std::cout << "**************************************"<<std::endl;

  LatticeFermionD pm_src(FGrid);
  pm_src = ComplexD(1.0);
  PowerMethod<LatticeFermionD>       fPM;
  fPM(HermOp,pm_src);

  
  std::cout << "**************************************"<<std::endl;
  std::cout << "         Fine Lanczos  (poly, low)    "<<std::endl;
  std::cout << "**************************************"<<std::endl;
  
  int Nk=80;
  int Nm=Nk*3;
  int Nstop=8;
  int Nconv_test_interval=1;
  
  //  Chebyshev<LatticeFermionD>      IRLChebyLo(0.2,64.0,201);  // 1 iter
  Chebyshev<LatticeFermionD>      IRLChebyLo(0.0,55.0,101);  // 1 iter
  FunctionHermOp<LatticeFermionD>    PolyOp(IRLChebyLo,HermOp);
  PlainHermOp<LatticeFermionD>          Op(HermOp);

  ImplicitlyRestartedLanczos IRL(PolyOp,
				 Op,
				 Nk, // sought vecs
				 Nk, // sought vecs
				 Nm, // spare vecs
				 1.0e-8,
				 10 // Max iterations
				 );

  int Nconv;
  std::vector<RealD>            eval(Nm);
  std::vector<LatticeFermionD>     evec(Nm,FGrid);
  LatticeFermionD     irl_src(FGrid);

  IRL.calc(eval,evec,irl_src,Nconv);

  Grid_finalize();
  return 0;
}
