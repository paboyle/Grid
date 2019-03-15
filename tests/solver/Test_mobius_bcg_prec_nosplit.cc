   /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_mrhs_cg.cc

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

#include <Grid/algorithms/iterative/BlockConjugateGradient.h>
using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  typedef typename DomainWallFermionR::FermionField FermionField; 
  typedef typename DomainWallFermionR::ComplexField ComplexField; 
  typename DomainWallFermionR::ImplParams params; 

  const int Ls=16;

  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();

  std::vector<ComplexD> boundary_phases(Nd,1.);
  boundary_phases[Nd-1]=-1.;
  params.boundary_phases = boundary_phases;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), 
								   GridDefaultSimd(Nd,vComplex::Nsimd()),
								   GridDefaultMpi());
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * rbGrid  = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);
 
  double stp = 1.e-8;
  int nrhs = 4;

  ///////////////////////////////////////////////
  // Set up the problem as a 4d spreadout job
  ///////////////////////////////////////////////
  std::vector<int> seeds({1,2,3,4});

  std::vector<FermionField> src(nrhs,FGrid);
  std::vector<FermionField> src_chk(nrhs,FGrid);
  std::vector<FermionField> result(nrhs,FGrid);
  FermionField tmp(FGrid);
  std::cout << GridLogMessage << "Made the Fermion Fields"<<std::endl;

  for(int s=0;s<nrhs;s++) result[s]=zero;
  GridParallelRNG pRNG5(FGrid);  pRNG5.SeedFixedIntegers(seeds);
  for(int s=0;s<nrhs;s++) {
    random(pRNG5,src[s]);
    std::cout << GridLogMessage << " src ["<<s<<"] "<<norm2(src[s])<<std::endl;
  }

  std::cout << GridLogMessage << "Intialised the Fermion Fields"<<std::endl;

  LatticeGaugeField Umu(UGrid); 

  int conf = 2;
  if(conf==0) { 
    FieldMetaData header;
    std::string file("./lat.in");
    NerscIO::readConfiguration(Umu,header,file);
    std::cout << GridLogMessage << " Config "<<file<<" successfully read" <<std::endl;
  } else if (conf==1){
    GridParallelRNG pRNG(UGrid );  

    pRNG.SeedFixedIntegers(seeds);
    SU3::HotConfiguration(pRNG,Umu);
    std::cout << GridLogMessage << "Intialised the HOT Gauge Field"<<std::endl;
  } else {
    SU3::ColdConfiguration(Umu);
    std::cout << GridLogMessage << "Intialised the COLD Gauge Field"<<std::endl;
  }

  ///////////////////////////////////////////////////////////////
  // Set up N-solvers as trivially parallel
  ///////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << " Building the solvers"<<std::endl;
  RealD mass=0.01;
  RealD M5=1.8;
  DomainWallFermionR Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*rbGrid,mass,M5,params);

  std::cout << GridLogMessage << "****************************************************************** "<<std::endl;
  std::cout << GridLogMessage << " Calling DWF CG "<<std::endl;
  std::cout << GridLogMessage << "****************************************************************** "<<std::endl;

  MdagMLinearOperator<DomainWallFermionR,FermionField> HermOp(Ddwf);
  ConjugateGradient<FermionField> CG((stp),100000);

  for(int rhs=0;rhs<1;rhs++){
    result[rhs] = zero;
    CG(HermOp,src[rhs],result[rhs]);
  }

  for(int rhs=0;rhs<1;rhs++){
    std::cout << " Result["<<rhs<<"] norm = "<<norm2(result[rhs])<<std::endl;
  }

  /////////////////////////////////////////////////////////////
  // Try block CG
  /////////////////////////////////////////////////////////////
  int blockDim = 0;//not used for BlockCGVec
  for(int s=0;s<nrhs;s++){
    result[s]=zero;
  }


  {
    BlockConjugateGradient<FermionField>    BCGV  (BlockCGrQVec,blockDim,stp,100000);
    SchurRedBlackDiagTwoSolve<FermionField> SchurSolver(BCGV);
    SchurSolver(Ddwf,src,result);
  }
  
  for(int rhs=0;rhs<nrhs;rhs++){
    std::cout << " Result["<<rhs<<"] norm = "<<norm2(result[rhs])<<std::endl;
  }

  Grid_finalize();
}
