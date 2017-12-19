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

  const int Ls=4;

  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  std::vector<int> mpi_split (mpi_layout.size(),1);

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * rbGrid  = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  /////////////////////////////////////////////
  // Split into 1^4 mpi communicators
  /////////////////////////////////////////////

  for(int i=0;i<argc;i++){
    if(std::string(argv[i]) == "--split"){
      for(int k=0;k<mpi_layout.size();k++){
	std::stringstream ss; 
	ss << argv[i+1+k]; 
	ss >> mpi_split[k];
      }
      break;
    }
  }

  int nrhs = 1;
  for(int i=0;i<mpi_layout.size();i++) nrhs *= (mpi_layout[i]/mpi_split[i]);

  GridCartesian         * SGrid = new GridCartesian(GridDefaultLatt(),
						    GridDefaultSimd(Nd,vComplex::Nsimd()),
						    mpi_split,
						    *UGrid); 

  GridCartesian         * SFGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,SGrid);
  GridRedBlackCartesian * SrbGrid  = SpaceTimeGrid::makeFourDimRedBlackGrid(SGrid);
  GridRedBlackCartesian * SFrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,SGrid);

  ///////////////////////////////////////////////
  // Set up the problem as a 4d spreadout job
  ///////////////////////////////////////////////
  std::vector<int> seeds({1,2,3,4});

  GridParallelRNG pRNG(UGrid );  pRNG.SeedFixedIntegers(seeds);
  GridParallelRNG pRNG5(FGrid);  pRNG5.SeedFixedIntegers(seeds);
  std::vector<FermionField>    src(nrhs,FGrid);
  std::vector<FermionField> src_chk(nrhs,FGrid);
  std::vector<FermionField> result(nrhs,FGrid);
  FermionField tmp(FGrid);

  for(int s=0;s<nrhs;s++) random(pRNG5,src[s]);
  for(int s=0;s<nrhs;s++) result[s]=zero;

  LatticeGaugeField Umu(UGrid); SU3::HotConfiguration(pRNG,Umu);

  /////////////////
  // MPI only sends
  /////////////////
  int me = UGrid->ThisRank();

  LatticeGaugeField s_Umu(SGrid);
  FermionField s_src(SFGrid);
  FermionField s_tmp(SFGrid);
  FermionField s_res(SFGrid);

  ///////////////////////////////////////////////////////////////
  // split the source out using MPI instead of I/O
  ///////////////////////////////////////////////////////////////
  Grid_split  (Umu,s_Umu);
  Grid_split  (src,s_src);

  ///////////////////////////////////////////////////////////////
  // Set up N-solvers as trivially parallel
  ///////////////////////////////////////////////////////////////
  RealD mass=0.01;
  RealD M5=1.8;
  DomainWallFermionR Dchk(Umu,*FGrid,*FrbGrid,*UGrid,*rbGrid,mass,M5);
  DomainWallFermionR Ddwf(s_Umu,*SFGrid,*SFrbGrid,*SGrid,*SrbGrid,mass,M5);

  std::cout << GridLogMessage << "****************************************************************** "<<std::endl;
  std::cout << GridLogMessage << " Calling DWF CG "<<std::endl;
  std::cout << GridLogMessage << "****************************************************************** "<<std::endl;

  MdagMLinearOperator<DomainWallFermionR,FermionField> HermOp(Ddwf);
  MdagMLinearOperator<DomainWallFermionR,FermionField> HermOpCk(Dchk);
  ConjugateGradient<FermionField> CG((1.0e-8/(me+1)),10000);
  s_res = zero;
  CG(HermOp,s_src,s_res);

  /////////////////////////////////////////////////////////////
  // Report how long they all took
  /////////////////////////////////////////////////////////////
  std::vector<uint32_t> iterations(nrhs,0);
  iterations[me] = CG.IterationsToComplete;

  for(int n=0;n<nrhs;n++){
    UGrid->GlobalSum(iterations[n]);
    std::cout << GridLogMessage<<" Rank "<<n<<" "<< iterations[n]<<" CG iterations"<<std::endl;
  }

  /////////////////////////////////////////////////////////////
  // Gather and residual check on the results
  /////////////////////////////////////////////////////////////
  std::cout << GridLogMessage<< "Unsplitting the result"<<std::endl;
  Grid_unsplit(result,s_res);

  std::cout << GridLogMessage<< "Checking the residuals"<<std::endl;
  for(int n=0;n<nrhs;n++){
    HermOpCk.HermOp(result[n],tmp); tmp = tmp - src[n];
    std::cout << GridLogMessage<<" resid["<<n<<"]  "<< norm2(tmp)<<std::endl;
  }

  Grid_finalize();
}
