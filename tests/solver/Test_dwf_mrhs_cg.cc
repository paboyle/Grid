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

  const int Ls=8;

  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  std::vector<int> mpi_split (mpi_layout.size(),1);

  std::cout << "UGrid (world root)"<<std::endl;
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());

  std::cout << "FGrid (child of UGrid)"<<std::endl;
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);

  int nrhs = UGrid->RankCount() ;

  /////////////////////////////////////////////
  // Split into 1^4 mpi communicators
  /////////////////////////////////////////////
  std::cout << "SGrid (world root)"<<std::endl;
  GridCartesian         * SGrid = new GridCartesian(GridDefaultLatt(),
						    GridDefaultSimd(Nd,vComplex::Nsimd()),
						    mpi_split,
						    *UGrid); 

  GridCartesian         * SFGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,SGrid);
  std::cout << "SFGrid"<<std::endl;
  GridRedBlackCartesian * SrbGrid  = SpaceTimeGrid::makeFourDimRedBlackGrid(SGrid);
  std::cout << "SrbGrid"<<std::endl;
  GridRedBlackCartesian * SFrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,SGrid);
  std::cout << "SFrbGrid"<<std::endl;

  ///////////////////////////////////////////////
  // Set up the problem as a 4d spreadout job
  ///////////////////////////////////////////////
  std::vector<int> seeds({1,2,3,4});

  GridParallelRNG pRNG(UGrid );  pRNG.SeedFixedIntegers(seeds);
  GridParallelRNG pRNG5(FGrid);  pRNG5.SeedFixedIntegers(seeds);
  std::vector<FermionField>    src(nrhs,FGrid);
  std::vector<FermionField> result(nrhs,FGrid);

  for(int s=0;s<nrhs;s++) random(pRNG5,src[s]);
  for(int s=0;s<nrhs;s++) result[s] = zero;

  LatticeGaugeField Umu(UGrid); SU3::HotConfiguration(pRNG,Umu);

  ///////////////////////////////////////////////////////////////
  // Bounce these fields to disk
  ///////////////////////////////////////////////////////////////

  std::cout << GridLogMessage << "****************************************************************** "<<std::endl;
  std::cout << GridLogMessage << " Writing out in parallel view "<<std::endl;
  std::cout << GridLogMessage << "****************************************************************** "<<std::endl;
  emptyUserRecord record;
  std::string file("./scratch.scidac");
  {
    ScidacWriter _ScidacWriter;
    _ScidacWriter.open(file);
    _ScidacWriter.writeScidacFieldRecord(Umu,record);
    for(int n=0;n<nrhs;n++){
      _ScidacWriter.writeScidacFieldRecord(src[n],record);
      
    }
    _ScidacWriter.close();
  }



  //////////////////////////////////////////
  // Read back into single rank fields
  //////////////////////////////////////////

  std::cout << GridLogMessage << "****************************************************************** "<<std::endl;
  std::cout << GridLogMessage << " Reading back in the single process view "<<std::endl;
  std::cout << GridLogMessage << "****************************************************************** "<<std::endl;

  int me = UGrid->ThisRank();
  LatticeGaugeField s_Umu(SGrid);
  FermionField s_src(SFGrid);
  FermionField s_res(SFGrid);
  {
    ScidacReader  _ScidacReader;
    _ScidacReader.open(file);
    std::cout << GridLogMessage << " Opened file "<<std::endl;
    _ScidacReader.readScidacFieldRecord(s_Umu,record);
    std::cout << GridLogMessage << " Read gauge field "<<std::endl;
    for(int n=0;n<nrhs;n++){
      if ( n==me ) { 
	std::cout << GridLogMessage << " Read record  "<<n<<std::endl;
	_ScidacReader.readScidacFieldRecord(s_src,record);
      } else { 
	std::cout << GridLogMessage << " Skip record  "<<n<<std::endl;
	_ScidacReader.skipScidacFieldRecord();
      }      
    }
    _ScidacReader.close();
  }

  ///////////////////////////////////////////////////////////////
  // Set up N-solvers as trivially parallel
  ///////////////////////////////////////////////////////////////

  RealD mass=0.01;
  RealD M5=1.8;
  DomainWallFermionR Ddwf(s_Umu,*SFGrid,*SFrbGrid,*SGrid,*SrbGrid,mass,M5);

  std::cout << GridLogMessage << "****************************************************************** "<<std::endl;
  std::cout << GridLogMessage << " Calling DWF CG "<<std::endl;
  std::cout << GridLogMessage << "****************************************************************** "<<std::endl;

  MdagMLinearOperator<DomainWallFermionR,FermionField> HermOp(Ddwf);
  ConjugateGradient<FermionField> CG((1.0e-8/(me+1)),10000);
  s_res = zero;
  CG(HermOp,s_src,s_res);

  ///////////////////////////////////////
  // Share the information
  ///////////////////////////////////////
  std::vector<uint32_t> iterations(nrhs,0);
  iterations[me] = CG.IterationsToComplete;

  for(int n=0;n<nrhs;n++){
    UGrid->GlobalSum(iterations[n]);
  }

  /////////////////////////////////////////////////////////////
  // Report how long they all took
  /////////////////////////////////////////////////////////////
  for(int r=0;r<nrhs;r++){
    std::cout << " Rank "<<r<<" "<< iterations[r]<<" CG iterations"<<std::endl;
  }
  Grid_finalize();
}
