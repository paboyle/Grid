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
  typedef typename MobiusFermionR::FermionField FermionField; 
  typedef typename MobiusFermionR::ComplexField ComplexField; 
  typename MobiusFermionR::ImplParams params; 

  const int Ls=12;

  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  std::vector<int> mpi_split (mpi_layout.size(),1);
  std::vector<int> split_coor (mpi_layout.size(),1);
  std::vector<int> split_dim (mpi_layout.size(),1);

  std::vector<ComplexD> boundary_phases(Nd,1.);
  boundary_phases[Nd-1]=-1.;
  params.boundary_phases = boundary_phases;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), 
								   GridDefaultSimd(Nd,vComplex::Nsimd()),
								   GridDefaultMpi());
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

 
  double stp = 1.e-8;
  int nrhs = 1;
  int me;
  for(int i=0;i<mpi_layout.size();i++){
//	split_dim[i] = (mpi_layout[i]/mpi_split[i]);
	nrhs *= (mpi_layout[i]/mpi_split[i]);
//	split_coor[i] = FGrid._processor_coor[i]/mpi_split[i];
  }
  std::cout << GridLogMessage << "Creating split grids " <<std::endl;
  GridCartesian         * SGrid = new GridCartesian(GridDefaultLatt(),
						    GridDefaultSimd(Nd,vComplex::Nsimd()),
						    mpi_split,
						    *UGrid,me); 
  std::cout << GridLogMessage <<"Creating split ferm grids " <<std::endl;

  GridCartesian         * SFGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,SGrid);
  std::cout << GridLogMessage <<"Creating split rb grids " <<std::endl;
  GridRedBlackCartesian * SrbGrid  = SpaceTimeGrid::makeFourDimRedBlackGrid(SGrid);
  std::cout << GridLogMessage <<"Creating split ferm rb grids " <<std::endl;
  GridRedBlackCartesian * SFrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,SGrid);
  std::cout << GridLogMessage << "Made the grids"<<std::endl;
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

  if(0) { 
    FieldMetaData header;
    std::string file("./lat.in");
    NerscIO::readConfiguration(Umu,header,file);
    std::cout << GridLogMessage << " "<<file<<" successfully read" <<std::endl;
  } else {
    GridParallelRNG pRNG(UGrid );  
    std::cout << GridLogMessage << "Intialising 4D RNG "<<std::endl;
    pRNG.SeedFixedIntegers(seeds);
    std::cout << GridLogMessage << "Intialised 4D RNG "<<std::endl;
    SU3::HotConfiguration(pRNG,Umu);
    std::cout << GridLogMessage << "Intialised the HOT Gauge Field"<<std::endl;
    std::cout << " Site zero "<< Umu._odata[0]   <<std::endl;
  } 

  /////////////////
  // MPI only sends
  /////////////////
  LatticeGaugeField s_Umu(SGrid);
  FermionField s_src(SFGrid);
  FermionField s_tmp(SFGrid);
  FermionField s_res(SFGrid);

  std::cout << GridLogMessage << "Made the split grid fields"<<std::endl;
  ///////////////////////////////////////////////////////////////
  // split the source out using MPI instead of I/O
  ///////////////////////////////////////////////////////////////
  Grid_split  (Umu,s_Umu);
  Grid_split  (src,s_src);
  std::cout << GridLogMessage << " split rank  " <<me << " s_src "<<norm2(s_src)<<std::endl;

  ///////////////////////////////////////////////////////////////
  // Set up N-solvers as trivially parallel
  ///////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << " Building the solvers"<<std::endl;
  //  RealD mass=0.00107;
  RealD mass=0.1;
  RealD M5=1.8;
  RealD mobius_factor=32./12.;
  RealD mobius_b=0.5*(mobius_factor+1.);
  RealD mobius_c=0.5*(mobius_factor-1.);
  MobiusFermionR Dchk(Umu,*FGrid,*FrbGrid,*UGrid,*rbGrid,mass,M5,mobius_b,mobius_c,params);
  MobiusFermionR Ddwf(s_Umu,*SFGrid,*SFrbGrid,*SGrid,*SrbGrid,mass,M5,mobius_b,mobius_c,params);

  std::cout << GridLogMessage << "****************************************************************** "<<std::endl;
  std::cout << GridLogMessage << " Calling DWF CG "<<std::endl;
  std::cout << GridLogMessage << "****************************************************************** "<<std::endl;

  MdagMLinearOperator<MobiusFermionR,FermionField> HermOp(Ddwf);
  MdagMLinearOperator<MobiusFermionR,FermionField> HermOpCk(Dchk);
  ConjugateGradient<FermionField> CG((stp),100000);
  s_res = zero;

  CG(HermOp,s_src,s_res);

  std::cout << GridLogMessage << " split residual norm "<<norm2(s_res)<<std::endl;
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
    std::cout << GridLogMessage<< " res["<<n<<"] norm "<<norm2(result[n])<<std::endl;
    HermOpCk.HermOp(result[n],tmp); tmp = tmp - src[n];
    std::cout << GridLogMessage<<" resid["<<n<<"]  "<< std::sqrt(norm2(tmp)/norm2(src[n]))<<std::endl;
  }


  for(int s=0;s<nrhs;s++){
    result[s]=zero;
  }

  /////////////////////////////////////////////////////////////
  // Try block CG
  /////////////////////////////////////////////////////////////
  int blockDim = 0;//not used for BlockCGVec
  BlockConjugateGradient<FermionField>    BCGV  (BlockCGrQVec,blockDim,stp,100000);
  {
    BCGV(HermOpCk,src,result);
  }

  
  
  Grid_finalize();
}
