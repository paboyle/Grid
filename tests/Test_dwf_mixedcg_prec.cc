    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_cg_prec.cc

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

int main (int argc, char ** argv)
{
  char hostname[HOST_NAME_MAX+1];
  gethostname(hostname, HOST_NAME_MAX+1);
  std::string host(hostname);
  
  Grid_init(&argc,&argv);

  const int Ls=12;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexD::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  GridCartesian         * UGrid_f   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid_f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_f);
  GridCartesian         * FGrid_f   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_f);
  GridRedBlackCartesian * FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_f);
  
  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeFermionD    src(FGrid); random(RNG5,src);
  LatticeFermionD result(FGrid); result=Zero();
  LatticeGaugeFieldD Umu(UGrid);
  LatticeGaugeFieldF Umu_f(UGrid_f); 
  
  SU<Nc>::HotConfiguration(RNG4,Umu);

  precisionChange(Umu_f,Umu);
  
  RealD mass=0.1;
  RealD M5=1.8;
  DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  DomainWallFermionF Ddwf_f(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);

  LatticeFermionD    src_o(FrbGrid);
  LatticeFermionD result_o(FrbGrid);
  LatticeFermionD result_o_2(FrbGrid);
  pickCheckerboard(Odd,src_o,src);
  result_o.Checkerboard() = Odd;
  result_o = Zero();
  result_o_2.Checkerboard() = Odd;
  result_o_2 = Zero();

  SchurDiagMooeeOperator<DomainWallFermionD,LatticeFermionD> HermOpEO(Ddwf);
  SchurDiagMooeeOperator<DomainWallFermionF,LatticeFermionF> HermOpEO_f(Ddwf_f);

  int nsecs=600;
  if( GridCmdOptionExists(argv,argv+argc,"--seconds") ){
    std::string arg = GridCmdOptionPayload(argv,argv+argc,"--seconds");
    GridCmdOptionInt(arg,nsecs);
  }
  
  std::cout << GridLogMessage << "::::::::::::: Starting mixed CG for "<<nsecs <<" seconds" << std::endl;

  MixedPrecisionConjugateGradient<LatticeFermionD,LatticeFermionF> mCG(1.0e-8, 10000, 50, FrbGrid_f, HermOpEO_f, HermOpEO);
  double t1,t2,flops;
  double MdagMsiteflops = 1452; // Mobius (real coeffs)
  // CG overhead: 8 inner product, 4+8 axpy_norm, 4+4 linear comb (2 of)
  double CGsiteflops = (8+4+8+4+4)*Nc*Ns ;
  std:: cout << " MdagM site flops = "<< 4*MdagMsiteflops<<std::endl;
  std:: cout << " CG    site flops = "<< CGsiteflops <<std::endl;
  int iters;

  time_t start = time(NULL);

  uint32_t csum, csumref;
  csumref=0;
  int iter=0;
  do {
    std::cerr << "******************* SINGLE PRECISION SOLVE "<<iter<<std::endl;
    result_o = Zero();
    t1=usecond();
    mCG(src_o,result_o);
    t2=usecond();
    iters = mCG.TotalInnerIterations; //Number of inner CG iterations
    flops = MdagMsiteflops*4*FrbGrid->gSites()*iters;
    flops+= CGsiteflops*FrbGrid->gSites()*iters;
    std::cout << " SinglePrecision iterations/sec "<< iters/(t2-t1)*1000.*1000.<<std::endl;
    std::cout << " SinglePrecision GF/s "<< flops/(t2-t1)/1000.<<std::endl;

    csum = crc(result_o);

    if ( csumref == 0 ) {
      csumref = csum;
    } else {
      if ( csum != csumref ) { 
	std::cerr << host<<" FAILURE " <<iter <<" csum "<<std::hex<<csum<< " != "<<csumref <<std::dec<<std::endl;
	assert(0);
      } else {
	std::cout << host <<" OK " <<iter <<" csum "<<std::hex<<csum<<std::dec<<" -- OK! "<<std::endl;
      }
    }
    iter ++;
  } while (time(NULL) < (start + nsecs/2) );
    
  std::cout << GridLogMessage << "::::::::::::: Starting double precision CG" << std::endl;
  ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
  csumref=0;
  int i=0;
  do { 
    std::cerr << "******************* DOUBLE PRECISION SOLVE "<<i<<std::endl;
    result_o_2 = Zero();
    t1=usecond();
    CG(HermOpEO,src_o,result_o_2);
    t2=usecond();
    iters = CG.IterationsToComplete;
    flops = MdagMsiteflops*4*FrbGrid->gSites()*iters; 
    flops+= CGsiteflops*FrbGrid->gSites()*iters;

    std::cout << " DoublePrecision iterations/sec "<< iters/(t2-t1)*1000.*1000.<<std::endl;
    std::cout << " DoublePrecision GF/s "<< flops/(t2-t1)/1000.<<std::endl;

    csum = crc(result_o);

    if ( csumref == 0 ) {
      csumref = csum;
    } else {
      if ( csum != csumref ) { 
	std::cerr << i <<" csum "<<std::hex<<csum<< " != "<<csumref <<std::dec<<std::endl;
	assert(0);
      } else {
	std::cout << i <<" csum "<<std::hex<<csum<<std::dec<<" -- OK! "<<std::endl;
      }
    }
    i++;
  } while (time(NULL) < (start + nsecs) );

  LatticeFermionD diff_o(FrbGrid);
  RealD diff = axpy_norm(diff_o, -1.0, result_o, result_o_2);

  std::cout << GridLogMessage << "::::::::::::: Diff between mixed and regular CG: " << diff << std::endl;
  assert(diff < 1e-4);
  
  Grid_finalize();
}
