    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/solver/Test_dwf_relupcg_prec.cc

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@bnl.gov>
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
  Grid_init(&argc,&argv);

  double relup_delta = 0.2;
  for(int i=1;i<argc-1;i++){
    std::string sarg = argv[i];
    if(sarg == "--relup_delta"){
      std::stringstream ss; ss << argv[i+1]; ss >> relup_delta;
      std::cout << GridLogMessage << "Set reliable update Delta to " << relup_delta << std::endl;
    }
  }   
  
  const int Ls=12;

  { 
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

  std::cout << GridLogMessage << "::::::::::::: Starting mixed CG" << std::endl;
  ConjugateGradientReliableUpdate<LatticeFermionD,LatticeFermionF> mCG(1e-8, 10000, relup_delta, FrbGrid_f, HermOpEO_f, HermOpEO);
  double t1,t2,flops;
  double MdagMsiteflops = 1452; // Mobius (real coeffs)
  // CG overhead: 8 inner product, 4+8 axpy_norm, 4+4 linear comb (2 of)
  double CGsiteflops = (8+4+8+4+4)*Nc*Ns ;
  std:: cout << " MdagM site flops = "<< 4*MdagMsiteflops<<std::endl;
  std:: cout << " CG    site flops = "<< CGsiteflops <<std::endl;
  int iters, iters_cleanup, relups, tot_iters;
  for(int i=0;i<10;i++){
    result_o = Zero();
    t1=usecond();
    mCG(src_o,result_o);
    t2=usecond();
    iters = mCG.IterationsToComplete; //Number of single prec CG iterations
    iters_cleanup = mCG.IterationsToCleanup;
    relups = mCG.ReliableUpdatesPerformed;
    tot_iters  = iters + iters_cleanup + relups; //relup cost MdagM application in double
    
    flops = MdagMsiteflops*4*FrbGrid->gSites()*tot_iters;
    flops+= CGsiteflops*FrbGrid->gSites()*tot_iters;
    std::cout << " SinglePrecision single prec iterations/sec "<< iters/(t2-t1)*1000.*1000.<<std::endl;
    std::cout << " SinglePrecision double prec cleanup iterations/sec "<< iters_cleanup/(t2-t1)*1000.*1000.<<std::endl;
    std::cout << " SinglePrecision reliable updates/sec "<< relups/(t2-t1)*1000.*1000.<<std::endl;
    std::cout << " SinglePrecision GF/s "<< flops/(t2-t1)/1000.<<std::endl;
  }
  std::cout << GridLogMessage << "::::::::::::: Starting regular CG" << std::endl;
  ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
  for(int i=0;i<1;i++){
    result_o_2 = Zero();
    t1=usecond();
    CG(HermOpEO,src_o,result_o_2);
    t2=usecond();
    iters = CG.IterationsToComplete;
    flops = MdagMsiteflops*4*FrbGrid->gSites()*iters; 
    flops+= CGsiteflops*FrbGrid->gSites()*iters;
    
    std::cout << " DoublePrecision iterations/sec "<< iters/(t2-t1)*1000.*1000.<<std::endl;
    std::cout << " DoublePrecision GF/s "<< flops/(t2-t1)/1000.<<std::endl;
  }
  
  //  MemoryManager::Print();

  LatticeFermionD diff_o(FrbGrid);
  RealD diff = axpy_norm(diff_o, -1.0, result_o, result_o_2);

  std::cout << GridLogMessage << "::::::::::::: Diff between mixed and regular CG: " << diff << std::endl;
  }
  
  MemoryManager::Print();

  Grid_finalize();
}
