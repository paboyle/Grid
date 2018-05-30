    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_pcg.cc

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
#include<bitset>
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

//Preconditioning:   M psi = chi
//               =  M P^-1 P psi = chi
//               =  M P^-1 psi' = chi

//Solve for psi' using M P^-1 as operator, then apply  P^-1 psi' = psi 

//Inexact preconditioned CG requires slight modification because we want to avoid computing P^-1 exactly


/////////////////////////////////////////////////////////////
// Base classes for iterative processes based on operators
// single input vec, single output vec.
/////////////////////////////////////////////////////////////

template<typename T>
T parse(const std::string &name, std::istream &in){
  std::string p;
  in >> p;
  assert(p==name);
  char eq;
  in >> eq;
  assert(eq == '=');
  T out;
  in >> out;
  return out;
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int Ls=12;
  RealD mass=0.01;
  RealD outer_tol = 1e-8;
  RealD inner_tol_full = 1e-5;
  RealD inner_tol_half = 1e-5;
  RealD inner_tol_16c = 1e-5;
  RealD inner_tol_8c = 1e-5;

  RealD relup_delta_full = 0.1;
  RealD relup_delta_half = 0.1;
  RealD relup_delta_16c = 0.1;
  RealD relup_delta_8c = 0.1;

  std::string config_file = "";
  RealD lo = 1.0;
  RealD hi = 64.0;
  int   order=10;
  for(int i=1;i<argc;i++){
    if(std::string(argv[i]) == "--params"){
      std::ifstream f(argv[i+1]);
      f.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
      Ls = parse<int>("Ls", f);
#define PARSEIT(NM) NM = parse<RealD>(#NM, f)
      PARSEIT(mass);
      PARSEIT(outer_tol);
      PARSEIT(inner_tol_full);
      PARSEIT(inner_tol_half);
      PARSEIT(inner_tol_16c);
      PARSEIT(inner_tol_8c);
      PARSEIT(relup_delta_full);
      PARSEIT(relup_delta_half);
      PARSEIT(relup_delta_16c);
      PARSEIT(relup_delta_8c);
#undef PARSEIT

      //f >> outer_tol >> inner_tol_full >> inner_tol_half >> inner_tol_16c >> inner_tol_8c;      
    }else if(std::string(argv[i]) == "--config"){
      config_file = argv[i+1];
    }else if(std::string(argv[i]) == "--order"){
      std::string ss(argv[i+1]);
      std::stringstream f(ss);
      f>>order; std::cout << " Order poly set to " <<order<<std::endl;
    }else if(std::string(argv[i]) == "--lo"){
      std::string ss(argv[i+1]);
      std::stringstream f(ss);
      f>>lo; std::cout << " Lo poly set to " <<lo<<std::endl;
    }else if(std::string(argv[i]) == "--hi"){
      std::string ss(argv[i+1]);
      std::stringstream f(ss);
      f>>hi; std::cout << " Hi poly set to " <<hi<<std::endl;
    }
  }
  
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
  LatticeFermionD result(FGrid); result=zero;
  LatticeGaugeFieldD Umu(UGrid);
  LatticeGaugeFieldF Umu_f(UGrid_f); 

  if(config_file.size() > 0){
    FieldMetaData header;
    NerscIO::readConfiguration(Umu,header,config_file);
  }else{
    SU3::HotConfiguration(RNG4,Umu);
  }

  precisionChange(Umu_f,Umu);
  
  RealD M5=1.8;

  LatticeFermionD    src_o(FrbGrid);
  pickCheckerboard(Odd,src_o,src);

  //if(0){ //Test preconditioned CG
    std::cout << "Test preconditioned CG" << std::endl;
    LatticeFermionD result_o(FrbGrid);
    LatticeFermionD result_o_2(FrbGrid);
    result_o.checkerboard = Odd;
    result_o = zero;
    result_o_2.checkerboard = Odd;
    result_o_2 = zero;

    DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
    SchurDiagMooeeOperator<DomainWallFermionD,LatticeFermionD> HermOpEO(Ddwf);
    //DoNothingLinearOperator<LatticeFermionD> Prec;
    //FixedIterConjugateGradientPreconditioner<LatticeFermionD> Prec(HermOpEO, 20);
    //    SloppyConjugateGradientPreconditioner<LatticeFermionD> Prec(HermOpEO, 1e-2, 1000);
    PolynomialPreconditioner<LatticeFermionD> Prec(HermOpEO,lo,hi,order) ;

    std::cout << "Preconditioned CG" << std::endl;
    InexactPreconditionedConjugateGradient<LatticeFermionD> pCG(Prec,1.0e-8,10000);
    pCG(HermOpEO,src_o,result_o);

    std::cout << "Starting regular CG" << std::endl;
    ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
    CG(HermOpEO,src_o,result_o_2);

    LatticeFermionD diff_o(FrbGrid);
    RealD diff = axpy_norm(diff_o, -1.0, result_o, result_o_2);

    std::cout << "pCG HermOp applications " << " Lo " << lo << " Hi " << hi << " Order " << order << " " <<  pCG.IterationsToComplete << "(outer) + " << Prec.InnerIterations << "(inner) = " << pCG.IterationsToComplete + Prec.InnerIterations << std::endl;
    std::cout << "CG HermOp applications " << CG.IterationsToComplete << std::endl;
    std::cout << "Diff between results: " << diff << std::endl;
  //}

  if(0){ //Test compressor
    LatticeFermionD result_o(FrbGrid);
    LatticeFermionD result_o_2(FrbGrid);
    result_o.checkerboard = Odd;
    result_o = zero;
    result_o_2.checkerboard = Odd;
    result_o_2 = zero;

    DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
    SchurDiagMooeeOperator<DomainWallFermionD,LatticeFermionD> HermOpEO(Ddwf);

    DomainWallFermionDF DdwfC(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
    SchurDiagMooeeOperator<DomainWallFermionDF,LatticeFermionD> HermOpEOC(DdwfC);

    std::cout << "Starting regular CG with compressed operator" << std::endl;
    Integer iter1;
    {
      ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
      CG.ErrorOnNoConverge = false;
      CG(HermOpEOC,src_o,result_o);
      iter1 = CG.IterationsToComplete;
    }
    Integer iter2;
    {
      std::cout << "Starting regular CG" << std::endl;
      ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
      CG(HermOpEO,src_o,result_o_2);
      iter2 = CG.IterationsToComplete;
    }

    LatticeFermionD diff_o(FrbGrid);
    RealD diff = axpy_norm(diff_o, -1.0, result_o, result_o_2);

    std::cout << "CG HermOp CC applications " << iter1 << std::endl;
    std::cout << "CG HermOp applications " << iter2 << std::endl;
    std::cout << "Diff between results: " << diff << std::endl;
  }
  
  if(1){ //Compare mixed prec restarted single/single internal with same but with single/compressed
    LatticeFermionD result_o_full(FrbGrid);
    LatticeFermionD result_o_half(FrbGrid);
    LatticeFermionD result_o_16(FrbGrid);
    LatticeFermionD result_o_8(FrbGrid);
    result_o_full.checkerboard = Odd;
    result_o_full = zero;
    result_o_16 = result_o_8 = result_o_half = result_o_full;

    //Std
    DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
    SchurDiagMooeeOperator<DomainWallFermionD,LatticeFermionD> HermOpEO(Ddwf);

    DomainWallFermionF Ddwf_f(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);
    SchurDiagMooeeOperator<DomainWallFermionF,LatticeFermionF> HermOpEO_f(Ddwf_f);

    //1/2 prec
    DomainWallFermionFH Ddwfhalf_f(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);
    SchurDiagMooeeOperator<DomainWallFermionFH,LatticeFermionF> HermOpEOhalf_f(Ddwfhalf_f);

  }


  Grid_finalize();
}
