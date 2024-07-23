    /*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_cayley_cg.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/qcd/action/fermion/Reconstruct5Dprop.h>

using namespace std;
using namespace Grid;




template<class What> 
void  TestCGinversions(What & Ddwf, 
		       LatticeGaugeField &Umu,
		       GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
		       GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
		       RealD mass, RealD M5,
		       GridParallelRNG *RNG4,
		       GridParallelRNG *RNG5);
template<class What> 
void  TestCGschur(What & Ddwf, 
		  GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
		  GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
		  RealD mass, RealD M5,
		  GridParallelRNG *RNG4,
		  GridParallelRNG *RNG5);

template<class What> 
void  TestCGunprec(What & Ddwf, 
		   GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
		   GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
		   RealD mass, RealD M5,
		   GridParallelRNG *RNG4,
		   GridParallelRNG *RNG5);

template<class What> 
void  TestCGprec(What & Ddwf, 
		 GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
		 GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
		 RealD mass, RealD M5,
		 GridParallelRNG *RNG4,
		 GridParallelRNG *RNG5);

template<class What> 
void  TestReconstruct5D(What & Ddwf, 
			LatticeGaugeField &Umu,
			GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
			GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
			RealD mass, RealD M5,
			GridParallelRNG *RNG4,
			GridParallelRNG *RNG5);

template<class What,class WhatF> 
void  TestReconstruct5DFA(What & Ddwf, 
			  WhatF & DdwfF, 
			  LatticeGaugeField &Umu,
			  GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
			  GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
			  RealD mass, RealD M5,
			  GridParallelRNG *RNG4,
			  GridParallelRNG *RNG5);

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  const int Ls=8;
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), 
								   GridDefaultSimd(Nd,vComplex::Nsimd()),
								   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);


  GridCartesian         * UGridF   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), 
								    GridDefaultSimd(Nd,vComplexF::Nsimd()),
								    GridDefaultMpi());
  GridRedBlackCartesian * UrbGridF = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridF);
  GridCartesian         * FGridF   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGridF);
  GridRedBlackCartesian * FrbGridF = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGridF);


  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid);
  LatticeGaugeFieldF UmuF(UGridF);
  SU<Nc>::HotConfiguration(RNG4,Umu);
  precisionChange(UmuF,Umu);
  std::vector<LatticeColourMatrix> U(4,UGrid);

  RealD mass=0.1;
  RealD M5  =1.8;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"DomainWallFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  DomainWallFermionF DdwfF(UmuF,*FGridF,*FrbGridF,*UGridF,*UrbGridF,mass,M5);
  TestCGinversions<DomainWallFermionD>(Ddwf,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);
  TestReconstruct5DFA<DomainWallFermionD,DomainWallFermionF>(Ddwf,DdwfF,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  RealD b=1.5;// Scale factor b+c=2, b-c=1
  RealD c=0.5;
  std::vector<ComplexD> gamma(Ls,ComplexD(1.0,0.0));

  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"MobiusFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  MobiusFermionD Dmob(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);
  MobiusFermionF DmobF(UmuF,*FGridF,*FrbGridF,*UGridF,*UrbGridF,mass,M5,b,c);
  TestCGinversions<MobiusFermionD>(Dmob,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);
  TestReconstruct5DFA<MobiusFermionD,MobiusFermionF>(Dmob,DmobF,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"ZMobiusFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  ZMobiusFermionD ZDmob(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,gamma,b,c);
  TestCGinversions<ZMobiusFermionD>(ZDmob,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);
  TestReconstruct5D<ZMobiusFermionD>(ZDmob,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"MobiusZolotarevFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  MobiusZolotarevFermionD Dzolo(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c,0.1,2.0);
  TestCGinversions<MobiusZolotarevFermionD>(Dzolo,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);
  TestReconstruct5D<MobiusZolotarevFermionD>(Dzolo,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"ScaledShamirFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  ScaledShamirFermionD Dsham(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,2.0);
  ScaledShamirFermionF DshamF(UmuF,*FGridF,*FrbGridF,*UGridF,*UrbGridF,mass,M5,2.0);
  TestCGinversions<ScaledShamirFermionD>(Dsham,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);
  TestReconstruct5DFA<ScaledShamirFermionD,ScaledShamirFermionF>(Dsham,DshamF,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"ShamirZolotarevFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  ShamirZolotarevFermionD Dshamz(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,0.1,2.0);
  TestCGinversions<ShamirZolotarevFermionD>(Dshamz,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);
  TestReconstruct5D<ShamirZolotarevFermionD>(Dshamz,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"OverlapWilsonCayleyTanhFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  OverlapWilsonCayleyTanhFermionD Dov(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,1.0);
  OverlapWilsonCayleyTanhFermionF DovF(UmuF,*FGridF,*FrbGridF,*UGridF,*UrbGridF,mass,M5,1.0);
  TestCGinversions<OverlapWilsonCayleyTanhFermionD>(Dov,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);
  TestReconstruct5DFA<OverlapWilsonCayleyTanhFermionD,OverlapWilsonCayleyTanhFermionF>(Dov,DovF,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"OverlapWilsonCayleyZolotarevFermion test"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  OverlapWilsonCayleyZolotarevFermionD Dovz(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,0.1,2.0);
  TestCGinversions<OverlapWilsonCayleyZolotarevFermionD>(Dovz,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);
  TestReconstruct5D<OverlapWilsonCayleyZolotarevFermionD>(Dovz,Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  Grid_finalize();
}
template<class What> 
void  TestCGinversions(What & Ddwf, 
		       LatticeGaugeField &Umu,
		       GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
		       GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
		       RealD mass, RealD M5,
		       GridParallelRNG *RNG4,
		       GridParallelRNG *RNG5)
{
  std::cout<<GridLogMessage << "Testing unpreconditioned inverter"<<std::endl;
  TestCGunprec<What>(Ddwf,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,RNG4,RNG5);
  std::cout<<GridLogMessage << "Testing red black preconditioned inverter"<<std::endl;
  TestCGprec<What>(Ddwf,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,RNG4,RNG5);
  std::cout<<GridLogMessage << "Testing red black Schur inverter"<<std::endl;
  TestCGschur<What>(Ddwf,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,RNG4,RNG5);
}


template<class What> 
void  TestCGunprec(What & Ddwf, 
		   GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
		   GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
		   RealD mass, RealD M5,
		   GridParallelRNG *RNG4,
		   GridParallelRNG *RNG5)
{
  LatticeFermion src   (FGrid); random(*RNG5,src);
  LatticeFermion result(FGrid); result=Zero();

  MdagMLinearOperator<What,LatticeFermion> HermOp(Ddwf);
  ConjugateGradient<LatticeFermion> CG(1.0e-8,10000);
  CG(HermOp,src,result);

}
template<class What> 
void  TestCGprec(What & Ddwf, 
		 GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
		 GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
		 RealD mass, RealD M5,
		 GridParallelRNG *RNG4,
		 GridParallelRNG *RNG5)
{
  LatticeFermion src   (FGrid); random(*RNG5,src);
  LatticeFermion    src_o(FrbGrid);
  LatticeFermion result_o(FrbGrid);
  pickCheckerboard(Odd,src_o,src);
  result_o=Zero();

  SchurDiagMooeeOperator<What,LatticeFermion> HermOpEO(Ddwf);
  ConjugateGradient<LatticeFermion> CG(1.0e-8,10000);
  CG(HermOpEO,src_o,result_o);
}

template<class What> 
void  TestReconstruct5D(What & Ddwf, 
			LatticeGaugeField & Umu,
			GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
			GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
			RealD mass, RealD M5,
			GridParallelRNG *RNG4,
			GridParallelRNG *RNG5)
{
  LatticeFermion src4   (UGrid); random(*RNG4,src4);
  LatticeFermion res4   (UGrid); res4 = Zero();

  LatticeFermion src   (FGrid);
  LatticeFermion src_NE(FGrid);
  LatticeFermion result(FGrid);
  LatticeFermion result_rec(FGrid);
  LatticeFermion result_madwf(FGrid);

  MdagMLinearOperator<What,LatticeFermion> HermOp(Ddwf);
  double Resid = 1.0e-12;
  double Residi = 1.0e-6;
  ConjugateGradient<LatticeFermion> CG(Resid,10000);
  ConjugateGradient<LatticeFermion> CGi(Residi,10000);

  Ddwf.ImportPhysicalFermionSource(src4,src);
  Ddwf.Mdag(src,src_NE);
  CG(HermOp,src_NE,result);

  Ddwf.ExportPhysicalFermionSolution(result, res4);

  Ddwf.M(result,src_NE);
  src_NE = src_NE - src;
  std::cout <<GridLogMessage<< " True residual is " << norm2(src_NE)<<std::endl;

  std::cout <<GridLogMessage<< " Reconstructing " <<std::endl;

  ////////////////////////////
  // RBprec PV inverse
  ////////////////////////////
  typedef LatticeFermion Field;
  typedef SchurRedBlackDiagTwoSolve<Field> SchurSolverType; 
  typedef SchurRedBlackDiagTwoSolve<Field> SchurSolverTypei; 
  typedef PauliVillarsSolverRBprec<Field,SchurSolverType> PVinverter;
  SchurSolverType SchurSolver(CG);
  PVinverter      PVinverse(SchurSolver);

  Reconstruct5DfromPhysical<LatticeFermion,PVinverter> reconstructor(PVinverse);

  reconstructor(Ddwf,res4,src4,result_rec);

  std::cout <<GridLogMessage << "Result     "<<norm2(result)<<std::endl;
  std::cout <<GridLogMessage << "Result_rec "<<norm2(result_rec)<<std::endl;

  result_rec = result_rec - result;
  std::cout <<GridLogMessage << "Difference "<<norm2(result_rec)<<std::endl;

  //////////////////////////////
  // Now try MADWF
  //////////////////////////////
  SchurSolverTypei SchurSolveri(CGi);
  ZeroGuesser<LatticeFermion> Guess;
  MADWF<What,What,PVinverter,SchurSolverTypei,ZeroGuesser<LatticeFermion> > 
    madwf(Ddwf,Ddwf,PVinverse,SchurSolveri,Guess,Resid,10);
  
  madwf(src4,result_madwf);
  result_madwf = result_madwf - result;
  std::cout <<GridLogMessage << "Difference "<<norm2(result_madwf)<<std::endl;


}
template<class What,class WhatF> 
void  TestReconstruct5DFA(What & Ddwf, 
			  WhatF & DdwfF, 
			  LatticeGaugeField & Umu,
			  GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
			  GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
			  RealD mass, RealD M5,
			  GridParallelRNG *RNG4,
			  GridParallelRNG *RNG5)
{
  LatticeFermion src4   (UGrid); random(*RNG4,src4);
  LatticeFermion res4   (UGrid); res4 = Zero();

  LatticeFermion src   (FGrid);
  LatticeFermion src_NE(FGrid);
  LatticeFermion result(FGrid);
  LatticeFermion result_rec(FGrid);
  LatticeFermion result_madwf(FGrid);

  MdagMLinearOperator<What,LatticeFermion> HermOp(Ddwf);
  double Resid = 1.0e-12;
  double Residi = 1.0e-5;
  ConjugateGradient<LatticeFermion> CG(Resid,10000);
  ConjugateGradient<LatticeFermionF> CGi(Residi,10000);

  Ddwf.ImportPhysicalFermionSource(src4,src);
  Ddwf.Mdag(src,src_NE);
  CG(HermOp,src_NE,result);

  Ddwf.ExportPhysicalFermionSolution(result, res4);

  Ddwf.M(result,src_NE);
  src_NE = src_NE - src;
  std::cout <<GridLogMessage<< " True residual is " << norm2(src_NE)<<std::endl;

  std::cout <<GridLogMessage<< " Reconstructing " <<std::endl;

  ////////////////////////////
  // Fourier accel PV inverse
  ////////////////////////////
  typedef LatticeFermion Field;
  typedef LatticeFermionF FieldF;
  typedef SchurRedBlackDiagTwoSolve<FieldF> SchurSolverTypei; 
  typedef PauliVillarsSolverFourierAccel<LatticeFermion,LatticeGaugeField> PVinverter;
  PVinverter PVinverse(Umu,CG);

  Reconstruct5DfromPhysical<LatticeFermion,PVinverter> reconstructor(PVinverse);

  reconstructor(Ddwf,res4,src4,result_rec);

  std::cout <<GridLogMessage << "Result     "<<norm2(result)<<std::endl;
  std::cout <<GridLogMessage << "Result_rec "<<norm2(result_rec)<<std::endl;

  result_rec = result_rec - result;
  std::cout <<GridLogMessage << "Difference "<<norm2(result_rec)<<std::endl;

  //////////////////////////////
  // Now try MADWF
  //////////////////////////////
  SchurSolverTypei SchurSolver(CGi);
  ZeroGuesser<LatticeFermionF> Guess;
  MADWF<What,WhatF,PVinverter,SchurSolverTypei,ZeroGuesser<LatticeFermionF> > 
    madwf(Ddwf,DdwfF,PVinverse,SchurSolver,Guess,Resid,10);
  
  madwf(src4,result_madwf);
  result_madwf = result_madwf - result;
  std::cout <<GridLogMessage << "Difference "<<norm2(result_madwf)<<std::endl;

}



template<class What> 
void  TestCGschur(What & Ddwf, 
		   GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
		   GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
		   RealD mass, RealD M5,
		   GridParallelRNG *RNG4,
		   GridParallelRNG *RNG5)
{
  LatticeFermion src   (FGrid); random(*RNG5,src);
  LatticeFermion result1(FGrid); result1=Zero();
  LatticeFermion result2(FGrid); result2=Zero();
  LatticeFermion result3(FGrid); result3=Zero();

  ConjugateGradient<LatticeFermion> CG(1.0e-8,10000);
  SchurRedBlackDiagMooeeSolve<LatticeFermion> SchurSolver(CG);
  SchurSolver(Ddwf,src,result1);

  SchurRedBlackDiagOneSolve<LatticeFermion> SchurSolverSymm1(CG);
  SchurSolverSymm1(Ddwf,src,result2);

  SchurRedBlackDiagTwoSolve<LatticeFermion> SchurSolverSymm2(CG);
  SchurSolverSymm2(Ddwf,src,result3);

  std::cout << GridLogMessage << " Standard " <<norm2(result1)<<std::endl;

  std::cout << GridLogMessage << " Symm1    " <<norm2(result2)<<std::endl; 
  result2=result2-result1;
  std::cout << GridLogMessage << " diff " <<norm2(result2) <<std::endl; 

  std::cout << GridLogMessage << " Symm2    " <<norm2(result3)<<std::endl; 
  result3=result3-result1;
  std::cout << GridLogMessage << " diff " <<norm2(result3) <<std::endl; 
}
