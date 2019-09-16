    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_contfrac_cg.cc

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
 ;

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


template<class What> 
void  TestCGinversions(What & Ddwf, 
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

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  const int Ls=9;
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);


  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid); SU3::HotConfiguration(RNG4,Umu);
  std::vector<LatticeColourMatrix> U(4,UGrid);

  RealD mass=0.1;
  RealD M5  =1.8;


  std::cout<<GridLogMessage <<"OverlapWilsonContFracTanhFermion  test"<<std::endl;
  OverlapWilsonContFracTanhFermionR Dcf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,1.0);
  TestCGinversions<OverlapWilsonContFracTanhFermionR>(Dcf,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  std::cout<<GridLogMessage <<"OverlapWilsonContFracZolotarevFermion  test"<<std::endl;
  OverlapWilsonContFracZolotarevFermionR Dcfz(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,0.1,6.0);
  TestCGinversions<OverlapWilsonContFracZolotarevFermionR>(Dcfz,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);


  std::cout<<GridLogMessage <<"OverlapWilsonPartialFractionTanhFermion  test"<<std::endl;
  OverlapWilsonPartialFractionTanhFermionR Dpf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,1.0);
  TestCGinversions<OverlapWilsonPartialFractionTanhFermionR>(Dpf,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  std::cout<<GridLogMessage <<"OverlapWilsonPartialFractionZolotarevFermion  test"<<std::endl;
  OverlapWilsonPartialFractionZolotarevFermionR Dpfz(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,0.1,6.0);
  TestCGinversions<OverlapWilsonPartialFractionZolotarevFermionR>(Dpfz,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);


  Grid_finalize();
}
template<class What> 
void  TestCGinversions(What & Ddwf, 
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
void  TestCGschur(What & Ddwf, 
		   GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
		   GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
		   RealD mass, RealD M5,
		   GridParallelRNG *RNG4,
		   GridParallelRNG *RNG5)
{
  LatticeFermion src   (FGrid); random(*RNG5,src);
  LatticeFermion result(FGrid); result=Zero();

  ConjugateGradient<LatticeFermion> CG(1.0e-8,10000);
  SchurRedBlackDiagMooeeSolve<LatticeFermion> SchurSolver(CG);
  SchurSolver(Ddwf,src,result);
}
