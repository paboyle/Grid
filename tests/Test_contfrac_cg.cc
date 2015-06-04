#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

template<class d>
struct scal {
  d internal;
};

  Gamma::GammaMatrix Gmu [] = {
    Gamma::GammaX,
    Gamma::GammaY,
    Gamma::GammaZ,
    Gamma::GammaT
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
  std::cout << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  const int Ls=9;
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);


  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid); random(RNG4,Umu);
  std::vector<LatticeColourMatrix> U(4,UGrid);

  RealD mass=0.1;
  RealD M5  =1.8;


  std::cout <<"OverlapWilsonContFracTanhFermion  test"<<std::endl;
  OverlapWilsonContFracTanhFermion Dcf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,1.0);
  TestCGinversions<OverlapWilsonContFracTanhFermion>(Dcf,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  std::cout <<"OverlapWilsonContFracZolotarevFermion  test"<<std::endl;
  OverlapWilsonContFracZolotarevFermion Dcfz(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,0.1,6.0);
  TestCGinversions<OverlapWilsonContFracZolotarevFermion>(Dcfz,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);


  std::cout <<"OverlapWilsonPartialFractionTanhFermion  test"<<std::endl;
  OverlapWilsonPartialFractionTanhFermion Dpf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,1.0);
  TestCGinversions<OverlapWilsonPartialFractionTanhFermion>(Dpf,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  std::cout <<"OverlapWilsonPartialFractionZolotarevFermion  test"<<std::endl;
  OverlapWilsonPartialFractionZolotarevFermion Dpfz(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,0.1,6.0);
  TestCGinversions<OverlapWilsonPartialFractionZolotarevFermion>(Dpfz,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);


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
  std::cout << "Testing unpreconditioned inverter"<<std::endl;
  TestCGunprec<What>(Ddwf,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,RNG4,RNG5);
  std::cout << "Testing red black preconditioned inverter"<<std::endl;
  TestCGprec<What>(Ddwf,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,RNG4,RNG5);
  std::cout << "Testing red black Schur inverter"<<std::endl;
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
  LatticeFermion result(FGrid); result=zero;

  HermitianOperator<What,LatticeFermion> HermOp(Ddwf);
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
  result_o=zero;

  HermitianCheckerBoardedOperator<What,LatticeFermion> HermOpEO(Ddwf);
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
  LatticeFermion result(FGrid); result=zero;

  ConjugateGradient<LatticeFermion> CG(1.0e-8,10000);
  SchurRedBlackSolve<LatticeFermion> SchurSolver(CG);
  SchurSolver(Ddwf,src,result);
}
