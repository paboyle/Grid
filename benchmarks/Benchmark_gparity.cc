#include <Grid/Grid.h>
#include <sstream>
using namespace std;
using namespace Grid;
using namespace Grid::QCD;

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

typedef typename GparityDomainWallFermionF::FermionField GparityLatticeFermionF;
typedef typename GparityDomainWallFermionD::FermionField GparityLatticeFermionD;



int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int Ls=16;
  for(int i=0;i<argc;i++)
    if(std::string(argv[i]) == "-Ls"){
      std::stringstream ss(argv[i+1]); ss >> Ls;
    }


  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
  std::cout<<GridLogMessage << "Ls = " << Ls << std::endl;

  std::vector<int> latt4 = GridDefaultLatt();

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  
  std::cout << GridLogMessage << "Initialising 4d RNG" << std::endl;
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  std::cout << GridLogMessage << "Initialising 5d RNG" << std::endl;
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  std::cout << GridLogMessage << "Initialised RNGs" << std::endl;

  GparityLatticeFermionF src   (FGrid); random(RNG5,src);
  RealD N2 = 1.0/::sqrt(norm2(src));
  src = src*N2;

  GparityLatticeFermionF result(FGrid); result=zero;
  GparityLatticeFermionF    ref(FGrid);    ref=zero;
  GparityLatticeFermionF    tmp(FGrid);
  GparityLatticeFermionF    err(FGrid);

  std::cout << GridLogMessage << "Drawing gauge field" << std::endl;
  LatticeGaugeFieldF Umu(UGrid); 
  SU3::HotConfiguration(RNG4,Umu); 
  std::cout << GridLogMessage << "Random gauge initialised " << std::endl;

  RealD mass=0.1;
  RealD M5  =1.8;

  RealD NP = UGrid->_Nprocessors;
  RealD NN = UGrid->NodeCount();

  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Kernel options --dslash-generic, --dslash-unroll, --dslash-asm" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Benchmarking DomainWallFermion::Dhop                  "<<std::endl;
  std::cout << GridLogMessage<< "* Vectorising space-time by "<<vComplexF::Nsimd()<<std::endl;
#ifdef GRID_OMP
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute ) std::cout << GridLogMessage<< "* Using Overlapped Comms/Compute" <<std::endl;
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsThenCompute) std::cout << GridLogMessage<< "* Using sequential comms compute" <<std::endl;
#endif
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;



  std::cout << GridLogMessage<< "* SINGLE/SINGLE"<<std::endl;
  GparityDomainWallFermionF Dw(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  int ncall =1000;
  if (1) {
    FGrid->Barrier();
    Dw.ZeroCounters();
    Dw.Dhop(src,result,0);
    std::cout<<GridLogMessage<<"Called warmup"<<std::endl;
    double t0=usecond();
    for(int i=0;i<ncall;i++){
      __SSC_START;
      Dw.Dhop(src,result,0);
      __SSC_STOP;
    }
    double t1=usecond();
    FGrid->Barrier();
    
    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=2*1344*volume*ncall;

    std::cout<<GridLogMessage << "Called Dw "<<ncall<<" times in "<<t1-t0<<" us"<<std::endl;
    //    std::cout<<GridLogMessage << "norm result "<< norm2(result)<<std::endl;
    //    std::cout<<GridLogMessage << "norm ref    "<< norm2(ref)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per rank =  "<< flops/(t1-t0)/NP<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per node =  "<< flops/(t1-t0)/NN<<std::endl;
    Dw.Report();
  }

  std::cout << GridLogMessage<< "* SINGLE/HALF"<<std::endl;
  GparityDomainWallFermionFH DwH(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  if (1) {
    FGrid->Barrier();
    DwH.ZeroCounters();
    DwH.Dhop(src,result,0);
    double t0=usecond();
    for(int i=0;i<ncall;i++){
      __SSC_START;
      DwH.Dhop(src,result,0);
      __SSC_STOP;
    }
    double t1=usecond();
    FGrid->Barrier();
    
    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=2*1344*volume*ncall;

    std::cout<<GridLogMessage << "Called half prec comms Dw "<<ncall<<" times in "<<t1-t0<<" us"<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per rank =  "<< flops/(t1-t0)/NP<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per node =  "<< flops/(t1-t0)/NN<<std::endl;
    DwH.Report();
  }

  GridCartesian         * UGrid_d   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexD::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid_d = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_d);
  GridCartesian         * FGrid_d   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_d);
  GridRedBlackCartesian * FrbGrid_d = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_d);

  
  std::cout << GridLogMessage<< "* DOUBLE/DOUBLE"<<std::endl;
  GparityLatticeFermionD src_d(FGrid_d);
  precisionChange(src_d,src);

  LatticeGaugeFieldD Umu_d(UGrid_d); 
  precisionChange(Umu_d,Umu);

  GparityLatticeFermionD result_d(FGrid_d);

  GparityDomainWallFermionD DwD(Umu_d,*FGrid_d,*FrbGrid_d,*UGrid_d,*UrbGrid_d,mass,M5);
  if (1) {
    FGrid_d->Barrier();
    DwD.ZeroCounters();
    DwD.Dhop(src_d,result_d,0);
    std::cout<<GridLogMessage<<"Called warmup"<<std::endl;
    double t0=usecond();
    for(int i=0;i<ncall;i++){
      __SSC_START;
      DwD.Dhop(src_d,result_d,0);
      __SSC_STOP;
    }
    double t1=usecond();
    FGrid_d->Barrier();
    
    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=2*1344*volume*ncall;

    std::cout<<GridLogMessage << "Called Dw "<<ncall<<" times in "<<t1-t0<<" us"<<std::endl;
    //    std::cout<<GridLogMessage << "norm result "<< norm2(result)<<std::endl;
    //    std::cout<<GridLogMessage << "norm ref    "<< norm2(ref)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per rank =  "<< flops/(t1-t0)/NP<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per node =  "<< flops/(t1-t0)/NN<<std::endl;
    DwD.Report();
  }

  Grid_finalize();
}

