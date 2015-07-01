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

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);


  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(latt_size,simd_layout,mpi_layout);

  int threads = GridThread::GetThreads();
  std::cout << "Grid is setup to use "<<threads<<" threads"<<std::endl;
  std::cout << "Grid floating point word size is REALF"<< sizeof(RealF)<<std::endl;
  std::cout << "Grid floating point word size is REALD"<< sizeof(RealD)<<std::endl;
  std::cout << "Grid floating point word size is REAL"<< sizeof(Real)<<std::endl;

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&Grid);
  pRNG.SeedFixedIntegers(seeds);
  //  pRNG.SeedRandomDevice();

  LatticeFermion src   (&Grid); random(pRNG,src);
  LatticeFermion result(&Grid); result=zero;
  LatticeFermion    ref(&Grid);    ref=zero;
  LatticeFermion    tmp(&Grid);    tmp=zero;
  LatticeFermion    err(&Grid);    tmp=zero;
  LatticeGaugeField Umu(&Grid); random(pRNG,Umu);
  std::vector<LatticeColourMatrix> U(4,&Grid);

  double volume=1;
  for(int mu=0;mu<Nd;mu++){
    volume=volume*latt_size[mu];
  }  

  // Only one non-zero (y)
#if 0
  Umu=zero;
  Complex cone(1.0,0.0);
  for(int nn=0;nn<Nd;nn++){
    random(pRNG,U[nn]);
    if(0) {
      if (nn==-1) { U[nn]=zero; std::cout << "zeroing gauge field in dir "<<nn<<std::endl; }
      else       { U[nn] = cone;std::cout << "unit gauge field in dir "<<nn<<std::endl; }
    }
    pokeIndex<LorentzIndex>(Umu,U[nn],nn);
  }
#endif

  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }
  
  { // Naive wilson implementation
    ref = zero;
    for(int mu=0;mu<Nd;mu++){
      //    ref =  src + Gamma(Gamma::GammaX)* src ; // 1-gamma_x
      tmp = U[mu]*Cshift(src,mu,1);
      for(int i=0;i<ref._odata.size();i++){
	ref._odata[i]+= tmp._odata[i] + Gamma(Gmu[mu])*tmp._odata[i]; ;
      }

      tmp =adj(U[mu])*src;
      tmp =Cshift(tmp,mu,-1);
      for(int i=0;i<ref._odata.size();i++){
	ref._odata[i]+= tmp._odata[i] - Gamma(Gmu[mu])*tmp._odata[i]; ;
      }
    }
  }
  ref = -0.5*ref;
  RealD mass=0.1;
  WilsonFermion Dw(Umu,Grid,RBGrid,mass);
  
  std::cout << "Calling Dw"<<std::endl;
  int ncall=10000;
  double t0=usecond();
  for(int i=0;i<ncall;i++){
    Dw.Dhop(src,result,0);
  }
  double t1=usecond();
  double flops=1344*volume*ncall;
  
  std::cout << "Called Dw"<<std::endl;
  std::cout << "norm result "<< norm2(result)<<std::endl;
  std::cout << "norm ref    "<< norm2(ref)<<std::endl;
  std::cout << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
  err = ref-result; 
  std::cout << "norm diff   "<< norm2(err)<<std::endl;


  //  for(int ss=0;ss<10;ss++ ){
  for(int ss=0;ss<0;ss++ ){
    for(int i=0;i<Ns;i++){
      for(int j=0;j<Nc;j++){
	ComplexF * ref_p = (ComplexF *)&ref._odata[ss]()(i)(j);
	ComplexF * res_p = (ComplexF *)&result._odata[ss]()(i)(j);
	std::cout << ss<< " "<<i<<" "<<j<<" "<< (*ref_p)<<" " <<(*res_p)<<std::endl;
      }
    }
  }

  { // Naive wilson dag implementation
    ref = zero;
    for(int mu=0;mu<Nd;mu++){

      //    ref =  src - Gamma(Gamma::GammaX)* src ; // 1+gamma_x
      tmp = U[mu]*Cshift(src,mu,1);
      for(int i=0;i<ref._odata.size();i++){
	ref._odata[i]+= tmp._odata[i] - Gamma(Gmu[mu])*tmp._odata[i]; ;
      }

      tmp =adj(U[mu])*src;
      tmp =Cshift(tmp,mu,-1);
      for(int i=0;i<ref._odata.size();i++){
	ref._odata[i]+= tmp._odata[i] + Gamma(Gmu[mu])*tmp._odata[i]; ;
      }
    }
  }
  ref = -0.5*ref;
  Dw.Dhop(src,result,1);
  std::cout << "Called DwDag"<<std::endl;
  std::cout << "norm result "<< norm2(result)<<std::endl;
  std::cout << "norm ref    "<< norm2(ref)<<std::endl;
  err = ref-result; 
  std::cout << "norm diff   "<< norm2(err)<<std::endl;

  Grid_finalize();
}
