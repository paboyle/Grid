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

  std::vector<int> simd_layout({1,1,2,2});
  std::vector<int> mpi_layout ({1,1,1,1});
  std::vector<int> latt_size  ({8,8,8,8});
    
  GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
  std::vector<int> seeds({1,2,3,4});

  GridParallelRNG          pRNG(&Grid);
  //  std::vector<int> seeds({1,2,3,4});
  //  pRNG.SeedFixedIntegers(seeds);
  pRNG.SeedRandomDevice();

  LatticeFermion src(&Grid); random(pRNG,src);
  LatticeFermion result(&Grid); result=zero;
  LatticeFermion    ref(&Grid);    ref=zero;
  LatticeFermion    err(&Grid);    
  LatticeFermion    tmp(&Grid);    tmp=zero;
  LatticeGaugeField Umu(&Grid); random(pRNG,Umu);
  std::vector<LatticeColourMatrix> U(4,&Grid);

  double volume=1;
  for(int mu=0;mu<Nd;mu++){
    volume=volume*latt_size[mu];
  }  

  for(int mu=0;mu<Nd;mu++){
    U[mu] = peekIndex<LorentzIndex>(Umu,mu);
  }
  
  std::vector<int> mask({1,1,1,1,1,1,1,1});
  { // Naive wilson implementation
    ref = zero;
    for(int mu=0;mu<Nd;mu++){

      //    ref =  src + Gamma(Gamma::GammaX)* src ; // 1-gamma_x
      if( mask[mu] ) {
	tmp = U[mu]*Cshift(src,mu,1);
	for(int i=0;i<ref._odata.size();i++){
	  ref._odata[i]+= tmp._odata[i] + Gamma(Gmu[mu])*tmp._odata[i]; ;
	}
      }

      if( mask[mu+4] ){
	tmp =adj(U[mu])*src;
	tmp =Cshift(tmp,mu,-1);
	for(int i=0;i<ref._odata.size();i++){
	  ref._odata[i]+= tmp._odata[i] - Gamma(Gmu[mu])*tmp._odata[i]; ;
	}
      }
    }
  }

  RealD mass=0.1;
  WilsonMatrix Dw(Umu,mass);
  
  std::cout << "Calling Dw"<<std::endl;
  int ncall=100;
  double t0=usecond();
  for(int i=0;i<ncall;i++){
    Dw.multiply(src,result);
  }
  double t1=usecond();
  double flops=1320*volume*ncall;
  
  std::cout << "Called Dw"<<std::endl;
  std::cout << "norm result "<< norm2(result)<<std::endl;
  std::cout << "norm ref    "<< norm2(ref)<<std::endl;
  std::cout << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
  err = ref -result;
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

  Grid_finalize();
}
