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

  int threads = GridThread::GetThreads();
  std::cout << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  std::vector<int> latt4 = GridDefaultLatt();
  std::vector<int> simd4 = GridDefaultSimd(Nd,vComplexF::Nsimd());
  std::vector<int> mpi4  = GridDefaultMpi();

  assert(latt4.size()==4 ); 
  assert(simd4.size()==4 );
  assert(mpi4.size() ==4 );

  const int Ls=1;
  std::vector<int> latt5({Ls,latt4[0],latt4[1],latt4[2],latt4[3]});
  std::vector<int> simd5({1 ,simd4[0],simd4[1],simd4[2],simd4[3]}); 
  std::vector<int>  mpi5({1 , mpi4[0], mpi4[1], mpi4[2], mpi4[3]}); 
  std::vector<int>   cb5({0,1,1,1,1}); // Checkerboard 4d only
  int                cbd=1;            // use dim-1 to reduce

  // Four dim grid for gauge field U
  GridCartesian               UGrid(latt4,simd4,mpi4); 
  GridRedBlackCartesian     UrbGrid(&UGrid);

  // Five dim grid for fermions F
  GridCartesian               FGrid(latt5,simd5,mpi5); 
  GridRedBlackCartesian     FrbGrid(latt5,simd5,mpi5,cb5,cbd); 

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});

  GridParallelRNG          RNG5(&FGrid);  RNG5.SeedFixedIntegers(seeds5);
  LatticeFermion src   (&FGrid); random(RNG5,src);
  LatticeFermion result(&FGrid); result=zero;
  LatticeFermion    ref(&FGrid);    ref=zero;
  LatticeFermion    tmp(&FGrid);
  LatticeFermion    err(&FGrid);

  ColourMatrix cm = Complex(1.0,0.0);

  GridParallelRNG          RNG4(&UGrid);  RNG4.SeedFixedIntegers(seeds4);
  LatticeGaugeField Umu(&UGrid); random(RNG4,Umu);
  LatticeGaugeField Umu5d(&FGrid); 

  // replicate across fifth dimension
  for(int ss=0;ss<Umu._grid->oSites();ss++){
    for(int s=0;s<Ls;s++){
      Umu5d._odata[Ls*ss+s] = Umu._odata[ss];
    }
  }

  ////////////////////////////////////
  // Naive wilson implementation
  ////////////////////////////////////
  std::vector<LatticeColourMatrix> U(4,&FGrid);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = peekIndex<LorentzIndex>(Umu5d,mu);
  }

  if (1)
  {
    ref = zero;
    for(int mu=0;mu<Nd;mu++){

      tmp = U[mu]*Cshift(src,mu+1,1);
      ref=ref + tmp + Gamma(Gmu[mu])*tmp;

      tmp =adj(U[mu])*src;
      tmp =Cshift(tmp,mu+1,-1);
      ref=ref + tmp - Gamma(Gmu[mu])*tmp;
    }
    ref = -0.5*ref;
  }

  RealD mass=0.1;
  FiveDimWilsonFermion Dw(Umu,FGrid,FrbGrid,UGrid,UrbGrid,mass);
  
  std::cout << "Calling Dw"<<std::endl;
  int ncall=1000;
  double t0=usecond();
  for(int i=0;i<ncall;i++){
    Dw.Dhop(src,result,0);
  }
  double t1=usecond();


  double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
  double flops=1344*volume*ncall;
  
  std::cout << "Called Dw"<<std::endl;
  std::cout << "norm result "<< norm2(result)<<std::endl;
  std::cout << "norm ref    "<< norm2(ref)<<std::endl;
  std::cout << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
  err = ref-result; 
  std::cout << "norm diff   "<< norm2(err)<<std::endl;


  if (1)
  { // Naive wilson dag implementation
    ref = zero;
    for(int mu=0;mu<Nd;mu++){

      //    ref =  src - Gamma(Gamma::GammaX)* src ; // 1+gamma_x
      tmp = U[mu]*Cshift(src,mu+1,1);
      for(int i=0;i<ref._odata.size();i++){
	ref._odata[i]+= tmp._odata[i] - Gamma(Gmu[mu])*tmp._odata[i]; ;
      }

      tmp =adj(U[mu])*src;
      tmp =Cshift(tmp,mu+1,-1);
      for(int i=0;i<ref._odata.size();i++){
	ref._odata[i]+= tmp._odata[i] + Gamma(Gmu[mu])*tmp._odata[i]; ;
      }
    }
    ref = -0.5*ref;
  }
  Dw.Dhop(src,result,1);
  std::cout << "Called DwDag"<<std::endl;
  std::cout << "norm result "<< norm2(result)<<std::endl;
  std::cout << "norm ref    "<< norm2(ref)<<std::endl;
  err = ref-result; 
  std::cout << "norm diff   "<< norm2(err)<<std::endl;

  Grid_finalize();
}
