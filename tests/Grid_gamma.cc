#include <Grid.h>
#include <parallelIO/GridNerscIO.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> simd_layout({1,1,2,2});
  std::vector<int> mpi_layout ({1,1,1,1});
  std::vector<int> latt_size  ({8,8,8,8});
    
  GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

  GridParallelRNG          pRNG(&Grid);
  pRNG.SeedRandomDevice();

  GridSerialRNG            sRNG;
  sRNG.SeedRandomDevice();

  SpinMatrix ident=zero;
  SpinMatrix rnd  ; random(sRNG,rnd);

  SpinMatrix ll=zero;
  SpinMatrix rr=zero;
  SpinMatrix result;

  SpinVector lv=zero;
  SpinVector rv=zero;

  for(int a=0;a<Ns;a++){
    ident()(a,a) = 1.0;
  }

  Gamma::GammaMatrix g [] = {
    Gamma::Identity,
    Gamma::GammaX,
    Gamma::GammaY,
    Gamma::GammaZ,
    Gamma::GammaT,
    Gamma::Gamma5,
    Gamma::MinusIdentity,
    Gamma::MinusGammaX,
    Gamma::MinusGammaY,
    Gamma::MinusGammaZ,
    Gamma::MinusGammaT,
    Gamma::MinusGamma5
  };
  const char *list[] = { 
    "Identity ",
    "GammaX   ",
    "GammaY   ",
    "GammaZ   ",
    "GammaT   ",
    "Gamma5   ",
    "-Identity",
    "-GammaX  ",
    "-GammaY  ",
    "-GammaZ  ",
    "-GammaT  ",
    "-Gamma5  ",
    "         "
  };
  result =ll*Gamma(g[0])*rr;
  result =ll*Gamma(g[0]);
  rv = Gamma(g[0])*lv;

  for(int mu=0;mu<12;mu++){

    result = Gamma(g[mu])* ident;

    for(int i=0;i<Ns;i++){

      if(i==0) std::cout << list[mu];
      else     std::cout << list[12];

      std::cout<<"(";
      for(int j=0;j<Ns;j++){
	if ( abs(result()(i,j)())==0 ) {
	  std::cout<< " 0";
	} else if ( abs(result()(i,j)() - ComplexF(0,1))==0){
	  std::cout<< " i";
	} else if ( abs(result()(i,j)() + ComplexF(0,1))==0){
	  std::cout<< "-i";
	} else if ( abs(result()(i,j)() - ComplexF(1,0))==0){
	  std::cout<< " 1";
	} else if ( abs(result()(i,j)() + ComplexF(1,0))==0){
	  std::cout<< "-1";
	}
	std::cout<< ((j==Ns-1) ? ")" : "," );
      }
      std::cout << std::endl;
    }

    std::cout << std::endl;

  }

  std::cout << "Testing Gamma^2 - 1 = 0"<<std::endl;
  for(int mu=0;mu<6;mu++){
    result =  Gamma(g[mu])* ident * Gamma(g[mu]);
    result = result - ident;
    double mag = TensorRemove(norm2l(result));
    std::cout << list[mu]<<" " << mag<<std::endl;
  }

  std::cout << "Testing (MinusGamma + G )M = 0"<<std::endl;
  for(int mu=0;mu<6;mu++){
    result =          rnd * Gamma(g[mu]);
    result = result + rnd * Gamma(g[mu+6]);
    double mag = TensorRemove(norm2l(result));
    std::cout << list[mu]<<" " << mag<<std::endl;
  }

  std::cout << "Testing M(MinusGamma + G )  = 0"<<std::endl;
  for(int mu=0;mu<6;mu++){
    result =           Gamma(g[mu])  *rnd;
    result = result +  Gamma(g[mu+6])*rnd;
    double mag = TensorRemove(norm2l(result));
    std::cout << list[mu]<<" " << mag<<std::endl;
  }

  Grid_finalize();
}
