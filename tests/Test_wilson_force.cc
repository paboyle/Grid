#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

#define parallel_for PARALLEL_FOR_LOOP for

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();

  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(latt_size,simd_layout,mpi_layout);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  std::vector<int> seeds({1,2,3,4});

  GridParallelRNG          pRNG(&Grid);
  pRNG.SeedRandomDevice();

  LatticeFermion phi        (&Grid); gaussian(pRNG,phi);
  LatticeFermion Mphi       (&Grid); 
  LatticeFermion MphiPrime  (&Grid); 

  LatticeGaugeField U(&Grid);

  SU3::HotConfiguration(pRNG,U);
  //  SU3::ColdConfiguration(pRNG,U);
  
  ////////////////////////////////////
  // Unmodified matrix element
  ////////////////////////////////////
  RealD mass=-4.0; //kills the diagonal term
  WilsonFermion Dw     (U,     Grid,RBGrid,mass);
  Dw.M   (phi,Mphi);

  ComplexD S    = innerProduct(Mphi,Mphi); // pdag MdagM p

  // get the deriv of phidag MdagM phi with respect to "U"
  LatticeGaugeField UdSdU(&Grid);
  LatticeGaugeField tmp(&Grid);

  Dw.MDeriv(tmp , Mphi,  phi,DaggerNo );  UdSdU=tmp;
  Dw.MDeriv(tmp , phi,  Mphi,DaggerYes ); UdSdU=UdSdU+tmp;


  LatticeFermion Ftmp      (&Grid);

  ////////////////////////////////////
  // Modify the gauge field a little 
  ////////////////////////////////////
  RealD dt = 1.0e-6;

  LatticeColourMatrix mommu(&Grid); 
  LatticeGaugeField mom(&Grid); 
  LatticeGaugeField Uprime(&Grid); 

  for(int mu=0;mu<Nd;mu++){

    SU3::GaussianLieAlgebraMatrix(pRNG, mommu); // Traceless antihermitian momentum; gaussian in lie alg

    PokeIndex<LorentzIndex>(mom,mommu,mu);
    parallel_for(auto i=mom.begin();i<mom.end();i++){
      Uprime[i](mu) =U[i](mu)+ mom[i](mu)*U[i](mu)*dt;
    }

  }

  Dw.DoubleStore(Dw.Umu,Uprime);
  Dw.M          (phi,MphiPrime);

  ComplexD Sprime    = innerProduct(MphiPrime   ,MphiPrime);

  //////////////////////////////////////////////
  // Use derivative to estimate dS
  //////////////////////////////////////////////
  LatticeComplex dS(&Grid); dS = zero;

  parallel_for(auto i=mom.begin();i<mom.end();i++){
    for(int mu=0;mu<Nd;mu++){
      //      dS[i]() = dS[i]()+trace(mom[i](mu) * UdSdU[i](mu) - mom[i](mu)* adj( UdSdU[i](mu)) )*dt;
      dS[i]() =    dS[i]()+trace(mom[i](mu) * (UdSdU[i](mu)))*dt;
      dS[i]() =    dS[i]()-trace(mom[i](mu) * adj(UdSdU[i](mu)))*dt;
    }
  }
  Complex dSpred    = sum(dS);

  std::cout << GridLogMessage << " S      "<<S<<std::endl;
  std::cout << GridLogMessage << " Sprime "<<Sprime<<std::endl;
  std::cout << GridLogMessage << "dS      "<<Sprime-S<<std::endl;
  std::cout << GridLogMessage << "predict dS    "<< dSpred <<std::endl;

  std::cout<< GridLogMessage << "Done" <<std::endl;
  Grid_finalize();
}
