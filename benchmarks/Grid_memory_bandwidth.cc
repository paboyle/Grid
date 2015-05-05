#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> simd_layout({1,2,2,2});
  std::vector<int> mpi_layout ({1,1,1,1});

  const int Nvec=8;
  typedef Lattice< iVector< vReal,Nvec> > LatticeVec;

  int Nloop=1000;

  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "= Benchmarking fused AXPY bandwidth"<<std::endl;
  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "  L  "<<"\t\t"<<"bytes"<<"\t\t"<<"GB/s"<<std::endl;
  std::cout << "----------------------------------------------------------"<<std::endl;

  for(int lat=4;lat<=24;lat+=4){

      std::vector<int> latt_size  ({lat,lat,lat,lat});

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      //GridParallelRNG          pRNG(&Grid);      pRNG.SeedRandomDevice();

      LatticeVec z(&Grid); //random(pRNG,z);
      LatticeVec x(&Grid); //random(pRNG,x);
      LatticeVec y(&Grid); //random(pRNG,y);
      double a=2.0;


      double start=usecond();
      for(int i=0;i<Nloop;i++){
	//	z=a*x+y;
	//   inline void axpy(Lattice<vobj> &ret,double a,const Lattice<vobj> &lhs,const Lattice<vobj> &rhs){
	axpy(z,a,x,y);
      }
      double stop=usecond();
      double time = (stop-start)/Nloop/1000;
      
      double bytes=3*lat*lat*lat*lat*Nvec*sizeof(Real);
      std::cout<<std::setprecision(3) << lat<<"\t\t"<<bytes<<"\t\t"<<bytes/time<<std::endl;

    }

  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "= Benchmarking a*x + y bandwidth"<<std::endl;
  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "  L  "<<"\t\t"<<"bytes"<<"\t\t"<<"GB/s"<<"\t\tGB/s (eff)"<<std::endl;
  std::cout << "----------------------------------------------------------"<<std::endl;

  
  for(int lat=4;lat<=24;lat+=4){

      std::vector<int> latt_size  ({lat,lat,lat,lat});

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      //GridParallelRNG          pRNG(&Grid);      pRNG.SeedRandomDevice();

      LatticeVec z(&Grid); //random(pRNG,z);
      LatticeVec x(&Grid); //random(pRNG,x);
      LatticeVec y(&Grid); //random(pRNG,y);
      double a=2.0;


      double start=usecond();
      for(int i=0;i<Nloop;i++){
	z=a*x+y;
      }
      double stop=usecond();
      double time = (stop-start)/Nloop/1000;
      
      double effbytes=3*lat*lat*lat*lat*Nvec*sizeof(Real);
      double bytes=5*lat*lat*lat*lat*Nvec*sizeof(Real);
      std::cout<<std::setprecision(3) << lat<<"\t\t"<<bytes<<"\t\t"<<bytes/time<<"\t\t"<<effbytes/time<<std::endl;

    }

  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "= Benchmarking SCALE bandwidth"<<std::endl;
  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "  L  "<<"\t\t"<<"bytes"<<"\t\t"<<"GB/s"<<std::endl;

  for(int lat=4;lat<=24;lat+=4){

      std::vector<int> latt_size  ({lat,lat,lat,lat});

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      //GridParallelRNG          pRNG(&Grid);      pRNG.SeedRandomDevice();

      LatticeVec z(&Grid); //random(pRNG,z);
      LatticeVec x(&Grid); //random(pRNG,x);
      LatticeVec y(&Grid); //random(pRNG,y);
      RealD a=2.0;


      double start=usecond();
      for(int i=0;i<Nloop;i++){
	x=a*z;
      }
      double stop=usecond();
      double time = (stop-start)/Nloop/1000;
      
      double bytes=2*lat*lat*lat*lat*Nvec*sizeof(Real);
      std::cout <<std::setprecision(3) << lat<<"\t\t"<<bytes<<"\t\t"<<bytes/time<<std::endl;

  }

  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "= Benchmarking READ bandwidth"<<std::endl;
  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "  L  "<<"\t\t"<<"bytes"<<"\t\t"<<"GB/s"<<std::endl;
  std::cout << "----------------------------------------------------------"<<std::endl;

  for(int lat=4;lat<=24;lat+=4){

      std::vector<int> latt_size  ({lat,lat,lat,lat});

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      //GridParallelRNG          pRNG(&Grid);      pRNG.SeedRandomDevice();

      LatticeVec z(&Grid); //random(pRNG,z);
      LatticeVec x(&Grid); //random(pRNG,x);
      LatticeVec y(&Grid); //random(pRNG,y);
      RealD a=2.0;
      ComplexD nn;

      double start=usecond();
      for(int i=0;i<Nloop;i++){
	nn=norm2(x);
      }
      double stop=usecond();
      double time = (stop-start)/Nloop/1000;
      
      double bytes=lat*lat*lat*lat*Nvec*sizeof(Real);
      std::cout<<std::setprecision(3) << lat<<"\t\t"<<bytes<<"\t\t"<<bytes/time<<std::endl;

  }    

  Grid_finalize();
}
