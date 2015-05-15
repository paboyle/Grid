#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int Nloop=1000;

  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();

  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "= Benchmarking SU3xSU3  x= x*y"<<std::endl;
  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "  L  "<<"\t\t"<<"bytes"<<"\t\t"<<"GB/s\t\t GFlop/s"<<std::endl;
  std::cout << "----------------------------------------------------------"<<std::endl;

  for(int lat=2;lat<=24;lat+=2){

      std::vector<int> latt_size  ({lat,lat,lat,lat});

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      //      GridParallelRNG          pRNG(&Grid);      pRNG.SeedRandomDevice();

      LatticeColourMatrix z(&Grid);// random(pRNG,z);
      LatticeColourMatrix x(&Grid);// random(pRNG,x);
      LatticeColourMatrix y(&Grid);// random(pRNG,y);

      double start=usecond();
      for(int i=0;i<Nloop;i++){
	x=x*y;
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000.0;
      
      double bytes=3.0*lat*lat*lat*lat*Nc*Nc*sizeof(Complex);
      double footprint=2.0*lat*lat*lat*lat*Nc*Nc*sizeof(Complex);
      double flops=Nc*Nc*(6.0+8.0+8.0)*lat*lat*lat*lat;
      std::cout<<std::setprecision(2) << lat<<"\t\t"<<footprint<<"\t\t"<<bytes/time<<"\t\t" << flops/time<<std::endl;

    }


  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "= Benchmarking SU3xSU3  z= x*y"<<std::endl;
  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "  L  "<<"\t\t"<<"bytes"<<"\t\t"<<"GB/s\t\t GFlop/s"<<std::endl;
  std::cout << "----------------------------------------------------------"<<std::endl;

  for(int lat=2;lat<=24;lat+=2){

      std::vector<int> latt_size  ({lat,lat,lat,lat});

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      //      GridParallelRNG          pRNG(&Grid);      pRNG.SeedRandomDevice();

      LatticeColourMatrix z(&Grid); //random(pRNG,z);
      LatticeColourMatrix x(&Grid); //random(pRNG,x);
      LatticeColourMatrix y(&Grid); //random(pRNG,y);

      double start=usecond();
      for(int i=0;i<Nloop;i++){
	z=x*y;
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000.0;
      
      double bytes=3*lat*lat*lat*lat*Nc*Nc*sizeof(Complex);
      double flops=Nc*Nc*(6+8+8)*lat*lat*lat*lat;
      std::cout<<std::setprecision(2) << lat<<"\t\t"<<bytes<<"\t\t"<<bytes/time<<"\t\t" << flops/time<<std::endl;

    }

  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "= Benchmarking SU3xSU3  mult(z,x,y)"<<std::endl;
  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "  L  "<<"\t\t"<<"bytes"<<"\t\t"<<"GB/s\t\t GFlop/s"<<std::endl;
  std::cout << "----------------------------------------------------------"<<std::endl;

  for(int lat=2;lat<=24;lat+=2){

      std::vector<int> latt_size  ({lat,lat,lat,lat});

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      //      GridParallelRNG          pRNG(&Grid);      pRNG.SeedRandomDevice();

      LatticeColourMatrix z(&Grid); //random(pRNG,z);
      LatticeColourMatrix x(&Grid); //random(pRNG,x);
      LatticeColourMatrix y(&Grid); //random(pRNG,y);

      double start=usecond();
      for(int i=0;i<Nloop;i++){
	mult(z,x,y);
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000.0;
      
      double bytes=3*lat*lat*lat*lat*Nc*Nc*sizeof(Complex);
      double flops=Nc*Nc*(6+8+8)*lat*lat*lat*lat;
      std::cout<<std::setprecision(2) << lat<<"\t\t"<<bytes<<"\t\t"<<bytes/time<<"\t\t" << flops/time<<std::endl;

    }

  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "= Benchmarking SU3xSU3  mac(z,x,y)"<<std::endl;
  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "  L  "<<"\t\t"<<"bytes"<<"\t\t"<<"GB/s\t\t GFlop/s"<<std::endl;
  std::cout << "----------------------------------------------------------"<<std::endl;

  for(int lat=2;lat<=24;lat+=2){

      std::vector<int> latt_size  ({lat,lat,lat,lat});

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
      //      GridParallelRNG          pRNG(&Grid);      pRNG.SeedRandomDevice();

      LatticeColourMatrix z(&Grid); //random(pRNG,z);
      LatticeColourMatrix x(&Grid); //random(pRNG,x);
      LatticeColourMatrix y(&Grid); //random(pRNG,y);

      double start=usecond();
      for(int i=0;i<Nloop;i++){
	mac(z,x,y);
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000.0;
      
      double bytes=3*lat*lat*lat*lat*Nc*Nc*sizeof(Complex);
      double flops=Nc*Nc*(8+8+8)*lat*lat*lat*lat;
      std::cout<<std::setprecision(2) << lat<<"\t\t"<<bytes<<"\t\t"<<bytes/time<<"\t\t" << flops/time<<std::endl;

    }

  Grid_finalize();
}
