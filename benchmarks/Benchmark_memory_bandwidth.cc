#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Nvec=8;
  typedef Lattice< iVector< vReal,Nvec> > LatticeVec;

  int Nloop=1000;

  std::vector<int> simd_layout = GridDefaultSimd(Nd,vReal::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();

  int threads = GridThread::GetThreads();
  std::cout << "Grid is setup to use "<<threads<<" threads"<<std::endl;
  
  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "= Benchmarking fused AXPY bandwidth ; sizeof(Real) "<<sizeof(Real)<<std::endl;
  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s"<<"\t\t"<<"Gflop/s"<<std::endl;
  std::cout << "----------------------------------------------------------"<<std::endl;

  for(int lat=4;lat<=32;lat+=4){

      std::vector<int> latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      //GridParallelRNG          pRNG(&Grid);      pRNG.SeedRandomDevice();

      LatticeVec z(&Grid); //random(pRNG,z);
      LatticeVec x(&Grid); //random(pRNG,x);
      LatticeVec y(&Grid); //random(pRNG,y);
      double a=2.0;


      double start=usecond();
      for(int i=0;i<Nloop;i++){
	//   inline void axpy(Lattice<vobj> &ret,double a,const Lattice<vobj> &lhs,const Lattice<vobj> &rhs){
	axpy(z,a,x,y);
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000;
      
      double flops=vol*Nvec*2;// mul,add
      double bytes=3*vol*Nvec*sizeof(Real);
      std::cout<<std::setprecision(3) << lat<<"\t\t"<<bytes<<"   \t\t"<<bytes/time<<"\t\t"<<flops/time<<std::endl;

    }

  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "= Benchmarking a*x + y bandwidth"<<std::endl;
  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s"<<"\t\t"<<"Gflop/s"<<std::endl;
  std::cout << "----------------------------------------------------------"<<std::endl;
  
  for(int lat=4;lat<=32;lat+=4){

      std::vector<int> latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      //GridParallelRNG          pRNG(&Grid);      pRNG.SeedRandomDevice();

      LatticeVec z(&Grid); //random(pRNG,z);
      LatticeVec x(&Grid); //random(pRNG,x);
      LatticeVec y(&Grid); //random(pRNG,y);
      double a=2.0;

      double start=usecond();
      for(int i=0;i<Nloop;i++){
	z=a*x-y;
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000;
     
      double flops=vol*Nvec*2;// mul,add
      double bytes=3*vol*Nvec*sizeof(Real);
      std::cout<<std::setprecision(3) << lat<<"\t\t"<<bytes<<"   \t\t"<<bytes/time<<"\t\t"<<flops/time<<std::endl;

    }

  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "= Benchmarking SCALE bandwidth"<<std::endl;
  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s"<<"\t\t"<<"Gflop/s"<<std::endl;

  for(int lat=4;lat<=32;lat+=4){

      std::vector<int> latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      //GridParallelRNG          pRNG(&Grid);      pRNG.SeedRandomDevice();

      LatticeVec z(&Grid); //random(pRNG,z);
      LatticeVec x(&Grid); //random(pRNG,x);
      LatticeVec y(&Grid); //random(pRNG,y);
      RealD a=2.0;


      double start=usecond();
      for(int i=0;i<Nloop;i++){
	z=a*x;
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000;
      
      double bytes=2*vol*Nvec*sizeof(Real);
      double flops=vol*Nvec*1;// mul
      std::cout <<std::setprecision(3) << lat<<"\t\t"<<bytes<<"   \t\t"<<bytes/time<<"\t\t"<<flops/time<<std::endl;

  }

  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "= Benchmarking READ bandwidth"<<std::endl;
  std::cout << "===================================================================================================="<<std::endl;
  std::cout << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s"<<"\t\t"<<"Gflop/s"<<std::endl;
  std::cout << "----------------------------------------------------------"<<std::endl;

  for(int lat=4;lat<=32;lat+=4){

      std::vector<int> latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int vol = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
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
      double time = (stop-start)/Nloop*1000;
      
      double bytes=vol*Nvec*sizeof(Real);
      double flops=vol*Nvec*2;// mul,add
      std::cout<<std::setprecision(3) << lat<<"\t\t"<<bytes<<"  \t\t"<<bytes/time<<"\t\t"<<flops/time<<std::endl;

  }    

  Grid_finalize();
}
