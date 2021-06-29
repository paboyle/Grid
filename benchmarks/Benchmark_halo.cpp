#include <Grid/Grid.h>

#define NCOL 2

using namespace Grid;

constexpr int Ndim = 3;
typedef ScalarAdjMatrixImplTypes<vComplex,NCOL>::Field SUNField;
typedef typename SUNField::vector_object vobj;
typedef CartesianStencil<vobj, vobj,int> Stencil;

int main(int argc, char **argv){

  // Initialise grid //////////////////////////////////////////////
  Grid_init(&argc,&argv);
  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  // Module ///////////////////////////////////////////////////////
  GridModule GridMod;
  if (GridDefaultLatt().size() != Ndim){
    std::cout << GridLogError << "Incorrect dimension of the grid\n. Expected dim=" << Ndim << std::endl;
    return EXIT_FAILURE;
  }
  if (GridDefaultMpi().size() != Ndim){
    std::cout << GridLogError << "Incorrect dimension of the mpi grid\n. Expected dim=" << Ndim << std::endl;
    return EXIT_FAILURE;
  }
  GridMod.set_full(new GridCartesian(GridDefaultLatt(),
  GridDefaultSimd(Ndim, vComplex::Nsimd()),
  GridDefaultMpi()));
  GridMod.set_rb(new GridRedBlackCartesian(GridMod.get_full()));
  auto grid = GridMod.get_full();

  GridParallelRNG pRNG(grid);
  pRNG.SeedFixedIntegers({11,84,79,47,90});

  // Stencil //////////////////////////////////////////////////////
  int npoint = 2 * Ndim;
  std::vector<int> directions(npoint);
  std::vector<int> displacements(npoint);

  for (int mu = 0; mu < Ndim; mu++){
    directions[mu] = mu;
    directions[mu + Ndim] = mu;
    displacements[mu] = 1;
    displacements[mu + Ndim] = -1;
  }

  Stencil Stencil_phi(grid, npoint, 0, directions, displacements,0);
  SimpleCompressor<vobj> compressor;

  // Field /////////////////////////////////////////////////////////
  SUNField phi(grid);

  // MPI sublattice surface area ///////////////////////////////////

  int mpi_area = 0;
  int mpi_face;
  // Calculates the total surface area of an MPI hypercube
  for (int mu_ex=0;mu_ex<Ndim;++mu_ex){
    mpi_face = 1;

    for (int mu=0; mu<Ndim; ++mu){
        if (mu != mu_ex) mpi_face *= GridDefaultLatt()[mu]/GridDefaultMpi()[mu];
    }

    mpi_area += 2*mpi_face;
  }

  std::cout << GridLogMessage << "Total MPI surface area = " << mpi_area << std::endl;

  // Benchmarking //////////////////////////////////////////////////

  int nloops = 100;
  double start;
  double time;
  double avgtime = 0;
  double bytes = sizeof(Complex)*NCOL*NCOL*mpi_area*4.;
  // 4 is for the two reads and writes in receiving and sending data
  // I don't know if I am to consider all data being sent and received in all mpi processes across all gpus 
  double avgbandwidth = 0;

  for (int i=0;i<nloops;++i){

    start = usecond();
    Stencil_phi.HaloExchange(phi, compressor);
    time = usecond();
    std::cout << GridLogMessage << "Exchange " << i << " time (us) = " << time-start << " | " << "Bandwidth (GB/s) = " << (bytes/1e9)/(time/1e6) << std::endl;
    avgtime += time/double(nloops);
    avgbandwidth += (bytes/1e9)/(time/1e6)/double(nloops);
    
  }

  std::cout << GridLogMessage << "Average time (us) = " << avgtime << " | Average bandwidth (GB/s) = " << avgbandwidth << std::endl;

  Grid_finalize();

}