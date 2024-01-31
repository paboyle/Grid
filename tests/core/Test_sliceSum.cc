#include <Grid/Grid.h>


int main (int argc, char ** argv) {
    
    using namespace Grid;

    Grid_init(&argc,&argv);


    Coordinate latt_size({64,64,64,16});
    auto simd_layout = GridDefaultSimd(Nd, vComplexD::Nsimd());
    auto mpi_layout = GridDefaultMpi();
    GridCartesian Grid(latt_size, simd_layout, mpi_layout);

    std::vector<int> seeds({1, 2, 3, 4});

    GridParallelRNG pRNG(&Grid);
    pRNG.SeedFixedIntegers(seeds);

    LatticeComplexD test_data(&Grid);
    gaussian(pRNG,test_data);

    std::vector<TComplexD> reduction_reference;
    std::vector<TComplexD> reduction_result;

    //warmup
    for (int sweeps = 0; sweeps < 5; sweeps++) {
      sliceSumGpu(test_data,reduction_result,0);
    }


    for (int i = 0; i < Nd; i++) {
        RealD t=-usecond();
        sliceSum(test_data,reduction_reference,i);
        t+=usecond();
        std::cout << " sliceSum took "<<t<<" usecs"<<std::endl;
        
        RealD tgpu=-usecond();
        tracePush("sliceSumGpu");
        sliceSumGpu(test_data,reduction_result,i);
        tracePop("sliceSumGpu");
        tgpu+=usecond();
        std::cout <<" sliceSumGpu took "<<tgpu<<" usecs"<<std::endl;

    for(int t=0;t<reduction_reference.size();t++){

      auto diff = reduction_reference[t]-reduction_result[t];
      assert(abs(TensorRemove(diff)) < 1e-8 );
    }

    
    }
    Grid_finalize();
    return 0;
}