#include <Grid/Grid.h>

template<class vobj> inline void sliceSumCPU(const Grid::Lattice<vobj> &Data,std::vector<typename vobj::scalar_object> &result,int orthogdim)
{
  using namespace Grid;
  ///////////////////////////////////////////////////////
  // FIXME precision promoted summation
  // may be important for correlation functions
  // But easily avoided by using double precision fields
  ///////////////////////////////////////////////////////
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_object::scalar_type scalar_type;
  GridBase  *grid = Data.Grid();
  assert(grid!=NULL);

  const int    Nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();

  assert(orthogdim >= 0);
  assert(orthogdim < Nd);

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  std::vector<vobj> lvSum(rd); // will locally sum vectors first
  std::vector<sobj> lsSum(ld,Zero());                    // sum across these down to scalars
  ExtractBuffer<sobj> extracted(Nsimd);                  // splitting the SIMD

  result.resize(fd); // And then global sum to return the same vector to every node 
  for(int r=0;r<rd;r++){
    lvSum[r]=Zero();
  }

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];
  int ostride=grid->_ostride[orthogdim];
  
  //Reduce Data down to lvSum
  sliceSumReduction_cpu(Data,lvSum,rd, e1,e2,stride,ostride,Nsimd);

  // Sum across simd lanes in the plane, breaking out orthog dir.
  Coordinate icoor(Nd);

  for(int rt=0;rt<rd;rt++){

    extract(lvSum[rt],extracted);

    for(int idx=0;idx<Nsimd;idx++){

      grid->iCoorFromIindex(icoor,idx);

      int ldx =rt+icoor[orthogdim]*rd;

      lsSum[ldx]=lsSum[ldx]+extracted[idx];

    }
  }
  
  // sum over nodes.
  for(int t=0;t<fd;t++){
    int pt = t/ld; // processor plane
    int lt = t%ld;
    if ( pt == grid->_processor_coor[orthogdim] ) {
      result[t]=lsSum[lt];
    } else {
      result[t]=Zero();
    }

  }
  scalar_type * ptr = (scalar_type *) &result[0];
  int words = fd*sizeof(sobj)/sizeof(scalar_type);
  grid->GlobalSumVector(ptr, words);
}


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
      reduction_result = sliceSum(test_data,0);
    }

    int trace_id = traceStart("sliceSum benchmark - ComplexD");
    std::cout << GridLogMessage << "Testing ComplexD" << std::endl;
    std::cout << GridLogMessage << "sizeof(ComplexD) = " << sizeof(ComplexD) << std::endl;
    std::cout << GridLogMessage << "sizeof(vComplexD) = " << sizeof(vComplexD) << std::endl;
    for (int i = 0; i < Nd; i++) {

      RealD t=-usecond();

      tracePush("sliceSum");
      sliceSumCPU(test_data,reduction_reference,i);
      tracePop("sliceSum");

      t+=usecond();
      std::cout << GridLogMessage << "Orthog. dir. = " << i << std::endl;
      std::cout << GridLogMessage << "CPU sliceSum took "<<t<<" usecs"<<std::endl;
      
      
      RealD tgpu=-usecond();

      tracePush("sliceSumGpu");
      reduction_result = sliceSum(test_data,i);
      tracePop("sliceSumGpu");

      tgpu+=usecond();

      std::cout << GridLogMessage <<"GPU sliceSum took "<<tgpu<<" usecs"<<std::endl<<std::endl;;


      for(int t=0;t<reduction_reference.size();t++) {

        auto diff = reduction_reference[t]-reduction_result[t];
        assert(abs(TensorRemove(diff)) < 1e-8 );

      }

    
    }
    traceStop(trace_id);

    LatticeSpinVectorD test_data_cv(&Grid);
    gaussian(pRNG,test_data_cv);

    std::vector<SpinVectorD> reduction_reference_cv;
    std::vector<SpinVectorD> reduction_result_cv;

    //warmup
    for (int sweeps = 0; sweeps < 5; sweeps++) {
      reduction_result_cv = sliceSum(test_data_cv,0);
    }
    trace_id = traceStart("sliceSum benchmark - SpinVectorD");

    std::cout << GridLogMessage << "Testing SpinVectorD" << std::endl;
    std::cout << GridLogMessage << "sizeof(SpinVectorD) = " << sizeof(SpinVectorD) << std::endl;
    std::cout << GridLogMessage << "sizeof(vSpinVectorD) = " << sizeof(vSpinVectorD) << std::endl;
    for (int i = 0; i < Nd; i++) {

      RealD t=-usecond();

      tracePush("sliceSum");
      sliceSumCPU(test_data_cv,reduction_reference_cv,i);
      tracePop("sliceSum");

      t+=usecond();
      std::cout << GridLogMessage << "Orthog. dir. = " << i << std::endl;
      std::cout << GridLogMessage << "CPU sliceSum took "<<t<<" usecs"<<std::endl;
      
      
      RealD tgpu=-usecond();

      tracePush("sliceSumGpu");
      reduction_result_cv = sliceSum(test_data_cv,i);
      tracePop("sliceSumGpu");

      tgpu+=usecond();

      std::cout << GridLogMessage <<"GPU sliceSum took "<<tgpu<<" usecs"<<std::endl<<std::endl;;


      for(int t=0;t<reduction_reference_cv.size();t++) {

        auto diff = reduction_reference_cv[t]-reduction_result_cv[t];
        assert(abs(diff()(0)()) < 1e-8 );
        assert(abs(diff()(1)()) < 1e-8 );
        assert(abs(diff()(2)()) < 1e-8 );
        assert(abs(diff()(3)()) < 1e-8 );

      }

    
    }
    traceStop(trace_id);

    LatticeSpinColourVectorD test_data_scv(&Grid);
    gaussian(pRNG,test_data_scv);

    std::vector<SpinColourVectorD> reduction_reference_scv;
    std::vector<SpinColourVectorD> reduction_result_scv;

    //warmup
    for (int sweeps = 0; sweeps < 5; sweeps++) {
      reduction_result_scv = sliceSum(test_data_scv,0);
    }
    trace_id = traceStart("sliceSum benchmark - SpinColourVectorD");

    std::cout << GridLogMessage << "Testing SpinColourVectorD" << std::endl;
    std::cout << GridLogMessage << "sizeof(SpinColourVectorD) = " << sizeof(SpinColourVectorD) << std::endl;
    std::cout << GridLogMessage << "sizeof(vSpinColourVectorD) = " << sizeof(vSpinColourVectorD) << std::endl;
    for (int i = 0; i < Nd; i++) {

      RealD t=-usecond();

      tracePush("sliceSum");
      sliceSumCPU(test_data_scv,reduction_reference_scv,i);
      tracePop("sliceSum");

      t+=usecond();
      std::cout << GridLogMessage << "Orthog. dir. = " << i << std::endl;
      std::cout << GridLogMessage << "CPU sliceSum took "<<t<<" usecs"<<std::endl;
      
      
      RealD tgpu=-usecond();

      tracePush("sliceSumGpu");
      reduction_result_scv = sliceSum(test_data_scv,i);
      tracePop("sliceSumGpu");

      tgpu+=usecond();

      std::cout << GridLogMessage <<"GPU sliceSum took "<<tgpu<<" usecs"<<std::endl<<std::endl;;


      for(int t=0;t<reduction_reference_scv.size();t++) {

        auto diff = reduction_reference_scv[t]-reduction_result_scv[t];
        // std::cout << diff <<std::endl;
        assert(abs(diff()(0)(0)) < 1e-8 );
        assert(abs(diff()(0)(1)) < 1e-8 );
        assert(abs(diff()(0)(2)) < 1e-8 );
        assert(abs(diff()(1)(0)) < 1e-8 );
        assert(abs(diff()(1)(1)) < 1e-8 );
        assert(abs(diff()(1)(2)) < 1e-8 );    
        assert(abs(diff()(2)(0)) < 1e-8 );
        assert(abs(diff()(2)(1)) < 1e-8 );
        assert(abs(diff()(2)(2)) < 1e-8 );    
        assert(abs(diff()(3)(0)) < 1e-8 );
        assert(abs(diff()(3)(1)) < 1e-8 );
        assert(abs(diff()(3)(2)) < 1e-8 );

      }

    
    }
    traceStop(trace_id);

    LatticeSpinColourMatrixD test_data_scm(&Grid);
    gaussian(pRNG,test_data_scm);

    std::vector<SpinColourMatrixD> reduction_reference_scm;
    std::vector<SpinColourMatrixD> reduction_result_scm;

    //warmup
    for (int sweeps = 0; sweeps < 5; sweeps++) {
      reduction_result_scm = sliceSum(test_data_scm,0);
    }
    trace_id = traceStart("sliceSum benchmark - SpinColourMatrixD");

    std::cout << GridLogMessage << "Testing SpinColourMatrixD" << std::endl;
    std::cout << GridLogMessage << "sizeof(SpinColourMatrixD) = " << sizeof(SpinColourMatrixD) << std::endl;
    std::cout << GridLogMessage << "sizeof(vSpinColourMatrixD) = " << sizeof(vSpinColourMatrixD) << std::endl;
    for (int i = 0; i < Nd; i++) {

      RealD t=-usecond();

      tracePush("sliceSum");
      sliceSumCPU(test_data_scm,reduction_reference_scm,i);
      tracePop("sliceSum");

      t+=usecond();
      std::cout << GridLogMessage << "Orthog. dir. = " << i << std::endl;
      std::cout << GridLogMessage << "CPU sliceSum took "<<t<<" usecs"<<std::endl;
      
      
      RealD tgpu=-usecond();

      tracePush("sliceSumGpu");
      reduction_result_scm = sliceSum(test_data_scm,i);
      tracePop("sliceSumGpu");

      tgpu+=usecond();

      std::cout << GridLogMessage <<"GPU sliceSum took "<<tgpu<<" usecs"<<std::endl<<std::endl;;


      for(int t=0;t<reduction_reference_scm.size();t++) {

        auto diff = reduction_reference_scm[t]-reduction_result_scm[t];
        // std::cout << diff <<std::endl;
        for (int is = 0; is < Ns; is++) {
          for (int js = 0; js < Ns; js++) {
            for (int ic = 0; ic < Nc; ic++) {
              for (int jc = 0; jc < Nc; jc++) {
                assert(abs(diff()(is,js)(ic,jc)) < 1e-8);
              }
            }
          }
        }

      }

    
    }
    traceStop(trace_id);

    Grid_finalize();
    return 0;
}
