#pragma once

#if defined(GRID_CUDA)

#include <cub/cub.cuh>
#define gpucub cub
#define gpuError_t cudaError_t
#define gpuSuccess cudaSuccess

#elif defined(GRID_HIP)

#include <hipcub/hipcub.hpp>
#define gpucub hipcub
#define gpuError_t hipError_t
#define gpuSuccess hipSuccess

#endif


NAMESPACE_BEGIN(Grid);


#if defined(GRID_CUDA) || defined(GRID_HIP)
template<class vobj> inline void sliceSumReduction_cub_small(const vobj *Data, Vector<vobj> &lvSum, const int rd, const int e1, const int e2, const int stride, const int ostride, const int Nsimd) {
  size_t subvol_size = e1*e2;
  commVector<vobj> reduction_buffer(rd*subvol_size);
  auto rb_p = &reduction_buffer[0];
  vobj zero_init;
  zeroit(zero_init);

  
  void *temp_storage_array = NULL;
  size_t temp_storage_bytes = 0;
  vobj *d_out;
  int* d_offsets;

  std::vector<int> offsets(rd+1,0);

  for (int i = 0; i < offsets.size(); i++) {
    offsets[i] = i*subvol_size;
  }
  
  //Allocate memory for output and offset arrays on device
  d_out = static_cast<vobj*>(acceleratorAllocDevice(rd*sizeof(vobj)));
  
  d_offsets = static_cast<int*>(acceleratorAllocDevice((rd+1)*sizeof(int)));
  
  //copy offsets to device
  acceleratorCopyToDeviceAsync(&offsets[0],d_offsets,sizeof(int)*(rd+1),computeStream);
  
  
  gpuError_t gpuErr = gpucub::DeviceSegmentedReduce::Reduce(temp_storage_array, temp_storage_bytes, rb_p,d_out, rd, d_offsets, d_offsets+1, ::gpucub::Sum(), zero_init, computeStream);
  if (gpuErr!=gpuSuccess) {
    std::cout << GridLogError << "Lattice_slicesum_gpu.h: Encountered error during gpucub::DeviceSegmentedReduce::Reduce (setup)! Error: " << gpuErr <<std::endl;
    exit(EXIT_FAILURE);
  }

  //allocate memory for temp_storage_array  
  temp_storage_array = acceleratorAllocDevice(temp_storage_bytes);
  
  //prepare buffer for reduction
  //use non-blocking accelerator_for to avoid syncs (ok because we submit to same computeStream)
  //use 2d accelerator_for to avoid launch latencies found when serially looping over rd 
  accelerator_for2dNB( s,subvol_size, r,rd, Nsimd,{ 
  
    int n = s / e2;
    int b = s % e2;
    int so=r*ostride; // base offset for start of plane 
    int ss= so+n*stride+b;

    coalescedWrite(rb_p[r*subvol_size+s], coalescedRead(Data[ss]));

  });
  
  //issue segmented reductions in computeStream
  gpuErr = gpucub::DeviceSegmentedReduce::Reduce(temp_storage_array, temp_storage_bytes, rb_p, d_out, rd, d_offsets, d_offsets+1,::gpucub::Sum(), zero_init, computeStream);
  if (gpuErr!=gpuSuccess) {
    std::cout << GridLogError << "Lattice_slicesum_gpu.h: Encountered error during gpucub::DeviceSegmentedReduce::Reduce! Error: " << gpuErr <<std::endl;
    exit(EXIT_FAILURE);
  }
  
  acceleratorCopyFromDeviceAsync(d_out,&lvSum[0],rd*sizeof(vobj),computeStream);
  
  //sync after copy
  accelerator_barrier();
 
  acceleratorFreeDevice(temp_storage_array);
  acceleratorFreeDevice(d_out);
  acceleratorFreeDevice(d_offsets);
  

}
#endif 


#if defined(GRID_SYCL)
template<class vobj> inline void sliceSumReduction_sycl_small(const vobj *Data, Vector <vobj> &lvSum, const int  &rd, const int &e1, const int &e2, const int &stride, const int &ostride, const int &Nsimd)
{
  size_t subvol_size = e1*e2;

  vobj *mysum = (vobj *) malloc_shared(rd*sizeof(vobj),*theGridAccelerator);
  vobj vobj_zero;
  zeroit(vobj_zero);
  for (int r = 0; r<rd; r++) { 
    mysum[r] = vobj_zero; 
  }

  commVector<vobj> reduction_buffer(rd*subvol_size);    

  auto rb_p = &reduction_buffer[0];

  // autoView(Data_v, Data, AcceleratorRead);

  //prepare reduction buffer 
  accelerator_for2d( s,subvol_size, r,rd, (size_t)Nsimd,{ 
  
      int n = s / e2;
      int b = s % e2;
      int so=r*ostride; // base offset for start of plane 
      int ss= so+n*stride+b;

      coalescedWrite(rb_p[r*subvol_size+s], coalescedRead(Data[ss]));

  });

  for (int r = 0; r < rd; r++) {
      theGridAccelerator->submit([&](cl::sycl::handler &cgh) {
          auto Reduction = cl::sycl::reduction(&mysum[r],std::plus<>());
          cgh.parallel_for(cl::sycl::range<1>{subvol_size},
          Reduction,
          [=](cl::sycl::id<1> item, auto &sum) {
              auto s = item[0];
              sum += rb_p[r*subvol_size+s];
          });
      });
      
     
  }
  theGridAccelerator->wait();
  for (int r = 0; r < rd; r++) {
    lvSum[r] = mysum[r];
  }
  free(mysum,*theGridAccelerator);
}
#endif

template<class vobj> inline void sliceSumReduction_large(const vobj *Data, Vector<vobj> &lvSum, const int rd, const int e1, const int e2, const int stride, const int ostride, const int Nsimd) {
  typedef typename vobj::vector_type vector;
  const int words = sizeof(vobj)/sizeof(vector);
  const int osites = rd*e1*e2;
  commVector<vector>buffer(osites);
  vector *dat = (vector *)Data;
  vector *buf = &buffer[0];
  Vector<vector> lvSum_small(rd);
  vector *lvSum_ptr = (vector *)&lvSum[0];

  for (int w = 0; w < words; w++) {
    accelerator_for(ss,osites,1,{
	    buf[ss] = dat[ss*words+w];
    });

    #if defined(GRID_CUDA) || defined(GRID_HIP)
      sliceSumReduction_cub_small(buf,lvSum_small,rd,e1,e2,stride, ostride,Nsimd);
    #elif defined(GRID_SYCL)
      sliceSumReduction_sycl_small(buf,lvSum_small,rd,e1,e2,stride, ostride,Nsimd);
    #endif

    for (int r = 0; r < rd; r++) {
      lvSum_ptr[w+words*r]=lvSum_small[r];
    }

  }

  
}

template<class vobj> inline void sliceSumReduction_gpu(const Lattice<vobj> &Data, Vector<vobj> &lvSum, const int rd, const int e1, const int e2, const int stride, const int ostride, const int Nsimd)
{
  autoView(Data_v, Data, AcceleratorRead); //reduction libraries cannot deal with large vobjs so we split into small/large case.
    if constexpr (sizeof(vobj) <= 256) { 

      #if defined(GRID_CUDA) || defined(GRID_HIP)
        sliceSumReduction_cub_small(&Data_v[0], lvSum, rd, e1, e2, stride, ostride, Nsimd);
      #elif defined (GRID_SYCL)
        sliceSumReduction_sycl_small(&Data_v[0], lvSum, rd, e1, e2, stride, ostride, Nsimd);
      #endif

    }
    else {
      sliceSumReduction_large(&Data_v[0], lvSum, rd, e1, e2, stride, ostride, Nsimd);
    }
}


template<class vobj> inline void sliceSumReduction_cpu(const Lattice<vobj> &Data, Vector<vobj> &lvSum, const int &rd, const int &e1, const int &e2, const int &stride, const int &ostride, const int &Nsimd)
{
  // sum over reduced dimension planes, breaking out orthog dir
  // Parallel over orthog direction
  autoView( Data_v, Data, CpuRead);
  thread_for( r,rd, {
    int so=r*ostride; // base offset for start of plane 
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
        int ss= so+n*stride+b;
        lvSum[r]=lvSum[r]+Data_v[ss];
      }
    }
  });
}

template<class vobj> inline void sliceSumReduction(const Lattice<vobj> &Data, Vector<vobj> &lvSum, const int &rd, const int &e1, const int &e2, const int &stride, const int &ostride, const int &Nsimd) 
{
  #if defined(GRID_CUDA) || defined(GRID_HIP) || defined(GRID_SYCL)
  
  sliceSumReduction_gpu(Data, lvSum, rd, e1, e2, stride, ostride, Nsimd);
  
  #else
  sliceSumReduction_cpu(Data, lvSum, rd, e1, e2, stride, ostride, Nsimd);

  #endif
}


NAMESPACE_END(Grid);
