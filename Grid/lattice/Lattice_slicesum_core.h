#pragma once
#if defined(GRID_CUDA)

#include <cub/cub.cuh>
#define gpucub cub
#define gpuMalloc cudaMalloc
#define gpuMemcpyAsync cudaMemcpyAsync
#define gpuMemcpyDeviceToHost cudaMemcpyDeviceToHost
#define gpuMemcpyHostToDevice cudaMemcpyHostToDevice
#define gpuError_t cudaError_t
#define gpuSuccess cudaSuccess

#elif defined(GRID_HIP)

#include <hipcub/hipcub.hpp>
#define gpucub hipcub
#define gpuMalloc hipMalloc
#define gpuMemcpyAsync hipMemcpyAsync
#define gpuMemcpyDeviceToHost hipMemcpyDeviceToHost
#define gpuMemcpyHostToDevice hipMemcpyHostToDevice
#define gpuError_t hipError_t
#define gpuSuccess hipSuccess

#endif


NAMESPACE_BEGIN(Grid);

#if defined(GRID_CUDA) || defined(GRID_HIP)
template<class vobj> inline void sliceSumReduction_cub(const Lattice<vobj> &Data, Vector<vobj> &lvSum, const int rd, const int e1, const int e2, const int stride, const int ostride, const int Nsimd)
{
  typedef typename vobj::scalar_object sobj;

  size_t subvol_size = e1*e2;

  commVector<vobj> reduction_buffer(rd*subvol_size);
  auto rb_p = &reduction_buffer[0];

  vobj vobj_zero; //Need to provide initial value for reduction operation
  zeroit(vobj_zero);
  

  void *temp_storage_array = NULL;
  size_t temp_storage_bytes = 0;
  vobj *d_out;
  int* d_offsets;

  std::vector<int> offsets(rd+1,0);

  for (int i = 0; i < offsets.size(); i++) {
    offsets[i] = i*subvol_size;
  }
  
  //Allocate memory for output and offset arrays on device
  gpuError_t gpuErr = gpuMalloc(&d_out,rd*sizeof(vobj));
  if (gpuErr != gpuSuccess) {
    std::cout << GridLogError << "Lattice_slicesum_gpu.h: Encountered error during gpuMalloc (d_out)! Error: " << gpuErr <<std::endl;
    exit(EXIT_FAILURE);
  }

  gpuErr = gpuMalloc(&d_offsets,sizeof(int)*(rd+1));
  if (gpuErr != gpuSuccess) {
    std::cout << GridLogError << "Lattice_slicesum_gpu.h: Encountered error during gpuMalloc (d_offsets)! Error: " << gpuErr <<std::endl;
    exit(EXIT_FAILURE);
  }

  //copy offsets to device
  gpuErr = gpuMemcpyAsync(d_offsets,&offsets[0],sizeof(int)*(rd+1),gpuMemcpyHostToDevice,computeStream);
  if (gpuErr != gpuSuccess) {
    std::cout << GridLogError << "Lattice_slicesum_gpu.h: Encountered error during gpuMemcpy (d_offsets)! Error: " << gpuErr <<std::endl;
    exit(EXIT_FAILURE);
  }

  //determine temp_storage_array size
  gpuErr = gpucub::DeviceSegmentedReduce::Reduce(temp_storage_array, temp_storage_bytes, rb_p,d_out, rd, d_offsets, d_offsets+1, ::gpucub::Sum(), vobj_zero, computeStream);
  if (gpuErr!=gpuSuccess) {
    std::cout << GridLogError << "Lattice_slicesum_gpu.h: Encountered error during gpucub::DeviceSegmentedReduce::Reduce (setup)! Error: " << gpuErr <<std::endl;
    exit(EXIT_FAILURE);
  }

  //allocate memory for temp_storage_array  
  gpuErr = gpuMalloc(&temp_storage_array,temp_storage_bytes);
  if (gpuErr!=gpuSuccess) {
    std::cout << GridLogError << "Lattice_slicesum_gpu.h: Encountered error during gpuMalloc (temp_storage_array)! Error: " << gpuErr <<std::endl;
    exit(EXIT_FAILURE);
  }
  
  autoView( Data_v, Data, AcceleratorRead);
  //prepare buffer for reduction
  //use non-blocking accelerator_for to avoid syncs (ok because we submit to same computeStream)
  //use 2d accelerator_for to avoid launch latencies found when serially looping over rd 
  
  accelerator_for2dNB( s,subvol_size, r,rd, Nsimd,{ 
  
    int n = s / e2;
    int b = s % e2;
    int so=r*ostride; // base offset for start of plane 
    int ss= so+n*stride+b;

    coalescedWrite(rb_p[r*subvol_size+s], coalescedRead(Data_v[ss]));

  });
  
  //issue segmented reductions in computeStream
  gpuErr = gpucub::DeviceSegmentedReduce::Reduce(temp_storage_array, temp_storage_bytes, rb_p, d_out, rd, d_offsets, d_offsets+1,::gpucub::Sum(), vobj_zero, computeStream);
  if (gpuErr!=gpuSuccess) {
    std::cout << GridLogError << "Lattice_slicesum_gpu.h: Encountered error during gpucub::DeviceSegmentedReduce::Reduce! Error: " << gpuErr <<std::endl;
    exit(EXIT_FAILURE);
  }
  
  gpuErr = gpuMemcpyAsync(&lvSum[0],d_out,rd*sizeof(vobj),gpuMemcpyDeviceToHost,computeStream);
  if (gpuErr!=gpuSuccess) {
    std::cout << GridLogError << "Lattice_slicesum_gpu.h: Encountered error during gpuMemcpy (d_out)! Error: " << gpuErr <<std::endl;
    exit(EXIT_FAILURE);
  }

  //sync after copy
  accelerator_barrier();
 

}
#endif


#if defined(GRID_SYCL)
template<class vobj> inline void sliceSumReduction_sycl(const Lattice<vobj> &Data, Vector <vobj> &lvSum, const int  &rd, const int &e1, const int &e2, const int &stride, const int &ostride, const int &Nsimd)
{
  typedef typename vobj::scalar_object sobj;
  size_t subvol_size = e1*e2;

  vobj *mysum = (vobj *) malloc_shared(sizeof(vobj),*theGridAccelerator);
  vobj vobj_zero;
  zeroit(vobj_zero);
    
  commVector<vobj> reduction_buffer(rd*subvol_size);    

  auto rb_p = &reduction_buffer[0];

  autoView(Data_v, Data, AcceleratorRead);

  //prepare reduction buffer 
  accelerator_for2d( s,subvol_size, r,rd, (size_t)Nsimd,{ 
  
      int n = s / e2;
      int b = s % e2;
      int so=r*ostride; // base offset for start of plane 
      int ss= so+n*stride+b;

      coalescedWrite(rb_p[r*subvol_size+s], coalescedRead(Data_v[ss]));

  });

  for (int r = 0; r < rd; r++) {
      mysum[0] = vobj_zero; //dirty hack: cannot pass vobj_zero as identity to sycl::reduction as its not device_copyable
      theGridAccelerator->submit([&](cl::sycl::handler &cgh) {
          auto Reduction = cl::sycl::reduction(mysum,std::plus<>());
          cgh.parallel_for(cl::sycl::range<1>{subvol_size},
          Reduction,
          [=](cl::sycl::id<1> item, auto &sum) {
              auto s = item[0];
              sum += rb_p[r*subvol_size+s];
          });
      });
      theGridAccelerator->wait();
      lvSum[r] = mysum[0];
  }

}
#endif

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
  #if defined(GRID_CUDA) || defined(GRID_HIP)
  
  sliceSumReduction_cub(Data, lvSum, rd, e1, e2, stride, ostride, Nsimd);
  
  #elif defined(GRID_SYCL)
  
  sliceSumReduction_sycl(Data, lvSum, rd, e1, e2, stride, ostride, Nsimd);
  
  #else
  sliceSumReduction_cpu(Data, lvSum, rd, e1, e2, stride, ostride, Nsimd);

  #endif
}


NAMESPACE_END(Grid);