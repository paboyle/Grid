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

template<class vobj> inline void sliceSumGpu(const Lattice<vobj> &Data,std::vector<typename vobj::scalar_object> &result,int orthogdim)
{
 
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

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];
  int ostride=grid->_ostride[orthogdim];
  size_t subvol_size = e1*e2;

  Vector<vobj> lvSum(rd); 
  Vector<sobj> lsSum(ld,Zero());                    
  commVector<vobj> reduction_buffer(rd*e1*e2);
  ExtractBuffer<sobj> extracted(Nsimd);                  
  
  result.resize(fd);

  for(int r=0;r<rd;r++){
    lvSum[r]=Zero();
  }

  vobj vobj_zero; //Need to provide initial value for reduction operation
  zeroit(vobj_zero);
  
  autoView( Data_v, Data, AcceleratorRead);

  auto rb_p = &reduction_buffer[0];
  void *helperArray = NULL;
  vobj *d_out;
  size_t temp_storage_bytes = 0;
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

  //determine helperArray size
  gpuErr = gpucub::DeviceSegmentedReduce::Reduce(helperArray, temp_storage_bytes, rb_p,d_out, rd, d_offsets, d_offsets+1, ::gpucub::Sum(), vobj_zero, computeStream);
  if (gpuErr!=gpuSuccess) {
    std::cout << GridLogError << "Lattice_slicesum_gpu.h: Encountered error during gpucub::DeviceSegmentedReduce::Reduce (setup)! Error: " << gpuErr <<std::endl;
    exit(EXIT_FAILURE);
  }

  //allocate memory for helperArray  
  gpuErr = gpuMalloc(&helperArray,temp_storage_bytes);
  if (gpuErr!=gpuSuccess) {
    std::cout << GridLogError << "Lattice_slicesum_gpu.h: Encountered error during gpuMalloc (helperArray)! Error: " << gpuErr <<std::endl;
    exit(EXIT_FAILURE);
  }

  //prepare buffer for reduction
  //use non-blocking accelerator_for to avoid syncs (ok because we submit to same computeStream)
  //use 2d accelerator_for to avoid launch latencies found when serially looping over rd 
  
  accelerator_for2dNB( s,subvol_size, r,rd, grid->Nsimd(),{ 
  
    int n = s / e2;
    int b = s % e2;
    int so=r*ostride; // base offset for start of plane 
    int ss= so+n*stride+b;

    coalescedWrite(rb_p[r*subvol_size+s], coalescedRead(Data_v[ss]));

  });
  
  //issue segmented reductions in computeStream
  gpuErr = gpucub::DeviceSegmentedReduce::Reduce(helperArray, temp_storage_bytes, rb_p, d_out, rd, d_offsets, d_offsets+1,::gpucub::Sum(), vobj_zero, computeStream);
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

NAMESPACE_END(Grid);