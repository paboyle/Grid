#pragma once
#if defined(GRID_CUDA)

#include <cub/cub.cuh>
#define gpucub cub
#define gpuMalloc cudaMalloc
#define gpuMemcpy cudaMemcpy
#define gpuMemcpyDeviceToHost cudaMemcpyDeviceToHost
#define gpuMemcpyHostToDevice cudaMemcpyHostToDevice
#define gpuError_t cudaError_t
#define gpuSuccess cudaSuccess

#elif defined(GRID_HIP)

#include <hipcub/hipcub.hpp>
#define gpucub hipcub
#define gpuMalloc hipMalloc
#define gpuMemcpy hipMemcpy
#define gpuMemcpyDeviceToHost hipMemcpyDeviceToHost
#define gpuMemcpyHostToDevice hipMemcpyHostToDevice
#define gpuError_t hipError_t
#define gpuSuccess hipSuccess

// extern hipStream_t computeStream;
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
  Vector<vobj> lvSum(rd); 
  Vector<sobj> lsSum(ld,Zero());                    
  commVector<vobj> reduction_buffer(rd*e1*e2);
  ExtractBuffer<sobj> extracted(Nsimd);                  
  
  result.resize(fd);
  for(int r=0;r<rd;r++){
    lvSum[r]=Zero();
  }
  vobj identity;
  zeroit(identity);
  
  autoView( Data_v, Data, AcceleratorRead);
  auto rb_p = &reduction_buffer[0];

  void *helperArray = NULL;
  vobj *d_out;
  size_t temp_storage_bytes = 0;
  size_t size = e1*e2;
  std::vector<int> offsets(rd+1,0);
  for (int i = 0; i < offsets.size(); i++) {
    offsets[i] = i*size;
  }
  int* d_offsets;

  gpuError_t gpuErr = gpuMalloc(&d_out,rd*sizeof(vobj));
   if (gpuErr != gpuSuccess) {
    std::cout << "Lattice_slicesum_gpu.h: Encountered error during gpuMalloc(1) Error: " << gpuErr <<std::endl;
  }
  gpuErr = gpuMalloc(&d_offsets,sizeof(int)*(rd+1));
   if (gpuErr != gpuSuccess) {
    std::cout << "Lattice_slicesum_gpu.h: Encountered error during gpuMalloc(2) Error: " << gpuErr <<std::endl;
  }
  gpuErr = gpuMemcpy(d_offsets,&offsets[0],sizeof(int)*(rd+1),gpuMemcpyHostToDevice);
  if (gpuErr != gpuSuccess) {
    std::cout << "Lattice_slicesum_gpu.h: Encountered error during gpuMemcpy(1) Error: " << gpuErr <<std::endl;
  }

  gpuErr = gpucub::DeviceSegmentedReduce::Reduce(helperArray, temp_storage_bytes, rb_p,d_out, rd, d_offsets, d_offsets+1, ::gpucub::Sum(), identity, computeStream);
  if (gpuErr!=gpuSuccess) {
    std::cout << "Lattice_slicesum_gpu.h: Encountered error during cub::DeviceReduce::Sum(1)! Error: " << gpuErr <<std::endl;
  }
  
  gpuErr = gpuMalloc(&helperArray,temp_storage_bytes);
  if (gpuErr!=gpuSuccess) {
    std::cout << "Lattice_slicesum_gpu.h: Encountered error during gpuMalloc Error: " << gpuErr <<std::endl;
  }

  //prepare buffer for reduction
  accelerator_for2dNB( s,e1*e2, r,rd, grid->Nsimd(),{ //use non-blocking accelerator_for to avoid syncs (ok because we submit to same computeStream)
                                                      //use 2d accelerator_for to avoid launch latencies found when looping over rd 
    int n = s / e2;
    int b = s % e2;
    int so=r*ostride; // base offset for start of plane 
    int ss= so+n*stride+b;

    coalescedWrite(rb_p[r*e1*e2+s], coalescedRead(Data_v[ss]));

  });
  
  //issue reductions in computeStream
  gpuErr =gpucub::DeviceSegmentedReduce::Reduce(helperArray, temp_storage_bytes, rb_p, d_out, rd, d_offsets, d_offsets+1,::gpucub::Sum(), identity, computeStream);
  if (gpuErr!=gpuSuccess) {
    std::cout << "Lattice_slicesum_gpu.h: Encountered error during cub::DeviceReduce::Sum(2)! Error: " << gpuErr <<std::endl;
  }
  
  //sync before copy
  accelerator_barrier();
  gpuErr = gpuMemcpy(&lvSum[0],d_out,rd*sizeof(vobj),gpuMemcpyDeviceToHost);
  if (gpuErr!=gpuSuccess) {
    std::cout << "Lattice_slicesum_gpu.h: Encountered error during gpuMemcpy(2) Error: " << gpuErr <<std::endl;
  }
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

template<class vobj> inline
std::vector<typename vobj::scalar_object> 
sliceSumGpu(const Lattice<vobj> &Data,int orthogdim)
{
  std::vector<typename vobj::scalar_object> result;
  sliceSumGpu(Data,result,orthogdim);
  return result;
}

NAMESPACE_END(Grid);