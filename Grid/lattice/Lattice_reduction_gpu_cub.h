/*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid
    Source file: ./Grid/lattice/Lattice_reduction_gpu_cub.h
    Copyright (C) 2015-2024
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Unified lattice reduction using CUB (CUDA/HIP) and sycl::reduction (SYCL).
//
// Strategy: one accelerator_for pass per site to extract SIMD lanes and promote to sobjD,
// then a single library reduce over the sobjD array.  No small/large split is needed:
// CUB DeviceReduce::Reduce and sycl::reduction both handle arbitrary object sizes by
// tuning block occupancy internally.
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(GRID_CUDA) || defined(GRID_HIP)

template<class vobj>
inline typename vobj::scalar_objectD sumD_gpu(const vobj *lat, Integer osites)
{
  typedef typename vobj::scalar_object  sobj;
  typedef typename vobj::scalar_objectD sobjD;

  // Per-site: sum SIMD lanes (Reduce) and promote to double precision.
  deviceVector<sobjD> per_site(osites);
  sobjD *per_site_p = &per_site[0];

  accelerator_for(ss, osites, 1, {
    sobj tmp = Reduce(lat[ss]);
    sobjD tmpD; tmpD = tmp;
    per_site_p[ss] = tmpD;
  });

  // CUB global reduction over the sobjD array.
  sobjD zero; zeroit(zero);
  sobjD *d_out = static_cast<sobjD *>(acceleratorAllocDevice(sizeof(sobjD)));
  void  *d_temp = nullptr;
  size_t temp_bytes = 0;

  gpuError_t gpuErr;
  gpuErr = gpucub::DeviceReduce::Reduce(d_temp, temp_bytes, per_site_p, d_out,
                                        (int)osites, gpucub::Sum(), zero, computeStream);
  if (gpuErr != gpuSuccess) {
    std::cout << GridLogError << "Lattice_reduction_gpu_cub.h: DeviceReduce size query failed: "
              << gpuErr << std::endl;
    exit(EXIT_FAILURE);
  }

  d_temp = acceleratorAllocDevice(temp_bytes);

  gpuErr = gpucub::DeviceReduce::Reduce(d_temp, temp_bytes, per_site_p, d_out,
                                        (int)osites, gpucub::Sum(), zero, computeStream);
  if (gpuErr != gpuSuccess) {
    std::cout << GridLogError << "Lattice_reduction_gpu_cub.h: DeviceReduce failed: "
              << gpuErr << std::endl;
    exit(EXIT_FAILURE);
  }

  accelerator_barrier();

  sobjD result;
  acceleratorCopyFromDevice(d_out, &result, sizeof(sobjD));
  acceleratorFreeDevice(d_temp);
  acceleratorFreeDevice(d_out);
  return result;
}

// sumD_gpu_small and sumD_gpu_large are preserved as aliases for API compatibility.
template<class vobj>
inline typename vobj::scalar_objectD sumD_gpu_small(const vobj *lat, Integer osites)
{
  return sumD_gpu(lat, osites);
}

template<class vobj>
inline typename vobj::scalar_objectD sumD_gpu_large(const vobj *lat, Integer osites)
{
  return sumD_gpu(lat, osites);
}

template<class vobj>
inline typename vobj::scalar_object sum_gpu(const vobj *lat, Integer osites)
{
  typedef typename vobj::scalar_object sobj;
  sobj result;
  result = sumD_gpu(lat, osites);
  return result;
}

template<class vobj>
inline typename vobj::scalar_object sum_gpu_large(const vobj *lat, Integer osites)
{
  return sum_gpu(lat, osites);
}

#endif // GRID_CUDA || GRID_HIP

#if defined(GRID_SYCL)

// Accumulates in sobjD throughout, fixing the precision bug in the original
// Lattice_reduction_sycl.h which accumulated in sobj then converted at the end.
template<class vobj>
inline typename vobj::scalar_objectD sumD_gpu(const vobj *lat, Integer osites)
{
  typedef typename vobj::scalar_object  sobj;
  typedef typename vobj::scalar_objectD sobjD;

  sobjD identity; zeroit(identity);
  sobjD ret;      zeroit(ret);
  {
    sycl::buffer<sobjD, 1> abuff(&ret, {1});
    theGridAccelerator->submit([&](sycl::handler &cgh) {
      auto Reduction = sycl::reduction(abuff, cgh, identity, std::plus<>());
      cgh.parallel_for(sycl::range<1>{(size_t)osites},
                       Reduction,
                       [=](sycl::id<1> item, auto &sum) {
                         sobj s = Reduce(lat[item[0]]);
                         sobjD sd; sd = s;
                         sum += sd;
                       });
    });
  }
  return ret;
}

template<class vobj>
inline typename vobj::scalar_objectD sumD_gpu_small(const vobj *lat, Integer osites)
{
  return sumD_gpu(lat, osites);
}

template<class vobj>
inline typename vobj::scalar_objectD sumD_gpu_large(const vobj *lat, Integer osites)
{
  return sumD_gpu(lat, osites);
}

template<class vobj>
inline typename vobj::scalar_object sum_gpu(const vobj *lat, Integer osites)
{
  typedef typename vobj::scalar_object sobj;
  sobj result;
  result = sumD_gpu(lat, osites);
  return result;
}

template<class vobj>
inline typename vobj::scalar_object sum_gpu_large(const vobj *lat, Integer osites)
{
  return sum_gpu(lat, osites);
}

#endif // GRID_SYCL

NAMESPACE_END(Grid);
