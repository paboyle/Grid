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
// CUDA/HIP: one accelerator_for pass per site to extract SIMD lanes and promote to sobjD,
// then CUB/hipCUB DeviceReduce::Reduce over the resulting array.
//
// rocPRIM's DeviceReduce requires warpSize(64) threads per block, each holding one element
// in shared memory: sizeof(T)*64 must fit in sharedMemPerBlock.  Large QCD objects such as
// LatticePropagator (sobjD = 2304 bytes, 64*2304 = 147 KB) exceed this budget.
//
// For those types sumD_gpu_large groups the vobj's vector_type words in bundles of 4,
// reducing each bundle as a WordBundle4<scalarD> (64 bytes, 64*64 = 4 KB — always safe).
// Words that do not fill a complete bundle are zero-padded.
//
// SYCL: sycl::reduction handles any type size through the runtime, so one path suffices.
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(GRID_CUDA) || defined(GRID_HIP)

// Bundles 4 scalar_typeD values for the radix-4 large-type reduction path.
// sizeof = 4 * sizeof(scalarD) <= 64 bytes; 64 * 64 = 4096 bytes, safely within
// rocPRIM's shared-memory budget on all supported devices.
template<class scalarD>
struct WordBundle4 {
  scalarD w[4];
  accelerator_inline WordBundle4 operator+(const WordBundle4 &rhs) const {
    WordBundle4 r;
    r.w[0] = w[0] + rhs.w[0];
    r.w[1] = w[1] + rhs.w[1];
    r.w[2] = w[2] + rhs.w[2];
    r.w[3] = w[3] + rhs.w[3];
    return r;
  }
};

// Direct CUB reduction on the full scalar_objectD.
// Only safe when sizeof(sobjD)*64 <= device sharedMemPerBlock.
// Do not call directly for large composite types (e.g. LatticePropagator).
template<class vobj>
inline typename vobj::scalar_objectD sumD_gpu_direct(const vobj *lat, Integer osites)
{
  typedef typename vobj::scalar_object  sobj;
  typedef typename vobj::scalar_objectD sobjD;

  deviceVector<sobjD> per_site(osites);
  sobjD *per_site_p = &per_site[0];

  accelerator_for(ss, osites, 1, {
    sobj tmp = Reduce(lat[ss]);
    sobjD tmpD; tmpD = tmp;
    per_site_p[ss] = tmpD;
  });

  sobjD zero; zeroit(zero);
  sobjD *d_out = static_cast<sobjD *>(acceleratorAllocDevice(sizeof(sobjD)));
  void  *d_temp = nullptr;
  size_t temp_bytes = 0;

  gpuError_t gpuErr;
  gpuErr = gpucub::DeviceReduce::Reduce(d_temp, temp_bytes, per_site_p, d_out,
                                        (int)osites, gpucub::Sum(), zero, computeStream);
  if (gpuErr != gpuSuccess) {
    std::cout << GridLogError << "sumD_gpu_direct: DeviceReduce size query failed: "
              << gpuErr << std::endl;
    exit(EXIT_FAILURE);
  }

  d_temp = acceleratorAllocDevice(temp_bytes);

  gpuErr = gpucub::DeviceReduce::Reduce(d_temp, temp_bytes, per_site_p, d_out,
                                        (int)osites, gpucub::Sum(), zero, computeStream);
  if (gpuErr != gpuSuccess) {
    std::cout << GridLogError << "sumD_gpu_direct: DeviceReduce failed: "
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

// Radix-4 word-bundle path for types too large for the direct CUB path.
// Treats vobj as words of vector_type; groups them in bundles of 4 and reduces
// each bundle as a WordBundle4<scalarD>.  If words % 4 != 0, the final partial
// bundle is zero-padded so all unused slots contribute zero to the sum.
template<class vobj>
inline typename vobj::scalar_objectD sumD_gpu_large(const vobj *lat, Integer osites)
{
  typedef typename vobj::vector_type    vector;
  typedef typename vobj::scalar_typeD   scalarD;
  typedef typename vobj::scalar_objectD sobjD;
  using R4 = WordBundle4<scalarD>;

  const int words = sizeof(vobj) / sizeof(vector);
  const int nfull = words / 4;
  const int rem   = words % 4;

  sobjD ret; zeroit(ret);
  scalarD *ret_p = (scalarD *)&ret;

  iScalar<vector> *idat = (iScalar<vector> *)lat;
  deviceVector<R4> buf(osites);
  R4 *buf_p = &buf[0];

  R4 zero4;
  zero4.w[0] = zero4.w[1] = zero4.w[2] = zero4.w[3] = scalarD(0);

  R4 *d_out = static_cast<R4 *>(acceleratorAllocDevice(sizeof(R4)));
  void  *d_temp = nullptr;
  size_t temp_bytes = 0;

  // Probe workspace size once — type R4 and count osites are fixed across all groups.
  gpuError_t gpuErr;
  gpuErr = gpucub::DeviceReduce::Reduce(d_temp, temp_bytes, buf_p, d_out,
                                        (int)osites, gpucub::Sum(), zero4, computeStream);
  if (gpuErr != gpuSuccess) {
    std::cout << GridLogError << "sumD_gpu_large: DeviceReduce size query failed: "
              << gpuErr << std::endl;
    exit(EXIT_FAILURE);
  }
  d_temp = acceleratorAllocDevice(temp_bytes);

  // Full groups of 4 words.
  for (int g = 0; g < nfull; g++) {
    int base = 4 * g;
    accelerator_for(ss, osites, 1, {
      R4 r4;
      r4.w[0] = TensorRemove(Reduce(idat[ss * words + base    ]));
      r4.w[1] = TensorRemove(Reduce(idat[ss * words + base + 1]));
      r4.w[2] = TensorRemove(Reduce(idat[ss * words + base + 2]));
      r4.w[3] = TensorRemove(Reduce(idat[ss * words + base + 3]));
      buf_p[ss] = r4;
    });
    gpuErr = gpucub::DeviceReduce::Reduce(d_temp, temp_bytes, buf_p, d_out,
                                          (int)osites, gpucub::Sum(), zero4, computeStream);
    if (gpuErr != gpuSuccess) {
      std::cout << GridLogError << "sumD_gpu_large: DeviceReduce failed (group "
                << g << "): " << gpuErr << std::endl;
      exit(EXIT_FAILURE);
    }
    accelerator_barrier();
    R4 group_result;
    acceleratorCopyFromDevice(d_out, &group_result, sizeof(R4));
    ret_p[base    ] = group_result.w[0];
    ret_p[base + 1] = group_result.w[1];
    ret_p[base + 2] = group_result.w[2];
    ret_p[base + 3] = group_result.w[3];
  }

  // Partial last group: zero-pad unused slots so they contribute nothing to the sum.
  if (rem > 0) {
    int base = 4 * nfull;
    accelerator_for(ss, osites, 1, {
      R4 r4;
      r4.w[0] = r4.w[1] = r4.w[2] = r4.w[3] = scalarD(0);
      for (int k = 0; k < rem; k++)
        r4.w[k] = TensorRemove(Reduce(idat[ss * words + base + k]));
      buf_p[ss] = r4;
    });
    gpuErr = gpucub::DeviceReduce::Reduce(d_temp, temp_bytes, buf_p, d_out,
                                          (int)osites, gpucub::Sum(), zero4, computeStream);
    if (gpuErr != gpuSuccess) {
      std::cout << GridLogError << "sumD_gpu_large: DeviceReduce failed (partial group): "
                << gpuErr << std::endl;
      exit(EXIT_FAILURE);
    }
    accelerator_barrier();
    R4 partial_result;
    acceleratorCopyFromDevice(d_out, &partial_result, sizeof(R4));
    for (int k = 0; k < rem; k++)
      ret_p[4 * nfull + k] = partial_result.w[k];
  }

  acceleratorFreeDevice(d_temp);
  acceleratorFreeDevice(d_out);
  return ret;
}

// Dispatch: direct CUB path for types that fit in the shared-memory budget,
// radix-4 word-bundle path for larger types.
// Threshold 512 bytes: 64 * 512 = 32768 bytes, within rocPRIM's
// ROCPRIM_SHARED_MEMORY_MAX on all supported devices.
template<class vobj>
inline typename vobj::scalar_objectD sumD_gpu(const vobj *lat, Integer osites)
{
  typedef typename vobj::scalar_objectD sobjD;
  if constexpr (sizeof(sobjD) > 512) {
    return sumD_gpu_large(lat, osites);
  } else {
    return sumD_gpu_direct(lat, osites);
  }
}

template<class vobj>
inline typename vobj::scalar_objectD sumD_gpu_small(const vobj *lat, Integer osites)
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
  typedef typename vobj::scalar_object sobj;
  sobj result;
  result = sumD_gpu_large(lat, osites);
  return result;
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
