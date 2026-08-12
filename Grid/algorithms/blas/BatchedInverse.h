/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: BatchedInverse.h

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

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

#include <Grid/algorithms/blas/BatchedBlas.h>

#ifdef GRID_HIP
#include <rocsolver/rocsolver.h>
#endif
// GRID_CUDA: batched LU inversion lives in cuBLAS (getrfBatched/getriBatched);
//            cublas_v2.h already included via BatchedBlas.h.
// GRID_SYCL: oneapi/mkl.hpp already included via BatchedBlas.h (lapack::getrf/getri).

NAMESPACE_BEGIN(Grid);

///////////////////////////////////////////////////////////////////////////////
// GridBLASInverse: cross-platform batched dense matrix inversion.
//
// HIGH LEVEL contract (deliberately NOT a getrf/getrs interface): invert a
// batch of dense N x N matrices IN PLACE,
//
//     A[i] <- A[i]^{-1}       i = 0 .. batchCount-1
//
// Layout: column major, lda = N, contiguous per batch element; pointer list
// exactly as GridBLAS::gemmBatched (deviceVector<T*> of device pointers).
// Each backend chooses HOW:
//   HIP  : rocSOLVER getrf_batched + getri_batched
//   CUDA : cuBLAS getrfBatched + getriBatched (out-of-place getri; workspace
//          hidden here, result copied back so the surface stays in-place)
//   SYCL : oneMKL LAPACK getrf + getri per batch element (USM, in-order queue)
//   CPU  : Eigen PartialPivLU (the correctness oracle for all of the above)
//
// The int32 vendor-batched entry points bound N < 2^31 (asserted); the huge
// single-matrix ILP64 path (getrf_64 + blocked identity-getrs harvest, proven
// in the dense coarse-coarse setup at N=69120) migrates here as a batch==1
// large-N dispatch in a follow-up -- the recursive Schur leaves are the
// batched consumers this surface is shaped for.
//
// NB GPU-backend call signatures are written to vendor documentation but the
// air-gapped development loop compiles only the CPU/Eigen path; verify the
// rocSOLVER/cuBLAS/oneMKL calls against headers on first device compile.
// Semantics are locked by the CPU unit test (Test_batched_blas).
///////////////////////////////////////////////////////////////////////////////
class GridBLASInverse {
public:

#ifdef GRID_HIP
  // rocSOLVER runs on a rocblas_handle (distinct type from hipblasHandle_t)
  static rocblas_handle & Handle(void) {
    static rocblas_handle h;
    static int init = 0;
    if ( !init ) {
      auto st = rocblas_create_handle(&h);
      GRID_ASSERT(st == rocblas_status_success);
      init = 1;
    }
    return h;
  }
#endif
#ifdef GRID_CUDA
  // cuBLAS batched LU shares the GridBLAS handle
  static cublasHandle_t & Handle(void) {
    GridBLAS::Init();
    return GridBLAS::gridblasHandle;
  }
#endif
#ifdef GRID_SYCL
  static sycl::queue * & Handle(void) {
    GridBLAS::Init();
    return GridBLAS::gridblasHandle;
  }
#endif

  GridBLASInverse() {};
  ~GridBLASInverse() {};

  void inverseBatched(int64_t N, deviceVector<ComplexF*> &Amat)
  {
    int32_t batchCount = Amat.size();
    GRID_ASSERT(batchCount > 0);

#ifdef GRID_HIP
    GRID_ASSERT( N < 2147483647L );
    rocblas_int n   = (rocblas_int)N;
    rocblas_int lda = (rocblas_int)N;

    deviceVector<rocblas_int> ipiv((uint64_t)batchCount*N);
    deviceVector<rocblas_int> info(batchCount);

    auto st1 = rocsolver_cgetrf_batched(Handle(), n, n,
					(rocblas_float_complex *const *)&Amat[0], lda,
					&ipiv[0], (rocblas_stride)N,
					&info[0], batchCount);
    GRID_ASSERT(st1 == rocblas_status_success);
    auto st2 = rocsolver_cgetri_batched(Handle(), n,
					(rocblas_float_complex *const *)&Amat[0], lda,
					&ipiv[0], (rocblas_stride)N,
					&info[0], batchCount);
    GRID_ASSERT(st2 == rocblas_status_success);
    accelerator_barrier();
    std::vector<rocblas_int> info_h(batchCount);
    acceleratorCopyFromDevice(&info[0],&info_h[0],batchCount*sizeof(rocblas_int));
    for(int i=0;i<batchCount;i++) GRID_ASSERT(info_h[i]==0); // singular pivot => abort loudly
#endif
#ifdef GRID_CUDA
    GRID_ASSERT( N < 2147483647L );
    int n = (int)N;

    deviceVector<int> ipiv((uint64_t)batchCount*N);
    deviceVector<int> info(batchCount);

    auto st1 = cublasCgetrfBatched(Handle(), n,
				   (cuComplex **)&Amat[0], n,
				   &ipiv[0], &info[0], batchCount);
    GRID_ASSERT(st1 == CUBLAS_STATUS_SUCCESS);

    // getri is OUT of place: hidden workspace keeps the surface in-place
    deviceVector<ComplexF>  work((uint64_t)batchCount*N*N);
    deviceVector<ComplexF*> Cptr(batchCount);
    std::vector<ComplexF*>  Cptr_h(batchCount);
    std::vector<ComplexF*>  Aptr_h(batchCount);
    for(int i=0;i<batchCount;i++) Cptr_h[i] = &work[(uint64_t)i*N*N];
    acceleratorCopyToDevice(&Cptr_h[0],&Cptr[0],batchCount*sizeof(ComplexF*));
    acceleratorCopyFromDevice(&Amat[0],&Aptr_h[0],batchCount*sizeof(ComplexF*));

    auto st2 = cublasCgetriBatched(Handle(), n,
				   (const cuComplex *const *)&Amat[0], n,
				   &ipiv[0],
				   (cuComplex **)&Cptr[0], n,
				   &info[0], batchCount);
    GRID_ASSERT(st2 == CUBLAS_STATUS_SUCCESS);
    accelerator_barrier();
    std::vector<int> info_h(batchCount);
    acceleratorCopyFromDevice(&info[0],&info_h[0],batchCount*sizeof(int));
    for(int i=0;i<batchCount;i++) GRID_ASSERT(info_h[i]==0);
    for(int i=0;i<batchCount;i++)
      acceleratorCopyDeviceToDevice(Cptr_h[i],Aptr_h[i],(uint64_t)N*N*sizeof(ComplexF));
#endif
#ifdef GRID_SYCL
    // Per-element oneMKL LAPACK on the in-order queue; group API optimisation later.
    sycl::queue *q = Handle();
    std::vector<ComplexF*> Aptr_h(batchCount);
    acceleratorCopyFromDevice(&Amat[0],&Aptr_h[0],batchCount*sizeof(ComplexF*));

    int64_t lwf = oneapi::mkl::lapack::getrf_scratchpad_size<std::complex<float> >(*q,N,N,N);
    int64_t lwi = oneapi::mkl::lapack::getri_scratchpad_size<std::complex<float> >(*q,N,N);
    deviceVector<ComplexF> scratchf(lwf);
    deviceVector<ComplexF> scratchi(lwi);
    deviceVector<int64_t>  ipiv(N);
    for(int i=0;i<batchCount;i++){
      oneapi::mkl::lapack::getrf(*q,N,N,(std::complex<float>*)Aptr_h[i],N,&ipiv[0],
				 (std::complex<float>*)&scratchf[0],lwf);
      oneapi::mkl::lapack::getri(*q,N,  (std::complex<float>*)Aptr_h[i],N,&ipiv[0],
				 (std::complex<float>*)&scratchi[0],lwi);
    }
    q->wait();
#endif
#if !defined(GRID_SYCL) && !defined(GRID_CUDA) && !defined(GRID_HIP)
    // Reference implementation; the oracle the unit test locks semantics with.
    thread_for (p, batchCount, {
      Eigen::Map<Eigen::MatrixXcf> eA(Amat[p],N,N);
      Eigen::PartialPivLU<Eigen::MatrixXcf> lu(eA);
      eA = lu.inverse();
      });
#endif
  }

  void inverseBatched(int64_t N, deviceVector<ComplexD*> &Amat)
  {
    int32_t batchCount = Amat.size();
    GRID_ASSERT(batchCount > 0);

#ifdef GRID_HIP
    GRID_ASSERT( N < 2147483647L );
    rocblas_int n   = (rocblas_int)N;
    rocblas_int lda = (rocblas_int)N;

    deviceVector<rocblas_int> ipiv((uint64_t)batchCount*N);
    deviceVector<rocblas_int> info(batchCount);

    auto st1 = rocsolver_zgetrf_batched(Handle(), n, n,
					(rocblas_double_complex *const *)&Amat[0], lda,
					&ipiv[0], (rocblas_stride)N,
					&info[0], batchCount);
    GRID_ASSERT(st1 == rocblas_status_success);
    auto st2 = rocsolver_zgetri_batched(Handle(), n,
					(rocblas_double_complex *const *)&Amat[0], lda,
					&ipiv[0], (rocblas_stride)N,
					&info[0], batchCount);
    GRID_ASSERT(st2 == rocblas_status_success);
    accelerator_barrier();
    std::vector<rocblas_int> info_h(batchCount);
    acceleratorCopyFromDevice(&info[0],&info_h[0],batchCount*sizeof(rocblas_int));
    for(int i=0;i<batchCount;i++) GRID_ASSERT(info_h[i]==0);
#endif
#ifdef GRID_CUDA
    GRID_ASSERT( N < 2147483647L );
    int n = (int)N;

    deviceVector<int> ipiv((uint64_t)batchCount*N);
    deviceVector<int> info(batchCount);

    auto st1 = cublasZgetrfBatched(Handle(), n,
				   (cuDoubleComplex **)&Amat[0], n,
				   &ipiv[0], &info[0], batchCount);
    GRID_ASSERT(st1 == CUBLAS_STATUS_SUCCESS);

    deviceVector<ComplexD>  work((uint64_t)batchCount*N*N);
    deviceVector<ComplexD*> Cptr(batchCount);
    std::vector<ComplexD*>  Cptr_h(batchCount);
    std::vector<ComplexD*>  Aptr_h(batchCount);
    for(int i=0;i<batchCount;i++) Cptr_h[i] = &work[(uint64_t)i*N*N];
    acceleratorCopyToDevice(&Cptr_h[0],&Cptr[0],batchCount*sizeof(ComplexD*));
    acceleratorCopyFromDevice(&Amat[0],&Aptr_h[0],batchCount*sizeof(ComplexD*));

    auto st2 = cublasZgetriBatched(Handle(), n,
				   (const cuDoubleComplex *const *)&Amat[0], n,
				   &ipiv[0],
				   (cuDoubleComplex **)&Cptr[0], n,
				   &info[0], batchCount);
    GRID_ASSERT(st2 == CUBLAS_STATUS_SUCCESS);
    accelerator_barrier();
    std::vector<int> info_h(batchCount);
    acceleratorCopyFromDevice(&info[0],&info_h[0],batchCount*sizeof(int));
    for(int i=0;i<batchCount;i++) GRID_ASSERT(info_h[i]==0);
    for(int i=0;i<batchCount;i++)
      acceleratorCopyDeviceToDevice(Cptr_h[i],Aptr_h[i],(uint64_t)N*N*sizeof(ComplexD));
#endif
#ifdef GRID_SYCL
    sycl::queue *q = Handle();
    std::vector<ComplexD*> Aptr_h(batchCount);
    acceleratorCopyFromDevice(&Amat[0],&Aptr_h[0],batchCount*sizeof(ComplexD*));

    int64_t lwf = oneapi::mkl::lapack::getrf_scratchpad_size<std::complex<double> >(*q,N,N,N);
    int64_t lwi = oneapi::mkl::lapack::getri_scratchpad_size<std::complex<double> >(*q,N,N);
    deviceVector<ComplexD> scratchf(lwf);
    deviceVector<ComplexD> scratchi(lwi);
    deviceVector<int64_t>  ipiv(N);
    for(int i=0;i<batchCount;i++){
      oneapi::mkl::lapack::getrf(*q,N,N,(std::complex<double>*)Aptr_h[i],N,&ipiv[0],
				 (std::complex<double>*)&scratchf[0],lwf);
      oneapi::mkl::lapack::getri(*q,N,  (std::complex<double>*)Aptr_h[i],N,&ipiv[0],
				 (std::complex<double>*)&scratchi[0],lwi);
    }
    q->wait();
#endif
#if !defined(GRID_SYCL) && !defined(GRID_CUDA) && !defined(GRID_HIP)
    thread_for (p, batchCount, {
      Eigen::Map<Eigen::MatrixXcd> eA(Amat[p],N,N);
      Eigen::PartialPivLU<Eigen::MatrixXcd> lu(eA);
      eA = lu.inverse();
      });
#endif
  }
};

NAMESPACE_END(Grid);
