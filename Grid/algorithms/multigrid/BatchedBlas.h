/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: BatchedBlas.h

    Copyright (C) 2023

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

#ifdef GRID_HIP
#include <hipblas/hipblas.h>
#endif
#ifdef GRID_CUDA
#include <hipblas/hipblas.h>
#endif
#ifdef GRID_SYCL
#error // need oneMKL version
#endif

///////////////////////////////////////////////////////////////////////	  
// Need to rearrange lattice data to be in the right format for a
// batched multiply. Might as well make these static, dense packed
///////////////////////////////////////////////////////////////////////
NAMESPACE_BEGIN(Grid);
#ifdef GRID_HIP
  typedef hipblasHandle_t gridblasHandle_t;
#endif
#ifdef GRID_CUDA
  typedef cudablasHandle_t gridblasHandle_t;
#endif
#ifdef GRID_SYCL
  typedef int32_t gridblasHandle_t;
#endif
#if !defined(GRID_SYCL) && !defined(GRID_CUDA) && !defined(GRID_HIP)
  typedef int32_t gridblasHandle_t;
#endif

class GridBLAS {
public:

     static gridblasHandle_t gridblasHandle;
     static int            gridblasInit;
  
     static void Init(void)
     {
       if ( ! gridblasInit ) {
#ifdef GRID_CUDA
	 std::cout << "cublasCreate"<<std::endl;
         cublasCreate(&gridblasHandle);
#endif
#ifdef GRID_HIP
	 std::cout << "hipblasCreate"<<std::endl;
         hipblasCreate(&gridblasHandle);
#endif
#ifdef GRID_SYCL
	 #error
#endif
       }
     }

  // Force construct once
  GridBLAS() { Init(); };
  ~GridBLAS() { };
  
/////////////////////////////////////////////////////////////////////////////////////
// BLAS GEMM conventions:
/////////////////////////////////////////////////////////////////////////////////////
// - C = alpha A * B + beta C
// Dimensions:
// - C_m.n
// - A_m.k
// - B_k.n
// - Flops = 8 M N K
// - Bytes = 2*sizeof(word) * (MN+MK+KN)
// M=60, N=12
// Flop/Byte = 8 . 60.60.12 / (60.12+60.60+60.12)/16 = 4 so expect about 4 TF/s on a GCD
/////////////////////////////////////////////////////////////////////////////////////
  void synchronise(void)
  {
#ifdef GRID_HIP
    auto err = hipDeviceSynchronize();
    assert(err==hipSuccess);
#endif
#ifdef GRID_CUDA
    auto err = cudaDeviceSynchronize();
    assert(err==cudaSuccess);
#endif
#ifdef GRID_SYCL
    accelerator_barrier();
#endif
  }
  void benchmark(int nbasis, int nrhs, int coarseVol, int nstencil)
     {
       int32_t N_A = nbasis*nbasis*coarseVol*nstencil;
       int32_t N_B = nbasis*nrhs*coarseVol*nstencil; // One leg of stencil at a time
       int32_t N_C = nbasis*nrhs*coarseVol*nstencil; 
       deviceVector<ComplexD> A(N_A); acceleratorMemSet(&A[0],0,N_A*sizeof(ComplexD));
       deviceVector<ComplexD> B(N_B); acceleratorMemSet(&B[0],0,N_B*sizeof(ComplexD));
       deviceVector<ComplexD> C(N_C); acceleratorMemSet(&C[0],0,N_C*sizeof(ComplexD));
       ComplexD alpha(1.0);
       ComplexD beta (1.0);
       for(int i=0;i<10;i++){
	 RealD t0 = usecond();
	 for(int s=0;s<nstencil;s++){
	   gemmStridedBatched(nbasis,nrhs,nbasis,
		     alpha,
		     &A[0], // m x k 
		     &B[0], // k x n
		     beta, 
		     &C[0], // m x n
		     coarseVol);
	 }
	 synchronise();
	 RealD t1 = usecond();
	 RealD flops = 8.0*nbasis*nbasis*nrhs*coarseVol*nstencil;
	 RealD bytes = 1.0*sizeof(ComplexD)*(nbasis*nbasis+nbasis*nrhs*3)*coarseVol*nstencil;
	 std::cout << " batched Blas call "<<i<<" "<< flops/(t1-t0)/1.e3 <<" GF/s "<<(t1-t0)/1.e3<<" ms "<<std::endl;
	 std::cout << " batched Blas call "<<i<<" "<< bytes/(t1-t0)/1.e3 <<" GB/s "<<(t1-t0)/1.e3<<" ms "<<std::endl;
       }
     }

     void gemmBatched(int m,int n, int k,
		      ComplexD alpha,
		      deviceVector<ComplexD*> &Amk,  // pointer list to matrices
		      deviceVector<ComplexD*> &Bkn,
		      ComplexD beta,
		      deviceVector<ComplexD*> &Cmn)
     {
       RealD t2=usecond();
       int32_t batchCount = Amk.size();
       // Use C-row major storage, so transpose calls
       int lda = m; // m x k column major
       int ldb = k; // k x n column major
       int ldc = m; // m x b column major
       static deviceVector<ComplexD> alpha_p(1);
       static deviceVector<ComplexD> beta_p(1);
       // can prestore the 1 and the zero on device
       acceleratorCopyToDevice((void *)&alpha,(void *)&alpha_p[0],sizeof(ComplexD));
       acceleratorCopyToDevice((void *)&beta ,(void *)&beta_p[0],sizeof(ComplexD));
       RealD t0=usecond();
#ifdef GRID_HIP
       std::cout << "hipblasZgemmBatched mnk  "<<m<<","<<n<<","<<k<<" count "<<batchCount<<std::endl;
       assert(Bkn.size()==batchCount);
       assert(Cmn.size()==batchCount);
       auto err = hipblasZgemmBatched(gridblasHandle,
				      HIPBLAS_OP_N,
				      HIPBLAS_OP_N,
				      m,n,k,
				      (hipblasDoubleComplex *) &alpha_p[0],
				      (hipblasDoubleComplex **)&Amk[0], lda,
				      (hipblasDoubleComplex **)&Bkn[0], ldb,
				      (hipblasDoubleComplex *) &beta_p[0],
				      (hipblasDoubleComplex **)&Cmn[0], ldc,
				      batchCount);
       //	 std::cout << " hipblas return code " <<(int)err<<std::endl;
       assert(err==HIPBLAS_STATUS_SUCCESS);
       synchronise();
#endif
#ifdef GRID_CUDA
     #error "CUDA implemenetation "
#endif
#ifdef GRID_SYCL
     #error "oneMKL implemenetation "
#endif
#if !defined(GRID_SYCL) && !defined(GRID_CUDA) && !defined(GRID_HIP)
     // Need a default/reference implementation
     for (int p = 0; p < batchCount; ++p) {
       for (int mm = 0; mm < m; ++mm) {
	 for (int nn = 0; nn < n; ++nn) {
	   ComplexD c_mn(0.0);
	   for (int kk = 0; kk < k, ++kk)
	     c_mn += Amk[mm + kk*lda + p*sda] * Bkn[kk + nn*ldb + p*sdb];
	   Cmn[mm + nn*ldc + p*sdc] =  (*alpha_p)*c_mn + (*beta_p)*Cmn[mm + nn*ldc + p*sdc];
	 }
       }
     }
#endif
     RealD t1=usecond();
	 //	 std::cout << " hipblas synchronised " <<std::endl;
     RealD flops = 8.0*m*n*k*batchCount;
     RealD bytes = 1.0*sizeof(ComplexD)*(m*k+k*n+m*n)*batchCount;
     std::cout << " batched Blas copy "<<(t0-t2)/1.e3 <<" ms "<<std::endl;
     std::cout << " batched Blas call "<<m<<","<<n<<","<<k<<" "<< flops/(t1-t0)/1.e3 <<" GF/s "<<(t1-t0)/1.e3<<" ms "<<std::endl;
     std::cout << " batched Blas call "<<m<<","<<n<<","<<k<<" "<< bytes/(t1-t0)/1.e3 <<" GB/s "<<(t1-t0)/1.e3<<" ms "<<std::endl;
     
  }

  
     void gemmStridedBatched(int m,int n, int k,
		      ComplexD alpha,
		      ComplexD* Amk,  // pointer list to matrices
		      ComplexD* Bkn,
		      ComplexD beta,
		      ComplexD* Cmn,
		      int batchCount)
     {
       // Use C-row major storage, so transpose calls
       int lda = m; // m x k column major
       int ldb = k; // k x n column major
       int ldc = m; // m x b column major
       int sda = m*k;
       int sdb = k*n;
       int sdc = m*n;
       deviceVector<ComplexD> alpha_p(1);
       deviceVector<ComplexD> beta_p(1);
       acceleratorCopyToDevice((void *)&alpha,(void *)&alpha_p[0],sizeof(ComplexD));
       acceleratorCopyToDevice((void *)&beta ,(void *)&beta_p[0],sizeof(ComplexD));
#ifdef GRID_HIP
       std::cout << "hipblasZgemmStridedBatched mnk  "<<m<<","<<n<<","<<k<<" count "<<batchCount<<std::endl;
       std::cout << "hipblasZgemmStridedBatched ld   "<<lda<<","<<ldb<<","<<ldc<<std::endl;
       std::cout << "hipblasZgemmStridedBatched sd   "<<sda<<","<<sdb<<","<<sdc<<std::endl;
       {
	 auto err = hipblasZgemmStridedBatched(gridblasHandle,
					       HIPBLAS_OP_N,
					       HIPBLAS_OP_N,
					       m,n,k,
					       (hipblasDoubleComplex *) &alpha_p[0],
					       (hipblasDoubleComplex *) Amk, lda, sda,
					       (hipblasDoubleComplex *) Bkn, ldb, sdb,
					       (hipblasDoubleComplex *) &beta_p[0],
					       (hipblasDoubleComplex *) Cmn, ldc, sdc,
					       batchCount);
	 std::cout << " hipblas return code " <<(int)err<<std::endl;
	 assert(err==HIPBLAS_STATUS_SUCCESS);
     }
#endif
#ifdef GRID_CUDA
     cublasZgemmStridedBatched(gridblasHandle,
			CUBLAS_OP_T,
			CUBLAS_OP_T,
			m,n,k,
			(cuDoubleComplex *)&alpha_p[0],
		        (cuDoubleComplex *) Amk, lda, sda,
			(cuDoubleComplex *) Bkn, ldb, sdb,
			(cuDoubleComplex *)&beta_p[],
			(cuDoubleComplex *) Cmn, ldc, sdc,
			batchCount);
#endif
#ifdef GRID_SYCL
     #error "oneMKL implemenetation "
#endif
#if !defined(GRID_SYCL) && !defined(GRID_CUDA) && !defined(GRID_HIP)
     // Need a default/reference implementation
     for (int p = 0; p < batchCount; ++p) {
       for (int mm = 0; mm < m; ++mm) {
	 for (int nn = 0; nn < n; ++nn) {
	   ComplexD c_mn(0.0);
	   for (int kk = 0; kk < k, ++kk)
	     c_mn += Amk[mm + kk*lda + p*sda] * Bkn[kk + nn*ldb + p*sdb];
	   Cmn[mm + nn*ldc + p*sdc] =  (*alpha_p)*c_mn + (*beta_p)*Cmn[mm + nn*ldc + p*sdc];
	 }
       }
     }
#endif
  }
};

NAMESPACE_END(Grid);
