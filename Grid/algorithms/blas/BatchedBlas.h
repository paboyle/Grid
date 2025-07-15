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
#include <cublas_v2.h>
#endif
#ifdef GRID_SYCL
#include <oneapi/mkl.hpp>
#endif
#if 0
#define GRID_ONE_MKL
#endif
#ifdef GRID_ONE_MKL
#include <oneapi/mkl.hpp>
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
  typedef cublasHandle_t gridblasHandle_t;
#endif
#ifdef GRID_SYCL
  typedef sycl::queue *gridblasHandle_t;
#endif
#ifdef GRID_ONE_MKL
  typedef sycl::queue *gridblasHandle_t;
#endif
#if !defined(GRID_SYCL) && !defined(GRID_CUDA) && !defined(GRID_HIP) && !defined(GRID_ONE_MKL)
  typedef int32_t gridblasHandle_t;
#endif

enum GridBLASOperation_t { GridBLAS_OP_N, GridBLAS_OP_T, GridBLAS_OP_C } ;
enum GridBLASPrecision_t { GridBLAS_PRECISION_DEFAULT, GridBLAS_PRECISION_16F, GridBLAS_PRECISION_16BF, GridBLAS_PRECISION_TF32 };

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
      cublasSetPointerMode(gridblasHandle, CUBLAS_POINTER_MODE_DEVICE);
#endif
#ifdef GRID_HIP
      std::cout << "hipblasCreate"<<std::endl;
      hipblasCreate(&gridblasHandle);
#endif
#ifdef GRID_SYCL
      gridblasHandle = theGridAccelerator;
#endif
#ifdef GRID_ONE_MKL
      sycl::gpu_selector selector;
      sycl::device selectedDevice { selector };
      sycl::property_list q_prop{sycl::property::queue::in_order()};
      gridblasHandle =new sycl::queue (selectedDevice,q_prop);
#endif
      gridblasInit=1;
    }
  }

#ifdef GRID_CUDA
  cublasComputeType_t toDataType(GridBLASPrecision_t p) {
    switch (p) {
    case GridBLAS_PRECISION_16F:
      return CUBLAS_COMPUTE_32F_FAST_16F;
    case GridBLAS_PRECISION_16BF:
      return CUBLAS_COMPUTE_32F_FAST_16BF;
    case GridBLAS_PRECISION_TF32:
      return CUBLAS_COMPUTE_32F_FAST_TF32;
    default:
      assert(0);
    }
  }
#endif
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
#ifdef GRID_ONE_MKL
    gridblasHandle->wait();
#endif
  }
  
  void gemmBatched(int m,int n, int k,
		   ComplexD alpha,
		   deviceVector<ComplexD*> &Amk,  // pointer list to matrices
		   deviceVector<ComplexD*> &Bkn,
		   ComplexD beta,
		   deviceVector<ComplexD*> &Cmn,
		   GridBLASPrecision_t precision = GridBLAS_PRECISION_DEFAULT)
  {
    assert(precision == GridBLAS_PRECISION_DEFAULT);
    gemmBatched(GridBLAS_OP_N,GridBLAS_OP_N,
		m,n,k,
		alpha,
		Amk,
		Bkn,
		beta,
		Cmn);
  }
  void gemmBatched(int m,int n, int k,
		   ComplexF alpha,
		   deviceVector<ComplexF*> &Amk,  // pointer list to matrices
		   deviceVector<ComplexF*> &Bkn,
		   ComplexF beta,
		   deviceVector<ComplexF*> &Cmn)
  {
    gemmBatched(GridBLAS_OP_N,GridBLAS_OP_N,
		m,n,k,
		alpha,
		Amk,
		Bkn,
		beta,
		Cmn);
  }
  void gemmBatched(int m,int n, int k,
		   RealD alpha,
		   deviceVector<RealD*> &Amk,  // pointer list to matrices
		   deviceVector<RealD*> &Bkn,
		   RealD beta,
		   deviceVector<RealD*> &Cmn)
  {
    gemmBatched(GridBLAS_OP_N,GridBLAS_OP_N,
		m,n,k,
		alpha,
		Amk,
		Bkn,
		beta,
		Cmn);
  }
  void gemmBatched(int m,int n, int k,
		   RealF alpha,
		   deviceVector<RealF*> &Amk,  // pointer list to matrices
		   deviceVector<RealF*> &Bkn,
		   RealF beta,
		   deviceVector<RealF*> &Cmn)
  {
    gemmBatched(GridBLAS_OP_N,GridBLAS_OP_N,
		m,n,k,
		alpha,
		Amk,
		Bkn,
		beta,
		Cmn);
  }

  void gemmBatched(GridBLASOperation_t OpA,
		   GridBLASOperation_t OpB,
		   int m,int n, int k,
		   ComplexD alpha,
		   deviceVector<ComplexD*> &Amk,  // pointer list to matrices
		   deviceVector<ComplexD*> &Bkn,
		   ComplexD beta,
		   deviceVector<ComplexD*> &Cmn,
		   GridBLASPrecision_t precision = GridBLAS_PRECISION_DEFAULT)
  {
    assert(precision == GridBLAS_PRECISION_DEFAULT);
    RealD t2=usecond();
    int32_t batchCount = Amk.size();
    assert(Bkn.size()==batchCount);
    assert(Cmn.size()==batchCount);

    //assert(OpA!=GridBLAS_OP_T); // Complex case expect no transpose
    //assert(OpB!=GridBLAS_OP_T);

    int lda = m; // m x k column major
    int ldb = k; // k x n column major
    int ldc = m; // m x b column major
    if(OpA!=GridBLAS_OP_N)
      lda = k;
    if(OpB!=GridBLAS_OP_N)
      ldb = n;
    
    static deviceVector<ComplexD> alpha_p(1);
    static deviceVector<ComplexD> beta_p(1);
    // can prestore the 1 and the zero on device
    acceleratorCopyToDevice((void *)&alpha,(void *)&alpha_p[0],sizeof(ComplexD));
    acceleratorCopyToDevice((void *)&beta ,(void *)&beta_p[0],sizeof(ComplexD));
    RealD t0=usecond();
    //    std::cout << "ZgemmBatched mnk  "<<m<<","<<n<<","<<k<<" count "<<batchCount<<std::endl;
#ifdef GRID_HIP
    hipblasOperation_t hOpA;
    hipblasOperation_t hOpB;
    if ( OpA == GridBLAS_OP_N ) hOpA = HIPBLAS_OP_N;
    if ( OpA == GridBLAS_OP_T ) hOpA = HIPBLAS_OP_T;
    if ( OpA == GridBLAS_OP_C ) hOpA = HIPBLAS_OP_C;
    if ( OpB == GridBLAS_OP_N ) hOpB = HIPBLAS_OP_N;
    if ( OpB == GridBLAS_OP_T ) hOpB = HIPBLAS_OP_T;
    if ( OpB == GridBLAS_OP_C ) hOpB = HIPBLAS_OP_C;
    auto err = hipblasZgemmBatched(gridblasHandle,
				   hOpA,
				   hOpB,
				   m,n,k,
				   (hipblasDoubleComplex *) &alpha_p[0],
				   (hipblasDoubleComplex **)&Amk[0], lda,
				   (hipblasDoubleComplex **)&Bkn[0], ldb,
				   (hipblasDoubleComplex *) &beta_p[0],
				   (hipblasDoubleComplex **)&Cmn[0], ldc,
				   batchCount);
    //	 std::cout << " hipblas return code " <<(int)err<<std::endl;
    assert(err==HIPBLAS_STATUS_SUCCESS);
#endif
#ifdef GRID_CUDA
    cublasOperation_t hOpA;
    cublasOperation_t hOpB;
    if ( OpA == GridBLAS_OP_N ) hOpA = CUBLAS_OP_N;
    if ( OpA == GridBLAS_OP_T ) hOpA = CUBLAS_OP_T;
    if ( OpA == GridBLAS_OP_C ) hOpA = CUBLAS_OP_C;
    if ( OpB == GridBLAS_OP_N ) hOpB = CUBLAS_OP_N;
    if ( OpB == GridBLAS_OP_T ) hOpB = CUBLAS_OP_T;
    if ( OpB == GridBLAS_OP_C ) hOpB = CUBLAS_OP_C;
    auto err = cublasZgemmBatched(gridblasHandle,
				  hOpA,
				  hOpB,
				  m,n,k,
				  (cuDoubleComplex *) &alpha_p[0],
				  (cuDoubleComplex **)&Amk[0], lda,
				  (cuDoubleComplex **)&Bkn[0], ldb,
				  (cuDoubleComplex *) &beta_p[0],
				  (cuDoubleComplex **)&Cmn[0], ldc,
				  batchCount);
    assert(err==CUBLAS_STATUS_SUCCESS);
#endif
#ifdef GRID_SYCL
      int64_t m64=m;
      int64_t n64=n;
      int64_t k64=k;
      int64_t lda64=lda;
      int64_t ldb64=ldb;
      int64_t ldc64=ldc;
      int64_t batchCount64=batchCount;

      oneapi::mkl::transpose iOpA;
      oneapi::mkl::transpose iOpB;
      
      if ( OpA == GridBLAS_OP_N ) iOpA = oneapi::mkl::transpose::N;
      if ( OpA == GridBLAS_OP_T ) iOpA = oneapi::mkl::transpose::T;
      if ( OpA == GridBLAS_OP_C ) iOpA = oneapi::mkl::transpose::C;
      if ( OpB == GridBLAS_OP_N ) iOpB = oneapi::mkl::transpose::N;
      if ( OpB == GridBLAS_OP_T ) iOpB = oneapi::mkl::transpose::T;
      if ( OpB == GridBLAS_OP_C ) iOpB = oneapi::mkl::transpose::C;

      oneapi::mkl::blas::column_major::gemm_batch(*gridblasHandle,
						  &iOpA,
						  &iOpB,
						  &m64,&n64,&k64,
						  (ComplexD *) &alpha_p[0],
						  (const ComplexD **)&Amk[0], (const int64_t *)&lda64,
						  (const ComplexD **)&Bkn[0], (const int64_t *)&ldb64,
						  (ComplexD *) &beta_p[0],
						  (ComplexD **)&Cmn[0], (const int64_t *)&ldc64,
						  (int64_t)1,&batchCount64,std::vector<sycl::event>());
      synchronise();
#if 0
      // This code was used to check the mat mul on Sunspot/OneMKL
      std::cerr << " Called SYCL batched ZGEMM OpA "<< OpA << " OpB "<<OpB <<std::endl;
      std::vector<ComplexD> A(m*k);  // pointer list to matrices
      std::vector<ComplexD> B(k*n);
      std::vector<ComplexD> C(m*n);
      //      int sda = lda*k;
      //      int sdb = ldb*k;
      //      int sdc = ldc*n;
      std::cerr << " Checking the GEMM results "<<std::endl;
      for (int p = 0; p < 1; ++p) {
	ComplexD * Amk_p;  // pointer list to matrices
	ComplexD * Bkn_p;  // pointer list to matrices
	ComplexD * Cmn_p;  // pointer list to matrices
	acceleratorCopyFromDevice((void *)&Amk[p],(void *)&Amk_p,sizeof(ComplexD*));
	acceleratorCopyFromDevice((void *)&Bkn[p],(void *)&Bkn_p,sizeof(ComplexD*));
	acceleratorCopyFromDevice((void *)&Cmn[p],(void *)&Cmn_p,sizeof(ComplexD*));
	std::cerr << " p " << p << " copied pointers "<<std::endl;
	acceleratorCopyFromDevice((void *)Amk_p,(void *)&A[0],m*k*sizeof(ComplexD));
	acceleratorCopyFromDevice((void *)Bkn_p,(void *)&B[0],k*n*sizeof(ComplexD));
	acceleratorCopyFromDevice((void *)Cmn_p,(void *)&C[0],m*n*sizeof(ComplexD));
	std::cerr << " p " << p << " copied matrices "<<std::endl;
	std::cerr << " C[0] "<<C[0]<<std::endl;
	std::cerr << " A[0] "<<A[0]<<std::endl;
	std::cerr << " B[0] "<<B[0]<<std::endl;
	std::cerr << " m "<<m<<std::endl;
	std::cerr << " n "<<n<<std::endl;
	std::cerr << " k "<<k<<std::endl;
	for (int mm = 0; mm < m; ++mm) {
	  for (int nn = 0; nn < n; ++nn) {
	    ComplexD c_mn(0.0);
	    for (int kk = 0; kk < k; ++kk) {
	      int idx_a, idx_b;
	      //    int lda = m; // m x k column major
	      //    int ldb = k; // k x n column major
	      //    int ldc = m; // m x b column major
	      if(OpA!=GridBLAS_OP_N) {
		idx_a =kk + mm*lda;
	      } else {
		idx_a =mm + kk*lda;
	      }
	      if(OpB!=GridBLAS_OP_N) {
		idx_b =nn + kk*ldb;
	      } else {
		idx_b =kk + nn*ldb;
	      }
	      //	      std::cerr << " idx_a "<<idx_a<<" idx_b "<<idx_b<<std::endl;

	      ComplexD Ac = A[idx_a];
	      ComplexD Bc = B[idx_b];
	      if(OpA==GridBLAS_OP_C) Ac = conjugate(Ac);
	      if(OpB==GridBLAS_OP_C) Bc = conjugate(Bc);
	      
	      c_mn += Ac*Bc;
	    }
	    std::cerr << " beta "<<beta<<" alpha "<<alpha<<" C_"<<mm<<","<<nn<<" "<<c_mn<<" "<<C[mm + nn*ldc]<<std::endl;
	  }
	}
      }
#endif
#endif
#if !defined(GRID_SYCL) && !defined(GRID_CUDA) && !defined(GRID_HIP)
    // Need a default/reference implementation; use Eigen
      if ( (OpA == GridBLAS_OP_N ) && (OpB == GridBLAS_OP_N) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcd> eAmk(Amk[p],m,k);
	  Eigen::Map<Eigen::MatrixXcd> eBkn(Bkn[p],k,n);
	  Eigen::Map<Eigen::MatrixXcd> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk * eBkn ;
	  else
	    eCmn = alpha * eAmk * eBkn ;
        });
      } else if ( (OpA == GridBLAS_OP_C ) && (OpB == GridBLAS_OP_N) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcd> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXcd> eBkn(Bkn[p],k,n);
	  Eigen::Map<Eigen::MatrixXcd> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk.adjoint() * eBkn ;
	  else
	    eCmn = alpha * eAmk.adjoint() * eBkn ;
	  });
      } else if ( (OpA == GridBLAS_OP_T ) && (OpB == GridBLAS_OP_N) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcd> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXcd> eBkn(Bkn[p],k,n);
	  Eigen::Map<Eigen::MatrixXcd> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk.transpose() * eBkn ;
	  else
	    eCmn = alpha * eAmk.transpose() * eBkn ;
	  });
      } else if ( (OpA == GridBLAS_OP_N ) && (OpB == GridBLAS_OP_C) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcd> eAmk(Amk[p],m,k);
	  Eigen::Map<Eigen::MatrixXcd> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXcd> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk * eBkn.adjoint() ;
	  else
	    eCmn = alpha * eAmk * eBkn.adjoint() ;
	  });
      } else if ( (OpA == GridBLAS_OP_N ) && (OpB == GridBLAS_OP_T) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcd> eAmk(Amk[p],m,k);
	  Eigen::Map<Eigen::MatrixXcd> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXcd> eCmn(Cmn[p],m,n);
	  eCmn = beta * eCmn + alpha * eAmk * eBkn.transpose() ;
	  });
      } else if ( (OpA == GridBLAS_OP_C ) && (OpB == GridBLAS_OP_C) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcd> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXcd> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXcd> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk.adjoint() * eBkn.adjoint() ;
	  else
	    eCmn = alpha * eAmk.adjoint() * eBkn.adjoint() ;
	  } );
      } else if ( (OpA == GridBLAS_OP_T ) && (OpB == GridBLAS_OP_T) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcd> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXcd> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXcd> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk.transpose() * eBkn.transpose() ;
	  else
	    eCmn = alpha * eAmk.transpose() * eBkn.transpose() ;
	  } );
      } else { 
	assert(0);
      }
#endif
     RealD t1=usecond();
     RealD flops = 8.0*m*n*k*batchCount;
     RealD bytes = 1.0*sizeof(ComplexD)*(m*k+k*n+m*n)*batchCount;
     //     std::cout <<GridLogMessage<< " batched Blas copy "<<(t0-t2)/1.e3 <<" ms "<<std::endl;
     //     std::cout <<GridLogMessage<< " batched Blas zGemm call "<<m<<","<<n<<","<<k<<" "<< flops/(t1-t0)/1.e3 <<" GF/s "<<(t1-t0)/1.e3<<" ms "<<std::endl;
     //     std::cout <<GridLogMessage<< " batched Blas zGemm call "<<m<<","<<n<<","<<k<<" "<< bytes/(t1-t0)/1.e3 <<" GB/s "<<(t1-t0)/1.e3<<" ms "<<std::endl;
  }

  void gemmBatched(GridBLASOperation_t OpA,
		   GridBLASOperation_t OpB,
		   int m,int n, int k,
		   ComplexF alpha,
		   deviceVector<ComplexF*> &Amk,  // pointer list to matrices
		   deviceVector<ComplexF*> &Bkn,
		   ComplexF beta,
		   deviceVector<ComplexF*> &Cmn,
		   GridBLASPrecision_t precision = GridBLAS_PRECISION_DEFAULT)
  {
    RealD t2=usecond();
    int32_t batchCount = Amk.size();

    //assert(OpA!=GridBLAS_OP_T); // Complex case expect no transpose
    //assert(OpB!=GridBLAS_OP_T);

    int lda = m; // m x k column major
    int ldb = k; // k x n column major
    int ldc = m; // m x b column major
    if(OpA!=GridBLAS_OP_N)
      lda = k;
    if(OpB!=GridBLAS_OP_N)
      ldb = n;
    static deviceVector<ComplexF> alpha_p(1);
    static deviceVector<ComplexF> beta_p(1);
    // can prestore the 1 and the zero on device
    acceleratorCopyToDevice((void *)&alpha,(void *)&alpha_p[0],sizeof(ComplexF));
    acceleratorCopyToDevice((void *)&beta ,(void *)&beta_p[0],sizeof(ComplexF));
    RealD t0=usecond();

    assert(Bkn.size()==batchCount);
    assert(Cmn.size()==batchCount);
#ifdef GRID_HIP
    assert(precision == GridBLAS_PRECISION_DEFAULT);
    hipblasOperation_t hOpA;
    hipblasOperation_t hOpB;
    if ( OpA == GridBLAS_OP_N ) hOpA = HIPBLAS_OP_N;
    if ( OpA == GridBLAS_OP_T ) hOpA = HIPBLAS_OP_T;
    if ( OpA == GridBLAS_OP_C ) hOpA = HIPBLAS_OP_C;
    if ( OpB == GridBLAS_OP_N ) hOpB = HIPBLAS_OP_N;
    if ( OpB == GridBLAS_OP_T ) hOpB = HIPBLAS_OP_T;
    if ( OpB == GridBLAS_OP_C ) hOpB = HIPBLAS_OP_C;
    auto err = hipblasCgemmBatched(gridblasHandle,
				   hOpA,
				   hOpB,
				   m,n,k,
				   (hipblasComplex *) &alpha_p[0],
				   (hipblasComplex **)&Amk[0], lda,
				   (hipblasComplex **)&Bkn[0], ldb,
				   (hipblasComplex *) &beta_p[0],
				   (hipblasComplex **)&Cmn[0], ldc,
				   batchCount);

    assert(err==HIPBLAS_STATUS_SUCCESS);
#endif
#ifdef GRID_CUDA
    cublasOperation_t hOpA;
    cublasOperation_t hOpB;
    if ( OpA == GridBLAS_OP_N ) hOpA = CUBLAS_OP_N;
    if ( OpA == GridBLAS_OP_T ) hOpA = CUBLAS_OP_T;
    if ( OpA == GridBLAS_OP_C ) hOpA = CUBLAS_OP_C;
    if ( OpB == GridBLAS_OP_N ) hOpB = CUBLAS_OP_N;
    if ( OpB == GridBLAS_OP_T ) hOpB = CUBLAS_OP_T;
    if ( OpB == GridBLAS_OP_C ) hOpB = CUBLAS_OP_C;
    cublasStatus_t err;
    if (precision == GridBLAS_PRECISION_DEFAULT) {
      err = cublasCgemmBatched(gridblasHandle,
			       hOpA,
			       hOpB,
			       m,n,k,
			       (cuComplex *) &alpha_p[0],
			       (cuComplex **)&Amk[0], lda,
			       (cuComplex **)&Bkn[0], ldb,
			       (cuComplex *) &beta_p[0],
			       (cuComplex **)&Cmn[0], ldc,
			       batchCount);
    } else {
      cublasComputeType_t compute_precision = toDataType(precision);
      err = cublasGemmBatchedEx(gridblasHandle,
				hOpA,
				hOpB,
				m,n,k,
				(void *) &alpha_p[0],
				(void **)&Amk[0], CUDA_C_32F, lda,
				(void **)&Bkn[0], CUDA_C_32F, ldb,
				(void *) &beta_p[0],
				(void **)&Cmn[0], CUDA_C_32F, ldc,
				batchCount, compute_precision, CUBLAS_GEMM_DEFAULT);
    }
    assert(err==CUBLAS_STATUS_SUCCESS);
#endif
#ifdef GRID_SYCL
    assert(precision == GridBLAS_PRECISION_DEFAULT);
    int64_t m64=m;
    int64_t n64=n;
    int64_t k64=k;
    int64_t lda64=lda;
    int64_t ldb64=ldb;
    int64_t ldc64=ldc;
    int64_t batchCount64=batchCount;
    
    oneapi::mkl::transpose iOpA;
    oneapi::mkl::transpose iOpB;
    
    if ( OpA == GridBLAS_OP_N ) iOpA = oneapi::mkl::transpose::N;
    if ( OpA == GridBLAS_OP_T ) iOpA = oneapi::mkl::transpose::T;
    if ( OpA == GridBLAS_OP_C ) iOpA = oneapi::mkl::transpose::C;
    if ( OpB == GridBLAS_OP_N ) iOpB = oneapi::mkl::transpose::N;
    if ( OpB == GridBLAS_OP_T ) iOpB = oneapi::mkl::transpose::T;
    if ( OpB == GridBLAS_OP_C ) iOpB = oneapi::mkl::transpose::C;
    
    oneapi::mkl::blas::column_major::gemm_batch(*gridblasHandle,
						&iOpA,
						&iOpB,
						&m64,&n64,&k64,
						(ComplexF *) &alpha_p[0],
						(const ComplexF **)&Amk[0], (const int64_t *)&lda64,
						(const ComplexF **)&Bkn[0], (const int64_t *)&ldb64,
						(ComplexF *) &beta_p[0],
						(ComplexF **)&Cmn[0], (const int64_t *)&ldc64,
						(int64_t)1,&batchCount64,std::vector<sycl::event>());
    synchronise();
#endif
#if !defined(GRID_SYCL) && !defined(GRID_CUDA) && !defined(GRID_HIP)
    assert(precision == GridBLAS_PRECISION_DEFAULT);
    // Need a default/reference implementation; use Eigen
      if ( (OpA == GridBLAS_OP_N ) && (OpB == GridBLAS_OP_N) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcf> eAmk(Amk[p],m,k);
	  Eigen::Map<Eigen::MatrixXcf> eBkn(Bkn[p],k,n);
	  Eigen::Map<Eigen::MatrixXcf> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk * eBkn ;
	  else
	    eCmn = alpha * eAmk * eBkn ;
	  });
      } else if ( (OpA == GridBLAS_OP_C ) && (OpB == GridBLAS_OP_N) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcf> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXcf> eBkn(Bkn[p],k,n);
	  Eigen::Map<Eigen::MatrixXcf> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk.adjoint() * eBkn ;
	  else
	    eCmn = alpha * eAmk.adjoint() * eBkn ;
	  });
      } else if ( (OpA == GridBLAS_OP_T ) && (OpB == GridBLAS_OP_N) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcf> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXcf> eBkn(Bkn[p],k,n);
	  Eigen::Map<Eigen::MatrixXcf> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk.transpose() * eBkn ;
	  else
	    eCmn = alpha * eAmk.transpose() * eBkn ;
	  });
      } else if ( (OpA == GridBLAS_OP_N ) && (OpB == GridBLAS_OP_C) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcf> eAmk(Amk[p],m,k);
	  Eigen::Map<Eigen::MatrixXcf> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXcf> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk * eBkn.adjoint() ;
	  else
	    eCmn = alpha * eAmk * eBkn.adjoint() ;
	  });
      } else if ( (OpA == GridBLAS_OP_N ) && (OpB == GridBLAS_OP_T) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcf> eAmk(Amk[p],m,k);
	  Eigen::Map<Eigen::MatrixXcf> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXcf> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk * eBkn.transpose() ;
	  else
	    eCmn = alpha * eAmk * eBkn.transpose() ;
	  });
      } else if ( (OpA == GridBLAS_OP_C ) && (OpB == GridBLAS_OP_C) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcf> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXcf> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXcf> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk.adjoint() * eBkn.adjoint() ;
	  else
	    eCmn = alpha * eAmk.adjoint() * eBkn.adjoint() ;
	  } );
      } else if ( (OpA == GridBLAS_OP_T ) && (OpB == GridBLAS_OP_T) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcf> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXcf> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXcf> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk.transpose() * eBkn.transpose() ;
	  else
	    eCmn = alpha * eAmk.transpose() * eBkn.transpose() ;
	  } );
      } else { 
	assert(0);
      }
#endif
     RealD t1=usecond();
     RealD flops = 8.0*m*n*k*batchCount;
     RealD bytes = 1.0*sizeof(ComplexF)*(m*k+k*n+m*n)*batchCount;
  }
  
  ///////////////////////////////////////////////////////////////////////////
  // Single precision real GEMM
  ///////////////////////////////////////////////////////////////////////////

  void gemmBatched(GridBLASOperation_t OpA,
		   GridBLASOperation_t OpB,
		   int m,int n, int k,
		   RealF alpha,
		   deviceVector<RealF*> &Amk,  // pointer list to matrices
		   deviceVector<RealF*> &Bkn,
		   RealF beta,
		   deviceVector<RealF*> &Cmn)
  {
    RealD t2=usecond();
    int32_t batchCount = Amk.size();

    assert(OpA!=GridBLAS_OP_C); // Real case no conjugate
    assert(OpB!=GridBLAS_OP_C);

    int lda = m; // m x k column major
    int ldb = k; // k x n column major
    int ldc = m; // m x b column major
    if(OpA!=GridBLAS_OP_N)
      lda = k;
    if(OpB!=GridBLAS_OP_N)
      ldb = n;
    static deviceVector<RealF> alpha_p(1);
    static deviceVector<RealF> beta_p(1);
    // can prestore the 1 and the zero on device
    acceleratorCopyToDevice((void *)&alpha,(void *)&alpha_p[0],sizeof(RealF));
    acceleratorCopyToDevice((void *)&beta ,(void *)&beta_p[0],sizeof(RealF));
    RealD t0=usecond();

    assert(Bkn.size()==batchCount);
    assert(Cmn.size()==batchCount);
#ifdef GRID_HIP
    hipblasOperation_t hOpA;
    hipblasOperation_t hOpB;
    if ( OpA == GridBLAS_OP_N ) hOpA = HIPBLAS_OP_N;
    if ( OpA == GridBLAS_OP_T ) hOpA = HIPBLAS_OP_T;
    if ( OpA == GridBLAS_OP_C ) hOpA = HIPBLAS_OP_C;
    if ( OpB == GridBLAS_OP_N ) hOpB = HIPBLAS_OP_N;
    if ( OpB == GridBLAS_OP_T ) hOpB = HIPBLAS_OP_T;
    if ( OpB == GridBLAS_OP_C ) hOpB = HIPBLAS_OP_C;
    auto err = hipblasSgemmBatched(gridblasHandle,
				   hOpA,
				   hOpB,
				   m,n,k,
				   (float *) &alpha_p[0],
				   (float **)&Amk[0], lda,
				   (float **)&Bkn[0], ldb,
				   (float *) &beta_p[0],
				   (float **)&Cmn[0], ldc,
				   batchCount);
    assert(err==HIPBLAS_STATUS_SUCCESS);
#endif
#ifdef GRID_CUDA
    cublasOperation_t hOpA;
    cublasOperation_t hOpB;
    if ( OpA == GridBLAS_OP_N ) hOpA = CUBLAS_OP_N;
    if ( OpA == GridBLAS_OP_T ) hOpA = CUBLAS_OP_T;
    if ( OpA == GridBLAS_OP_C ) hOpA = CUBLAS_OP_C;
    if ( OpB == GridBLAS_OP_N ) hOpB = CUBLAS_OP_N;
    if ( OpB == GridBLAS_OP_T ) hOpB = CUBLAS_OP_T;
    if ( OpB == GridBLAS_OP_C ) hOpB = CUBLAS_OP_C;
    auto err = cublasSgemmBatched(gridblasHandle,
				  hOpA,
				  hOpB,
				  m,n,k,
				  (float *) &alpha_p[0],
				  (float **)&Amk[0], lda,
				  (float **)&Bkn[0], ldb,
				  (float *) &beta_p[0],
				  (float **)&Cmn[0], ldc,
				  batchCount);
    assert(err==CUBLAS_STATUS_SUCCESS);
#endif
#ifdef GRID_SYCL
      int64_t m64=m;
      int64_t n64=n;
      int64_t k64=k;
      int64_t lda64=lda;
      int64_t ldb64=ldb;
      int64_t ldc64=ldc;
      int64_t batchCount64=batchCount;

      oneapi::mkl::transpose iOpA;
      oneapi::mkl::transpose iOpB;
      
      if ( OpA == GridBLAS_OP_N ) iOpA = oneapi::mkl::transpose::N;
      if ( OpA == GridBLAS_OP_T ) iOpA = oneapi::mkl::transpose::T;
      if ( OpA == GridBLAS_OP_C ) iOpA = oneapi::mkl::transpose::C;
      if ( OpB == GridBLAS_OP_N ) iOpB = oneapi::mkl::transpose::N;
      if ( OpB == GridBLAS_OP_T ) iOpB = oneapi::mkl::transpose::T;
      if ( OpB == GridBLAS_OP_C ) iOpB = oneapi::mkl::transpose::C;

      oneapi::mkl::blas::column_major::gemm_batch(*gridblasHandle,
						  &iOpA,
						  &iOpB,
						  &m64,&n64,&k64,
						  (float *) &alpha_p[0],
						  (const float **)&Amk[0], (const int64_t *)&lda64,
						  (const float **)&Bkn[0], (const int64_t *)&ldb64,
						  (float *) &beta_p[0],
						  (float **)&Cmn[0], (const int64_t *)&ldc64,
						  (int64_t)1,&batchCount64,std::vector<sycl::event>());
      synchronise();
#endif
#if !defined(GRID_SYCL) && !defined(GRID_CUDA) && !defined(GRID_HIP)
    // Need a default/reference implementation; use Eigen
      if ( (OpA == GridBLAS_OP_N ) && (OpB == GridBLAS_OP_N) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXf> eAmk(Amk[p],m,k);
	  Eigen::Map<Eigen::MatrixXf> eBkn(Bkn[p],k,n);
	  Eigen::Map<Eigen::MatrixXf> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk * eBkn ;
	  else
	    eCmn = alpha * eAmk * eBkn ;
	  });
      } else if ( (OpA == GridBLAS_OP_T ) && (OpB == GridBLAS_OP_N) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXf> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXf> eBkn(Bkn[p],k,n);
	  Eigen::Map<Eigen::MatrixXf> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk.transpose() * eBkn ;
	  else
	    eCmn = alpha * eAmk.transpose() * eBkn ;
	  });
      } else if ( (OpA == GridBLAS_OP_N ) && (OpB == GridBLAS_OP_T) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXf> eAmk(Amk[p],m,k);
	  Eigen::Map<Eigen::MatrixXf> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXf> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk * eBkn.transpose() ;
	  else
	    eCmn = alpha * eAmk * eBkn.transpose() ;	  
	  });
      } else if ( (OpA == GridBLAS_OP_T ) && (OpB == GridBLAS_OP_T) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXf> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXf> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXf> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk.transpose() * eBkn.transpose() ;
	  else
	    eCmn = alpha * eAmk.transpose() * eBkn.transpose() ;
	  });
      } else { 
	assert(0);
      }
#endif
     RealD t1=usecond();
     RealD flops = 2.0*m*n*k*batchCount;
     RealD bytes = 1.0*sizeof(RealF)*(m*k+k*n+m*n)*batchCount;
  }
  
  
  ///////////////////////////////////////////////////////////////////////////
  // Double precision real GEMM
  ///////////////////////////////////////////////////////////////////////////
  void gemmBatched(GridBLASOperation_t OpA,
		   GridBLASOperation_t OpB,
		   int m,int n, int k,
		   RealD alpha,
		   deviceVector<RealD*> &Amk,  // pointer list to matrices
		   deviceVector<RealD*> &Bkn,
		   RealD beta,
		   deviceVector<RealD*> &Cmn)
  {
    RealD t2=usecond();
    int32_t batchCount = Amk.size();

    assert(OpA!=GridBLAS_OP_C); // Real case no conjugate
    assert(OpB!=GridBLAS_OP_C);

    int lda = m; // m x k column major
    int ldb = k; // k x n column major
    int ldc = m; // m x b column major
    if(OpA!=GridBLAS_OP_N)
      lda = k;
    if(OpB!=GridBLAS_OP_N)
      ldb = n;
    
    static deviceVector<RealD> alpha_p(1);
    static deviceVector<RealD> beta_p(1);
    // can prestore the 1 and the zero on device
    acceleratorCopyToDevice((void *)&alpha,(void *)&alpha_p[0],sizeof(RealD));
    acceleratorCopyToDevice((void *)&beta ,(void *)&beta_p[0],sizeof(RealD));
    RealD t0=usecond();

    assert(Bkn.size()==batchCount);
    assert(Cmn.size()==batchCount);
#ifdef GRID_HIP
    hipblasOperation_t hOpA;
    hipblasOperation_t hOpB;
    if ( OpA == GridBLAS_OP_N ) hOpA = HIPBLAS_OP_N;
    if ( OpA == GridBLAS_OP_T ) hOpA = HIPBLAS_OP_T;
    if ( OpA == GridBLAS_OP_C ) hOpA = HIPBLAS_OP_C;
    if ( OpB == GridBLAS_OP_N ) hOpB = HIPBLAS_OP_N;
    if ( OpB == GridBLAS_OP_T ) hOpB = HIPBLAS_OP_T;
    if ( OpB == GridBLAS_OP_C ) hOpB = HIPBLAS_OP_C;
    auto err = hipblasDgemmBatched(gridblasHandle,
				   HIPBLAS_OP_N,
				   HIPBLAS_OP_N,
				   m,n,k,
				   (double *) &alpha_p[0],
				   (double **)&Amk[0], lda,
				   (double **)&Bkn[0], ldb,
				   (double *) &beta_p[0],
				   (double **)&Cmn[0], ldc,
				   batchCount);
    assert(err==HIPBLAS_STATUS_SUCCESS);
#endif
#ifdef GRID_CUDA
    cublasOperation_t hOpA;
    cublasOperation_t hOpB;
    if ( OpA == GridBLAS_OP_N ) hOpA = CUBLAS_OP_N;
    if ( OpA == GridBLAS_OP_T ) hOpA = CUBLAS_OP_T;
    if ( OpA == GridBLAS_OP_C ) hOpA = CUBLAS_OP_C;
    if ( OpB == GridBLAS_OP_N ) hOpB = CUBLAS_OP_N;
    if ( OpB == GridBLAS_OP_T ) hOpB = CUBLAS_OP_T;
    if ( OpB == GridBLAS_OP_C ) hOpB = CUBLAS_OP_C;
    auto err = cublasDgemmBatched(gridblasHandle,
				  hOpA,
				  hOpB,
				  m,n,k,
				  (double *) &alpha_p[0],
				  (double **)&Amk[0], lda,
				  (double **)&Bkn[0], ldb,
				  (double *) &beta_p[0],
				  (double **)&Cmn[0], ldc,
				  batchCount);
    assert(err==CUBLAS_STATUS_SUCCESS);
#endif
#ifdef GRID_SYCL
      int64_t m64=m;
      int64_t n64=n;
      int64_t k64=k;
      int64_t lda64=lda;
      int64_t ldb64=ldb;
      int64_t ldc64=ldc;
      int64_t batchCount64=batchCount;

      oneapi::mkl::transpose iOpA;
      oneapi::mkl::transpose iOpB;
      
      if ( OpA == GridBLAS_OP_N ) iOpA = oneapi::mkl::transpose::N;
      if ( OpA == GridBLAS_OP_T ) iOpA = oneapi::mkl::transpose::T;
      if ( OpA == GridBLAS_OP_C ) iOpA = oneapi::mkl::transpose::C;
      if ( OpB == GridBLAS_OP_N ) iOpB = oneapi::mkl::transpose::N;
      if ( OpB == GridBLAS_OP_T ) iOpB = oneapi::mkl::transpose::T;
      if ( OpB == GridBLAS_OP_C ) iOpB = oneapi::mkl::transpose::C;

      oneapi::mkl::blas::column_major::gemm_batch(*gridblasHandle,
						  &iOpA,
						  &iOpB,
						  &m64,&n64,&k64,
						  (double *) &alpha_p[0],
						  (const double **)&Amk[0], (const int64_t *)&lda64,
						  (const double **)&Bkn[0], (const int64_t *)&ldb64,
						  (double *) &beta_p[0],
						  (double **)&Cmn[0], (const int64_t *)&ldc64,
						  (int64_t)1,&batchCount64,std::vector<sycl::event>());
      synchronise();
#endif
#if !defined(GRID_SYCL) && !defined(GRID_CUDA) && !defined(GRID_HIP)
    // Need a default/reference implementation; use Eigen
      if ( (OpA == GridBLAS_OP_N ) && (OpB == GridBLAS_OP_N) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXd> eAmk(Amk[p],m,k);
	  Eigen::Map<Eigen::MatrixXd> eBkn(Bkn[p],k,n);
	  Eigen::Map<Eigen::MatrixXd> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk * eBkn ;
	  else
	    eCmn = alpha * eAmk * eBkn ;
	  });
      } else if ( (OpA == GridBLAS_OP_T ) && (OpB == GridBLAS_OP_N) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXd> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXd> eBkn(Bkn[p],k,n);
	  Eigen::Map<Eigen::MatrixXd> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk.transpose() * eBkn ;
	  else
	    eCmn = alpha * eAmk.transpose() * eBkn ;
	  });
      } else if ( (OpA == GridBLAS_OP_N ) && (OpB == GridBLAS_OP_T) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXd> eAmk(Amk[p],m,k);
	  Eigen::Map<Eigen::MatrixXd> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXd> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk * eBkn.transpose() ;
	  else
	    eCmn = alpha * eAmk * eBkn.transpose() ;
	  });
      } else if ( (OpA == GridBLAS_OP_T ) && (OpB == GridBLAS_OP_T) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXd> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXd> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXd> eCmn(Cmn[p],m,n);
	  if (std::abs(beta) != 0.0)
	    eCmn = beta * eCmn + alpha * eAmk.transpose() * eBkn.transpose() ;
	  else
	    eCmn = alpha * eAmk.transpose() * eBkn.transpose() ;
	  });
      } else { 
	assert(0);
      }
#endif
     RealD t1=usecond();
     RealD flops = 2.0*m*n*k*batchCount;
     RealD bytes = 1.0*sizeof(RealD)*(m*k+k*n+m*n)*batchCount;
  }

  /*
    Inverse and Determinant

    - CPU version uses Eigen
    - GPU version uses LAPACK-compatible getrf / getri

    Design comment: Eigen does not expose getrf / getri in a LAPACK compatible manner.
                    Overhead to go through getrf / getri for CPU version too large.
		    Current interface therefore only guarantees the inverse and determinant
		    functions on all platforms but not the getrf / getri ones.
  */
#if !defined(GRID_SYCL) && !defined(GRID_CUDA) && !defined(GRID_HIP)

  void inverseBatched(int64_t n,
		      deviceVector<ComplexD*> &Ann,
		      deviceVector<ComplexD*> &Cnn) {

    int64_t batchCount = Ann.size();
    assert(batchCount == Cnn.size());
    thread_for(p,batchCount, {
	Eigen::Map<Eigen::MatrixXcd> eAnn(Ann[p],n,n);
	Eigen::Map<Eigen::MatrixXcd> eCnn(Cnn[p],n,n);
	eCnn = eAnn.inverse();
      });
  }

  void inverseBatched(int64_t n,
		      deviceVector<ComplexF*> &Ann,
		      deviceVector<ComplexF*> &Cnn) {

    int64_t batchCount = Ann.size();
    assert(batchCount == Cnn.size());
    thread_for(p,batchCount, {
	Eigen::Map<Eigen::MatrixXcf> eAnn(Ann[p],n,n);
	Eigen::Map<Eigen::MatrixXcf> eCnn(Cnn[p],n,n);
	eCnn = eAnn.inverse();
      });
  }

  void determinantBatched(int64_t n,
			  deviceVector<ComplexD*> &Ann,
			  deviceVector<ComplexD*> &C) {

    int64_t batchCount = Ann.size();
    assert(batchCount == C.size());
    thread_for(p,batchCount, {
	Eigen::Map<Eigen::MatrixXcd> eAnn(Ann[p],n,n);
	*C[p] = eAnn.determinant();
      });
  }

  void determinantBatched(int64_t n,
			  deviceVector<ComplexF*> &Ann,
			  deviceVector<ComplexF*> &C) {

    int64_t batchCount = Ann.size();
    assert(batchCount == C.size());
    thread_for(p,batchCount, {
	Eigen::Map<Eigen::MatrixXcf> eAnn(Ann[p],n,n);
	*C[p] = eAnn.determinant();
      });
  }
  
#else

#ifdef GRID_SYCL
  template<typename T>
  void getrfBatchedSYCL(int64_t n,
			deviceVector<T*> &Ann,
			deviceVector<int64_t> &ipiv,
			deviceVector<int64_t> &info) {
    
    int64_t batchCount = Ann.size();

    static deviceVector<T> scratchpad;
    int64_t sp_size = oneapi::mkl::lapack::getrf_batch_scratchpad_size<T>(*gridblasHandle, &n, &n, &n, (int64_t)1, &batchCount);
    if (sp_size > scratchpad.size())
      scratchpad.resize(sp_size);

    static deviceVector<int64_t*> _ipiv;
    if (batchCount > _ipiv.size())
      _ipiv.resize(batchCount);
    int64_t** p_ipiv = &_ipiv[0];
    int64_t* pipiv = &ipiv[0];

    accelerator_for(i, batchCount, 1, { p_ipiv[i] = &pipiv[i*n]; });

    oneapi::mkl::lapack::getrf_batch(*gridblasHandle,
				    &n, &n,
				    (T **)&Ann[0],
				    &n,
				    (int64_t**)&_ipiv[0],
				    (int64_t)1, &batchCount,
				    (T*)&scratchpad[0], (int64_t)scratchpad.size(),
				    std::vector<sycl::event>());
    synchronise();
  }
#endif

  void getrfBatched(int64_t n,
		    deviceVector<ComplexD*> &Ann,
		    deviceVector<int64_t> &ipiv,
		    deviceVector<int64_t> &info)
  {
    int64_t batchCount = Ann.size();
    assert(ipiv.size()==batchCount*n);
    assert(info.size()==batchCount);

#ifdef GRID_HIP
    auto err = hipblasZgetrfBatched(gridblasHandle,(int)n,
				    (hipblasDoubleComplex **)&Ann[0], (int)n,
				    (int*) &ipiv[0],
				    (int*) &info[0],
				    (int)batchCount);
    assert(err==HIPBLAS_STATUS_SUCCESS);
#endif
#ifdef GRID_CUDA
    auto err = cublasZgetrfBatched(gridblasHandle, (int)n,
				   (cuDoubleComplex **)&Ann[0], (int)n,
				   (int*) &ipiv[0],
				   (int*) &info[0],
				   (int)batchCount);
    assert(err==CUBLAS_STATUS_SUCCESS);
#endif
#ifdef GRID_SYCL
    getrfBatchedSYCL(n, Ann, ipiv, info);
#endif
  }

  void getrfBatched(int64_t n,
		    deviceVector<ComplexF*> &Ann,
		    deviceVector<int64_t> &ipiv,
		    deviceVector<int64_t> &info)
  {
    int64_t batchCount = Ann.size();
    assert(ipiv.size()==batchCount*n);
    assert(info.size()==batchCount);

#ifdef GRID_HIP
    auto err = hipblasCgetrfBatched(gridblasHandle,(int)n,
				    (hipblasComplex **)&Ann[0], (int)n,
				    (int*) &ipiv[0],
				    (int*) &info[0],
				    (int)batchCount);
    assert(err==HIPBLAS_STATUS_SUCCESS);
#endif
#ifdef GRID_CUDA
    auto err = cublasCgetrfBatched(gridblasHandle, (int)n,
				   (cuComplex **)&Ann[0], (int)n,
				   (int*) &ipiv[0],
				   (int*) &info[0],
				   (int)batchCount);
    assert(err==CUBLAS_STATUS_SUCCESS);
#endif
#ifdef GRID_SYCL
    getrfBatchedSYCL(n, Ann, ipiv, info);
#endif
  }

#ifdef GRID_SYCL
  template<typename T>
  void getriBatchedSYCL(int64_t n,
			deviceVector<T*> &Ann,
			deviceVector<int64_t> &ipiv,
			deviceVector<int64_t> &info,
			deviceVector<T*> &Cnn) {

    int64_t batchCount = Ann.size();

    static deviceVector<T> scratchpad;
    int64_t sp_size = oneapi::mkl::lapack::getri_batch_scratchpad_size<T>(*gridblasHandle, &n, &n, (int64_t)1, &batchCount);
    if (sp_size > scratchpad.size())
      scratchpad.resize(sp_size);

    static deviceVector<int64_t*> _ipiv;
    if (batchCount > _ipiv.size())
      _ipiv.resize(batchCount);
    int64_t** p_ipiv = &_ipiv[0];
    int64_t* pipiv = &ipiv[0];

    accelerator_for(i, batchCount, 1, { p_ipiv[i] = &pipiv[i*n]; });

    oneapi::mkl::lapack::getri_batch(*gridblasHandle,
				     &n,
				     (T **)&Ann[0],
				     &n,
				     (int64_t**)p_ipiv,
				     (int64_t)1, &batchCount,
				     (T *)&scratchpad[0], (int64_t)scratchpad.size(),
				     std::vector<sycl::event>());

    synchronise();

    T** pA = &Ann[0];
    T** pC = &Cnn[0];
    accelerator_for(i, batchCount*n*n, 1, {
	auto j = i / batchCount;
	auto k = i % batchCount;
	pC[k][j] = pA[k][j];
      });
  }

#endif

  void getriBatched(int64_t n,
		    deviceVector<ComplexD*> &Ann,
		    deviceVector<int64_t> &ipiv,
		    deviceVector<int64_t> &info,
		    deviceVector<ComplexD*> &Cnn)
  {
    int64_t batchCount = Ann.size();
    assert(ipiv.size()==batchCount*n);
    assert(info.size()==batchCount);
    assert(Cnn.size()==batchCount);

#ifdef GRID_HIP
    auto err = hipblasZgetriBatched(gridblasHandle,(int)n,
				    (hipblasDoubleComplex **)&Ann[0], (int)n,
				    (int*) &ipiv[0],
				    (hipblasDoubleComplex **)&Cnn[0], (int)n,
				    (int*) &info[0],
				    (int)batchCount);
    assert(err==HIPBLAS_STATUS_SUCCESS);
#endif
#ifdef GRID_CUDA
    auto err = cublasZgetriBatched(gridblasHandle, (int)n,
				   (cuDoubleComplex **)&Ann[0], (int)n,
				   (int*) &ipiv[0],
				   (cuDoubleComplex **)&Cnn[0], (int)n,
				   (int*) &info[0],
				   (int)batchCount);
    assert(err==CUBLAS_STATUS_SUCCESS);
#endif
#ifdef GRID_SYCL
    getriBatchedSYCL(n, Ann, ipiv, info, Cnn);
#endif
  }

  void getriBatched(int64_t n,
		    deviceVector<ComplexF*> &Ann,
		    deviceVector<int64_t> &ipiv,
		    deviceVector<int64_t> &info,
		    deviceVector<ComplexF*> &Cnn)
  {
    int64_t batchCount = Ann.size();
    assert(ipiv.size()==batchCount*n);
    assert(info.size()==batchCount);
    assert(Cnn.size()==batchCount);

#ifdef GRID_HIP
    auto err = hipblasCgetriBatched(gridblasHandle,(int)n,
				    (hipblasComplex **)&Ann[0], (int)n,
				    (int*) &ipiv[0],
				    (hipblasComplex **)&Cnn[0], (int)n,
				    (int*) &info[0],
				    (int)batchCount);
    assert(err==HIPBLAS_STATUS_SUCCESS);
#endif
#ifdef GRID_CUDA
    auto err = cublasCgetriBatched(gridblasHandle, (int)n,
				   (cuComplex **)&Ann[0], (int)n,
				   (int*) &ipiv[0],
				   (cuComplex **)&Cnn[0], (int)n,
				   (int*) &info[0],
				   (int)batchCount);
    assert(err==CUBLAS_STATUS_SUCCESS);
#endif
#ifdef GRID_SYCL
    getriBatchedSYCL(n, Ann, ipiv, info, Cnn);
#endif
  }

  template<typename dtype>
  void inverseBatched(int64_t n,
		      deviceVector<dtype*> &Ann, // this will be overwritten with LU decomposition
		      deviceVector<dtype*> &Cnn  // this will be overwritten with the inverse
		      ) {

    int64_t batchCount = Ann.size();
    RealD t0 = usecond();
    deviceVector<int64_t> ipiv(batchCount*n);
    deviceVector<int64_t> info(batchCount);

    //RealD t1 = usecond();
    getrfBatched(n, Ann, ipiv, info);
    // test info for non-invertibility?  set to nan if yes?
    getriBatched(n, Ann, ipiv, info, Cnn);
    //synchronise();
    //RealD t2 = usecond();
    //std::cout << GridLogMessage << "Temp " << t1-t0 << " rf/ri " << t2-t1  << std::endl;
  }

  template<typename dtype>
  void determinantBatched(int64_t n,
			  deviceVector<dtype*> &Ann, // this will be overwritten with LU decomposition
			  deviceVector<dtype*> &C    // this will be overwritten with determinant
			  ) {

    int64_t batchCount = Ann.size();
    //RealD t0 = usecond();
    deviceVector<int64_t> ipiv(batchCount*n);
    deviceVector<int64_t> info(batchCount);
    
    dtype** pAnn = (dtype**)&Ann[0];
    dtype** pC = (dtype**)&C[0];
#if defined(GRID_CUDA) || defined(GRID_HIP)
    int* pipiv = (int*)&ipiv[0];
#else
    int64_t* pipiv = (int64_t*)&ipiv[0];
#endif

    //RealD t1 = usecond();
    getrfBatched(n, Ann, ipiv, info);
    //RealD t2 = usecond();
    accelerator_for(i,batchCount,1,{
	dtype det = 1.0;
	for (int64_t j=0;j<n;j++) {
	  det *= pAnn[i][n*j + j];
	  // branchless signs
	  det *= (pipiv[i*n + j] == j+1) ? (1.0) : (-1.0);
	}
	*pC[i] = det;
      });
    
    //RealD t3 = usecond();
    //std::cout << GridLogMessage << "Temp " << t1 - t0 << " rf/ri " << t2-t1  << "final" << t3 - t2 << std::endl;
  }
#endif
  
  
  template<class CComplex>
  double benchmark(int M, int N, int K, int BATCH)
  {
    int32_t N_A = M*K*BATCH;
    int32_t N_B = K*N*BATCH;
    int32_t N_C = M*N*BATCH;
    deviceVector<CComplex> A(N_A); acceleratorMemSet(&A[0],0,N_A*sizeof(CComplex));
    deviceVector<CComplex> B(N_B); acceleratorMemSet(&B[0],0,N_B*sizeof(CComplex));
    deviceVector<CComplex> C(N_C); acceleratorMemSet(&C[0],0,N_C*sizeof(CComplex));
    CComplex alpha(1.0);
    CComplex beta (1.0);
    RealD flops = 8.0*M*N*K*BATCH;
    int ncall=1000;
    deviceVector<CComplex *> As(BATCH);
    deviceVector<CComplex *> Bs(BATCH);
    deviceVector<CComplex *> Cs(BATCH);
    for(int b = 0 ; b < BATCH;b++) {
      CComplex *ptr;
      ptr = &A[b*M*K];      acceleratorPut(As[b],ptr);
      ptr = &B[b*K*N];      acceleratorPut(Bs[b],ptr);
      ptr = &C[b*M*N];      acceleratorPut(Cs[b],ptr);
    }

    // Warm up call
    gemmBatched(M,N,K,
		alpha,
		As, // m x k 
		Bs, // k x n
		beta, 
		Cs);
    synchronise();

    RealD t0 = usecond();
    for(int i=0;i<ncall;i++){
      gemmBatched(M,N,K,
		  alpha,
		  As, // m x k 
		  Bs, // k x n
		  beta, 
		  Cs);
      synchronise();
    }
    RealD t1 = usecond();
    RealD bytes = 1.0*sizeof(CComplex)*(M*N*2+N*K+M*K)*BATCH;
    flops = 8.0*M*N*K*BATCH*ncall;
    flops = flops/(t1-t0)/1.e3;
    return flops; // Returns gigaflops
  }

};

NAMESPACE_END(Grid);
