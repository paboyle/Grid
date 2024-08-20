#include <cassert>
#include <complex>
#include <memory>
#include <vector>
#include <algorithm>
#include <array>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <ctime>
#include <iostream>
#include <sys/time.h>

#define GRID_SYCL
#undef  GRID_HIP
#undef  GRID_CUDA

#ifdef GRID_HIP
#include <hipblas/hipblas.h>
#endif
#ifdef GRID_CUDA
#include <cublas_v2.h>
#endif
#ifdef GRID_SYCL
#include <oneapi/mkl.hpp>
#endif

#ifdef GRID_SYCL
#include <sycl/CL/sycl.hpp>
#include <sycl/usm.hpp>
cl::sycl::queue *theAccelerator;
void acceleratorInit(void)
{
  int nDevices = 1;
  cl::sycl::gpu_selector selector;
  cl::sycl::device selectedDevice { selector };
  theAccelerator = new sycl::queue (selectedDevice);
  auto name = theAccelerator->get_device().get_info<sycl::info::device::name>();
  printf("AcceleratorSyclInit: Selected device is %s\n",name.c_str()); fflush(stdout);
}
inline void *acceleratorAllocDevice(size_t bytes){ return malloc_device(bytes,*theAccelerator);};
inline void acceleratorFreeDevice(void *ptr){free(ptr,*theAccelerator);};
inline void acceleratorFreeDevice(void *ptr,size_t bytes){free(ptr,*theAccelerator);};
inline void acceleratorMemSet(void *base,int value,size_t bytes) { theAccelerator->memset(base,value,bytes); theAccelerator->wait();}
inline void acceleratorCopyToDevice(void *from,void *to,size_t bytes)  { theAccelerator->memcpy(to,from,bytes); theAccelerator->wait();}
inline void acceleratorCopyFromDevice(void *from,void *to,size_t bytes){ theAccelerator->memcpy(to,from,bytes); theAccelerator->wait();}
template<class T> void acceleratorPut(T& dev,T&host)
{
  acceleratorCopyToDevice(&host,&dev,sizeof(T));
}
template<class T> T acceleratorGet(T& dev)
{
  T host;
  acceleratorCopyFromDevice(&dev,&host,sizeof(T));
  return host;
}
#define accelerator_barrier(dummy) { theAccelerator->wait(); }
#endif


/**************************************************************
 * Allocator
 **************************************************************
 */
template<typename _Tp>
class devAllocator {
public: 
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef _Tp*       pointer;
  typedef const _Tp* const_pointer;
  typedef _Tp&       reference;
  typedef const _Tp& const_reference;
  typedef _Tp        value_type;

  template<typename _Tp1>  struct rebind { typedef devAllocator<_Tp1> other; };
  devAllocator() throw() { }
  devAllocator(const devAllocator&) throw() { }
  template<typename _Tp1> devAllocator(const devAllocator<_Tp1>&) throw() { }
  ~devAllocator() throw() { }
  pointer       address(reference __x)       const { return &__x; }
  size_type  max_size() const throw() { return size_t(-1) / sizeof(_Tp); }

  pointer allocate(size_type __n, const void* _p= 0)
  { 
    size_type bytes = __n*sizeof(_Tp);
    _Tp *ptr = (_Tp*) acceleratorAllocDevice(bytes);
    if ( (_Tp*)ptr == (_Tp *) NULL ) {
      printf("Grid Device Allocator got NULL for %lu bytes\n",(unsigned long) bytes );
    }
    assert( ( (_Tp*)ptr != (_Tp *)NULL ) );
    return ptr;
  }

  void deallocate(pointer __p, size_type __n) 
  { 
    size_type bytes = __n * sizeof(_Tp);
    acceleratorFreeDevice((void *)__p,bytes);
  }
  void construct(pointer __p, const _Tp& __val) { };
  void construct(pointer __p) { };
  void destroy(pointer __p) { };
};

template<class T> using deviceVector  = std::vector<T,devAllocator<T> >;

/**************************************************************
 * Microsecond timer
 **************************************************************
 */
inline double usecond(void) {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return 1.0e6*tv.tv_sec + 1.0*tv.tv_usec;
}


typedef float  RealF;
typedef double RealD;
typedef std::complex<float>  ComplexF;
typedef std::complex<double> ComplexD;

///////////////////////////////////////////////////////////////////////	  
// Need to rearrange lattice data to be in the right format for a
// batched multiply. Might as well make these static, dense packed
///////////////////////////////////////////////////////////////////////

#ifdef GRID_HIP
  typedef hipblasHandle_t gridblasHandle_t;
#endif
#ifdef GRID_CUDA
  typedef cublasHandle_t gridblasHandle_t;
#endif
#ifdef GRID_SYCL
  typedef cl::sycl::queue *gridblasHandle_t;
#endif
#ifdef GRID_ONE_MKL
  typedef cl::sycl::queue *gridblasHandle_t;
#endif
#if !defined(GRID_SYCL) && !defined(GRID_CUDA) && !defined(GRID_HIP) && !defined(GRID_ONE_MKL)
  typedef int32_t gridblasHandle_t;
#endif

enum GridBLASOperation_t { GridBLAS_OP_N, GridBLAS_OP_T, GridBLAS_OP_C } ;

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
      gridblasHandle = theAccelerator;
#endif
#ifdef GRID_ONE_MKL
      cl::sycl::gpu_selector selector;
      cl::sycl::device selectedDevice { selector };
      cl::sycl::property_list q_prop{cl::sycl::property::queue::in_order()};
      gridblasHandle =new sycl::queue (selectedDevice,q_prop);
#endif
      gridblasInit=1;
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
#ifdef GRID_ONE_MKL
    gridblasHandle->wait();
#endif
  }
  
  void gemmBatched(int m,int n, int k,
		   ComplexD alpha,
		   deviceVector<ComplexD*> &Amk,  // pointer list to matrices
		   deviceVector<ComplexD*> &Bkn,
		   ComplexD beta,
		   deviceVector<ComplexD*> &Cmn)
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
		   deviceVector<ComplexD*> &Cmn)
  {
    RealD t2=usecond();
    int32_t batchCount = Amk.size();
    assert(Bkn.size()==batchCount);
    assert(Cmn.size()==batchCount);

    assert(OpA!=GridBLAS_OP_T); // Complex case expect no transpose
    assert(OpB!=GridBLAS_OP_T);

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
	  eCmn = beta * eCmn + alpha * eAmk * eBkn ;
        });
      } else if ( (OpA == GridBLAS_OP_C ) && (OpB == GridBLAS_OP_N) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcd> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXcd> eBkn(Bkn[p],k,n);
	  Eigen::Map<Eigen::MatrixXcd> eCmn(Cmn[p],m,n);
	  eCmn = beta * eCmn + alpha * eAmk.adjoint() * eBkn ;
	  });
      } else if ( (OpA == GridBLAS_OP_N ) && (OpB == GridBLAS_OP_C) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcd> eAmk(Amk[p],m,k);
	  Eigen::Map<Eigen::MatrixXcd> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXcd> eCmn(Cmn[p],m,n);
	  eCmn = beta * eCmn + alpha * eAmk * eBkn.adjoint() ;
	  });
      } else if ( (OpA == GridBLAS_OP_C ) && (OpB == GridBLAS_OP_C) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcd> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXcd> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXcd> eCmn(Cmn[p],m,n);
	  eCmn = beta * eCmn + alpha * eAmk.adjoint() * eBkn.adjoint() ;
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
		   deviceVector<ComplexF*> &Cmn)
  {
    RealD t2=usecond();
    int32_t batchCount = Amk.size();

    assert(OpA!=GridBLAS_OP_T); // Complex case expect no transpose
    assert(OpB!=GridBLAS_OP_T);

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
    auto err = cublasCgemmBatched(gridblasHandle,
				  hOpA,
				  hOpB,
				  m,n,k,
				  (cuComplex *) &alpha_p[0],
				  (cuComplex **)&Amk[0], lda,
				  (cuComplex **)&Bkn[0], ldb,
				  (cuComplex *) &beta_p[0],
				  (cuComplex **)&Cmn[0], ldc,
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
						  (ComplexF *) &alpha_p[0],
						  (const ComplexF **)&Amk[0], (const int64_t *)&lda64,
						  (const ComplexF **)&Bkn[0], (const int64_t *)&ldb64,
						  (ComplexF *) &beta_p[0],
						  (ComplexF **)&Cmn[0], (const int64_t *)&ldc64,
						  (int64_t)1,&batchCount64,std::vector<sycl::event>());
    synchronise();
#endif
#if !defined(GRID_SYCL) && !defined(GRID_CUDA) && !defined(GRID_HIP)
    // Need a default/reference implementation; use Eigen
      if ( (OpA == GridBLAS_OP_N ) && (OpB == GridBLAS_OP_N) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcf> eAmk(Amk[p],m,k);
	  Eigen::Map<Eigen::MatrixXcf> eBkn(Bkn[p],k,n);
	  Eigen::Map<Eigen::MatrixXcf> eCmn(Cmn[p],m,n);
	  eCmn = beta * eCmn + alpha * eAmk * eBkn ;
	  });
      } else if ( (OpA == GridBLAS_OP_C ) && (OpB == GridBLAS_OP_N) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcf> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXcf> eBkn(Bkn[p],k,n);
	  Eigen::Map<Eigen::MatrixXcf> eCmn(Cmn[p],m,n);
	  eCmn = beta * eCmn + alpha * eAmk.adjoint() * eBkn ;
	  });
      } else if ( (OpA == GridBLAS_OP_N ) && (OpB == GridBLAS_OP_C) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcf> eAmk(Amk[p],m,k);
	  Eigen::Map<Eigen::MatrixXcf> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXcf> eCmn(Cmn[p],m,n);
	  eCmn = beta * eCmn + alpha * eAmk * eBkn.adjoint() ;
	  });
      } else if ( (OpA == GridBLAS_OP_C ) && (OpB == GridBLAS_OP_C) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXcf> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXcf> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXcf> eCmn(Cmn[p],m,n);
	  eCmn = beta * eCmn + alpha * eAmk.adjoint() * eBkn.adjoint() ;
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
	  eCmn = beta * eCmn + alpha * eAmk * eBkn ;
	  });
      } else if ( (OpA == GridBLAS_OP_T ) && (OpB == GridBLAS_OP_N) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXf> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXf> eBkn(Bkn[p],k,n);
	  Eigen::Map<Eigen::MatrixXf> eCmn(Cmn[p],m,n);
	  eCmn = beta * eCmn + alpha * eAmk.transpose() * eBkn ;
	  });
      } else if ( (OpA == GridBLAS_OP_N ) && (OpB == GridBLAS_OP_T) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXf> eAmk(Amk[p],m,k);
	  Eigen::Map<Eigen::MatrixXf> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXf> eCmn(Cmn[p],m,n);
	  eCmn = beta * eCmn + alpha * eAmk * eBkn.transpose() ;
	  });
      } else if ( (OpA == GridBLAS_OP_T ) && (OpB == GridBLAS_OP_T) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXf> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXf> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXf> eCmn(Cmn[p],m,n);
	  eCmn = beta * eCmn + alpha * eAmk.transpose() * eBkn.transpose() ;
	  } );
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
	  eCmn = beta * eCmn + alpha * eAmk * eBkn ;
	  });
      } else if ( (OpA == GridBLAS_OP_T ) && (OpB == GridBLAS_OP_N) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXd> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXd> eBkn(Bkn[p],k,n);
	  Eigen::Map<Eigen::MatrixXd> eCmn(Cmn[p],m,n);
	  eCmn = beta * eCmn + alpha * eAmk.transpose() * eBkn ;
	  });
      } else if ( (OpA == GridBLAS_OP_N ) && (OpB == GridBLAS_OP_T) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXd> eAmk(Amk[p],m,k);
	  Eigen::Map<Eigen::MatrixXd> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXd> eCmn(Cmn[p],m,n);
	  eCmn = beta * eCmn + alpha * eAmk * eBkn.transpose() ;
	  });
      } else if ( (OpA == GridBLAS_OP_T ) && (OpB == GridBLAS_OP_T) ) {
	thread_for (p, batchCount, {
	  Eigen::Map<Eigen::MatrixXd> eAmk(Amk[p],k,m);
	  Eigen::Map<Eigen::MatrixXd> eBkn(Bkn[p],n,k);
	  Eigen::Map<Eigen::MatrixXd> eCmn(Cmn[p],m,n);
	  eCmn = beta * eCmn + alpha * eAmk.transpose() * eBkn.transpose() ;
	  });
      } else { 
	assert(0);
      }
#endif
     RealD t1=usecond();
     RealD flops = 2.0*m*n*k*batchCount;
     RealD bytes = 1.0*sizeof(RealD)*(m*k+k*n+m*n)*batchCount;
  }

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
    int ncall=10;
    deviceVector<CComplex *> As(BATCH);
    deviceVector<CComplex *> Bs(BATCH);
    deviceVector<CComplex *> Cs(BATCH);
    for(int b = 0 ; b < BATCH;b++) {
      CComplex *ptr;
      ptr = &A[b*M*K];      acceleratorPut(As[b],ptr);
      ptr = &B[b*K*N];      acceleratorPut(Bs[b],ptr);
      ptr = &C[b*M*N];      acceleratorPut(Cs[b],ptr);
    }

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


gridblasHandle_t GridBLAS::gridblasHandle;
int              GridBLAS::gridblasInit;
FILE * FP;

template<class CComplex>
static void BLAS(void)
{
  //int nbasis, int nrhs, int coarseVol
  int  basis[] = { 16,32,64 };
  int  rhs[]   = { 8,12,16 };
  int  vol  = 8*8*8*8;
  int  blk  = 4*4*4*4;
  
  GridBLAS blas;
  
  int fpbits = sizeof(CComplex)*4;
  std::cout<< "=================================================================================="<<std::endl;
  std::cout<< "= batched GEMM fp"<<fpbits<<std::endl;
  std::cout<< "=================================================================================="<<std::endl;
  std::cout << "  M  "<<"\t\t"<<"N"<<"\t\t\t"<<"K"<<"\t\t"<<"Gflop/s / rank (coarse mrhs)"<<std::endl;
  std::cout << "----------------------------------------------------------"<<std::endl;
  
  fprintf(FP,"GEMM\n\n M, N, K, BATCH, GF/s per rank fp%d\n",fpbits);
  
  for(int b=0;b<3;b++){
    for(int r=0;r<3;r++){
      int M=basis[b];
      int N=rhs[r];
      int K=basis[b];
      int BATCH=vol;
      double p=blas.benchmark<CComplex>(M,N,K,BATCH);
      
      fprintf(FP,"%d, %d, %d, %d, %f\n", M, N, K, BATCH, p);
      
      std::cout<< M<<"\t\t"<<N<<"\t\t"<<K<<"\t\t"<<BATCH<<"\t\t"<<p<<std::endl;
    }}
  std::cout << "----------------------------------------------------------"<<std::endl;
  std::cout << "  M  "<<"\t\t"<<"N"<<"\t\t\t"<<"K"<<"\t\t"<<"Gflop/s / rank (block project)"<<std::endl;
  std::cout << "----------------------------------------------------------"<<std::endl;
  for(int b=0;b<3;b++){
    for(int r=0;r<3;r++){
      int M=basis[b];
      int N=rhs[r];
      int K=blk;
      int BATCH=vol;
      double p=blas.benchmark<CComplex>(M,N,K,BATCH);
      
      fprintf(FP,"%d, %d, %d, %d, %f\n", M, N, K, BATCH, p);
      std::cout<< M<<"\t\t"<<N<<"\t\t"<<K<<"\t\t"<<BATCH<<"\t\t"<<p<<std::endl;
    }}
  std::cout << "----------------------------------------------------------"<<std::endl;
  std::cout << "  M  "<<"\t\t"<<"N"<<"\t\t\t"<<"K"<<"\t\t"<<"Gflop/s / rank (block promote)"<<std::endl;
  std::cout << "----------------------------------------------------------"<<std::endl;
  for(int b=0;b<3;b++){
    for(int r=0;r<3;r++){
      int M=rhs[r];
      int N=blk;
      int K=basis[b];
      int BATCH=vol;
      double p=blas.benchmark<CComplex>(M,N,K,BATCH);
      
      fprintf(FP,"%d, %d, %d, %d, %f\n", M, N, K, BATCH, p);
      std::cout<< M<<"\t\t"<<N<<"\t\t"<<K<<"\t\t"<<BATCH<<"\t\t"<<p<<std::endl;
    }}
  fprintf(FP,"\n\n\n");
  std::cout << "=================================================================================="<<std::endl;
};


int main (int argc, char ** argv)
{
  acceleratorInit();
  FP = fopen("Benchmark_usqcd.csv","w");
  std::cout << "=================================================================================="<<std::endl;
  std::cout << " Batched BLAS benchmark " <<std::endl;
  std::cout << "=================================================================================="<<std::endl;
  BLAS<ComplexD>();
  BLAS<ComplexF>();
  fclose(FP);
}
