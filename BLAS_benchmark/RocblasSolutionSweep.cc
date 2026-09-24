///////////////////////////////////////////////////////////////////////////////
// Which rocBLAS kernel does the coarse multigrid GEMM get, and is a better
// one available?
//
// The coarse operator issues a batched complex GEMM per stencil point with
// M = K = nbasis, N = Nrhs and batch = the local coarse volume.  Traces of a
// production run show the complex-single kernel taking as long as the
// complex-double one while reading half the bytes, and the kernel names say
// why: the library picks a 64x64 macro tile for single and 16x16 for double,
// while the output block is only nbasis x Nrhs.  At nbasis 16 single
// precision comes out absolutely slower than double.
//
// hipBLAS offers no way to override that choice, but rocBLAS's extended
// entry point takes a solution index.  This walks the indices, keeps the ones
// the library accepts for this problem, times each, and prints them against
// the default heuristic.  No Grid, no MPI: one GPU, a few seconds.
//
//   hipcc -O3 -std=c++17 RocblasSolutionSweep.cc -lrocblas -lamdhip64 \
//         -o RocblasSolutionSweep
//   ./RocblasSolutionSweep [nbasis] [nrhs] [batch]
//
// If the accepted solutions include one roughly twice the double-precision
// rate, the application can select it by index and the tuning gap is ours to
// close.  If none does, the gap is the library's, and this output is the
// report.
///////////////////////////////////////////////////////////////////////////////
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>
#include <sys/time.h>
#include <hip/hip_runtime.h>
#include <rocblas/rocblas.h>

static double usecond(void)
{
  struct timeval tv; gettimeofday(&tv,NULL);
  return 1.0e6*tv.tv_sec + 1.0*tv.tv_usec;
}
#define HIP_CHECK(x)  do { hipError_t e=(x); if(e!=hipSuccess){ printf("HIP error %s at line %d\n",hipGetErrorString(e),__LINE__); exit(1);} } while(0)

struct Result { int index; double gflops; double gbs; };

// One shape, one datatype, every solution index the library accepts.
static void sweep(rocblas_handle handle,
		  rocblas_datatype type, const char *name, size_t elembytes,
		  int M,int N,int K,int BATCH,int maxindex,int ncall)
{
  // Device data.  Contents are irrelevant to the timing; the shapes are not.
  void *A,*B,*C;
  HIP_CHECK(hipMalloc(&A,(size_t)M*K*BATCH*elembytes));
  HIP_CHECK(hipMalloc(&B,(size_t)K*N*BATCH*elembytes));
  HIP_CHECK(hipMalloc(&C,(size_t)M*N*BATCH*elembytes));
  HIP_CHECK(hipMemset(A,0,(size_t)M*K*BATCH*elembytes));
  HIP_CHECK(hipMemset(B,0,(size_t)K*N*BATCH*elembytes));
  HIP_CHECK(hipMemset(C,0,(size_t)M*N*BATCH*elembytes));

  // Pointer arrays, as the batched interface wants them, on the device.
  std::vector<void *> hA(BATCH),hB(BATCH),hC(BATCH);
  for(int b=0;b<BATCH;b++){
    hA[b] = (char *)A + (size_t)b*M*K*elembytes;
    hB[b] = (char *)B + (size_t)b*K*N*elembytes;
    hC[b] = (char *)C + (size_t)b*M*N*elembytes;
  }
  void **dA,**dB,**dC;
  HIP_CHECK(hipMalloc(&dA,BATCH*sizeof(void*)));
  HIP_CHECK(hipMalloc(&dB,BATCH*sizeof(void*)));
  HIP_CHECK(hipMalloc(&dC,BATCH*sizeof(void*)));
  HIP_CHECK(hipMemcpy(dA,&hA[0],BATCH*sizeof(void*),hipMemcpyHostToDevice));
  HIP_CHECK(hipMemcpy(dB,&hB[0],BATCH*sizeof(void*),hipMemcpyHostToDevice));
  HIP_CHECK(hipMemcpy(dC,&hC[0],BATCH*sizeof(void*),hipMemcpyHostToDevice));

  // alpha = beta = 1: the stencil accumulates into C, so C is read as well
  // as written and the traffic below counts it twice.
  double alpha_d[2] = {1.0,0.0}, beta_d[2] = {1.0,0.0};
  float  alpha_f[2] = {1.0f,0.0f}, beta_f[2] = {1.0f,0.0f};
  const void *alpha = (type==rocblas_datatype_f64_c) ? (const void *)alpha_d : (const void *)alpha_f;
  const void *beta  = (type==rocblas_datatype_f64_c) ? (const void *)beta_d  : (const void *)beta_f;

  double flops = 8.0*M*N*K*BATCH;                                  // complex madd
  double bytes = (double)elembytes*((double)M*K + (double)K*N + 2.0*(double)M*N)*BATCH;

  auto run = [&](rocblas_gemm_algo algo,int solution)->double
  {
    rocblas_status st = rocblas_gemm_batched_ex(handle,
			   rocblas_operation_none, rocblas_operation_none,
			   M,N,K, alpha,
			   (const void *)dA, type, M,
			   (const void *)dB, type, K, beta,
			   (const void *)dC, type, M,
			   (void *)dC,       type, M,
			   BATCH, type, algo, solution, 0);
    if ( st != rocblas_status_success ) return -1.0;
    HIP_CHECK(hipDeviceSynchronize());
    double t0=usecond();
    for(int i=0;i<ncall;i++){
      rocblas_gemm_batched_ex(handle,
			   rocblas_operation_none, rocblas_operation_none,
			   M,N,K, alpha,
			   (const void *)dA, type, M,
			   (const void *)dB, type, K, beta,
			   (const void *)dC, type, M,
			   (void *)dC,       type, M,
			   BATCH, type, algo, solution, 0);
    }
    HIP_CHECK(hipDeviceSynchronize());
    return (usecond()-t0)/ncall;
  };

  printf("\n=== %s : M=%d N=%d K=%d batch=%d ===\n",name,M,N,K,BATCH);

  double tdef = run(rocblas_gemm_algo_standard,0);
  printf("  default heuristic          : %8.1f Gflop/s  %8.1f GB/s  (%.3f ms)\n",
	 flops/tdef/1.0e3, bytes/tdef/1.0e3, tdef/1.0e3);

  // The extended call only accepts indices the library has itself listed for
  // this problem, so ask it for the list rather than guessing.  If this does
  // not compile, check the argument order in rocblas.h: the enumeration call
  // has moved an argument or two between releases.
  rocblas_int nsol = 0;
  rocblas_status qst = rocblas_gemm_batched_ex_get_solutions(handle,
			   rocblas_operation_none, rocblas_operation_none,
			   M,N,K, alpha,
			   (const void *)dA, type, M,
			   (const void *)dB, type, K, beta,
			   (const void *)dC, type, M,
			   (void *)dC,       type, M,
			   BATCH, type, rocblas_gemm_flags_none,
			   NULL, &nsol);
  std::vector<rocblas_int> sols;
  if ( qst == rocblas_status_success && nsol > 0 ) {
    sols.resize(nsol);
    rocblas_gemm_batched_ex_get_solutions(handle,
			   rocblas_operation_none, rocblas_operation_none,
			   M,N,K, alpha,
			   (const void *)dA, type, M,
			   (const void *)dB, type, K, beta,
			   (const void *)dC, type, M,
			   (void *)dC,       type, M,
			   BATCH, type, rocblas_gemm_flags_none,
			   &sols[0], &nsol);
  } else {
    // Older library: fall back to probing indices.
    for(int s=1;s<=maxindex;s++) sols.push_back(s);
  }

  std::vector<Result> ok;
  for(size_t i=0;i<sols.size();i++){
    double t = run(rocblas_gemm_algo_solution_index,sols[i]);
    if ( t < 0 ) continue;                      // library rejects it for this problem
    ok.push_back({(int)sols[i],flops/t/1.0e3,bytes/t/1.0e3});
  }
  printf("  library lists %d solution(s) for this problem\n",(int)sols.size());
  std::sort(ok.begin(),ok.end(),[](const Result &a,const Result &b){ return a.gflops>b.gflops; });

  printf("  %d solutions accepted of %d tried\n",(int)ok.size(),maxindex);
  for(int i=0;i<(int)ok.size() && i<8;i++){
    printf("    solution %5d            : %8.1f Gflop/s  %8.1f GB/s\n",
	   ok[i].index,ok[i].gflops,ok[i].gbs);
  }
  if ( ok.size() ) {
    printf("  best/default               : %.2f\n", ok[0].gflops*tdef*1.0e3/flops);
  }

  HIP_CHECK(hipFree(dA)); HIP_CHECK(hipFree(dB)); HIP_CHECK(hipFree(dC));
  HIP_CHECK(hipFree(A));  HIP_CHECK(hipFree(B));  HIP_CHECK(hipFree(C));
}

int main(int argc,char **argv)
{
  int nbasis = (argc>1) ? atoi(argv[1]) : 60;
  int nrhs   = (argc>2) ? atoi(argv[2]) : 12;
  int batch  = (argc>3) ? atoi(argv[3]) : 1024;
  int maxidx = (argc>4) ? atoi(argv[4]) : 400;

  rocblas_handle handle;
  rocblas_create_handle(&handle);
  rocblas_set_pointer_mode(handle,rocblas_pointer_mode_host);

  printf("rocBLAS solution sweep: the coarse stencil GEMM shape\n");
  printf("nbasis %d, Nrhs %d, batch %d, solution indices 1..%d\n",nbasis,nrhs,batch,maxidx);
  printf("(batch is the LOCAL coarse volume; x9 is the grouped call)\n");

  sweep(handle,rocblas_datatype_f32_c,"complex single",8, nbasis,nrhs,nbasis,batch,maxidx,20);
  sweep(handle,rocblas_datatype_f64_c,"complex double",16,nbasis,nrhs,nbasis,batch,maxidx,20);

  rocblas_destroy_handle(handle);
  return 0;
}
