/*
 * Minimal reproducer for hipfftMakePlanMany / hipfftPlanMany failures.
 *
 * Compile on Frontier (no Grid headers needed):
 *   hipcc -o Test_hipfft_minimal Test_hipfft_minimal.cc -lhipfft
 *
 * Run:
 *   ./Test_hipfft_minimal
 */

#include <cstdio>
#include <cstdlib>
#include <hipfft/hipfft.h>
#include <hip/hip_runtime.h>

static const char *hipfftResultString(hipfftResult r) {
  switch (r) {
  case HIPFFT_SUCCESS:                  return "HIPFFT_SUCCESS";
  case HIPFFT_INVALID_PLAN:             return "HIPFFT_INVALID_PLAN";
  case HIPFFT_ALLOC_FAILED:             return "HIPFFT_ALLOC_FAILED";
  case HIPFFT_INVALID_TYPE:             return "HIPFFT_INVALID_TYPE";
  case HIPFFT_INVALID_VALUE:            return "HIPFFT_INVALID_VALUE";
  case HIPFFT_INTERNAL_ERROR:           return "HIPFFT_INTERNAL_ERROR";
  case HIPFFT_EXEC_FAILED:              return "HIPFFT_EXEC_FAILED";
  case HIPFFT_SETUP_FAILED:             return "HIPFFT_SETUP_FAILED";
  case HIPFFT_INVALID_SIZE:             return "HIPFFT_INVALID_SIZE";
  case HIPFFT_UNALIGNED_DATA:           return "HIPFFT_UNALIGNED_DATA";
  case HIPFFT_INCOMPLETE_PARAMETER_LIST:return "HIPFFT_INCOMPLETE_PARAMETER_LIST";
  case HIPFFT_INVALID_DEVICE:           return "HIPFFT_INVALID_DEVICE";
  case HIPFFT_PARSE_ERROR:              return "HIPFFT_PARSE_ERROR";
  case HIPFFT_NO_WORKSPACE:             return "HIPFFT_NO_WORKSPACE";
  case HIPFFT_NOT_IMPLEMENTED:          return "HIPFFT_NOT_IMPLEMENTED";
  case HIPFFT_NOT_SUPPORTED:            return "HIPFFT_NOT_SUPPORTED";
  default:                              return "UNKNOWN";
  }
}

// Plan creation + execution for (G, howmany) using hipfftCreate+hipfftMakePlanMany.
// This is the path Grid's FFT.h now uses.
static void tryPlanAndExec(int G, long howmany) {
  int n[] = {G};
  long nelems = (long)G * howmany;

  printf("--- G=%-4d  howmany=%-10ld  total_elems=%-12ld ---\n",
         G, howmany, nelems);

  // Allocate device buffer (hipfftDoubleComplex = 16 bytes each)
  hipfftDoubleComplex *dbuf = nullptr;
  hipError_t herr = hipMalloc(&dbuf, nelems * sizeof(hipfftDoubleComplex));
  if (herr != hipSuccess) {
    printf("  hipMalloc failed (%d) for %ld elems — skipping\n\n", (int)herr, nelems);
    return;
  }
  hipMemset(dbuf, 0, nelems * sizeof(hipfftDoubleComplex));

  // 1. hipfftPlanMany (one-step, nullptr embed) — current Grid path
  {
    hipfftHandle p;
    hipfftResult rv = hipfftPlanMany(&p, 1, n,
                                     nullptr, 1, G,
                                     nullptr, 1, G,
                                     HIPFFT_Z2Z, (int)howmany);
    printf("  hipfftPlanMany   create : %d (%s)\n", (int)rv, hipfftResultString(rv));
    if (rv == HIPFFT_SUCCESS) {
      rv = hipfftExecZ2Z(p, dbuf, dbuf, HIPFFT_FORWARD);
      hipDeviceSynchronize();
      printf("  hipfftPlanMany   execFwd: %d (%s)\n", (int)rv, hipfftResultString(rv));
      hipfftDestroy(p);
    }
  }

  // 2. hipfftCreate + hipfftMakePlanMany (two-step) — also current Grid path
  {
    hipfftHandle p;
    size_t workSize = 0;
    hipfftResult rc = hipfftCreate(&p);
    if (rc == HIPFFT_SUCCESS) {
      hipfftResult rv = hipfftMakePlanMany(p, 1, n,
                                           nullptr, 1, G,
                                           nullptr, 1, G,
                                           HIPFFT_Z2Z, (int)howmany, &workSize);
      printf("  hipfftMakePlanMany      : %d (%s)  workSize=%zu\n",
             (int)rv, hipfftResultString(rv), workSize);
      if (rv == HIPFFT_SUCCESS) {
        rv = hipfftExecZ2Z(p, dbuf, dbuf, HIPFFT_FORWARD);
        hipDeviceSynchronize();
        printf("  hipfftMakePlanMany exec : %d (%s)\n", (int)rv, hipfftResultString(rv));
      }
      hipfftDestroy(p);
    } else {
      printf("  hipfftCreate            : %d (%s)\n", (int)rc, hipfftResultString(rc));
    }
  }

  // 3. hipfftPlan1d (simplest API, batch = howmany)
  {
    hipfftHandle p;
    hipfftResult rv = hipfftPlan1d(&p, G, HIPFFT_Z2Z, (int)howmany);
    printf("  hipfftPlan1d     create : %d (%s)\n", (int)rv, hipfftResultString(rv));
    if (rv == HIPFFT_SUCCESS) {
      rv = hipfftExecZ2Z(p, dbuf, dbuf, HIPFFT_FORWARD);
      hipDeviceSynchronize();
      printf("  hipfftPlan1d     execFwd: %d (%s)\n", (int)rv, hipfftResultString(rv));
      hipfftDestroy(p);
    }
  }

  hipFree(dbuf);
  printf("\n");
}

int main(void) {
  // Print HIP device info
  int device = 0;
  hipGetDevice(&device);
  hipDeviceProp_t prop;
  hipGetDeviceProperties(&prop, device);
  printf("Device %d: %s  warpSize=%d\n\n", device, prop.name, prop.warpSize);

#ifdef hipfftVersionMinor
  printf("hipFFT version: %d.%d.%d\n\n",
         hipfftVersionMajor, hipfftVersionMinor, hipfftVersionPatch);
#endif

  // Original sweep with small howmany (these passed first time)
  printf("=== Small howmany (original sweep) ===\n\n");
  for (int G : {4, 8, 12, 16, 24, 32, 48, 64})
    tryPlanAndExec(G, 512);

  // Grid-realistic howmany values derived from actual lattice geometries.
  // howmany = Ncomp * product(ldimensions[d] for d != dim)
  // For LatticeComplexD: Ncomp=1.
  printf("=== Grid-realistic parameters ===\n\n");

  // --grid 16.16.16.16  4D FFT  (KNOWN TO FAIL in Grid)
  // Each dim: G=16, Nperp=16^3=4096
  tryPlanAndExec(16, 4096);

  // --grid 32.32.32.32  4D FFT  (KNOWN TO SUCCEED in Grid)
  // Each dim: G=32, Nperp=32^3=32768
  tryPlanAndExec(32, 32768);

  // --grid 32.32.32.32 Ls=8  5D DWF FFT  (KNOWN TO FAIL on dim 0 in Grid)
  // dim 0: G=8,  Nperp=32^4=1048576
  tryPlanAndExec(8, 1048576);
  // dim 1-4: G=32, Nperp=8*32^3=262144
  tryPlanAndExec(32, 262144);

  // Extra intermediate cases to bracket the failure
  tryPlanAndExec(16, 1024);
  tryPlanAndExec(16, 2048);
  tryPlanAndExec(16, 8192);
  tryPlanAndExec(8,  4096);
  tryPlanAndExec(8,  65536);
  tryPlanAndExec(8,  262144);

  return 0;
}
