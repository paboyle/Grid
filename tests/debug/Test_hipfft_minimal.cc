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

// Try both hipfftPlanMany and hipfftCreate+hipfftMakePlanMany for given G and howmany.
static void tryPlans(int G, int howmany) {
  int n[] = {G};

  printf("--- G=%-4d  howmany=%-6d ---\n", G, howmany);

  // 1. hipfftPlanMany (one-step)
  {
    hipfftHandle p;
    hipfftResult rv = hipfftPlanMany(&p, 1, n,
                                     nullptr, 1, G,
                                     nullptr, 1, G,
                                     HIPFFT_Z2Z, howmany);
    printf("  hipfftPlanMany          : %d (%s)\n", (int)rv, hipfftResultString(rv));
    if (rv == HIPFFT_SUCCESS) hipfftDestroy(p);
  }

  // 2. hipfftPlanMany with inembed=n (old Grid behaviour)
  {
    hipfftHandle p;
    hipfftResult rv = hipfftPlanMany(&p, 1, n,
                                     n, 1, G,
                                     n, 1, G,
                                     HIPFFT_Z2Z, howmany);
    printf("  hipfftPlanMany(inembed=n): %d (%s)\n", (int)rv, hipfftResultString(rv));
    if (rv == HIPFFT_SUCCESS) hipfftDestroy(p);
  }

  // 3. hipfftCreate + hipfftMakePlanMany (two-step, nullptr embed)
  {
    hipfftHandle p;
    size_t workSize = 0;
    hipfftResult rc = hipfftCreate(&p);
    printf("  hipfftCreate            : %d (%s)\n", (int)rc, hipfftResultString(rc));
    if (rc == HIPFFT_SUCCESS) {
      hipfftResult rv = hipfftMakePlanMany(p, 1, n,
                                           nullptr, 1, G,
                                           nullptr, 1, G,
                                           HIPFFT_Z2Z, howmany, &workSize);
      printf("  hipfftMakePlanMany      : %d (%s)  workSize=%zu\n",
             (int)rv, hipfftResultString(rv), workSize);
      hipfftDestroy(p);
    }
  }

  // 4. hipfftPlan1d (simplest API)
  {
    hipfftHandle p;
    hipfftResult rv = hipfftPlan1d(&p, G, HIPFFT_Z2Z, howmany);
    printf("  hipfftPlan1d            : %d (%s)\n", (int)rv, hipfftResultString(rv));
    if (rv == HIPFFT_SUCCESS) hipfftDestroy(p);
  }

  printf("\n");
}

int main(void) {
  // Print HIP device info
  int device = 0;
  hipGetDevice(&device);
  hipDeviceProp_t prop;
  hipGetDeviceProperties(&prop, device);
  printf("Device %d: %s  warpSize=%d\n\n", device, prop.name, prop.warpSize);

  // Print hipFFT version if available
#ifdef hipfftVersionMinor
  printf("hipFFT version: %d.%d.%d\n\n",
         hipfftVersionMajor, hipfftVersionMinor, hipfftVersionPatch);
#endif

  // Sweep over transform sizes G with a fixed representative howmany
  // (512 = typical Nperp*Ncomp for a small lattice)
  const int howmany = 512;
  for (int G : {4, 8, 12, 16, 24, 32, 48, 64}) {
    tryPlans(G, howmany);
  }

  // Also try the exact parameters from the failing Grid FFT
  printf("=== Grid-specific parameters ===\n\n");
  tryPlans(8,  512);   // Ls=8,  small 4D lattice
  tryPlans(16, 512);   // 16^4
  tryPlans(32, 512);   // 32^4 (known to work)

  return 0;
}
