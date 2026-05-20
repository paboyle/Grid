/*
 * Minimal program demonstrating the workaround for the hipfft ROCm 7 bug.
 *
 * Workaround: create the hipfft plan BEFORE any hipMalloc.  Plan creation
 * for G < 32 then succeeds even with an empty rocFFT cache.
 *
 * Compile:
 *   hipcc -o Test_hipfft_bug_pass Test_hipfft_bug_pass.cc -lhipfft
 *
 * Run:
 *   rm -rf ~/.cache/rocfft
 *   ./Test_hipfft_bug_pass
 *
 * Expected: all G values succeed.
 * Compare with Test_hipfft_bug_fail.cc which uses the opposite ordering.
 */

#include <cstdio>
#include <hipfft/hipfft.h>
#include <hip/hip_runtime.h>

int main(void) {
  hipDeviceProp_t prop;
  hipGetDeviceProperties(&prop, 0);
  printf("Device: %s\n", prop.name);
#ifdef hipfftVersionMinor
  printf("hipFFT version: %d.%d.%d\n\n",
         hipfftVersionMajor, hipfftVersionMinor, hipfftVersionPatch);
#endif

  for (int G : {8, 16, 32}) {
    int howmany = 512;
    int n[] = {G};
    long nelems = (long)G * howmany;

    // Plan created BEFORE hipMalloc — succeeds for all G
    hipfftHandle p;
    size_t workSize = 0;
    hipfftCreate(&p);
    hipfftResult rv = hipfftMakePlanMany(p, 1, n,
                                         nullptr, 1, G, nullptr, 1, G,
                                         HIPFFT_Z2Z, howmany, &workSize);
    printf("G=%-4d  plan-then-hipMalloc: %d (%s)\n",
           G, (int)rv, rv == HIPFFT_SUCCESS ? "HIPFFT_SUCCESS" : "HIPFFT_PARSE_ERROR");

    if (rv == HIPFFT_SUCCESS) {
      hipfftDoubleComplex *buf = nullptr;
      hipMalloc(&buf, nelems * sizeof(hipfftDoubleComplex));
      hipMemset(buf, 0, nelems * sizeof(hipfftDoubleComplex));
      rv = hipfftExecZ2Z(p, buf, buf, HIPFFT_FORWARD);
      hipDeviceSynchronize();
      printf("G=%-4d  execFwd:             %d (%s)\n",
             G, (int)rv, rv == HIPFFT_SUCCESS ? "HIPFFT_SUCCESS" : "FAILED");
      hipFree(buf);
    }

    hipfftDestroy(p);
  }

  return 0;
}
