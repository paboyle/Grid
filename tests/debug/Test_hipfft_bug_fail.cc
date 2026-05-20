/*
 * Minimal program demonstrating hipfft HIPFFT_PARSE_ERROR on ROCm 7 / hipFFT 1.0.20.
 *
 * Bug: hipfftMakePlanMany returns HIPFFT_PARSE_ERROR (12) for transform sizes
 * smaller than 32 when a hipMalloc is issued before plan creation.
 *
 * Compile:
 *   hipcc -o Test_hipfft_bug_fail Test_hipfft_bug_fail.cc -lhipfft
 *
 * Run with empty rocFFT cache:
 *   rm -rf ~/.cache/rocfft
 *   ./Test_hipfft_bug_fail
 *
 * Expected (broken): hipfftMakePlanMany returns 12 (HIPFFT_PARSE_ERROR) for G < 32.
 * See Test_hipfft_bug_pass.cc for the workaround.
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

    // hipMalloc BEFORE plan creation — triggers HIPFFT_PARSE_ERROR for G < 32
    hipfftDoubleComplex *buf = nullptr;
    hipMalloc(&buf, nelems * sizeof(hipfftDoubleComplex));

    hipfftHandle p;
    size_t workSize = 0;
    hipfftCreate(&p);
    hipfftResult rv = hipfftMakePlanMany(p, 1, n,
                                         nullptr, 1, G, nullptr, 1, G,
                                         HIPFFT_Z2Z, howmany, &workSize);
    printf("G=%-4d  hipMalloc-then-plan: %d (%s)\n",
           G, (int)rv, rv == HIPFFT_SUCCESS ? "HIPFFT_SUCCESS" : "HIPFFT_PARSE_ERROR");
    hipfftDestroy(p);
    hipFree(buf);
  }

  return 0;
}
