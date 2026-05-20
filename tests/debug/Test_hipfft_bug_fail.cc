/*
 * Isolating the hipfft HIPFFT_PARSE_ERROR on ROCm 7 / hipFFT 1.0.20.
 *
 * Tests three orderings with an empty rocFFT cache to find which GPU
 * operation before plan creation triggers the failure:
 *   A) hipMalloc only            — hypothesis: passes (no async GPU work)
 *   B) hipMalloc + hipMemset     — hypothesis: fails  (async GPU work in flight)
 *   C) hipMalloc + hipMemset     — hypothesis: passes (work completed before plan)
 *      + hipDeviceSynchronize
 *
 * Compile:
 *   hipcc -o Test_hipfft_bug_fail Test_hipfft_bug_fail.cc -lhipfft
 *
 * Run with empty cache:
 *   rm -rf ~/.cache/
 *   ./Test_hipfft_bug_fail
 */

#include <cstdio>
#include <hipfft/hipfft.h>
#include <hip/hip_runtime.h>

static const char *res(hipfftResult rv) {
  return rv == HIPFFT_SUCCESS ? "SUCCESS" : "PARSE_ERROR";
}

static hipfftResult makePlan(int G, int howmany) {
  int n[] = {G};
  hipfftHandle p;
  size_t workSize = 0;
  hipfftCreate(&p);
  hipfftResult rv = hipfftMakePlanMany(p, 1, n,
                                       nullptr, 1, G, nullptr, 1, G,
                                       HIPFFT_Z2Z, howmany, &workSize);
  hipfftDestroy(p);
  return rv;
}

int main(void) {
  hipDeviceProp_t prop;
  hipGetDeviceProperties(&prop, 0);
  printf("Device: %s\n", prop.name);
#ifdef hipfftVersionMinor
  printf("hipFFT version: %d.%d.%d\n\n",
         hipfftVersionMajor, hipfftVersionMinor, hipfftVersionPatch);
#endif

  for (int G : {4, 8, 16, 32}) {
    int howmany = 512;
    long nelems = (long)G * howmany;
    hipfftDoubleComplex *buf = nullptr;
    hipMalloc(&buf, nelems * sizeof(hipfftDoubleComplex));

    // Tests ordered so each runs before a prior success can populate the cache.

    // B first: hipMalloc + hipMemset (async GPU work in flight)
    // If this fails, A (no hipMemset) will pass, confirming hipMemset is the trigger.
    hipMemset(buf, 0, nelems * sizeof(hipfftDoubleComplex));
    hipfftResult rvB = makePlan(G, howmany);
    printf("G=%-4d  B) hipMalloc + hipMemset       : %s\n", G, res(rvB));

    // C: hipMalloc + hipMemset + sync — does syncing before plan creation fix it?
    hipMemset(buf, 0, nelems * sizeof(hipfftDoubleComplex));
    hipDeviceSynchronize();
    hipfftResult rvC = makePlan(G, howmany);
    printf("G=%-4d  C) hipMalloc + hipMemset + sync: %s\n", G, res(rvC));

    // A last: hipMalloc only, no async GPU work — should always pass
    hipfftResult rvA = makePlan(G, howmany);
    printf("G=%-4d  A) hipMalloc only             : %s\n\n", G, res(rvA));

    hipFree(buf);
  }

  return 0;
}
