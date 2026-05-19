---
name: gpu-memory-performance
description: Diagnose and fix GPU memory bandwidth and occupancy problems in Grid HPC kernels — acceleratorThreads() pitfalls, LambdaApply thread mapping, coalescedRead/Write idiom, when to use accelerator_for vs a hand-rolled __global__ kernel, and fused vs staged HBM access patterns.
user-invocable: true
allowed-tools:
  - Read
  - Bash(grep -r)
---

# GPU Memory Performance in Grid

## The acceleratorThreads() Trap

`acceleratorThreads()` is a runtime-settable global (default **2**) that controls the `blockDim.y` of every `accelerator_for` launch.  It is NOT the SIMD width — it is the number of sites processed per block in the y-dimension.

```cpp
// Grid/threads/Accelerator.cc
uint32_t accelerator_threads = 2;   // <-- default
```

With `accelerator_for(ss, osites, nsimd, ...)`, the launch is:

```
dim3 threads(nsimd, acceleratorThreads(), 1)
dim3 blocks ((osites + acceleratorThreads() - 1) / acceleratorThreads(), 1, 1)
```

For `nsimd=1` and the default `acceleratorThreads()=2`:
- **2 threads per block** on a 64-thread AMD wavefront → **3% occupancy**
- Expected bandwidth ≈ peak × 3% ≈ 50 GB/s on MI250X

**Diagnostic**: observed bandwidth << peak, kernel time >> expected from data volume.  Check with `--accelerator-threads 16` or `--accelerator-threads 32` at runtime.  A large speedup confirms occupancy starvation.

**Fix options** (in order of preference):
1. Kernel needs its own thread count — use `getNumBlocksAndThreads` and launch a `__global__` kernel directly (see below).
2. Temporarily acceptable: set `--accelerator-threads 16` or 32 at the application level.  Note this affects every `accelerator_for` site in the binary.

## LambdaApply Thread Mapping

`accelerator_for` and `accelerator_for2d` go through `LambdaApply`:

```cpp
// HIP/CUDA LambdaApply kernel:
uint64_t x = threadIdx.y + blockDim.y * blockIdx.x;  // iter1 (site index)
uint64_t y = threadIdx.z + blockDim.z * blockIdx.y;  // iter2
uint64_t z = threadIdx.x;                             // lane (SIMD lane)
Lambda(x, y, z);
```

`threadIdx.x` is the **fast** (lane) dimension — consecutive thread IDs within a warp/wavefront correspond to consecutive lane values on the **same** site, not consecutive sites.

Consequence: for coalesced access from a `vobj` array (AoS layout, stride = `sizeof(vobj)` between adjacent sites), adjacent threads in a wavefront address the **same** site at different lanes, not adjacent sites.  With `Nsimd=1` (GPU scalar build), `threadIdx.x` is always 0 and provides no coalescing benefit at all.

## coalescedRead / coalescedWrite

These are Grid's canonical way to read/write one SIMD lane from a vector type inside a `GRID_SIMT` kernel:

```cpp
// accelerator_for(ss, osites, Nsimd, {
//   lane = acceleratorSIMTlane(Nsimd) = threadIdx.x
auto scalar_val = coalescedRead(field[ss]);       // extractLane(lane, field[ss])
coalescedWrite(field[ss], scalar_val);             // insertLane(lane, field[ss], scalar_val)
```

For `vobj` aggregate types, `coalescedRead` calls `extractLane(lane, vobj)` which recurses through the tensor hierarchy and returns `vobj::scalar_object`.

For `vsimd` (raw SIMD vector) types, it casts to `scalar_type*` and indexes with `lane`.

**When Nsimd=1** (GPU scalar build): `lane=0` always, so `coalescedRead`/`coalescedWrite` are effectively no-ops (direct read/write).  Coalescing must be achieved through the iteration structure instead.

## Coalescing the Iteration Structure

For an AoS input array where each site is `words` 16-byte elements, adjacent threads reading the same site's consecutive words achieve coalesced access:

```cpp
// Good: k varies across threads in a block → consecutive 16-byte reads
accelerator_for2d(k, R, ss, osites, Nsimd, {
    coalescedWrite(out[ss]._internal[k],
                   coalescedRead(idat[ss * words + base + k]));
});
// dim3(Nsimd, nt, 1): threadIdx.y = k (consecutive words, coalesced)
//                      threadIdx.x = lane (SIMD sub-lane, coalesced for Nsimd>1)
```

```cpp
// Bad: each thread reads all R words of its site serially
accelerator_for(ss, osites, 1, {
    Bundle b;
    for (int k = 0; k < R; k++)
        b._internal[k] = idat[ss * words + base + k];  // serial, not coalesced across threads
    out[ss] = b;                                         // bulk struct write
});
```

The bad pattern also accumulates a large struct in registers (192 bytes for R=12), increasing register pressure and reducing occupancy further.

## When to Use a __global__ Kernel Instead of accelerator_for

`accelerator_for` is correct for site-parallel work where `acceleratorThreads()` is tuned appropriately.  Use a direct `__global__` kernel when:

- The kernel requires a **specific thread count** for correctness or performance (reductions, shared-memory algorithms).
- The optimal thread count depends on `sizeof(sobj)` and `sharedMemPerBlock`, not on a runtime global.
- You need the retirement-count pattern for cross-block final reduction.

Pattern: use `getNumBlocksAndThreads` to pick `numThreads` and `numBlocks`:

```cpp
Integer numThreads, numBlocks;
int ok = getNumBlocksAndThreads(n, sizeof(sobj), numThreads, numBlocks);
// starts at warpSize (32/64), doubles while 2*threads*sizeof(sobj) < sharedMemPerBlock
// gives 64–256 threads/block → near-100% wavefront occupancy
Integer smemSize = numThreads * sizeof(sobj);
myKernel<<<numBlocks, numThreads, smemSize, computeStream>>>(args...);
```

This gives 64–256 threads/block regardless of `acceleratorThreads()`.  Grid's `reduceKernel` uses this pattern and achieves ~400 GB/s on MI250X.

## Fused vs Staged HBM Access

A staged pack+reduce reads the data **three times**:

```
pack kernel:   reads vobj array (N bytes), writes bundle buffer (N bytes)
reduce kernel: reads bundle buffer (N bytes), writes tiny result buffer
```

Total HBM: 3N bytes for N bytes of useful input.

A fused kernel reads the data **once**:

```
packReduceKernel: reads R words of vobj array (N bytes), reduces in-place
```

Total HBM: N bytes.  Register pressure increases (R words held per thread) but the 3× HBM saving dominates for large objects.

The fused pattern in Grid's `sumD_gpu_reduce_words<R>`:

```cpp
template <int R, class vobj, class sobj, class Iterator>
__device__ void packReduceBlocks(
    const iScalar<typename vobj::vector_type> *idat,
    sobj *g_odata, Iterator osites, int base, int words)
{
  // sobj = iVector<iScalar<scalarD>, R>  (R double-precision scalars per site)
  constexpr Iterator nsimd = vobj::Nsimd();
  ...
  while (i < osites * nsimd) {
    Iterator lane = i % nsimd;
    Iterator ss   = i / nsimd;
    sobj tmpD; zeroit(tmpD);
    for (int k = 0; k < R; k++) {
      auto w = extractLane(lane, idat[ss * words + base + k]);
      iScalar<typename vobj::scalar_typeD> wd; wd = w;   // float→double promotion
      tmpD._internal[k] = wd;
    }
    mySum += tmpD;
    ...
  }
  reduceBlock(sdata, mySum, tid);
}
```

Launched with `getNumBlocksAndThreads` → 128–256 threads/block → correct occupancy without depending on `acceleratorThreads()`.

## Observed Numbers on MI250X (32^4 LatticePropagatorD, Nsimd=1)

| Configuration | pack µs/group | reduce µs/group | total µs | GB/s |
|---|---|---|---|---|
| acceleratorThreads=2, staged | 10,080 | 470 | 126,909 | 50 |
| acceleratorThreads=16, staged | 342 | 310 | 8,251 | 297 |
| acceleratorThreads=16, fused | — | 349 | 4,584 | 546 |

The fused kernel at 349 µs/group reads 201 MB at 576 GB/s — 36% of MI250X HBM peak.  The remaining gap from peak is the in-kernel serial loop over R=12 words and the 12 serial kernel launches.

## Quick Checklist When a Kernel Is Slow

1. Check threads per block: `accelerator_for(ss, N, 1, ...)` with default `acceleratorThreads()=2` = 2 threads/block = 3% occupancy on AMD.  Try `--accelerator-threads 16` at runtime; if it helps a lot, occupancy is the problem.
2. Check for bulk struct accumulation in registers (`Bundle b; for(...) b._internal[k] = ...;`).  Replace with per-element writes via `coalescedWrite`.
3. Check for staged HBM access (pack → buffer → reduce).  Count the passes; fuse if ≥ 2 passes over the same data.
4. For reduction kernels, always use `getNumBlocksAndThreads` rather than `accelerator_for` so thread count is independent of `acceleratorThreads()`.
