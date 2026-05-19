---
name: gpu-memory-performance
description: Diagnose and fix GPU memory bandwidth and occupancy problems in Grid HPC kernels — acceleratorThreads() pitfalls, LambdaApply thread mapping, coalescedRead/Write idiom, when to use accelerator_for vs a hand-rolled __global__ kernel, and fused vs staged HBM access patterns.
user-invocable: true
allowed-tools:
  - Read
  - Bash(grep -r)
---

# GPU Memory Performance in Grid

## Nsimd on GPU builds

With `GEN_SIMD_WIDTH=64B` (the typical production setting), `Nsimd` is **not 1**:

| Scalar type | `sizeof` | `Nsimd = 64 / sizeof` |
|---|---|---|
| `ComplexD` | 16 B | **4** |
| `ComplexF` | 8 B  | **8** |
| `RealD`    | 8 B  | **8** |
| `RealF`    | 4 B  | **16** |

So for `LatticePropagatorD` (scalar type `ComplexD`), `Nsimd=4` and the SIMD lane runs `threadIdx.x ∈ {0,1,2,3}`.  True `Nsimd=1` scalar-GPU builds are the exception, not the rule.

## The acceleratorThreads() Trap

`acceleratorThreads()` is a runtime-settable global (default **8**) that controls the `blockDim.y` of every `accelerator_for` launch.  It is NOT the SIMD width — it is the number of sites processed per block in the y-dimension.

```cpp
// Grid/threads/Accelerator.cc
uint32_t accelerator_threads = 8;
```

With `accelerator_for(ss, osites, nsimd, ...)`, the launch is:

```
dim3 threads(nsimd, acceleratorThreads(), 1)
dim3 blocks ((osites + acceleratorThreads() - 1) / acceleratorThreads(), 1, 1)
```

Total threads per block = `Nsimd × acceleratorThreads()`.  With `GEN_SIMD_WIDTH=64B` and `Nsimd=8` (fp32 / ComplexF):

| `acceleratorThreads()` | threads/block | AMD wavefront | note |
|---|---|---|---|
| 2 (original default) | 16 | 25% | sub-wavefront, poor |
| 4 | 32 | 50% | half wavefront |
| **8 (current default)** | **64** | **100%** | **one full wavefront — Dslash sweet spot** |
| 16 | 128 | 200% | two wavefronts/block; register-pressure cliff for heavy kernels |
| 32 | 256 | 400% | severe register spill for stencil kernels |

**Why 8 and not higher?**  Compute-heavy kernels like the Domain Wall Dslash carry many live registers per thread (spinors + gauge links + projections).  Doubling `acceleratorThreads` from 8→16 doubles the register demand per block, which on AMD GFX90A triggers a hard occupancy cliff: `Benchmark_dwf_fp32` drops from 1.7 TF/s (nt=8) to ~300 GF/s (nt=16).  The sweet spot is one full wavefront per block, which with `Nsimd=8` (fp32) means `nt=8`.

For fp64 work (`Nsimd=4`), `nt=8` gives 32 threads = half a wavefront (AMD pads to 64 with idle lanes).  Kernels that are not register-limited (e.g. simple lattice arithmetic) can benefit from `--accelerator-threads 16` at runtime.  Reduction kernels bypass `acceleratorThreads()` entirely via `getNumBlocksAndThreads`.

**Why was the default ever 2?**  Before the `threadIdx.x`/`threadIdx.y` remap (see LambdaApply section below), the site index lived in `threadIdx.x` — the fast, coalescing dimension.  Increasing `acceleratorThreads` widened the block in the *site* direction, so adjacent threads in a warp hit adjacent sites, each stride `sizeof(vobj)` apart in AoS memory — breaking coalescing.  On early NVIDIA ports with `Nsimd≈8`, `nt=2` gave 16 threads = 50% of a 32-thread warp; NVIDIA recovers this via multiple concurrent blocks per SM, so occupancy was barely tolerable.  AMD has no such multiplier when blocks are already sub-wavefront.  After the remap put the site index in `threadIdx.y` and the SIMD lane in `threadIdx.x`, coalescing became independent of `acceleratorThreads`, removing the constraint.

**Diagnostic**: observed bandwidth << peak, kernel time >> expected from data volume.  Check with `--accelerator-threads 32` at runtime.  A large speedup confirms occupancy starvation.

**Fix options** (in order of preference):
1. Kernel needs its own thread count — use `getNumBlocksAndThreads` and launch a `__global__` kernel directly (see below).
2. Temporarily acceptable: set `--accelerator-threads 32` at the application level.  Note this affects every `accelerator_for` site in the binary.

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

With `GEN_SIMD_WIDTH=64B` and `Nsimd=4` (PropagatorD), a 64-thread AMD wavefront contains 64/4 = 16 sites, each processed by 4 lanes.  Adjacent threads within the wavefront read different lanes of the same site — this is a broadcast pattern (hardware handles this efficiently), not a stride.  The stride between consecutive *sites* (`sizeof(vobj)`) only appears between groups of `Nsimd` threads, spaced `Nsimd` apart in threadIdx.x — not between adjacent threads within a warp.

## coalescedRead / coalescedWrite

These are Grid's canonical way to read/write one SIMD lane from a vector type inside a `GRID_SIMT` kernel:

```cpp
// accelerator_for(ss, osites, Nsimd, {
//   lane = acceleratorSIMTlane(Nsimd) = threadIdx.x  ∈ {0..Nsimd-1}
auto scalar_val = coalescedRead(field[ss]);       // extractLane(lane, field[ss])
coalescedWrite(field[ss], scalar_val);             // insertLane(lane, field[ss], scalar_val)
```

For `vobj` aggregate types, `coalescedRead` calls `extractLane(lane, vobj)` which recurses through the tensor hierarchy and returns `vobj::scalar_object`.

For `vsimd` (raw SIMD vector) types, it casts to `scalar_type*` and indexes with `lane`.

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
// gives 64–256 threads/block → correct occupancy independent of acceleratorThreads()
Integer smemSize = numThreads * sizeof(sobj);
myKernel<<<numBlocks, numThreads, smemSize, computeStream>>>(args...);
```

Grid's `reduceKernel` uses this pattern and achieves ~400 GB/s on MI250X.

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

Launched with `getNumBlocksAndThreads` → 128 threads/block for R=12 (`BundleScalarD`=192 B, sharedMem=64 KB) → correct occupancy without depending on `acceleratorThreads()`.

## Observed Numbers on MI250X (32^4 LatticePropagatorD, Nsimd=4, GEN_SIMD_WIDTH=64B)

| Configuration | pack µs/group | reduce µs/group | total µs | GB/s |
|---|---|---|---|---|
| acceleratorThreads=2 (8 threads/block), staged | 10,080 | 470 | 126,909 | 50 |
| acceleratorThreads=16 (64 threads/block), staged | 342 | 310 | 8,251 | 297 |
| acceleratorThreads=16, fused (128 threads/block via getNumBlocksAndThreads) | — | 349 | 4,584 | 546 |

The fused kernel at 349 µs/group reads 201 MB at 576 GB/s — 36% of MI250X HBM peak.  The remaining gap from peak is the in-kernel serial loop over R=12 words and the 12 serial kernel launches.

## Quick Checklist When a Kernel Is Slow

1. **Check Nsimd**: `GEN_SIMD_WIDTH=64B` → Nsimd=4 (ComplexD), 8 (ComplexF).  Total threads/block = `Nsimd × acceleratorThreads()`.  With old default nt=2 and Nsimd=4: 8 threads = 12.5% of AMD wavefront.
2. Check threads per block: for `accelerator_for` kernels use `--accelerator-threads 32` and measure; a large speedup confirms occupancy starvation.
3. Check for bulk struct accumulation in registers (`Bundle b; for(...) b._internal[k] = ...;`).  Replace with per-element writes via `coalescedWrite`.
4. Check for staged HBM access (pack → buffer → reduce).  Count the passes; fuse if ≥ 2 passes over the same data.
5. For reduction kernels, always use `getNumBlocksAndThreads` rather than `accelerator_for` so thread count is independent of `acceleratorThreads()`.
