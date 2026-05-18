---
name: compiler-validation
description: Identify GPU compiler code generation bugs, distinguish them from hardware and runtime bugs, construct minimal reproducers, and validate correctness of generated assembly for performance-critical HPC kernels.
user-invocable: true
allowed-tools:
  - Read
  - Bash(grep -r)
  - Bash(objdump)
---

# Compiler Validation for GPU HPC Codes

## Why Compiler Bugs Are Distinct

Compiler bugs have a unique diagnostic signature: they produce *deterministically wrong* results. The same input always produces the same wrong output. This distinguishes them from:

- Hardware bugs: usually stochastic (wrong answer sometimes, correct answer other times)
- Runtime bugs (premature barrier, buffer aliasing): often stochastic or history-dependent
- Race conditions: non-deterministic

**The determinism test**: run the same kernel 100 times with the same input. If the wrong answer is always the same wrong answer, suspect the compiler.

## The Minimal Reproducer Protocol

When a kernel produces wrong results, isolate the compiler as quickly as possible:

**Step 1: Eliminate the physics**. Reduce the failing kernel to the smallest possible computation that still exhibits the bug. Replace QCD fields with `double` arrays. Replace lattice operations with scalar arithmetic. The goal is a 20-line CUDA/HIP/SYCL file that any compiler engineer can compile and run.

**Step 2: Binary search over optimisation levels**. Compile at `-O0` (or equivalent). If the answer becomes correct, the bug is in an optimisation pass. Then test `-O1`, `-O2`, `-O3` individually to find which optimisation level introduces the bug.

```bash
# HIP example
hipcc -O0 minimal_repro.cc -o test_O0 && ./test_O0   # should be correct
hipcc -O1 minimal_repro.cc -o test_O1 && ./test_O1   # compare
hipcc -O2 minimal_repro.cc -o test_O2 && ./test_O2   # compare
```

**Step 3: Identify the optimisation pass**. For LLVM-based compilers (clang, hipcc, dpcpp, nvcc via ptxas):
```bash
# Disable individual optimisation passes:
hipcc -O2 -mllvm -disable-loop-unrolling minimal_repro.cc -o test
hipcc -O2 -fno-vectorize minimal_repro.cc -o test
hipcc -O2 -fno-slp-vectorize minimal_repro.cc -o test
```

**Step 4: Inspect the generated code**. For CUDA/HIP, use `--generate-line-info` and `cuobjdump` or `roc-obj-extract` to get annotated assembly:
```bash
# CUDA
nvcc -O2 --generate-line-info --keep minimal_repro.cu
cuobjdump --dump-ptx minimal_repro.o
# HIP/ROCm
hipcc -O2 --save-temps minimal_repro.cc
llvm-objdump -d minimal_repro.o
# SYCL/DPC++
icpx -O2 -fsycl -Xclang -ast-dump minimal_repro.cc 2>&1 | grep -A5 "suspicious_expr"
```

Look for: incorrect register spill/fill sequences, loop trip count miscalculation, vectorisation across iteration boundaries, incorrect address arithmetic.

## Known Compiler Bug Patterns in GPU Code

### Register Pressure / Spill Bugs
High register usage forces spills to local memory. Some compiler versions generate incorrect spill/fill code — the value is written to local memory but a stale register value is read back instead of the spilled value.

**Signature**: Wrong answer with high-register-count kernels; becomes correct when `--maxrregcount=N` forces lower register count (more spilling) or higher (`--maxrregcount=256`, fewer spills).

**Diagnostic**: Check register usage:
```bash
nvcc -O2 --ptxas-options=-v minimal_repro.cu 2>&1 | grep "registers"
hipcc -O2 --offload-arch=gfx90a --save-temps minimal_repro.cc
llvm-mc --arch=amdgcn minimal_repro.s 2>&1 | grep "VGPRs"
```

### Vectorisation Across Loop Boundaries
The compiler vectorises two successive loop iterations as a SIMD unit when they have a data dependency that the compiler has incorrectly determined does not exist.

**Signature**: Wrong answer that becomes correct when the loop body is extracted to a non-inlined function (disabling auto-vectorisation across iterations).

### Incorrect Constant Propagation
The compiler evaluates a compile-time expression incorrectly, substituting a wrong constant. Common in template-heavy code where `sizeof(T)` or `alignof(T)` is used in arithmetic that the compiler folds at compile time.

**Signature**: Wrong array index or wrong stride. Inspecting the generated assembly shows a literal constant where you expect a computed value.

## Stress Patterns for Compiler Validation

These patterns exercise the compiler in ways that commonly expose bugs:

```cpp
// 1. Aliased pointer write followed by immediate read
// (tests correct handling of write-after-write in register allocation)
__global__ void alias_stress(double *a, double *b, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        a[i] = a[i] * 2.0;
        b[i] = a[i] + 1.0;  // must read the updated value, not the original
    }
}

// 2. Mixed-precision accumulation
// (tests correct type promotion in FMA sequences)
__global__ void precision_stress(float *in, double *out, int n) {
    double acc = 0.0;
    for (int i = 0; i < n; i++) acc += (double)in[i];
    *out = acc;
}

// 3. Large struct in shared memory
// (tests alignment and offset calculation for non-power-of-2-sized objects)
struct S { double x[3]; };  // sizeof = 24 bytes, not a power of 2
__global__ void struct_stress(S *in, S *out, int n) {
    extern __shared__ S smem[];
    int tid = threadIdx.x;
    smem[tid] = in[tid];
    __syncthreads();
    out[tid] = smem[(tid + 1) % blockDim.x];
}
```

## Separating Compiler from Runtime/Hardware

When results are deterministically wrong:

| Test | Compiler bug | Runtime/hardware bug |
|---|---|---|
| Recompile at -O0 | Fixes it | No effect |
| Run on CPU (host code equivalent) | Fixes it | No effect |
| Reorder loop iterations | Changes wrong answer | No effect or different pattern |
| Different compiler version | Fixes or changes wrong answer | No effect |
| Different GPU of same model | Same wrong answer | Different or no error |
| Different GPU model | Fixes it (ISA-specific codegen bug) | May or may not fix |

## Reporting to Compiler Teams

A compiler bug report needs:
1. Minimal reproducer (< 50 lines)
2. Compiler version (`hipcc --version`, `nvcc --version`, `icpx --version`)
3. GPU model and driver version
4. Exact wrong and correct answers (hexfloat for reproducibility)
5. Which compile flags change the behaviour
6. Generated assembly for the correct and incorrect variants

File with: LLVM Bugzilla (for hipcc/clang/dpcpp backends), NVIDIA bug portal (nvcc/ptxas), or vendor-specific developer forum. The minimal reproducer is the single most important element — without it, compiler teams cannot prioritise.

## Pragmatic In-Production Workaround

When a compiler bug is confirmed but the fix is not yet available, the lowest-risk workaround is to mark the affected function with reduced optimisation:

```cpp
#pragma clang optimize off  // clang/hipcc/dpcpp
void __attribute__((optimize("O0"))) affected_kernel_host_wrapper() { ... }
// For device code, use per-file compilation flags via CMake/Makefile
```

Document the workaround with a comment referencing the compiler bug report number so it can be removed when the compiler is updated.
