---
name: gpu-runtime-correctness
description: Detect and work around GPU runtime correctness failures — premature completion signalling, infinite poll hangs, stale completion flags, and the double-wait diagnostic pattern. Covers CUDA, HIP/ROCm, and SYCL/Level Zero runtimes.
user-invocable: true
allowed-tools:
  - Read
  - Bash(grep -r)
---

# GPU Runtime Correctness

## The Completion Signalling Problem

GPU runtimes expose a synchronisation primitive — `cudaDeviceSynchronize()`, `hipDeviceSynchronize()`, `q.wait()` — that is supposed to block until all previously submitted GPU work is complete. On several production systems, this guarantee has been violated in two distinct ways:

### Failure Mode A: Premature Return
The wait returns before the GPU work is done. The subsequent CPU code reads stale data from the output buffer. This is the most dangerous failure because it looks like a numerical instability, not a crash. Results are wrong but the program exits normally.

**Identifying Premature Return**: Insert a second, independent wait immediately after the first. If a second `q.wait()` "fixes" incorrect results that appeared with a single `q.wait()`, the first wait was returning prematurely.

```cpp
// Diagnostic version — if this stabilises results, you have premature return
accelerator_barrier();  // first wait
accelerator_barrier();  // second wait (diagnostic)
```

Production fix: submit a trivially cheap no-op kernel after the real work and wait for it. The no-op kernel cannot complete until all previous commands in the queue are done (command queue ordering guarantee), so waiting for the no-op is a stronger barrier than waiting for the queue itself:

```cpp
// Lightweight fence kernel
template<class T>
__global__ void noop_kernel(T *p) { if (threadIdx.x == 0) (void)(*p); }

void strong_barrier(T *device_ptr) {
    noop_kernel<<<1, 1, 0, computeStream>>>(device_ptr);
    cudaStreamSynchronize(computeStream);  // wait for the no-op
}
```

### Failure Mode B: Infinite Poll
The wait enters a polling loop that never terminates. The process consumes 100% CPU in a runtime library. The GPU has either stopped signalling progress entirely, or the completion flag is in a memory region that has become incoherent.

This is distinct from Failure Mode A: with premature return the CPU proceeds; with infinite poll the CPU is stuck.

**Identifying Infinite Poll**: `top` shows the MPI rank at 100% CPU. `perf top -p PID` or `strace -p PID` shows the process burning cycles inside the GPU runtime library (e.g. `libze_intel_gpu.so`, `libamdhip64.so`).

**Documented instances**: 
- Intel Level Zero on Pontevecchio (Aurora): both premature return *and* infinite poll have been observed as independent bugs on the same system.
- The two failure modes can co-exist and have overlapping symptoms at the application level.

## Completion Signalling Architecture

Understanding why these bugs happen requires knowing how completion signalling works:

```
GPU command processor
  → signals completion by writing to a host-visible memory address
  → CPU runtime polls that address (or uses OS event notification via ioctl)
```

A premature return means the memory write happened before the actual work completed (e.g. the signal is on a different command stream that has not been serialised with the work stream). An infinite poll means the memory write never happens (hardware or driver bug preventing the signal from being written).

**Implication**: `accelerator_barrier()` is not an unconditional correctness guarantee on all production systems. Application-level verification (double-run, checksums) is necessary as a second line of defence.

## The Double-Wait Pattern in Practice

The double-wait is a pragmatic workaround when premature return is suspected but not yet confirmed. It adds latency but does not change correctness if the barrier is working properly, so it is safe to enable in production:

```cpp
#ifdef WORKAROUND_PREMATURE_BARRIER
#define accelerator_barrier() do { \
    real_accelerator_barrier();    \
    real_accelerator_barrier();    \
} while(0)
#endif
```

Monitor whether this changes observed behaviour. If double-wait eliminates wrong answers, you have confirmed premature return. If it does not help but inserting a no-op kernel does, the issue is with the wait primitive specifically, not with the underlying completion signal.

## SYCL/Level Zero Specifics

Level Zero (the backend for Intel GPU runtimes) separates command submission from synchronisation. A `q.wait()` should wait for all previously submitted command lists to retire. Documented bugs include:

- `q.wait()` returning before the associated fence in Level Zero has been signalled.
- `q.wait()` entering an `ioctl(i915, I915_GEM_WAIT)` call that never returns (kernel driver bug, not runtime bug).

The latter requires a node reboot and cannot be worked around in application code. Detect it by checking process state (`D` in `ps aux`) and the kernel function via `/proc/PID/wchan`.

## Stream Ordering and Compute Streams

All GPU work must be submitted to the *same* stream/queue if you rely on in-order execution guarantees. Mixing default stream and non-default streams invalidates ordering assumptions on some backends.

Grid uses `computeStream` (CUDA/HIP) or `theGridAccelerator` (SYCL) consistently throughout. If mixing Grid with third-party GPU code, ensure the third-party code is directed to the same stream, or insert explicit inter-stream barriers.

## Checklist for New GPU Code

1. Every kernel launch is followed by an `accelerator_barrier()` before reading device-side output on the host.
2. All device-to-host copies use an explicit stream synchronisation after the copy, not before.
3. If results are non-deterministic across runs, insert a second barrier and observe whether reproducibility improves.
4. For correctness-critical operations (reductions that will be compared against reference values), add the double-run checksum test from `correctness-verification.md`.
5. If the process hangs at 100% CPU in a runtime library function, this is a driver/runtime bug — there is no application-level fix beyond scheduling a node reboot.

## Making an asynchronous GPU fault attributable

A "Memory access fault by GPU node-N ... Reason: Unknown" is raised by the
runtime, asynchronously, possibly long after the offending kernel or copy
was *queued*.  `--debug-signals` sees nothing useful: the host is elsewhere
by then.  Before exhaustive logging or bisection, make the fault land on the
call that caused it:

| runtime | serialise launches/copies | name kernels as they launch |
|---|---|---|
| HIP/ROCm | `AMD_SERIALIZE_KERNEL=3 AMD_SERIALIZE_COPY=3` (every launch and copy synchronous) | `AMD_LOG_LEVEL=3` (4 is very verbose) |
| CUDA | `CUDA_LAUNCH_BLOCKING=1` | `compute-sanitizer --tool memcheck ./bin ...` (names the kernel and the bad address) |
| SYCL/L0 | `SYCL_PI_TRACE=2` / `ZE_DEBUG=1` | (Level Zero validation layer: `ZE_ENABLE_VALIDATION_LAYER=1`) |

5-10x slower, so run on the smallest reproducer (one node, small lattice).
The last kernel name printed before the fault is the culprit; combined with
Grid's `--log ...,Memory` (every MemoryManager transfer with size, direction,
device and host pointers, in program order) it separates "stale pointer after
eviction" from "out-of-range index" in one run.  First used 2026-08-28 for the
NRHS>=6 fault in the 3-level example (`systems/Frontier/nrhs_fault.job`).

**GPU core dumps (ROCm):** on a fault the runtime writes `gpucore.<pid>` (a
GCD's memory, ~22 GB on MI250X) into the process's working directory; it obeys
`ulimit -c` exactly as the kernel core does, so `ulimit -c 0` in the batch
script suppresses it (Slurm propagates rlimits to `srun` tasks by default) and
`ulimit -c unlimited` + `HSA_ENABLE_DEBUG=1` produces one ROCgdb can load.
Deleted-while-open dumps on NFS persist as hidden `.nfs*` files against quota
until the holder exits.  Run fault hunts from a Lustre directory.
Source: ROCgdb documentation, "AMD GPU" chapter, core-dump section.
