---
name: mpi-heterogeneous
description: Diagnose and work around MPI correctness bugs on heterogeneous (CPU+GPU) systems — device buffer aliasing in MPI_Sendrecv, AARCH64 PLT corruption from libfabric, topology-dependent allreduce hangs, and deterministic point-to-point reduction trees as a replacement for MPI_Allreduce.
user-invocable: true
allowed-tools:
  - Read
  - Bash(grep -r)
---

# MPI Correctness on Heterogeneous HPC Systems

## The Core Problem

MPI libraries were designed for CPU-resident buffers. When GPU-resident buffers are passed directly (GPU-aware MPI / GPU direct RDMA), several correctness assumptions break:

- **Buffer aliasing**: The MPI library may internally alias send/receive buffer addresses for `MPI_Sendrecv` in ways that are safe for CPU memory but wrong for GPU memory with different cache coherency rules.
- **RDMA bandwidth**: GPU direct RDMA on some fabrics operates at a fraction of peak wirespeed (documented at ~30% on Pontevecchio/Aurora), making host-staging mandatory for performance even when correctness is not an issue.
- **Collective tree topology**: `MPI_Allreduce` implementations may select reduction trees based on process count or communicator topology that expose rank-ordering bugs, causing hangs on some configurations but not others.

## Bug Class 1: Device Buffer Aliasing in MPI_Sendrecv

**Symptom**: `MPI_Sendrecv` with GPU-resident send and receive buffers produces wrong results. The received data matches neither the expected payload nor a host-staged copy. The failure is *deterministic* for a given problem size and process count, but *history-dependent* — earlier sends affect which alias is selected.

**Root cause**: The MPI library internally reuses GPU buffer addresses for temporary staging without proper device memory ordering. When the same physical GPU memory pages appear in both the send and receive paths, writes from one path corrupt the other.

**Diagnosis**: 
1. Enable per-packet checksumming (see `correctness-verification.md`). If the checksum on the received packet does not match the sent checksum, the data was corrupted in transit.
2. Replace `MPI_Sendrecv` with separate `MPI_Isend` + `MPI_Irecv` + `MPI_Waitall`. If this fixes the problem, the bug is in the `MPI_Sendrecv` implementation's internal buffer handling.
3. Stage through host memory (`cudaMemcpy`/`hipMemcpy` to a host buffer, then `MPI_Sendrecv` on host buffers, then copy back). If this fixes the problem, confirms GPU-specific aliasing.

**Reported as**: MPICH issue #7302. Affects MPICH on Intel Pontevecchio (Aurora) with device-resident buffers.

**Workaround**: Do not use `MPI_Sendrecv` with GPU buffers. Use asynchronous send/receive pairs or host-staging. See `communication-overlap.md` for the full pipeline pattern.

## Bug Class 2: PLT Corruption on AARCH64 (libfabric)

**Symptom**: Application crashes or hangs on first `MPI_Comm_dup` call on AARCH64 systems (e.g. NVIDIA Grace/H200). Backtrace shows a bad instruction in the PLT (Procedure Linkage Table) for `MPI_Comm_dup` — specifically a `br x15` instruction that should instead be a proper trampoline.

**Root cause**: `libfabric`'s memory registration cache monitor patches PLT entries at runtime to intercept memory allocation calls. Its AARCH64 trampoline generation writes an incorrect instruction sequence, leaving `br x15` (branch to whatever happens to be in x15) in the PLT entry. The next call through that PLT entry executes garbage.

**Diagnosis**:
```bash
# Check if the PLT entry is corrupted
objdump -d /proc/PID/exe | grep -A5 "MPI_Comm_dup@plt"
# Look for "br x15" — this should be a proper stub, not a register branch
```

Or check the disassembly of the live process:
```bash
gdb -p PID -batch -ex "disassemble 'MPI_Comm_dup@plt'"
```

**Workaround**:
```bash
export FI_MR_CACHE_MONITOR=disabled
```
This prevents libfabric from patching PLT entries. It may reduce MR cache performance but restores correctness.

**Reported as**: libfabric issue #11451. Affects systems using AARCH64 + libfabric OFI provider (Cray Slingshot, AWS EFA) with memory registration cache enabled.

## Bug Class 3: Topology-Dependent Allreduce Hangs

**Symptom**: `MPI_Allreduce` hangs indefinitely on some node configurations but completes correctly on others. The failure correlates with process count (e.g. fails at 512 ranks, works at 256) or network topology (fails when crossing specific router boundaries).

**Root cause**: The MPI library's collective selection algorithm picks a reduction tree implementation that assumes symmetric participation from all ranks. A bug in one rank's contribution path (e.g. a GPU-side buffer not yet flushed when MPI reads it, due to premature barrier — see `gpu-runtime-correctness.md`) causes that rank to send wrong or incomplete data, and the tree-reduction protocol deadlocks waiting for data that never arrives correctly.

**Diagnosis**: Flight recorder step logging (see `hang-diagnosis.md`). SIGHUP broadcast to all ranks. Ranks that are hung will show step name `MPI_Allreduce::...`; ranks that completed will show the next step. The hung ranks are the *recipients* of the stale data, not necessarily the *cause*.

**Workaround — deterministic P2P reduction tree**:

Replace `MPI_Allreduce` with an explicit point-to-point binary tree reduction. This is slower for large communicators but:
1. Is immune to topology-dependent collective bugs.
2. Is deterministic in floating-point ordering (the tree is fixed, not chosen at runtime).
3. Makes the hang location explicit — each P2P operation is a named step in the flight recorder.

```cpp
// Binary tree reduction: rank 0 collects, then broadcasts
void GlobalSumP2P(double *data, int count, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &size);
    
    // Reduce phase: even ranks receive from odd neighbours
    for (int stride = 1; stride < size; stride *= 2) {
        if (rank % (2*stride) == 0) {
            int partner = rank + stride;
            if (partner < size) {
                std::vector<double> tmp(count);
                MPI_Recv(tmp.data(), count, MPI_DOUBLE, partner, 0, comm, MPI_STATUS_IGNORE);
                for (int i = 0; i < count; i++) data[i] += tmp[i];
            }
        } else if (rank % stride == 0) {
            int partner = rank - stride;
            MPI_Send(data, count, MPI_DOUBLE, partner, 0, comm);
            break;
        }
    }
    // Broadcast phase
    for (int stride = /* highest power of 2 <= size */; stride >= 1; stride /= 2) {
        if (rank % (2*stride) == 0) {
            int partner = rank + stride;
            if (partner < size)
                MPI_Send(data, count, MPI_DOUBLE, partner, 0, comm);
        } else if (rank % stride == 0) {
            int partner = rank - stride;
            MPI_Recv(data, count, MPI_DOUBLE, partner, 0, comm, MPI_STATUS_IGNORE);
        }
    }
}
```

Grid reference: `USE_GRID_REDUCTION` macro in `Grid/communicator/Communicator_mpi3.cc`.

## Compile-Time Guard Structure

Recommended macro structure to switch between the workaround paths:

```cpp
// In configure / CMake, expose as options:
// ACCELERATOR_AWARE_MPI  — use GPU direct (fast, potentially broken)
// GRID_CHECKSUM_COMMS    — per-packet checksums (overhead: ~5%)
// USE_GRID_REDUCTION     — P2P tree instead of MPI_Allreduce (slower, deterministic)
// FI_MR_CACHE_MONITOR    — libfabric PLT workaround (env var, not compile-time)
```

On a known-good system, enable `ACCELERATOR_AWARE_MPI` and disable the others. On a system with known bugs, disable `ACCELERATOR_AWARE_MPI` and enable `GRID_CHECKSUM_COMMS` + `USE_GRID_REDUCTION` as needed.

## Escalation Checklist

Before concluding a bug is in your code:

1. [ ] Can you reproduce with a minimal reproducer (two MPI ranks, no physics code)?
2. [ ] Does the failure rate correlate with buffer size, process count, or network route?
3. [ ] Does staging through host memory eliminate the failure?
4. [ ] Is the failure deterministic for a given input (same answer, always wrong) or stochastic?
5. [ ] Does the failure appear on a different MPI implementation (e.g. OpenMPI vs MPICH)?

Deterministic wrong answers that reproduce with minimal reproducers and disappear with host-staging are strong evidence of an MPI library bug. File with the MPI library issue tracker with the minimal reproducer.
