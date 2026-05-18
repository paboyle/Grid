---
name: hang-diagnosis
description: Diagnose and isolate process hangs on HPC systems — distinguishing kernel-level ioctl hangs, infinite poll loops, collective deadlocks, and GPU completion signalling failures using async-safe signal handlers and flight recorder step counters.
user-invocable: true
allowed-tools:
  - Read
  - Bash(grep -r)
  - Bash(strace)
  - Bash(gdb)
---

# Hang Diagnosis on HPC Systems

## Taxonomy of Hangs

Not all hangs are the same. Misidentifying the type leads to wrong mitigation. The four distinct classes encountered on production leadership systems:

### 1. Kernel-level ioctl hang (never returns)
The process is in `D` (uninterruptible sleep) state. `strace` shows it blocked in an `ioctl` syscall. The GPU device driver has entered an unrecoverable state.

**Diagnosis**: `ps aux | grep D` — the process shows `D` state. `cat /proc/PID/wchan` shows `i915_gem_wait_for_error` or similar.

**Resolution**: Only a driver reload or node reboot recovers it. Log the node identifier and request replacement from the facility scheduler.

### 2. Infinite poll loop (`q.wait()` or `cudaDeviceSynchronize()` never returns)
The process is in `R` (running) state, consuming 100% CPU. A polling loop inside the runtime is checking a completion flag that never becomes true, either because the hardware never sets it or because the flag is in a memory region not visible to the polling thread.

**Diagnosis**: `top` shows the rank at 100% CPU. `strace -p PID` shows repeated `futex` or `read` syscalls with zero-length results, or no syscalls at all (pure spinloop). `perf top -p PID` shows the process burning cycles in a single tight loop in a runtime library (e.g., `ze_intel_gpu.so`).

**Resolution**: The double-wait workaround — submit a trivially cheap kernel after the operation under test to act as a fence, then wait for the trivial kernel. See `gpu-runtime-correctness.md`.

### 3. Collective deadlock
One or more ranks are blocked in an MPI call, usually `MPI_Allreduce` or `MPI_Barrier`, while others are not. Root cause: a topology-dependent bug in the MPI library's collective algorithm where some ranks' contributions never arrive.

**Diagnosis**: Flight recorder step logs show some ranks at step N (inside the collective) while others are at step N+1 or stuck at step N with different `step_name` strings. The hung ranks will show `D` or `S` state in `ps`.

**Resolution**: Replace `MPI_Allreduce` with a deterministic point-to-point tree reduction. See `mpi-heterogeneous.md`.

### 4. Premature return from wait (silent wrong answer, not a hang)
The runtime returns from `q.wait()` before the GPU work is complete. The next operation reads stale data. This is not a hang — it manifests as a wrong answer or non-deterministic floating-point results. It is listed here because it is the most confusing failure mode: the code appears to run correctly and completes normally.

**Diagnosis**: Double-run with checksum (see `correctness-verification.md`). Insert a second `q.wait()` after the first and observe if results become reproducible. If inserting the second wait "fixes" wrong answers, the first wait was returning prematurely.

## Flight Recorder for Hang Localization

The most important diagnostic tool is knowing *which operation* a process is in when it hangs. Maintain a named step counter:

```cpp
// Call at the start of every major operation
FlightRecorder::StepLog("MPI_Allreduce::norm");
// ... do the operation ...
FlightRecorder::StepLog("MPI_Allreduce::done");
```

On SIGHUP, dump rank, step counter value, and step name to stderr in an async-safe manner:

```cpp
static void sighup_handler(int) {
    char buf[256];
    int n = snprintf(buf, sizeof(buf), "rank %d: step %llu '%s'\n",
        comm_rank,
        (unsigned long long)step_counter,
        step_name);
    write(2, buf, n);
    // backtrace_symbols_fd is async-safe on Linux glibc
    void *frames[32];
    backtrace_symbols_fd(frames, backtrace(frames, 32), 2);
}
signal(SIGHUP, sighup_handler);
```

Broadcast SIGHUP to all ranks from outside the job:
```bash
# In a separate shell while the job is hung
squeue --job $JOBID -o "%i %N" | awk '{print $2}' | \
  xargs -I{} ssh {} "pkill -SIGHUP -f my_application"
```

The step names from all ranks will reveal which collective operation has diverged.

## Distinguishing Driver Hang from MPI Hang

| Symptom | Driver hang | MPI hang |
|---|---|---|
| Process state | `D` (ioctl) or `R` (spinloop) | `S` (blocked in syscall) |
| `strace` | blocked `ioctl` or tight loop | blocked `recvmsg` / `read` |
| Scope | single rank / single node | subset of ranks, pattern-dependent |
| Recovery | reboot node | cancel job |
| Flight recorder | step name is a GPU operation | step name is a collective |

## Reducing Diagnostic Time

1. **Name every collective operation** in the flight recorder before calling it.
2. **Separate GPU work from MPI work** in the code so the step name unambiguously identifies which subsystem is hung.
3. **Log node identifiers** alongside step names so flaky nodes can be identified and blacklisted.
4. **Request flight recorder dumps from all ranks simultaneously** (SIGHUP broadcast) rather than attaching a debugger — attaching `gdb` to one rank of a hung MPI job usually deadlocks the debugger too.

## What Not to Do

- Do not `kill -9` a hung rank immediately — get the flight recorder dump first, otherwise diagnostic information is lost.
- Do not assume the first rank that prints an error is the faulty one — collective hangs are frequently caused by the *last* rank to arrive at the barrier.
- Do not use `MPI_Abort` in the hang handler — it may itself hang on some implementations. Use `_exit(1)` to force termination.
