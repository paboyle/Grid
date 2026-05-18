---
name: correctness-verification
description: Implement application-level correctness verification for HPC codes on unreliable hardware — double-run pattern, deterministic reductions, per-packet checksums, and flight recorder step logging.
user-invocable: true
allowed-tools:
  - Read
  - Bash(grep -r)
---

# Correctness Verification Infrastructure for HPC Codes

## The Problem

Leadership computing facilities sometimes have hardware or firmware bugs below the level visible to application code. The accelerator runtime can return from `q.wait()` or `cudaDeviceSynchronize()` before work is actually complete, or silently produce wrong answers in DMA transfers. Standard testing does not catch these because they are non-deterministic and often topology-dependent (fail only at specific process counts or on specific node configurations).

The symptoms look like numerical instabilities, random MPI hangs, or wrong physics results — not like crashes. Without deliberate infrastructure, diagnosing root cause takes months.

## The Double-Run Pattern

The most reliable correctness check for non-deterministic hardware bugs is to run every computation twice and compare bit-identical fingerprints.

**Key constraint**: the second run must use a *deterministic* code path. Non-deterministic floating-point ordering (e.g. from MPI_Allreduce with different reduction trees on retry) produces false mismatches. See `mpi-heterogeneous.md` for how to make reductions deterministic.

```cpp
// Pseudocode: double-run a step and compare CRC fingerprints
void run_step_verified(State &state) {
    state.save_checkpoint();
    
    uint64_t crc_a = run_step_and_fingerprint(state);
    state.restore_checkpoint();
    uint64_t crc_b = run_step_and_fingerprint(state);
    
    if (crc_a != crc_b) {
        report_mismatch("step", crc_a, crc_b);
        // Policy: abort, retry from checkpoint, or continue with alarm
    }
}
```

**Fingerprinting**: XOR-fold a CRC32 over all floating-point data after each step. XOR is order-independent, so it works across distributed nodes without communication. For field data:

```cpp
uint64_t fingerprint(const double *data, size_t n) {
    uint64_t acc = 0;
    for (size_t i = 0; i < n; i++) {
        uint64_t bits;
        memcpy(&bits, &data[i], sizeof(bits));
        acc ^= crc32(bits);
    }
    return acc;
}
```

On GPU, compute the XOR reduction on-device (avoids D2H transfer of the full field):

```cpp
// SYCL
uint64_t svm_xor(uint64_t *vec, uint64_t L) {
    uint64_t ret = 0;
    { sycl::buffer<uint64_t,1> abuff(&ret, {1});
      theGridAccelerator->submit([&](sycl::handler &cgh) {
        auto R = sycl::reduction(abuff, cgh, uint64_t(0), std::bit_xor<>());
        cgh.parallel_for(sycl::range<1>{L}, R,
          [=](sycl::id<1> i, auto &sum) { sum ^= vec[i]; });
      }); }
    theGridAccelerator->wait();
    return ret;
}
```

## Per-Packet Communication Checksums

Silent data corruption in MPI buffers (documented in MPICH with device-resident buffers; see `mpi-heterogeneous.md`) requires per-packet verification, not just end-to-end. The pattern:

1. Before packing a send buffer, compute a GPU-side checksum of the payload.
2. Append the checksum to the host staging buffer alongside the data.
3. After receiving and copying to device, recompute the checksum on-device and compare.

Salt each checksum with `packet_index + 1000 * mpi_tag` to detect transposition (packet A landing in packet B's slot):

```cpp
uint64_t salt = (uint64_t)packet_index + 1000ULL * mpi_tag;
checksum_send = checksum_gpu(payload_gpu, payload_words) ^ salt;
// ... transmit payload + checksum_send ...
checksum_recv = checksum_gpu(payload_gpu_recv, payload_words) ^ salt;
assert(checksum_recv == checksum_send);
```

Grid reference: `Grid/communicator/Communicator_mpi3.cc`, `#ifdef GRID_CHECKSUM_COMMS`.

## Flight Recorder: Step-Level Logging

Maintain a monotonic counter that names the current operation. On a hang, this is the only way to know *which* operation the process is stuck in without a debugger.

```cpp
struct FlightRecorder {
    std::atomic<uint64_t> step_counter{0};
    const char *step_name = "init";
    
    void step_log(const char *name) {
        step_name = name;
        step_counter.fetch_add(1, std::memory_order_relaxed);
    }
};
extern FlightRecorder gRecorder;
```

In Record mode, also store floating-point norms and communication checksums to vectors. In Verify mode, compare against stored values:

```cpp
void norm_log(double val) {
    if (mode == Record) norm_log_vec.push_back(val);
    if (mode == Verify) {
        double expected = norm_log_vec[norm_counter];
        if (val != expected) {  // bit-exact for deterministic paths
            std::cerr << "MISMATCH at step " << step_counter
                      << " (" << step_name << "): "
                      << std::hexfloat << val << " vs " << expected << "\n";
            print_backtrace();
        }
        norm_counter++;
    }
}
```

Grid reference: `Grid/util/FlightRecorder.h`, `Grid/util/FlightRecorder.cc`.

## Signal Handler for Hang Detection

Install a SIGHUP handler that dumps the current flight recorder state. This is async-safe only if the handler writes to a pre-allocated buffer using `write()` (not `printf`):

```cpp
static char hang_buf[4096];

static void sighup_handler(int) {
    int n = snprintf(hang_buf, sizeof(hang_buf),
        "rank=%d step=%llu name=%s\n",
        mpi_rank,
        (unsigned long long)gRecorder.step_counter.load(),
        gRecorder.step_name);
    write(STDERR_FILENO, hang_buf, n);
    // Optional: call backtrace_symbols_fd (async-safe on Linux)
    void *frames[64];
    int depth = backtrace(frames, 64);
    backtrace_symbols_fd(frames, depth, STDERR_FILENO);
}

// In main():
signal(SIGHUP, sighup_handler);
```

To diagnose a hang across all ranks: `kill -HUP $(pgrep my_app)` or via job scheduler.

## What to Verify at Each Step

| Data type | Fingerprint method | Frequency |
|---|---|---|
| Lattice fields | XOR of CRC32 over float64 words | Every algorithmic step |
| Communication buffers | GPU XOR reduction, salted | Every MPI operation |
| Scalar reductions | Bit-exact match of double | Every GlobalSum |
| Iteration counters | Exact integer match | Every solver iteration |

## When to Abort vs Continue

- **Abort immediately**: communication checksum mismatch (data is corrupt, continuing will silently propagate errors).
- **Log and continue**: norm mismatch in Verify mode if you need to map out which operations are unreliable.
- **Retry from checkpoint**: double-run mismatch when the underlying bug is non-deterministic (second retry will usually pass).

Track the mismatch rate over a production run. A rate above ~1/1000 steps indicates a systemic hardware issue that should be escalated to the facility.
