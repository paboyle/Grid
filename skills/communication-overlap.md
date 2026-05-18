---
name: communication-overlap
description: Design and implement communication/computation overlap pipelines for GPU+MPI codes — per-packet event tracking, host-staging through pinned memory, internode/intranode bandwidth separation, and the 7-phase pipeline pattern that replaces broken accelerator-aware MPI paths.
user-invocable: true
allowed-tools:
  - Read
  - Bash(grep -r)
---

# Communication/Computation Overlap Pipeline Design

## Why GPU-Direct MPI Is Often Not the Right Default

GPU-direct RDMA (passing GPU buffer pointers directly to MPI) is appealing because it eliminates explicit D2H/H2D copies. In practice on several leadership systems:

- **Bandwidth**: RDMA at 30% of wirespeed has been observed on Pontevecchio/Aurora. The overhead of staging through pinned host memory can be *lower* total latency than slow RDMA.
- **Correctness**: Device buffer aliasing in `MPI_Sendrecv` (see `mpi-heterogeneous.md`) makes direct GPU-to-GPU transfer unreliable.
- **Overlap**: Host-staging enables fine-grained overlap — each packet's D2H can be issued as a separate asynchronous event, and the corresponding MPI send can fire as soon as *that packet* arrives in host memory, not after all packets are ready.

The pipeline pattern below was developed to replace broken MPICH accelerator-aware paths. It achieves genuine computation/communication overlap by tracking per-packet GPU events.

## The 7-Phase Pipeline

Given a set of halo exchange operations (each identified by a `packet_index`):

### Phase 0: Prepare data on device
Pack halo data into contiguous GPU buffers. One buffer per direction/neighbour.

### Phase 1: Post receives + start D2H
Post all `MPI_Irecv` calls immediately (into pinned host buffers). Simultaneously, start asynchronous D2H copies for all send buffers:

```cpp
for (auto &pkt : send_packets) {
    MPI_Irecv(pkt.host_recv_buf, pkt.bytes, MPI_BYTE,
               pkt.src_rank, pkt.tag, comm, &pkt.recv_req);
    
    acceleratorCopyFromDeviceAsync(pkt.device_send_buf,
                                   pkt.host_send_buf,
                                   pkt.bytes, &pkt.d2h_event);
}
```

The key: `pkt.d2h_event` is a per-packet GPU event (e.g. `cudaEvent_t`, `hipEvent_t`, or SYCL event). We can poll individual packet completion rather than waiting for all.

### Phase 2: Fire sends as D2H completes (packet by packet)
Poll packet D2H events. As each packet becomes ready in host memory, immediately fire the corresponding `MPI_Isend`. Also start intranode D2D copies at this point — these are deferred until now to avoid competing with the internode D2H on PCIe bandwidth:

```cpp
bool all_sent = false;
while (!all_sent) {
    all_sent = true;
    for (auto &pkt : send_packets) {
        if (!pkt.sent && acceleratorEventIsComplete(pkt.d2h_event)) {
            MPI_Isend(pkt.host_send_buf, pkt.bytes, MPI_BYTE,
                       pkt.dst_rank, pkt.tag, comm, &pkt.send_req);
            pkt.sent = true;
            start_intranode_copy(pkt);  // now safe, D2H is done
        }
        if (!pkt.sent) all_sent = false;
    }
}
```

### Phase 3: Poll receives + start H2D as each arrives
`MPI_Test` individual receive requests. As each completes, immediately start the H2D copy into device-resident halo buffer:

```cpp
bool all_recvd = false;
while (!all_recvd) {
    all_recvd = true;
    for (auto &pkt : recv_packets) {
        if (!pkt.h2d_started) {
            int flag = 0;
            MPI_Test(&pkt.recv_req, &flag, MPI_STATUS_IGNORE);
            if (flag) {
                acceleratorCopyToDeviceAsync(pkt.host_recv_buf,
                                              pkt.device_recv_buf,
                                              pkt.bytes, &pkt.h2d_event);
                pkt.h2d_started = true;
            }
        }
        if (!pkt.h2d_started) all_recvd = false;
    }
}
```

### Phase 4: Wait for all sends
```cpp
std::vector<MPI_Request> send_reqs;
for (auto &pkt : send_packets) send_reqs.push_back(pkt.send_req);
MPI_Waitall(send_reqs.size(), send_reqs.data(), MPI_STATUSES_IGNORE);
```

### Phase 5: Wait for all H2D copies
```cpp
for (auto &pkt : recv_packets) acceleratorEventWait(pkt.h2d_event);
```

### Phase 6: Run interior computation
The interior (non-halo) computation can run from Phase 1 onwards, overlapped with all of the above:

```cpp
// Launched in Phase 1, runs in parallel with the pipeline
accelerator_for(ss, interior_sites, ...) { compute_interior(ss); }
```

Synchronise with interior before using the full field:
```cpp
accelerator_barrier();  // interior kernel done
// Halo H2D is also complete (Phase 5 above)
// Now safe to use full field
```

## Per-Packet Event Tracking Data Structure

```cpp
struct Packet {
    // Buffers
    void *device_send_buf;
    void *host_send_buf;    // pinned
    void *device_recv_buf;
    void *host_recv_buf;    // pinned
    size_t bytes;
    
    // MPI
    int src_rank, dst_rank, tag;
    MPI_Request send_req, recv_req;
    
    // GPU events (one per packet, not one global barrier)
    AcceleratorEvent d2h_event;
    AcceleratorEvent h2d_event;
    
    // State flags
    bool sent = false;
    bool h2d_started = false;
};
```

The critical design point: `d2h_event` and `h2d_event` are **per-packet**, not global. This allows the MPI send for packet 0 to fire while packet 1's D2H is still in progress.

## Internode vs Intranode Separation

PCIe (GPU-to-CPU) and NVLink/xGMI (GPU-to-GPU within a node) are separate bandwidth resources. They do not compete with each other, but they *do* compete with each other for transactions if both are active simultaneously.

Strategy: complete all internode D2H copies first (to maximise NIC injection bandwidth), then start intranode D2D copies (which use NVLink/xGMI and do not contend with PCIe for internode traffic):

```cpp
// In Phase 2: start intranode D2D only after D2H is confirmed complete
if (pkt.is_intranode && pkt.d2h_done) {
    // Use peer access (cudaMemcpyPeerAsync / hipMemcpyPeerAsync)
    // rather than staging through host for intranode
    cudaMemcpyPeerAsync(pkt.peer_recv_buf, pkt.dst_device,
                         pkt.device_send_buf, pkt.src_device,
                         pkt.bytes, computeStream);
}
```

Grid reference: `Grid/communicator/Communicator_mpi3.cc` — search for `NVLINK_GET` and `ACCELERATOR_AWARE_MPI` conditional blocks.

## Pinned Memory Allocation

All host staging buffers must be pinned (page-locked) for async D2H/H2D:

```cpp
// CUDA
cudaMallocHost(&host_buf, bytes);
cudaFreeHost(host_buf);

// HIP  
hipHostMalloc(&host_buf, bytes, hipHostMallocDefault);
hipHostFree(host_buf);

// SYCL
host_buf = sycl::malloc_host(bytes, *queue);
sycl::free(host_buf, *queue);
```

Pre-allocate at startup. Repeated `cudaMallocHost` in the hot path adds latency from the OS memory manager.

## Checksumming in the Pipeline

Insert checksum computation before D2H (on the GPU-resident data) and verification after H2D (on the received GPU-resident data). See `correctness-verification.md` for the checksum pattern. The salting (`packet_index + 1000 * tag`) detects packet transposition — critical for diagnosing MPI buffer aliasing bugs where two packets' contents are swapped.

## Smoke Test for a New System

Before running physics, validate the pipeline on a synthetic benchmark:

```cpp
// Send a buffer of known values, receive and check
// Run at multiple message sizes: 4KB, 64KB, 1MB, 16MB
// Run at multiple process counts: 2, 8, 64, 512
// Verify checksums on every packet
// Measure bandwidth: should be ≥ 80% of FDR/HDR/NDR peak for host-staged
```

Any bandwidth below 50% of theoretical, or any checksum failure, indicates a problem in the communication stack that must be resolved before production runs.
