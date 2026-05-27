---
name: ref_lattice_vs_vector
description: When to use Lattice<T> vs std::vector<T> for GPU-portable field storage in Grid
metadata: 
  node_type: memory
  type: reference
  originSessionId: 956e80aa-401d-481a-80bb-17f8abe1c131
---

## Rule

Use `Lattice<vobj>` (or `std::vector<Lattice<vobj>>`) for any field that will be read or written inside `accelerator_for`. `std::vector<vobj>` is host memory and is NOT device-accessible.

## Before vs after GPU offload

```cpp
// CPU-only (host memory, not GPU accessible)
std::vector<SpinColourVector_v> tloopv(oSites, Zero());
// accessed directly: tloopv[ss]

// GPU-portable
Lattice<SpinColourVector_v> tloop(grid);
// accessed via view: autoView(tloop_v, tloop, AcceleratorWrite);
//                    coalescedWrite(tloop_v[ss], val);
```

## Corollary: function signatures

CPU-only version:
```cpp
void PackLeft(const std::vector<std::vector<vobj>> &leftv);
```

GPU-portable version:
```cpp
void PackLeft(const std::vector<Lattice<vobj>> &leftv);
```

## deviceVector for raw device buffers

`deviceVector<T>` (defined in Grid) is like `std::vector<T>` but in device-accessible memory. Use for raw scalar scratch/pack buffers (e.g. GEMM input/output staging). Not for structured lattice data.

## Pointer arrays for batched BLAS

`deviceVector<scalar *>` holds batch pointer arrays. Populate with `acceleratorPut(ptrs[t], base + offset)` — sets device-side pointer from host. See `A2ASpatialSum::Allocate`.

## Related
[[ref_coalesced_views]] [[ref_accelerator_for]]
