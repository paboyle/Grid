---
name: ref_a2a_emf_work
description: "A2A Extended Meson Field GPU offload work — status, file locations, pending task"
metadata: 
  node_type: memory
  type: project
  originSessionId: 956e80aa-401d-481a-80bb-17f8abe1c131
---

## What was built

`Grid/algorithms/blas/A2ASpatialSum.h` — batched GEMM spatial sum replacing scalar SIMD accumulation. Included via `Grid/algorithms/Algorithms.h`.

`tests/Test_extended_meson_field.cc` — test with class `A2AExtendedMesonFieldRef` containing:
- CPU reference path (`use_blas=false`)
- BLAS path (`use_blas=true`) using `A2ASpatialSum`
- Per-phase timing with `[ref type=N]` / `[blas type=N]` labels
- 4 contraction types (0-3), all verified at machine precision (~4e-16 rel_err)

## Pending task: GPU offload class

**Goal**: Write `A2AExtendedMesonFieldGPU` in the same test file, replacing all `thread_for` loops with `accelerator_for`-based free function kernels.

The `thread_for` blocks to replace all have the form:
```cpp
thread_for(r, rd, {
  int so = r * grid->_ostride[orthogdim];
  for (int n = 0; n < e1; n++)
  for (int b = 0; b < e2; b++) {
    int ss = so + n * stride + b;
    // work
  }
});
```
Replace with `accelerator_for(ss, grid->oSites(), Nsimd, { ... })`.

**Free functions to write** (each takes `Lattice<T>` args, opens views internally):
- `A2ALoopPropagator` — outerProduct sum (loop build)
- `A2APackLeftConjugated` — conjugate left fermion fields into `Lattice<SpinColourVector_v>`
- `A2ALoopLeftContractionType0/1/2/3` — per-site loop × loop propagator → `tloop`
- `A2ALoopRightContractionType0/1/2/3` — per-site tloop × right → `loopRight[j]`

**Data structure changes required**:
- `tloopv`: `std::vector<SpinColourMatrix_v>` → `Lattice<SpinColourMatrix_v>` (PropagatorField)
- `leftv[i]`: `std::vector<SpinColourVector_v>` → `Lattice<SpinColourVector_v>`
- `loopRight[j]`: `std::vector<SpinColourVector_v>` → `Lattice<SpinColourVector_v>`

**Why**: `std::vector<vobj>` is host memory, not GPU accessible. See [[ref_lattice_vs_vector]].

**`A2ASpatialSum` impact**: `PackLeft`/`PackRight` currently take `std::vector<std::vector<vobj>>`. Once leftv/loopRight become `std::vector<Lattice<vobj>>`, those signatures must change to match.

## Timing on 8.8.8.16 (N_i=N_j=8, Nloop=4, 1 MPI rank)

Dominant costs:
- `loop_build`: 4-6 ms (outerProduct over 4 propagators)
- `pack_loopright`: 0.9-2.2 ms (type-dependent)
- `spatial_sum` (ref): ~1.5 ms
- `A2ASpatialSum TOTAL`: 2.5-4.3 ms (PackLeft+PackRight dominate GEMM on small volume)

## Related
[[ref_accelerator_for]] [[ref_coalesced_views]] [[ref_lattice_vs_vector]] [[ref_grid_simt_pattern]]
