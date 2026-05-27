---
name: ref_grid_simt_pattern
description: Grid GRID_SIMT
metadata: 
  node_type: memory
  type: reference
  originSessionId: 956e80aa-401d-481a-80bb-17f8abe1c131
---

## The problem

On CPU, `accelerator_for(sf, oSites, Nsimd, {...})` expands to `thread_for(sf, oSites, {...})` — one thread per osite. `acceleratorSIMTlane(Nsimd)` always returns **0** on CPU. If you need to iterate all Nsimd lanes (e.g. to extract SIMD-packed data), you must loop explicitly on CPU.

On GPU, `accelerator_for` launches one SIMT thread per (osite × lane). `acceleratorSIMTlane(Nsimd)` returns the actual lane index [0, Nsimd).

## Correct pattern (from MomentumProject::ImportVector)

```cpp
accelerator_for(sf, osites, Nsimd, {
#ifdef GRID_SIMT
  {
    int lane = acceleratorSIMTlane(Nsimd);
#else
    for (int lane = 0; lane < Nsimd; lane++) {
#endif
      // body using lane
    }
  });
```

- On GPU: `GRID_SIMT` is defined → single-lane body, lane from hardware
- On CPU: `GRID_SIMT` is not defined → explicit lane loop inside the osite thread

## When is this needed?

Only when you explicitly need the lane index, e.g.:
- Extracting scalar data from SIMD-packed `vobj` via `extractLane(lane, src[sf])`
- Computing full local coordinates from (osite, lane) → `Lexicographic::CoorFromIndex(icoor, lane, simd_layout)`

When using `coalescedRead`/`coalescedWrite`, this pattern is **not needed** — those handle lane selection transparently.

## Pitfall that caused a bug

`A2ASpatialSum::PackVectors` originally used `accelerator_for` without the `#ifdef GRID_SIMT` lane loop. On CPU, only lane=0 was extracted, giving wrong norms (~8× too small for `GEN_SIMD_WIDTH=64`, `Nsimd=4`). Fix: add the `#ifdef GRID_SIMT` pattern. See [[ref_accelerator_for]].

## Model file
`Grid/algorithms/blas/MomentumProject.h`, function `ImportVector`, lines ~166-207.
