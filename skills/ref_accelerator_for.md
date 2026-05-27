---
name: ref_accelerator_for
description: Grid accelerator_for usage — converting block-strided thread_for to GPU-portable oSites loops
metadata: 
  node_type: memory
  type: reference
  originSessionId: 956e80aa-401d-481a-80bb-17f8abe1c131
---

## Pattern: block-strided thread_for → accelerator_for over oSites

Old CPU-only pattern (block-strided over orthog dimension):
```cpp
thread_for(r, rd, {
  int so = r * grid->_ostride[orthogdim];
  for (int n = 0; n < e1; n++)
  for (int b = 0; b < e2; b++) {
    int ss = so + n * stride + b;
    // work on site ss
  }
});
```

GPU-portable replacement:
```cpp
accelerator_for(ss, grid->oSites(), Nsimd, {
  // work on site ss — one SIMT thread per (osite, lane) on GPU
  // one thread per osite (lane loop implicit via GRID_SIMT) on CPU
});
```

Key rules:
- `accelerator_for(iter, count, Nsimd, body)` — Nsimd is `vobj::Nsimd()` or `grid->Nsimd()`
- On CPU: expands to `thread_for` over count, `acceleratorSIMTlane` always returns 0 — must use `#ifdef GRID_SIMT` pattern if iterating lanes explicitly (see [[ref_grid_simt_pattern]])
- On GPU: one SIMT thread per (iter × lane), `acceleratorSIMTlane(Nsimd)` returns actual lane
- Loop body must capture only scalar/POD by value or via device-accessible pointers; no `std::vector` or host containers inside the body
- `Coordinate` inside `accelerator_for` must be `AcceleratorVector<int, MaxDims>` (stack-allocated, device-safe) — Grid's `Coordinate` typedef already satisfies this

## Where defined
`Grid/threads/Accelerator.h` — CPU path ~line 607; GPU paths in conditional blocks above.

## Model file
`Grid/algorithms/blas/MomentumProject.h` — `ImportVector` is the canonical example of correct `accelerator_for` + SIMD lane extraction.
