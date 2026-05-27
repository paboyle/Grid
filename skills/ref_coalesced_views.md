---
name: ref_coalesced_views
description: Grid coalescedRead/coalescedWrite and autoView — GPU-portable field access inside accelerator_for
metadata: 
  node_type: memory
  type: reference
  originSessionId: 956e80aa-401d-481a-80bb-17f8abe1c131
---

## View access modes

```cpp
autoView(v, field, AcceleratorRead);   // read-only, device-accessible
autoView(v, field, AcceleratorWrite);  // write-only, device-accessible
autoView(v, field, AcceleratorReadWrite); // read-write, device-accessible
autoView(v, field, CpuRead);           // CPU only (avoids GPU migration)
autoView(v, field, CpuWrite);          // CPU only
```

Views must be opened **before** `accelerator_for` and closed (go out of scope) **after**. Never open a view inside the accelerator_for body.

## coalescedRead / coalescedWrite

Inside `accelerator_for(ss, oSites, Nsimd, { ... })`:

```cpp
auto site = coalescedRead(v[ss]);          // reads SIMT lane; returns scalar_object on GPU, vobj on CPU
coalescedWrite(v[ss], site);               // writes SIMT lane
```

- `coalescedRead(v[ss])` calls `v.operator()(ss)` which on GPU returns `extractLane(lane, v[ss])` — one lane per SIMT thread, contiguous across threads → coalesced
- On CPU returns the full vobj (no lane extraction needed; handled transparently)
- The returned type is `decltype(coalescedRead(v[ss]))` — use `auto` or match with scalar_object

## Typical kernel pattern

```cpp
autoView(out_v, out, AcceleratorWrite);
autoView(in_v,  in,  AcceleratorRead);
accelerator_for(ss, grid->oSites(), vobj::Nsimd(), {
  auto x = coalescedRead(in_v[ss]);
  // modify x ...
  coalescedWrite(out_v[ss], x);
});
```

## Free function kernel signature

```cpp
template<class vobj>
void MyKernel(Lattice<vobj> &out, const Lattice<vobj> &in)
{
  GridBase *grid = in.Grid();
  autoView(out_v, out, AcceleratorWrite);
  autoView(in_v,  in,  AcceleratorRead);
  accelerator_for(ss, grid->oSites(), vobj::Nsimd(), {
    auto x = coalescedRead(in_v[ss]);
    coalescedWrite(out_v[ss], x);
  });
}
```

## What NOT to do
- Do not access `std::vector` elements inside `accelerator_for` — not device-accessible
- Do not use `CpuRead`/`CpuWrite` views inside `accelerator_for` — GPU will fault
- Do not assign to `v[ss]` directly inside `accelerator_for` — use `coalescedWrite`
- Do not open multiple write views on the same field simultaneously

## Related
[[ref_accelerator_for]] [[ref_lattice_vs_vector]]
