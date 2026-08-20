# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

Grid is a data-parallel C++ library for lattice QCD. It provides SIMD-vectorised lattice containers, MPI-based domain decomposition, GPU acceleration (CUDA/HIP/SYCL), and a full suite of QCD algorithms including HMC.

## Build

Uses GNU Autotools. The bootstrap step only needs to run once (or after `configure.ac` changes).

```bash
./bootstrap.sh                    # downloads Eigen 3.4.0, generates configure
mkdir build && cd build
../configure [options]
make -j$(nproc)
make check                        # run root-level tests
make install
```

Key configure options:

| Option | Common values |
|--------|---------------|
| `--enable-simd=` | `AVX2`, `AVX512`, `KNL`, `A64FX`, `NEONv8`, `GPU` |
| `--enable-comms=` | `mpi-auto`, `mpi3-auto`, `none` |
| `--enable-accelerator=` | `cuda`, `hip`, `sycl` |
| `--enable-shm=` | `shmopen`, `hugetlbfs`, `nvlink` |
| `--enable-Nc=` | `3` (default), `2`, `4`, `5` |
| `--with-gmp=`, `--with-mpfr=`, `--with-fftw=`, `--with-lime=` | paths to libs |
| `--enable-hdf5`, `--enable-mkl`, `--enable-lapack` | optional features |

GPU builds additionally need `--enable-gen-simd-width=64` (sets 512-bit SIMD width for GPU warp/wavefront sizing) and `--enable-unified=no --enable-shm=nvlink` for multi-GPU runs.

To speed up compilation, `--disable-fermion-reps --disable-gparity` skips instantiating G-parity and higher-representation fermion operators.

Platform recipes from `README.md`:
- **KNL**: `--enable-simd=KNL --enable-comms=mpi3-auto --enable-mkl`
- **Skylake/Haswell**: `--enable-simd=AVX512` or `AVX2` + `--enable-comms=mpi3-auto`
- **AMD EPYC**: `--enable-simd=AVX2 --enable-comms=mpi3`
- **A64FX (Fugaku)**: `--enable-simd=A64FX --enable-comms=mpi3 --enable-shm=shmget` (see `SVE_README.txt`)

Complete, working `configure` invocations for specific HPC systems (Frontier/ROCm, Perlmutter/CUDA, Summit, SDCC-A100, etc.) live in `systems/<platform>/config-command`. These are the canonical references for production builds.

Required external libs: GMP, MPFR, OpenSSL, zlib.

### Use `systems/` for real machines

`systems/<machine>/` holds the known-good build for each production platform (`Frontier`, `Aurora`, `Perlmutter`, `Summit`, `Tursa`, `Lumi`, `Booster`, `Crusher`, `SDCC-*`, `mac-arm`, …). Each contains a `config-command` (the exact `../../configure` invocation) and a `sourceme.sh` (module loads and env). **Prefer copying/adapting these over hand-rolling configure flags** — they encode compiler workarounds, `LDFLAGS`, and shared-memory settings that are easy to get wrong. `systems/WorkArounds.txt` records known vendor bugs.

Note the GPU builds use `--enable-simd=GPU --enable-gen-simd-width=64`, so `Nsimd` is *not* 1 on device (it is `64/sizeof(scalar)`).

### Regenerating `Make.inc` — required after adding or deleting source files

`Make.inc` files are generated, not tracked in git (`.gitignore`d). `scripts/filelist` walks `Grid/`, `tests/*`, `benchmarks/`, `examples/`, and `HMC/` and writes the file lists and per-test `bin_PROGRAMS` rules. Every new `.cc`/`.h` in `Grid/`, and every new `Test_*.cc` / `Benchmark_*.cc` / `Example_*.cc`, is invisible to the build until you run:

```bash
./scripts/filelist    # from the source root, then re-run configure/make
```

`bootstrap.sh` runs it for you on the first setup.

## Running Tests and Benchmarks

```bash
# From build directory
make check                          # root-level tests (Test_simd, Test_cshift, etc.)
make -C tests/<subdir> tests        # build tests in a subdirectory
make tests                          # build all tests across all subdirectories
./tests/core/Test_simd              # run a single test binary directly
mpirun -n 4 ./tests/core/Test_cshift --grid 16.16.16.16 --mpi 1.1.1.4
```

`make check` is a thin smoke test — building a subdirectory with `make -C tests/<subdir> tests` and running the relevant binaries directly is the normal development loop. Test binaries take Grid's standard command-line arguments (`--grid`, `--mpi`, `--accelerator-threads`, `--threads`, `--debug-signals`, `--log`); see `Grid/util/Init.cc`.

Test subdirectories and their focus: `core` (SIMD, stencil, comms), `solver` (CG, GMRES, eigensolvers), `hmc` (MD integrators), `forces` (fermion forces), `lanczos`, `IO`, `smearing`, `sp2n`, `debug`.

Tests and benchmarks that need optional fermion representations are guarded by `disable_tests_without_instantiations.h` / `disable_benchmarks_without_instantiations.h`, so a `--disable-fermion-reps --disable-gparity` build silently compiles them to no-ops.

## Architecture

### Layer stack (bottom to top)

1. **SIMD layer** (`Grid/simd/`) — platform-specific intrinsics wrapped into `vRealF`, `vComplexD`, etc. The SIMD width and layout are compile-time constants controlled by `--enable-simd`.

2. **Tensor layer** (`Grid/tensors/`) — Lorentz/colour/spin tensor algebra built on top of SIMD types. `iMatrix`, `iVector`, `iScalar` templates compose into QCD types like `ColourMatrix`, `SpinColourVector`.

3. **Lattice layer** (`Grid/lattice/`) — `Lattice<T>` container: a site-local tensor replicated across a distributed Cartesian grid. All arithmetic is site-parallel and expression-template-fused.

4. **Cartesian/comms layer** (`Grid/cartesian/`, `Grid/communicator/`) — `GridCartesian` holds the MPI topology and local/global geometry. `Grid/cshift/` implements nearest-neighbour halo exchange; `Grid/stencil/` is the optimised multi-hop stencil used by Dirac operators.

5. **Algorithm layer** (`Grid/algorithms/`) — iterative solvers (CG, GMRES, BiCGSTAB, mixed-precision), eigensolvers (Lanczos, LAPACK), FFT, smearing, and multigrid.

6. **QCD layer** (`Grid/qcd/`) — gauge and fermion actions, HMC integrators, observables.

### QCD subsystem (`Grid/qcd/`)

- `action/fermion/` — Wilson, Clover, DWF (Mobius), Staggered, twisted-mass, G-parity variants
- `action/gauge/` — Wilson gauge, Symanzik, Iwasaki, DBW2, plaquette+rect
- `representations/` — Fundamental, Adjoint, Two-index, Sp(2n)
- `hmc/` — Leapfrog, OMF2/OMF4 integrators; pseudofermion refreshment; Metropolis accept/reject
- `smearing/` — APE, Stout, HEX, gradient flow
- `observables/` — Polyakov loop, plaquette, topological charge

### GPU acceleration and the view/memory-manager discipline
### Multigrid (`Grid/algorithms/multigrid/`)

Aggregation-based algebraic multigrid for Wilson-type fermions. Key files: `CoarsenedMatrix.h` (coarse operator), `GeneralCoarsenedMatrix.h` and `GeneralCoarsenedMatrixMultiRHS.h` (general coarsening supporting multi-RHS solves), `Aggregates.h` (near-null vector construction), `Geometry.h` (coarse-grid geometry). `MultiGrid.h` is the top-level include.

### GPU acceleration

GPU support is injected via macros in `Grid/threads/Accelerator.h` — `accelerator_for(i, n, nsimd, {...})`, `accelerator_forNB` (non-blocking, must be followed by `accelerator_barrier()`), `accelerator_for2dNB`, and `accelerator_inline`. On a CPU build these degrade to `thread_for` (OpenMP). Unified virtual memory is on by default (`--enable-unified=yes`); device-aware MPI (`--enable-accelerator-aware-mpi`) avoids device→host copies on transfers.

Lattice data is **not** directly addressable inside a kernel. You must open a view with the correct access mode so `Grid/allocator/MemoryManager.h` can move/mark the data:

```cpp
autoView(out_v, out, AcceleratorWriteDiscard);   // RAII; closes at end of scope
autoView(in_v,  in,  AcceleratorRead);
accelerator_for(ss, grid->oSites(), Nsimd, {
    coalescedWrite(out_v[ss], coalescedRead(in_v[ss]));
});
```

Modes are `AcceleratorRead/Write/WriteDiscard` and `CpuRead/Write/WriteDiscard`. Getting the mode wrong (e.g. `AcceleratorRead` on a field you write) produces stale-data bugs that only appear on GPU builds. Inside kernels use `coalescedRead`/`coalescedWrite` rather than raw `operator[]` — they map the SIMD lane onto `threadIdx.x` so accesses stay coalesced.

### Repo-local debugging skills (`skills/`)

`skills/` contains hard-won, Grid-specific playbooks written as invocable skill files. Consult them before debugging in these areas rather than reasoning from first principles:

| File | Covers |
|---|---|
| `gpu-memory-performance.md` | `acceleratorThreads()`, LambdaApply thread mapping, `coalescedRead` idiom, fused vs staged HBM access |
| `gpu-runtime-correctness.md` | GPU runtime returning early from sync, silent wrong answers |
| `communication-overlap.md` | 7-phase halo pipeline, per-packet events, host-staging vs GPU-direct RDMA |
| `mpi-heterogeneous.md` | `MPI_Sendrecv` device-buffer aliasing, deterministic reductions |
| `compiler-validation.md` | Isolating GPU compiler codegen bugs, minimal reproducers |
| `correctness-verification.md` | Double-run fingerprinting, per-packet checksums, flight recorder |
| `hang-diagnosis.md` | Diagnosing MPI/accelerator hangs |

The key loop macros (defined in `Grid/threads/Accelerator.h`) are:
- `accelerator_for(iter, num, nsimd, {...})` — maps to CUDA/HIP kernel or OpenMP loop; `nsimd` is the innermost SIMD lane count
- `accelerator_forNB(...)` — non-blocking variant (no implicit barrier)
- `accelerator_for2dNB(iter1, num1, iter2, num2, nsimd, {...})` — 2D kernel launch
- `thread_for(iter, num, {...})` — CPU OpenMP loop (never dispatches to GPU)

On CPU builds, `accelerator_for` aliases to `thread_for`.

### Solver patterns

`SchurRedBlack` (`Grid/algorithms/iterative/SchurRedBlack.h`) implements red-black (even/odd) preconditioning for fermion operators. Most production fermion solves use `SchurRedBlackDiagMooeeSolve` or similar wrappers that internally call a `ConjugateGradient` on the Schur complement.

Mixed-precision solvers (`ConjugateGradientMixedPrec`, `BiCGSTABMixedPrec`) drive a double-precision outer loop with single-precision inner solves.

### Memory and I/O

- `Grid/allocator/` — aligned/NUMA-aware allocators; caching allocator via `--enable-alloc-cache`
- `Grid/parallelIO/` — distributed parallel reader/writer for ILDG (via LIME), SciDAC, and native binary formats
- `Grid/serialisation/` — text, binary, HDF5, XML/JSON serialisation of arbitrary Grid objects

### Executables

- `HMC/` — production HMC driver programmes (e.g. `Mobius2p1f.cc`, `DWF_plus_DSDR_nf2plus1_Shamir_Gparity.cc`)
- `benchmarks/` — `Benchmark_dwf`, `Benchmark_ITT`, `Benchmark_comms`, `Benchmark_memory_bandwidth`, … used to qualify a new machine
- `examples/` — small, readable programmes (`Example_plaquette.cc`, `Example_Mobius_spectrum.cc`) that are the best starting point for learning the API

Each of these directories auto-builds every top-level `.cc` as its own binary via `scripts/filelist`.

Every programme is wrapped in `Grid_init(&argc, &argv)` / `Grid_finalize()` (`Grid/util/Init.h`).

## Key Conventions

- **C++17** is required throughout.
- Template structure: most classes are templated on `<_FImpl>` (fermion impl) or `<Gimpl>` (gauge impl), which encode the representation and precision. Instantiation is controlled by `--enable-fermion-instantiations`.
- **Tensor indices are positional, not labelled.** The `Grid/tensors/` arithmetic recurses structurally over the `iScalar`/`iVector`/`iMatrix` nest: each level defines only the {scalar,vector,matrix}² products at its own level, with element types resolved by automatic type deduction, so every colour/spin/lorentz combination composes from ~200 lines (versus the pre-C++11 QDP++/PETE approach of machine-generating every case). An index's meaning derives entirely from its nesting depth counted from the outside; `iScalar` is the identity/broadcast case at every level. Never insert or remove a nesting level casually — the multiplication tables contract by position.
- **Multigrid coarsening deepens the tensor nest by one level.** A coarse site vector is `iVector<CComplex,nbasis>`, and `innerProduct` on it returns `iScalar<CComplex>` — one level deeper than the fine block scalar. So the block-inner-product scalar type gains one `iScalar` wrapper per MG level (fine: `vTComplex`; level 2: `iScalar<vTComplex>`; see `examples/Example_pvdagm_3level.cc`). When calling `blockInnerProduct`/`blockZAXPY`/`blockOrthogonalise` on coarse fields, the coarse scalar type must match `decltype(innerProduct(siteVector(),siteVector()))` exactly; a wrong depth fails to compile (no viable `operator=` deep in the instantiation chain) rather than mis-contracting.
- **Grids are borrowed, never owned.** `conformable` is pointer identity, so every object that interoperates must hold the *same* `GridCartesian *`; a class that minted its own grid internally could never conform with anything else. Ownership is therefore not available, and lifetime is managed by scope discipline instead of reference counting: whoever creates a grid retains it beyond every object it handed a reference to. Anything *derived* from a grid inherits this — `~PaddedCell` dereferences its `unpadded_grid`, so a `PaddedCell` cannot even be **destroyed** after its parent grid, only used. Where a consumer must let go early, it offers an explicit hand-back (`MultiGeneralCoarsenedOperatorV2::ReleaseGrid()`) to be called *before* the grid is destroyed.

- The `RealD`/`RealF`/`ComplexD`/`ComplexF` typedefs are used everywhere; avoid raw `double`/`float`.
- Use `GRID_ASSERT(cond)` (defined in `Grid/GridStd.h`), not bare `assert` — it prints a Grid-formatted message and aborts cleanly under MPI.
- Logging is stream-based, not macro-based: `std::cout << GridLogMessage << ... << std::endl;`. Channels declared in `Grid/log/Log.h` include `GridLogError`, `GridLogWarning`, `GridLogDebug`, `GridLogPerformance`, `GridLogIterative`, `GridLogSolver`, `GridLogHMC`, `GridLogComms`, `GridLogMemory`, `GridLogDslash`, `GridLogIRL`, `GridLogMG`. A subset is switched on at runtime with e.g. `--log Error,Warning,Message,Performance,Iterative,Integrator,Debug,Colours` (names given without the `GridLog` prefix).
- Performance-critical paths use `GRID_TRACE(name)` from `Grid/perfmon/Tracing.h` (compiled out unless `--enable-tracing` selects a backend) and the `GridStopWatch` timers in `Grid/perfmon/Timer.h`.
- Reductions across MPI ranks go through `GridBase::GlobalSum` / `GlobalMax`; never reduce with bare MPI calls inside library code.
- Everything lives in `NAMESPACE_BEGIN(Grid)` / `NAMESPACE_END(Grid)` macros; follow the surrounding file rather than writing `namespace Grid { }`.
- British spelling is used in identifiers and comments (`colour`, `neighbour`, `serialisation`).
