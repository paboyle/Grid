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
| `--enable-debug=yes` | adds `-g`, removes `-O3` |

Use `make V=1` for verbose compiler output (shows full flags; required for bug reports).

Platform recipes from `README.md`:
- **KNL**: `--enable-simd=KNL --enable-comms=mpi3-auto --enable-mkl`
- **Skylake/Haswell**: `--enable-simd=AVX512` or `AVX2` + `--enable-comms=mpi3-auto`
- **AMD EPYC**: `--enable-simd=AVX2 --enable-comms=mpi3`
- **A64FX (Fugaku)**: `--enable-simd=A64FX --enable-comms=mpi3 --enable-shm=shmget` (see `SVE_README.txt`)

Required external libs: GMP, MPFR, OpenSSL, zlib.

## Running Tests

```bash
# From build directory
make check                          # root-level tests (Test_simd, Test_cshift, etc.)
make -C tests/<subdir> tests        # build tests in a subdirectory
./tests/core/Test_simd              # run a single test binary directly
```

Test subdirectories and their focus: `core` (SIMD, stencil, comms), `solver` (CG, GMRES, eigensolvers), `hmc` (MD integrators), `forces` (fermion forces), `lanczos`, `IO`, `smearing`, `sp2n`, `debug`.

## Architecture

### Layer stack (bottom to top)

1. **SIMD layer** (`Grid/simd/`) — platform-specific intrinsics wrapped into `vRealF`, `vComplexD`, etc. The SIMD width and layout are compile-time constants controlled by `--enable-simd`.

2. **Tensor layer** (`Grid/tensors/`) — Lorentz/colour/spin tensor algebra built on top of SIMD types. `iMatrix`, `iVector`, `iScalar` templates compose into QCD types like `ColourMatrix`, `SpinColourVector`.

3. **Lattice layer** (`Grid/lattice/`) — `Lattice<T>` container: a site-local tensor replicated across a distributed Cartesian grid. All arithmetic is site-parallel and expression-template-fused.

4. **Cartesian/comms layer** (`Grid/cartesian/`, `Grid/communicator/`) — `GridCartesian` holds the MPI topology and local/global geometry. `Grid/cshift/` implements nearest-neighbour halo exchange; `Grid/stencil/` is the optimised multi-hop stencil used by Dirac operators.

5. **Algorithm layer** (`Grid/algorithms/`) — iterative solvers (CG, GMRES, BiCGSTAB, mixed-precision), eigensolvers (Lanczos, LAPACK), FFT, smearing.

6. **QCD layer** (`Grid/qcd/`) — gauge and fermion actions, HMC integrators, observables.

### QCD subsystem (`Grid/qcd/`)

- `action/fermion/` — Wilson, Clover, DWF (Mobius), Staggered, twisted-mass, G-parity variants
- `action/gauge/` — Wilson gauge, Symanzik, Iwasaki, DBW2, plaquette+rect
- `representations/` — Fundamental, Adjoint, Two-index, Sp(2n)
- `hmc/` — Leapfrog, OMF2/OMF4 integrators; pseudofermion refreshment; Metropolis accept/reject
- `smearing/` — APE, Stout, HEX, gradient flow
- `observables/` — Polyakov loop, plaquette, topological charge

### GPU acceleration

GPU support is injected via macros (`accelerator_for`, `accelerator_for2dNB`). The `Grid/simd/` SIMD types map to scalar on GPU device code; host code paths remain vectorised. Unified virtual memory is on by default (`--enable-unified=yes`); device-aware MPI (`--enable-accelerator-aware-mpi`) avoids device→host copies on transfers.

### Memory and I/O

- `Grid/allocator/` — aligned/NUMA-aware allocators; caching allocator via `--enable-alloc-cache`
- `Grid/parallelIO/` — distributed parallel reader/writer for ILDG (via LIME), SciDAC, and native binary formats
- `Grid/serialisation/` — text, binary, HDF5, XML/JSON serialisation of arbitrary Grid objects

### HMC applications

`HMC/` contains production-ready HMC driver programmes (e.g. `Mobius2p1f.cc`, `DWF_plus_DSDR_nf2plus1_Shamir_Gparity.cc`). These are built separately from the library tests.

## Key Conventions

- **C++17** is required throughout.
- Template structure: most classes are templated on `<_FImpl>` (fermion impl) or `<Gimpl>` (gauge impl), which encode the representation and precision. Instantiation is controlled by `--enable-fermion-instantiations`.
- The `RealD`/`RealF`/`ComplexD`/`ComplexF` typedefs are used everywhere; avoid raw `double`/`float`.
- Logging uses `Grid_log`, `Grid_error` macros (from `Grid/log/`); performance-critical paths use the `GRID_TRACE` / timer macros from `Grid/perfmon/`.
- Reductions across MPI ranks go through `GridBase::GlobalSum` / `GlobalMax`; never reduce with bare MPI calls inside library code.

## Skills

`skills/` contains seven user-invocable Claude Code skills encoding deep domain knowledge for HPC work in this repo. Invoke with `/skill-name` or ask Claude to use them by name:

| Skill | When to use |
|-------|-------------|
| `gpu-memory-performance` | Bandwidth/occupancy problems — `acceleratorThreads()` pitfalls, `coalescedRead/Write`, fused vs staged HBM patterns |
| `gpu-runtime-correctness` | Wrong answers, non-deterministic results, premature `q.wait()` returns |
| `communication-overlap` | Designing GPU+MPI overlap pipelines; replacing broken accelerator-aware MPI paths with host-staged 7-phase pipeline |
| `mpi-heterogeneous` | Collective hangs, buffer aliasing in `MPI_Sendrecv`, heterogeneous topology bugs |
| `hang-diagnosis` | Distinguishing ioctl hangs, infinite poll loops, collective deadlocks, and silent wrong-answer races |
| `correctness-verification` | Reproducibility checksums, double-wait testing, bisecting non-deterministic failures |
| `compiler-validation` | Confirming compiler/optimisation flags are safe before production runs |
