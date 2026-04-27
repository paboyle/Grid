# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

VTK-based visualisation and analysis tools for Grid lattice QCD eigenvector density and HMC force data. All programmes link against both Grid (for reading Scidac/ILDG lattice files) and VTK (for rendering).

## Build

```bash
cd /Users/peterboyle/QCD/AmSC/Grid/visualisation/build
cmake .. -DVTK_DIR=$HOME/QCD/vtk/VTK-9.4.2-install/lib/cmake/vtk-9.4
make <target>    # e.g. make ControlledVisualise5D
```

All executables are built as macOS bundles (`.app`) except `ForceAnalysis`, `FindPeak`, and `DumpField`.

## Programmes

### ControlledVisualise5D

Interactive VTK renderer for 5D DWF eigenvector density (`LatticeComplexD`). Driven via named pipe `/tmp/visualise_cmd`.

**Launch script**: `/Volumes/X9Pro/visualisation/Grid/visualisation/build/Hdwf_1_long/visualise_controlled.sh`

**Wire protocol** (one command per line to `/tmp/visualise_cmd`):

| Command | Effect |
|---------|--------|
| `file <N>` / `file +N` / `file -N` | Jump to file by index or relative |
| `slice <dim> <N>` / `+N` / `-N` | Set or shift slice coordinate in dimension dim |
| `spin <deg>` | Continuous azimuth rotation at deg/tick (100ms tick); `spin 0` stops |
| `azimuth <deg>` | Single azimuth rotation step |
| `elevation <deg>` | Single elevation rotation step |
| `zoom <factor>` | Camera dolly (>1 = in) |
| `iso <value>` | Isosurface threshold in RMS units |
| `status` | Print current state |
| `quit` | Exit |

**Dimension indices for 5D DWF grid** (`--grid 48.32.32.32.32`):

| dim | axis | size |
|-----|------|------|
| 0   | s (Ls) | 48 |
| 1   | x    | 32   |
| 2   | y    | 32   |
| 3   | z    | 32   |
| 4   | t    | 32   |

**MD time mapping** for trajectory 702 (241 files, τ=3.3–4.0):
- File index N → τ = 3.300000 + N × (1/480)
- τ → file index = round((τ − 3.3) × 480)

**Display axes**: `--xyz 0.3.4` shows s, z, t. The `--slice` argument sets initial values for all dims; dims not in `--xyz`, `--sum`, or `--loop` are the fixed slice dimensions (x=dim1, y=dim2 with `--xyz 0.3.4`).

**Spin**: Implemented via `g_spinDeg` global applied on every 100ms poll timer tick inside `CommandHandler::Execute()`. Does not flood the pipe.

### FindPeak

Reads a `LatticeComplexD` Scidac file, prints the top-N sites by real value to stderr.

```bash
./FindPeak --grid 48.32.32.32.32 --mpi 1.1.1.1.1 <file> 2>peaks.txt
```

Key result: At τ=3.670833 the tunneling hotsite on the s=0 wall is (x=21, y=24, z=2, t=23).

### ForceAnalysis

Reads 4D `LatticeComplexD` force snapshot files (Shuhei's snapshots at `/Volumes/X9Pro/visualisation/Shuhei/snapshots/`). Outputs TSV of RMS and hotsite value per file to stderr.

```bash
./ForceAnalysis --grid 32.32.32.32 --mpi 1.1.1 --hotsite 21.24.2.23 \
    <files...> 2>force.tsv 1>/dev/null
```

Force components: `Gauge_lat`, `Gauge_smr`, `Jacobian_smr`, `Ferm0047_lat`, `Ferm0047_smr`.

### DumpField

Reads a `LatticeComplexD` and dumps via Grid's `<<` operator to stdout for verification.

### TranscriptToVideo

Renders a conversation transcript to an MP4 video (1280×720, 10 fps) with a typewriter animation effect, scrolling history, and optional captions. Does **not** link against Grid — pure VTK only.

#### Transcript format

```
[USER] First question text, possibly
continuing on the next line.

A blank line within a turn creates a paragraph break (visual spacer).

[ASSISTANT] Response text.
Multiple continuation lines are preserved
as separate display lines, not merged.

[CAPTION] Caption text shown at bottom of screen in white italic.
[CAPTION]    (whitespace-only body clears the caption)
[USER] Next question...
```

- Lines beginning `[USER]`, `[ASSISTANT]`, `[CAPTION]` start a new turn.
- Continuation lines (no `[TAG]` prefix) are joined with `\n` — each becomes its own wrapped display line.
- Blank lines within a turn become paragraph-break spacers.
- Markdown emphasis markers (`**`, `*`, `` ` ``) are stripped automatically.
- UTF-8 smart quotes, em-dashes, ellipses, arrows are transliterated to ASCII.

#### Usage

```bash
cd /Users/peterboyle/QCD/AmSC/Grid/visualisation/build
# Set runtime library paths first (see Runtime Environment below)
./TranscriptToVideo <transcript_file> <output.mp4>
```

Transcript files live in `/Users/peterboyle/QCD/AmSC/Grid/visualisation/` (e.g. `transcript`, `transcript2`, `transcript3`).

#### Visual layout

| Element | Detail |
|---------|--------|
| Background | Near-black navy `(0.04, 0.04, 0.10)` |
| `[USER]` text | Gold `(1.00, 0.84, 0.00)` |
| `[ASSISTANT]` text | Steel blue `(0.68, 0.85, 0.90)` |
| History | Up to 18 lines; brightness fades linearly from 0.85 (newest) to 0.20 (oldest) |
| Caption | Arial italic 20pt white with shadow, centred at bottom |
| Progress bar | Blue, top of frame |
| Typewriter speed | 50 chars/sec (5 chars/frame at 10 fps) |
| Pause between lines | 3 frames (0.3 s) |
| Word-wrap column | 60 chars (body only, after prefix) |

#### Key implementation notes

- **Persistent render context**: a single `vtkRenderWindow` is created once and reused for all frames. Creating a new window per frame exhausts the macOS Metal GPU command buffer after ~33 frames (`MTLCommandBufferErrorDomain Code=8`).
- **`SanitiseASCII()`**: replaces multi-byte UTF-8 sequences before passing to VTK's font renderer (which crashes on non-ASCII input).
- Output format is MP4 via `vtkFFMPEGWriter`. `SetOffScreenRendering(1)` is required for headless rendering.

## Runtime Environment

All executables in `build/` require Spack-installed HDF5/FFTW/GMP/MPFR on the dynamic linker path:

```bash
SPACK=/Users/peterboyle/QCD/Spack/spack/opt/spack/darwin-m1
export DYLD_LIBRARY_PATH=\
$SPACK/hdf5-1.14.6-2265ms4kymgw6hcnwi6vqehslyfv74t4/lib:\
$SPACK/fftw-3.3.10-aznn6h3nac5cycidlhrhgjxvntpcbg57/lib:\
$SPACK/gmp-6.3.0-cwiz4n7ww33fnb3aban2iup4orcr6c7i/lib:\
$SPACK/mpfr-4.2.1-exgbz4qshmet6tmmuttdewdlunfvtrlb/lib:\
$DYLD_LIBRARY_PATH
```

(These paths are also set by the ControlledVisualise5D launch script.)

## Key Physics Context

See `/Volumes/X9Pro/visualisation/analysis_notes_20260407.md` for full analysis. Summary:

- Near-zero mode of H_DWF localises on the two walls (s=0 and s=47) of the 5D domain wall geometry
- Topology change transfers norm between walls, mediated by a near-zero mode of H_w (Hermitian Wilson at m=−1.8)
- Tunneling hotsite on s=0 wall: (x=21, y=24, z=2, t=23); s=47 wall: (x=4, y=8, z=0, t=20)
- Light fermion pseudofermion force (Ferm0047_smr) peaks at ~20× RMS at the hotsite during tunneling — this is the restoring force that causes topological bounces

## Grid/VTK interaction notes

- Grid log messages go to stdout; all data output in analysis programmes uses stderr to avoid interleaving
- `TensorRemove()` is required when extracting a scalar from `peekSite()` result: `real(TensorRemove(peekSite(field, site)))`
- For runtime-determined grid dimensionality use `GridDefaultSimd(latt_size.size(), vComplex::Nsimd())`
- DYLD_LIBRARY_PATH must include Spack HDF5/FFTW/GMP/MPFR paths (see launch script)
