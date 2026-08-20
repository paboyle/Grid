---
name: multigrid-design-notes
description: "Design, development and tuning of LQCD multigrid solvers for physical-mass Möbius DWF on Frontier (AMD MI250X); covers HDCG, PVdagM two-level solver, coarse operator performance, and Lüscher deflation of the coarse solve"
metadata: 
  node_type: memory
  type: project
  originSessionId: cc1844e3-ab6f-4425-bf7e-a091b9554290
---

# LQCD Multigrid: Design, Development and Tuning

## Physical problem

Physical-mass Möbius DWF: Ls=24, b=1.5, c=0.5, M5=1.8, mass=0.00078.
48³×96 lattice. MPI geometry 3.6.4.4 (288 ranks / GCDs on Frontier).
Target: accelerate HMC fermion force and CG solves.

## Two solver paths

### Path 1: HDCG (TwoLevelADEF2 on MdagM)

- File: `examples/Example_mdagm.cc`
- Operator: MdagM (Hermitian positive definite); outer solver is ADEF2 CG.
- Subspace: 60 near-null vectors via CG inverse iteration (`CreateSubspace`).
- Coarse geometry: block {4,4,3,4}, coarse lattice 12×12×16×24, Ls_coarse=1.
- Coarse operator: `GeneralCoarsenedMatrix` (Petrov-Galerkin, npoint=33, NextToNearestStencil).
- Smoother: fixed-iteration CG on shifted operator (M†M + lo), lo=hi/80, 20 iters.
- Coarse solve: CG with DeflatedGuesser using 60 chi deflation vectors.
- Chi vectors: extracted BEFORE block-GS by diagonalising W_ij=<ψ_i|M†M|ψ_j>, computing
  chi_k = Σ_i V[i,k] ψ_i. These are global near-null combinations; block-GS destroys this.
- Deflation effect: 1089 coarse CG iters (undeflated) → 329 with 60 chi vectors.
- Best result: 274 outer ADEF2 iters, ~400s on 288 ranks Frontier.
- Coarse MVM performance: 541.7 GFlop/s kernel, 1119.6 GB/s (70% HBM), 97% roofline.
  MPI latency = 1188 μs = 37% of 3.2 ms per coarse MVM call. Single-RHS is bandwidth-bound.

### Path 2: PVdagM two-level PGCR

- File: `examples/Example_pvdagm.cc` and `examples/Example_pvdagm_defl.cc`
- Operator: PVdagM = PV†M (non-Hermitian); outer solver is PGCR.
- PV is the Pauli-Villars (mass=1) Möbius operator. PVdagM has exact zero modes.
- Subspace: 60 near-null vectors via GCR inverse iteration (`CreateSubspaceGCR`).
  GCR setup is slow: each vector takes O(600) PGCR steps, total ~4100s setup on Frontier.
- Coarse geometry: block 2^4 (lattice halved each dim), Ls_coarse=1.
- Coarse operator: `GeneralCoarsenedMatrix` with non-Hermitian coarsening.
- Preconditioner: `MGPreconditioner` V-cycle (pre-smooth, project, coarse solve, promote, post-smooth).
- Smoother: PGCR on ShiftedPVdagM (shift=0.01).
- Baseline (no deflation, 5e-2 coarse tol): 59 outer iters, 300s solve time.
  Coarse solve: 4.58s/call, 250 PGCR steps, NEVER converges ("did not converge" every call).
- Outer iteration count vs coarse tolerance: 5e-2→59 iters, 1e-1→63 iters, 3e-2→34 iters.
  But at 3e-2 without deflation: 1000+ coarse PGCR steps per call (useless).

## Key implementation work: GeneralCoarsenedMatrix performance

`Grid/algorithms/multigrid/GeneralCoarsenedMatrix.h`:

1. **accelerator_barrier fix**: `acceleratorBarrier()` is not a Grid macro; correct call is
   `accelerator_barrier(dummy)` (takes a dummy argument). This caused SIGBUS on Frontier.

2. **Coalesced FT kernel**: In `CoarsenOperator`, the loop filling A_v[sss](i,j) was serialised
   over j. Changed to:
   ```cpp
   accelerator_for(sss, osites, nbasis, {
       int j = acceleratorSIMTlane(nbasis);
       A_v[sss](i,j) = FT_v[sss](j);
   });
   ```
   This gives coalesced HBM access (nbasis consecutive elements per warp lane).

3. **Batched CoarsenOperator**: Used `MultiRHSBlockProject` to batch all npoint=33 stencil
   directions in one GEMM call per basis vector, replacing serial blockProject calls.
   Reduced projection from dominant bottleneck to 12% of CoarsenOperator time.
   mat (linop applications) now dominates at 83%.

## Lüscher deflation of the coarse solve (Example_pvdagm_defl.cc)

### Theory (Lüscher arXiv:0706.2298, Section A.3)

For near-null vectors {ψ_s} of operator D, the Petrov-Galerkin initial guess is:
  guess = Ψ W⁻¹ Ψ† src
where W_st = <ψ_s|D|ψ_t> and Ψ is the matrix of ψ columns.

Condition <ψ_s | src - D*guess> = 0 gives W c = b, b_t = <ψ_t|src>.
No SVD needed — W is dense, invert directly (LU). The U,V from SVD are unitaries
within the ψ-basis and cancel in W⁻¹; direct inverse is cleaner.

### Diagnostic results (job 4948520, before deflation)

Fine projected matrix W (60×60):
- ||W|| = 0.02399 — near-null vectors are genuinely small.
- Singular values: range [0.00169, 0.00529], ratio ~3:1. Well-conditioned inverse.

Coarse null matrix C_kl = <P ψ_k | A_coarse | P ψ_l>:
- ||C|| = 0.02399 — **identical to ||W||**. Galerkin property is exact.
- ||C - C†|| / ||C|| = 2.59e-9 — **C is Hermitian to machine precision**.
  Despite PVdagM being non-Hermitian, the projected coarse matrix is numerically Hermitian.
- C singular values match W singular values exactly (coarsening faithful).

This means: the coarse near-null vectors (projections of fine ψ_k) are exact near-null
vectors of the coarse operator, AND the deflation can use a Hermitian eigensolver.

### Implementation

**CoarseDeflatedGuesser** (in Example_pvdagm_defl.cc, before MGPreconditioner):
```cpp
template<class Field>
class CoarseDeflatedGuesser : public LinearFunction<Field> {
    const std::vector<Field> &chi;   // coarse eigenvectors of C_sym
    const std::vector<RealD> &eval;  // eigenvalues
public:
    void operator()(const Field &src, Field &guess) {
        guess = Zero();
        for (int k = 0; k < chi.size(); k++)
            axpy(guess, TensorRemove(innerProduct(chi[k], src)) / eval[k], chi[k], guess);
    }
};
```

**Chi vector construction** (in runMG, after CoarsenOperator and C computation):
1. Project pre-GS fine subspace to coarse grid: psi_coarse[k] = P ψ_k
2. Compute C_sym = (C + C†)/2
3. SelfAdjointEigenSolver(C_sym) → eigenvalues lambda[k], eigenvectors V
4. chi_coarse[k] = Σ_i V[i,k] * psi_coarse[i]  (satisfies <chi_j|A_c|chi_k> = lambda_k δ_jk)
5. chi_eval[k] = lambda[k]

**MGPreconditioner** modified to hold `CoarseSolver &_CoarseGuesser`:
- Added as constructor parameter and member reference.
- In V-cycle: replaced `Csol = Zero()` with `_CoarseGuesser(Csrc, Csol)` before `_CoarseSolve`.

**Note on pre-GS subspace**: The subspace is copied at the top of runMG BEFORE CoarsenOperator
modifies it via block-GS. This pre-GS copy is essential for computing W (fine projected matrix)
and psi_coarse. After CoarsenOperator, AggregatesPD.subspace is the block-orthonormal basis.

### Expected benefit
Without deflation at 3e-2 coarse tolerance: 1000+ PGCR steps per call, never converges.
With deflation: near-null modes removed by guesser; PGCR sees well-conditioned complement.
Target: O(100-200) coarse PGCR steps per call, 34 outer iterations, large total time reduction.

## Memory considerations (Frontier, MI250X, 64 GB HBM per GCD)

Each fine fermion field: ~162 MB per rank (Ls=24, 48³×96/288 local vol, spincolour=12 complex).
60 fine fields (subspace): ~9.7 GB.
The pre-GS subspace COPY in runMG adds another 9.7 GB.
Chi/null vector construction: build one at a time into a temp field, project immediately.
Do NOT allocate 2×nbasis additional fine fields (would add 19 GB → OOM).

Coarse vectors (~2 MB each): 120 coarse vectors (psi_coarse + chi_coarse) = ~240 MB. Fine.

## Roofline analysis (coarse MVM, single-RHS)

Arithmetic intensity: ~0.5 flops/byte. Roofline crossover: 119.7 flops/byte.
Deeply memory-bandwidth bound. Peak HBM: 1600 GB/s. Achieved: 1119.6 GB/s (70%).
MPI latency: 1188 μs = 37% of 3.2 ms per call. Irreducible for single-RHS.
Multi-RHS is the fundamental solution for throughput, but cannot be used in HMC
(each trajectory has a new gauge field → new coarse operator).

## Three-level perspective

The chi deflation of the coarse solve can be viewed as a third level: coarsening all
the way to 1⁴ × nbasis. The 60 near-null vectors of the coarse operator span this
third level's null space. The Lüscher guesser is the exact inverse on this space.

## References

- HDCG paper (2014): arXiv:1409.xxxx (P. Boyle) — ADEF2 CG with coarse deflation.
- Lüscher (2007): arXiv:0706.2298 — non-Hermitian deflation, Section A.3.
- Physical DWF multigrid: arXiv:2409.03904.
- Grid library: github.com/paboyle/Grid.
