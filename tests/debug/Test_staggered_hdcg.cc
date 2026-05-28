/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: tests/debug/Test_staggered_hdcg.cc

    Authors: Thomas Blum, Peter Boyle

    HDCG (Hierarchical Deflation Conjugate Gradient) multigrid solver
    for naive staggered fermions, based on arXiv:2409.03904.

    Adapts the DWF HDCG infrastructure (Test_general_coarse_hdcg_phys48.cc) to:
      - NaiveStaggeredFermion (nearest-neighbour only, no Naik 3-hop term)
      - 4D SchurStaggeredOperator:  Mpc = m^2 - D_oe * D_eo  (hermitian, positive-definite)
      - vColourVector fine field type (staggered has colour but no spin)
      - NextToNearestStencilGeometry4D: 33-point coarse stencil

    Stencil count: D_oe*D_eo has 2-hop fine range.  With blocking B >= 2 the coarse
    shifts have L1-distance <= 2, giving 33 stencil points in 4D:
      1 (identity) + 8 (+-e_mu) + 24 (+-e_mu +- e_nu).
    NaiveStaggeredFermion has no Naik term, so any B >= 2 suffices.
    To extend to ImprovedStaggeredFermion later, use B >= 6.

    Reference: arXiv:2409.03904 (mrhs hermitian multigrid for DWF).

    Usage (after build):
      ./Test_staggered_hdcg --grid 16.16.16.16 --mpi 1.1.1.1

*************************************************************************************/
#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/AdefMrhs.h>
#include <Grid/algorithms/iterative/ImplicitlyRestartedBlockLanczos.h>
#include <Grid/algorithms/iterative/ImplicitlyRestartedBlockLanczosCoarse.h>

using namespace Grid;

// Non-converging CG used as a smoother (fixed number of iterations)
template<class Field>
class CGSmoother : public LinearFunction<Field>
{
public:
  typedef LinearOperatorBase<Field> FineOperator;
  FineOperator &_Op;
  int iters;
  CGSmoother(int _iters, FineOperator &Op) : _Op(Op), iters(_iters) {}
  void operator()(const Field &in, Field &out)
  {
    ConjugateGradient<Field> CG(0.0, iters, false);
    out = Zero();
    CG(_Op, in, out);
  }
};

int main(int argc, char **argv)
{
  fprintf(stderr, "TRACE: entering main\n"); fflush(stderr);
  Grid_init(&argc, &argv);
  fprintf(stderr, "TRACE: Grid_init done\n"); fflush(stderr);

  //--------------------------------------------------------------------
  // Parameters — tune for production
  //--------------------------------------------------------------------
  const int nbasis = 24;   // near-null space dimension
  const int cb     = 0;    // even checkerboard

  RealD mass = 0.05;

  // NaiveStaggeredFermion: nearest-neighbour hop only (no Naik term).
  // c1 = coefficient of the hopping term (1.0 = standard normalisation).
  // u0 = tadpole factor (1.0 = no tadpole improvement).
  RealD c1 = 1.0;
  RealD u0 = 1.0;

  //--------------------------------------------------------------------
  // Grids
  // Fine:   UGrid (4D full), UrbGrid (4D red-black)
  // Coarse: Coarse4d  with dimensions = GridDefaultLatt() / Block
  //
  // Recommended: GridDefaultLatt() >= 16^4, Block = {4,4,4,4}
  // NaiveStaggeredFermion works with any Block >= {2,2,2,2}
  //--------------------------------------------------------------------
  fprintf(stderr, "TRACE: making UGrid\n"); fflush(stderr);
  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  fprintf(stderr, "TRACE: making UrbGrid\n"); fflush(stderr);
  GridRedBlackCartesian *UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  Coordinate Block({4, 4, 4, 4});
  Coordinate clatt = GridDefaultLatt();
  for (int d = 0; d < clatt.size(); d++) clatt[d] /= Block[d];
  Coordinate csimd = GridDefaultSimd(Nd, vComplex::Nsimd());
  Coordinate cmpi  = GridDefaultMpi();
  fprintf(stderr, "TRACE: making Coarse4d clatt=%d %d %d %d simd=%d %d %d %d mpi=%d %d %d %d Nsimd=%d\n",
          clatt[0],clatt[1],clatt[2],clatt[3],
          csimd[0],csimd[1],csimd[2],csimd[3],
          cmpi[0],cmpi[1],cmpi[2],cmpi[3],
          (int)vComplex::Nsimd()); fflush(stderr);

  GridCartesian *Coarse4d = SpaceTimeGrid::makeFourDimGrid(clatt, csimd, cmpi);
  fprintf(stderr, "TRACE: Coarse4d made\n"); fflush(stderr);

  //--------------------------------------------------------------------
  // RNG + gauge field
  //--------------------------------------------------------------------
  fprintf(stderr, "TRACE: RNG4\n"); fflush(stderr);
  GridParallelRNG RNG4(UGrid);    RNG4.SeedFixedIntegers({1,2,3,4});
  fprintf(stderr, "TRACE: RNGrb\n"); fflush(stderr);
  GridParallelRNG RNGrb(UGrid); RNGrb.SeedFixedIntegers({5,6,7,8}); // must use full grid, not UrbGrid
  fprintf(stderr, "TRACE: Umu\n"); fflush(stderr);
  LatticeGaugeField Umu(UGrid);
  fprintf(stderr, "TRACE: HotConfig\n"); fflush(stderr);
  SU<Nc>::HotConfiguration(RNG4, Umu);
  fprintf(stderr, "TRACE: NaiveStaggeredFermionD\n"); fflush(stderr);
  NaiveStaggeredFermionD Ds(Umu, *UGrid, *UrbGrid, mass, c1, u0);
  fprintf(stderr, "TRACE: SchurStaggeredOperator\n"); fflush(stderr);
  SchurStaggeredOperator<NaiveStaggeredFermionD, LatticeStaggeredFermionD> HermOp(Ds);
  fprintf(stderr, "TRACE: HermOp done\n"); fflush(stderr);

  //--------------------------------------------------------------------
  // Subspace: inverse-iteration near-null vectors
  //
  // CreateSubspace applies CG (4 solves, tol=1e-4) to random noise vectors,
  // converging naturally to the low modes of HermOp without needing spectral
  // bound tuning.  Switch to CreateSubspaceChebyshevNew once the spectrum is
  // well characterised (hi ~ 5.0 for naive staggered SchurStaggeredOperator).
  //--------------------------------------------------------------------
  typedef Aggregation<vColourVector, vTComplex, nbasis> Subspace;
  Subspace Aggregates(Coarse4d, UrbGrid, cb);

  Aggregates.CreateSubspace(RNGrb, HermOp);
  Aggregates.Orthogonalise();

  //--------------------------------------------------------------------
  // Coarse geometry: NextToNearestStencilGeometry4D
  //   hops=2  ->  33 stencil points in 4D
  //--------------------------------------------------------------------
  NextToNearestStencilGeometry4D geom(Coarse4d);

  std::cout << GridLogMessage << "Coarse stencil: " << geom.npoint << " points" << std::endl;

  //--------------------------------------------------------------------
  // Single-RHS coarse operator (used for correctness check below)
  //--------------------------------------------------------------------
  typedef GeneralCoarsenedMatrix<vColourVector, vTComplex, nbasis> LittleDiracOp;
  typedef LittleDiracOp::CoarseVector CoarseVector;

  LittleDiracOp LDO(geom, UrbGrid, Coarse4d);
  LDO.CoarsenOperator(HermOp, Aggregates);

  //--------------------------------------------------------------------
  // Correctness check: P M_fine P^T c  ≈  M_coarse c
  //
  // Promote a random coarse vector into the fine subspace, apply the
  // fine operator, project back, and compare with the coarse operator
  // applied directly.  Error should be at the level of subspace
  // approximation quality (smaller = better basis vectors).
  //--------------------------------------------------------------------
  {
    GridParallelRNG RNGc(Coarse4d); RNGc.SeedFixedIntegers({9,10,11,12});
    CoarseVector c_src(Coarse4d), c_ldop(Coarse4d), c_proj(Coarse4d);
    random(RNGc, c_src);

    LatticeStaggeredFermionD f_v(UrbGrid), f_Mv(UrbGrid);
    Aggregates.PromoteFromSubspace(c_src, f_v);
    HermOp.Op(f_v, f_Mv);
    Aggregates.ProjectToSubspace(c_proj, f_Mv);

    LDO.M(c_src, c_ldop);

    c_proj -= c_ldop;
    RealD err = norm2(c_proj) / norm2(c_ldop);
    std::cout << GridLogMessage
              << "Coarsen check |P*M_fine - M_coarse| / |M_coarse| = " << err << std::endl;
  }

  //--------------------------------------------------------------------
  // Multi-RHS coarse grid
  //
  // The extra leading dimension holds nrhs right-hand sides packed into
  // SIMD lanes, matching the pattern of Test_general_coarse_hdcg_phys48.
  //--------------------------------------------------------------------
  const int nrhs = vComplex::Nsimd() * 2;

  Coordinate mpi   = GridDefaultMpi();
  Coordinate rhMpi ({1,    mpi[0],  mpi[1],  mpi[2],  mpi[3]});
  Coordinate rhLatt({nrhs, clatt[0], clatt[1], clatt[2], clatt[3]});
  Coordinate rhSimd({vComplex::Nsimd(), 1, 1, 1, 1});

  GridCartesian *CoarseMrhs = new GridCartesian(rhLatt, rhSimd, rhMpi);

  typedef MultiGeneralCoarsenedMatrix<vColourVector, vTComplex, nbasis> MultiCoarseOp;
  MultiCoarseOp mrhs(geom, CoarseMrhs);
  mrhs.CoarsenOperator(HermOp, Aggregates, Coarse4d);

  //--------------------------------------------------------------------
  // Coarse-grid Lanczos for deflation
  //--------------------------------------------------------------------
  typedef HermitianLinearOperator<MultiCoarseOp, CoarseVector> MrhsHermOp;
  MrhsHermOp MrhsCoarseOp(mrhs);

  // Estimate spectral bounds for Lanczos Chebyshev filter
  CoarseVector pm_src(CoarseMrhs); pm_src = ComplexD(1.0);
  PowerMethod<CoarseVector> cPM;
  RealD lambda_max = cPM(MrhsCoarseOp, pm_src);
  // Chebyshev filter window [lo, hi]:
  //   lo must sit in the spectral gap between the Nstop-th and (Nstop+1)-th
  //   coarse eigenvalues so that only the target modes receive cosh amplification.
  //
  // From a pilot run (16^4 fine, 4^4 coarse, mass=0.05, hot config):
  //   Group 1 (near-null, 24 modes): lambda in [0.002647, 0.002746]  ~= mass^2
  //   Spectral gap: factor 60 (lambda_24/lambda_23 = 0.165/0.00275)
  //   Group 2 (second group): lambda in [0.165, 0.179]
  //
  // lo = 0.02 sits in the spectral gap (factor 7x above lambda_23=0.00275,
  // factor 8x below lambda_24=0.165).
  //   hi = lambda_max_coarse * 1.1 ~= 2.121
  //   y(lambda_0=0.002647)  ~ -1.016 -> T_70 ~ 1.7e5  (cosh(70*0.182))
  //   y(lambda_23=0.002746) ~ -1.015 -> T_70 ~ 1.6e5
  //   Relative spread across near-null cluster: ~4.3%
  //   y(lambda_24=0.165)    ~ -0.862 -> inside [lo,hi] -> |T_70| <= 1
  //
  // order=71 (degree 70) is needed to give ~4% relative spread across the
  // near-null cluster of 24 nearly-degenerate eigenvalues; order=31 (tried)
  // gave only ~1.7% spread, insufficient for Nk=24/Nm=48 to converge.
  // Absolute amplification ~1e5; what matters for IRL convergence is the
  // relative spread, not the absolute value.
  // lo=0.005 failed (T_70~53, 0/24 modes in 10 restarts).
  // lo=0.01 worked but needed 2 restarts (13/24 then 24/24); lo=0.02 converges in 1.
  RealD lambda_lo  = 0.02;

  std::cout << GridLogMessage << "Chebyshev filter: lo=" << lambda_lo
            << " hi=" << lambda_max*1.1 << " order=71" << std::endl;

  Chebyshev<CoarseVector> IRLCheby(lambda_lo, lambda_max * 1.1, 71);

  // 24 near-null modes (eigenvalues ~mass^2) converge to resid^2~1e-28
  // in the first Lanczos restart.  The remaining modes (~0.165) are a
  // second spectral group that needs more Krylov vectors; handle them
  // separately once the basic HDCG solve is validated.
  int Nk    = 24;
  int Nm    = 48;
  int Nstop = Nk;

  GridParallelRNG CRNG(Coarse4d); CRNG.SeedFixedIntegers({13,14,15,16});

  ImplicitlyRestartedBlockLanczosCoarse<CoarseVector>
    IRL(MrhsCoarseOp, Coarse4d, CoarseMrhs, nrhs, IRLCheby,
        Nstop, /*conv_test_interval*/1, nrhs, Nk, Nm, 1.0e-5, 10);

  int Nconv;
  std::vector<RealD>        eval(Nm);
  std::vector<CoarseVector> evec(Nm,   Coarse4d);  // evec on f_grid (single-RHS coarse)
  std::vector<CoarseVector> c_srcs(nrhs, Coarse4d); // src on same grid as evec
  for (int r = 0; r < nrhs; r++) random(CRNG, c_srcs[r]);
  IRL.calc(eval, evec, c_srcs, Nconv, LanczosType::irbl);

  //--------------------------------------------------------------------
  // HDCG solver assembly
  //--------------------------------------------------------------------
  MultiRHSDeflation<CoarseVector> MrhsGuesser;
  MrhsGuesser.ImportEigenBasis(evec, eval);

  // MrhsProjector maps between fine (UrbGrid) and coarse (Coarse4d) spaces
  MultiRHSBlockProject<LatticeStaggeredFermionD> MrhsProjector;
  MrhsProjector.Allocate(nbasis, UrbGrid, Coarse4d);
  MrhsProjector.ImportBasis(Aggregates.subspace);

  ConjugateGradient<CoarseVector> CoarseCG(5.0e-2, 5000, false);
  DoNothingGuesser<CoarseVector>  DoNothing;
  HPDSolver<CoarseVector> HPDSolve(MrhsCoarseOp, CoarseCG, DoNothing);

  // Spectral radius of the fine operator, needed for the smoother shift.
  // Use a random checkerboard vector (UrbGrid) as starting guess for PowerMethod.
  LatticeStaggeredFermionD fine_pm_src(UrbGrid);
  random(RNGrb, fine_pm_src);
  PowerMethod<LatticeStaggeredFermionD> finePM;
  RealD fine_lambda_max = finePM(HermOp, fine_pm_src);

  // Shifted smoother: CG on (HermOp + shift*I) with shift = lambda_max / 100.
  //
  // The O(8) CG polynomial has 8 roots. With this shift all 8 roots lie in the
  // interval [shift, lambda_max + shift] ~ [0.046, 4.65], so the polynomial
  // focuses entirely on the HIGH-frequency part of the spectrum and leaves
  // near-null modes (lambda << shift) essentially untouched (polynomial ~ 1 there).
  //
  // This is the right target because the coarse-grid correction always introduces
  // high-frequency spectral leakage: the blocked coarse-grid degrees of freedom
  // are piecewise constant across coarse cells and therefore have sharp edges at
  // cell boundaries (like lego-block edges). Smoothness is measured by the
  // covariant Dirac derivative, so promoting the coarse solution back to the
  // fine grid inevitably excites high-frequency components — just as a step
  // function always carries high-frequency Fourier content.  The smoother must
  // repair exactly these high modes.
  //
  // The smoother and the coarse-grid correction are applied alternately: together
  // they both lift the low eigenvalues and pull down the upper eigenvalues of the
  // composite preconditioned operator, reducing the condition number seen by the
  // outer HDCG iterations.
  //
  // DWF HDCG convention; using mass^2 = 0.0025 was far too small: it scattered
  // the 8 roots over [0.005, 4.6] and diluted their effect on the high modes.
  RealD smootherShift = fine_lambda_max / 200.0;
  std::cout << GridLogMessage << "Smoother shift: lambda_max_fine/200 = "
            << fine_lambda_max << "/200 = " << smootherShift << std::endl;
  ShiftedHermOpLinearOperator<LatticeStaggeredFermionD> ShiftedOp(HermOp, smootherShift);
  CGSmoother<LatticeStaggeredFermionD> smoother(8, ShiftedOp);

  TwoLevelADEF2mrhs<LatticeStaggeredFermionD, CoarseVector>
    HDCG(1.0e-8, 500,
         HermOp,
         smoother,
         HPDSolve,   // M1 (coarse correction)
         HPDSolve,   // Vstart (initial guess projection)
         MrhsProjector,
         MrhsGuesser,
         CoarseMrhs);

  //--------------------------------------------------------------------
  // Solve: nrhs right-hand sides simultaneously
  //--------------------------------------------------------------------
  std::vector<LatticeStaggeredFermionD> src(nrhs, UrbGrid);
  std::vector<LatticeStaggeredFermionD> sol(nrhs, UrbGrid);

  GridParallelRNG RNGrb2(UGrid); RNGrb2.SeedFixedIntegers({17,18,19,20}); // must use full grid, not UrbGrid
  for (int r = 0; r < nrhs; r++) {
    random(RNGrb2, src[r]);
    sol[r] = Zero();
  }

  //--------------------------------------------------------------------
  // Baseline: standard single-RHS CG on HermOp (no preconditioning)
  // Run before HDCG to establish the unpreconditioned iteration count
  // and wall-clock time for direct comparison.
  //--------------------------------------------------------------------
  {
    ConjugateGradient<LatticeStaggeredFermionD> CG(1.0e-8, 50000, false);
    std::vector<LatticeStaggeredFermionD> cg_sol(nrhs, UrbGrid);
    for (int r = 0; r < nrhs; r++) cg_sol[r] = Zero();

    RealD t0 = usecond();
    int total_iters = 0;
    for (int r = 0; r < nrhs; r++) {
      std::cout << GridLogMessage << "====== CG baseline RHS " << r
                << " ======" << std::endl;
      CG(HermOp, src[r], cg_sol[r]);
      total_iters += CG.IterationsToComplete;
    }
    RealD t1 = usecond();
    std::cout << GridLogMessage << "CG baseline: " << nrhs << " RHS, "
              << total_iters << " total iterations, "
              << (t1 - t0) / 1.0e6 << " s total, "
              << (t1 - t0) / 1.0e6 / nrhs << " s/RHS" << std::endl;
  }

  //--------------------------------------------------------------------
  // HDCG solve
  //--------------------------------------------------------------------
  HDCG(src, sol);

  Grid_finalize();
  return 0;
}
