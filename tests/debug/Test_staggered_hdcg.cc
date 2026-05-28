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
  Grid_init(&argc, &argv);

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
  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  Coordinate Block({4, 4, 4, 4});
  Coordinate clatt = GridDefaultLatt();
  for (int d = 0; d < clatt.size(); d++) clatt[d] /= Block[d];

  GridCartesian *Coarse4d = SpaceTimeGrid::makeFourDimGrid(
      clatt, GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());

  //--------------------------------------------------------------------
  // RNG + gauge field
  //--------------------------------------------------------------------
  GridParallelRNG RNG4(UGrid);    RNG4.SeedFixedIntegers({1,2,3,4});
  GridParallelRNG RNGrb(UrbGrid); RNGrb.SeedFixedIntegers({5,6,7,8});

  LatticeGaugeField Umu(UGrid);
  // Hot gauge start for testing; replace with NerscIO::readConfiguration for production
  SU<Nc>::HotConfiguration(RNG4, Umu);

  //--------------------------------------------------------------------
  // Naive staggered fermion + Schur hermitian operator
  //
  // NaiveStaggeredFermion: nearest-neighbour hops only, no Naik term.
  // SchurStaggeredOperator::Op = HermOp = Mpc = m^2 - D_oe * D_eo
  // Positive-definite hermitian; CG and HDCG apply directly.
  // No HermOpAdaptor wrapper needed (Op == HermOp for SchurStaggeredOperator).
  //--------------------------------------------------------------------
  NaiveStaggeredFermionD Ds(Umu, *UGrid, *UrbGrid, mass, c1, u0);
  SchurStaggeredOperator<NaiveStaggeredFermionD, LatticeStaggeredFermionD> HermOp(Ds);

  //--------------------------------------------------------------------
  // Subspace: Chebyshev-filtered near-null vectors
  //
  // hi = spectral upper bound of HermOp = m^2 + lambda_max(D_oe*D_eo).
  // For naive staggered (c1=1), lambda_max(SchurStaggeredOperator) ~ 4.6.
  // Use 5.0 to give a small margin; refine with PowerMethod for each config.
  // PowerMethod for production.  Progressive Chebyshev filtering
  // (ChebyshevNew) uses two-pass polynomial filter.
  //--------------------------------------------------------------------
  typedef Aggregation<vColourVector, vTComplex, nbasis> Subspace;
  Subspace Aggregates(Coarse4d, UrbGrid, cb);

  RealD hi = 5.0;
  Aggregates.CreateSubspaceChebyshevNew(RNGrb, HermOp, hi);
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
  RealD lambda_lo  = mass * mass * 0.5;

  Chebyshev<CoarseVector> IRLCheby(lambda_lo, lambda_max * 1.1, 71);

  int Nk    = 48;
  int Nm    = 96;
  int Nstop = Nk;

  GridParallelRNG CRNG(CoarseMrhs); CRNG.SeedFixedIntegers({13,14,15,16});

  ImplicitlyRestartedBlockLanczosCoarse<CoarseVector>
    IRL(MrhsCoarseOp, Coarse4d, CoarseMrhs, nrhs, IRLCheby,
        Nstop, /*conv_test_interval*/1, nrhs, Nk, Nm, 1.0e-5, 10);

  int Nconv;
  std::vector<RealD>        eval(Nm);
  std::vector<CoarseVector> evec(Nm, Coarse4d);
  std::vector<CoarseVector> c_srcs(nrhs, CoarseMrhs);
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

  // Shifted smoother: CG on (HermOp + shift) damps low modes efficiently.
  // Shift ~ m^2 sits just above the spectral gap.
  RealD smootherShift = mass * mass;
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

  GridParallelRNG RNGrb2(UrbGrid); RNGrb2.SeedFixedIntegers({17,18,19,20});
  for (int r = 0; r < nrhs; r++) {
    random(RNGrb2, src[r]);
    sol[r] = Zero();
  }

  HDCG(src, sol);

  Grid_finalize();
  return 0;
}
