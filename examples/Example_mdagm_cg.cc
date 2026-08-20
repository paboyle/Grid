/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_mdagm_cg.cc

    Copyright (C) 2024

Author: Peter Boyle <pboyle@bnl.gov>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

// Two-level multigrid preconditioned CG for M†M (full Möbius DWF operator).
// Single RHS, targeting HMC use case.
//
// Operator hierarchy:
//   Fine:   M†M, acted on via HermOpAdaptor so Op = HermOp = M†M
//   Coarse: 33-point GeneralCoarsenedMatrix (NextToNearest, 2-hop M†M)
//
// Setup:
//   1. Chebyshev filter (T_600 x T_2500) on M†M to build near-null subspace {ψᵢ}
//   2. Before block-local Gram-Schmidt, compute W_ij = <ψᵢ|M†M|ψⱼ> (nbasis×nbasis)
//      and diagonalize it. The eigenvectors of W are approximate zero modes of the
//      coarse operator; projecting them onto the coarse grid gives coarse deflation
//      vectors. Block GS destroys the near-null property of individual basis vectors
//      (they become local linear combinations), so this step must precede Orthogonalise.
//   3. Block-local Gram-Schmidt (Orthogonalise), then CoarsenOperator.
//   4. Project W eigenvectors χₖ = Σᵢ V[i,k] ψᵢ to coarse grid → deflation vectors.
//
// Solvers:
//   Smoother:    CG(nIter) on (M†M + lo*I)
//   Coarse:      DeflatedGuesser (nbasis W-eigenvectors) + CG to tolerance 5e-2
//   Outer:       PrecGeneralisedConjugateResidual with two-level V-cycle preconditioner

#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidual.h>
#include <Grid/algorithms/iterative/NormalEquations.h>

using namespace std;
using namespace Grid;

// Wraps any LinearOperatorBase so that Op = AdjOp = HermOp.
// Required when coarsening an HPD operator whose Op != HermOp
// (e.g. MdagMLinearOperator where Op=M, HermOp=M†M).
template<class Field>
class HermOpAdaptor : public LinearOperatorBase<Field>
{
  LinearOperatorBase<Field> &wrapped;
public:
  HermOpAdaptor(LinearOperatorBase<Field> &wrapme) : wrapped(wrapme) {}
  void Op     (const Field &in, Field &out)   { wrapped.HermOp(in, out); }
  void HermOp (const Field &in, Field &out)   { wrapped.HermOp(in, out); }
  void AdjOp  (const Field &in, Field &out)   { wrapped.HermOp(in, out); }
  void OpDiag (const Field &in, Field &out)                    { GRID_ASSERT(0); }
  void OpDir  (const Field &in, Field &out, int dir, int disp) { GRID_ASSERT(0); }
  void OpDirAll(const Field &in, std::vector<Field> &out)      { GRID_ASSERT(0); }
  void HermOpAndNorm(const Field &in, Field &out, RealD &n1, RealD &n2) {
    wrapped.HermOp(in, out);
    ComplexD dot = innerProduct(in, out);
    n1 = real(dot); n2 = norm2(out);
  }
};

// Fixed-iteration CG as a smoother (LinearFunction).
// Used as the IR-shifted smoother: solves (M†M + lo*I) x = b approximately.
template<class Field>
class CGSmoother : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  LinearOperatorBase<Field> &_op;
  int iters;
  CGSmoother(int _iters, LinearOperatorBase<Field> &op) : _op(op), iters(_iters) {
    std::cout << GridLogMessage << "CGSmoother order " << iters << std::endl;
  }
  void operator()(const Field &in, Field &out) {
    ConjugateGradient<Field> CG(0.0, iters, false);
    out = Zero();
    CG(_op, in, out);
  }
};

// Two-level V-cycle preconditioner (LinearFunction).
template<class Fobj, class CComplex, int nbasis>
class MGPreconditioner : public LinearFunction<Lattice<Fobj>>
{
public:
  using LinearFunction<Lattice<Fobj>>::operator();
  typedef Aggregation<Fobj, CComplex, nbasis>  Aggregates;
  typedef typename Aggregates::FineField        FineField;
  typedef typename Aggregates::CoarseVector     CoarseVector;
  typedef LinearOperatorBase<FineField>         FineOperator;
  typedef LinearFunction<FineField>             FineSmoother;
  typedef LinearOperatorBase<CoarseVector>      CoarseOperator;
  typedef LinearFunction<CoarseVector>          CoarseSolver;

  Aggregates    &_Aggregates;
  FineOperator  &_FineOperator;
  FineSmoother  &_PreSmoother;
  FineSmoother  &_PostSmoother;
  CoarseOperator &_CoarseOperator;
  CoarseSolver  &_CoarseSolve;

  MGPreconditioner(Aggregates     &Agg,
                   FineOperator   &Fine,
                   FineSmoother   &Pre,
                   FineSmoother   &Post,
                   CoarseOperator &CoarseOp,
                   CoarseSolver   &CoarseSolve)
    : _Aggregates(Agg), _FineOperator(Fine),
      _PreSmoother(Pre), _PostSmoother(Post),
      _CoarseOperator(CoarseOp), _CoarseSolve(CoarseSolve) {}

  virtual void operator()(const FineField &in, FineField &out)
  {
    GridBase *CoarseGrid = _Aggregates.CoarseGrid;
    CoarseVector Csrc(CoarseGrid);
    CoarseVector Csol(CoarseGrid);
    FineField vec1(in.Grid());
    FineField vec2(in.Grid());

    double t;

    // Pre-smooth
    out = Zero();
    t = -usecond();
    _PreSmoother(in, out);
    t += usecond();
    std::cout << GridLogMessage << "PreSmoother " << t/1e3 << " ms" << std::endl;

    // Residual -> coarse projection
    _FineOperator.Op(out, vec1);
    sub(vec1, in, vec1);
    t = -usecond();
    _Aggregates.ProjectToSubspace(Csrc, vec1);
    t += usecond();
    std::cout << GridLogMessage << "ProjectToSubspace " << t/1e3 << " ms" << std::endl;

    // Coarse solve
    Csol = Zero();
    t = -usecond();
    _CoarseSolve(Csrc, Csol);
    t += usecond();
    std::cout << GridLogMessage << "CoarseSolve " << t/1e3 << " ms" << std::endl;

    // Promote correction and accumulate
    t = -usecond();
    _Aggregates.PromoteFromSubspace(Csol, vec1);
    t += usecond();
    std::cout << GridLogMessage << "PromoteFromSubspace " << t/1e3 << " ms" << std::endl;
    add(out, out, vec1);

    // Post-smooth on updated residual
    _FineOperator.Op(out, vec1);
    sub(vec1, in, vec1);
    vec2 = Zero();
    t = -usecond();
    _PostSmoother(vec1, vec2);
    t += usecond();
    std::cout << GridLogMessage << "PostSmoother " << t/1e3 << " ms" << std::endl;
    add(out, out, vec2);
  }
};

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  // Physical parameters matching the 48^3 x 96 test ensemble.
  // Target HMC ensemble is 128^3 x 288, Ls=12 at a^{-1}=3.5 GeV.
  const int Ls    = 24;
  const int nbasis = 62;
  RealD mass = 0.00078;
  RealD M5   = 1.8;
  RealD b    = 1.5;
  RealD c    = 0.5;
  { const char *e = getenv("MASS"); if (e && *e) mass = atof(e); }

  std::cout << GridLogMessage << "MdagM multigrid CG: mass=" << mass
            << " Ls=" << Ls << " b=" << b << " c=" << c << std::endl;

  GridCartesian         *UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                     GridDefaultSimd(Nd, vComplex::Nsimd()),
                                     GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         *FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
  GridRedBlackCartesian *FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);

  // Coarse grid: divide by 2 in each dimension, matching Example_mdagm.cc.
  // Coarse dims {24,24,24,48} -> local {8,8,6,6} with MPI 3.3.4.8: all even,
  // compatible with any power-of-2 SIMD layout.
  Coordinate clatt = GridDefaultLatt();
  for (int d = 0; d < (int)clatt.size(); d++) clatt[d] /= 2;

  GridCartesian *Coarse4d = SpaceTimeGrid::makeFourDimGrid(clatt,
                              GridDefaultSimd(Nd, vComplex::Nsimd()),
                              GridDefaultMpi());
  GridCartesian *Coarse5d = SpaceTimeGrid::makeFiveDimGrid(1, Coarse4d);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG RNG5(FGrid); RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG RNG4(UGrid); RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid);
  FieldMetaData header;
  std::string file("ckpoint_lat.1000");
  NerscIO::readConfiguration(Umu, header, file);

  ///////////////////////////////////////////////////////////
  // Fine operator: M†M on the full 5D lattice
  // MdagMLinearOperator::Op = M (not M†M), HermOp = M†M.
  // HermOpAdaptor makes Op = HermOp = M†M everywhere.
  ///////////////////////////////////////////////////////////
  MobiusFermionD Ddwf(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, b, c);
  MdagMLinearOperator<MobiusFermionD, LatticeFermionD> MdagMOp(Ddwf);
  HermOpAdaptor<LatticeFermionD> HermFineOp(MdagMOp);

  ///////////////////////////////////////////////////////////
  // Coarse operator: 33-point stencil for 2-hop M†M
  ///////////////////////////////////////////////////////////
  typedef GeneralCoarsenedMatrix<vSpinColourVector, vTComplex, nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector CoarseVector;

  NextToNearestStencilGeometry5D geom(Coarse5d);
  LittleDiracOperator LittleDiracOp(geom, FGrid, Coarse5d);

  ///////////////////////////////////////////////////////////
  // Near-null subspace via Chebyshev filter on M†M.
  // After CreateSubspaceChebyshevNew, subspace[i] = ψᵢ are the
  // original near-null vectors. Block-local Gram-Schmidt (Orthogonalise)
  // will rotate them into orthonormal block basis, destroying the near-null
  // property of individual vectors while preserving the subspace span.
  // We must compute the cheap coarse deflation BEFORE Orthogonalise.
  ///////////////////////////////////////////////////////////
  typedef Aggregation<vSpinColourVector, vTComplex, nbasis> Subspace;
  Subspace Aggregates(Coarse5d, FGrid, 0);

  RealD hi = 92.0;
  Aggregates.CreateSubspaceChebyshevNew(RNG5, HermFineOp, hi);

  // Compute W_ij = <ψᵢ|M†M|ψⱼ> and diagonalise BEFORE Orthogonalise.
  // Orthogonalise() applies block-local Gram-Schmidt which rotates subspace[i]
  // into an orthonormal block basis, destroying the near-null property of
  // individual vectors. The original span is preserved, but we need the
  // pre-GS ψᵢ to form χₖ = Σᵢ V[i,k] ψᵢ (approximate coarse zero modes).
  std::vector<LatticeFermionD> chi(nbasis, FGrid);
  {
    std::vector<LatticeFermionD> &psi = Aggregates.subspace; // pre-GS near-null vecs
    LatticeFermionD tmp(FGrid);
    Eigen::MatrixXcd W = Eigen::MatrixXcd::Zero(nbasis, nbasis);
    for (int i = 0; i < nbasis; i++) {
      HermFineOp.Op(psi[i], tmp);
      for (int j = 0; j < nbasis; j++)
        W(j, i) = TensorRemove(innerProduct(psi[j], tmp));
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> esolver(W);
    for (int k = 0; k < nbasis; k++) {
      chi[k] = Zero();
      for (int i = 0; i < nbasis; i++)
        chi[k] += ComplexD(esolver.eigenvectors()(i, k)) * psi[i];
    }
  }

  Aggregates.Orthogonalise();

  ///////////////////////////////////////////////////////////
  // Build coarse operator from subspace (post-GS basis)
  ///////////////////////////////////////////////////////////
  std::cout << GridLogMessage << "CoarsenOperator" << std::endl;
  LittleDiracOp.CoarsenOperator(HermFineOp, Aggregates);

  ///////////////////////////////////////////////////////////
  // Coarse linear operator (HPD: Op = HermOp = coarse M†M)
  ///////////////////////////////////////////////////////////
  HermitianLinearOperator<LittleDiracOperator, CoarseVector> LinOpCoarse(LittleDiracOp);

  ///////////////////////////////////////////////////////////
  // Project χₖ onto the coarse grid using the post-GS subspace.
  // ProjectToSubspace implicitly applies U†(x_c) at each block,
  // correctly mapping the pre-GS near-null combinations to coarse vectors.
  // Rayleigh quotients on the coarse operator give accurate eigenvalues.
  ///////////////////////////////////////////////////////////
  std::vector<CoarseVector> coarse_deflation_vecs(nbasis, Coarse5d);
  std::vector<RealD>        coarse_deflation_evals(nbasis);
  {
    CoarseVector Ac(Coarse5d);
    for (int k = 0; k < nbasis; k++) {
      Aggregates.ProjectToSubspace(coarse_deflation_vecs[k], chi[k]);
      RealD n = norm2(coarse_deflation_vecs[k]);
      coarse_deflation_vecs[k] *= 1.0 / std::sqrt(n);
      LinOpCoarse.HermOp(coarse_deflation_vecs[k], Ac);
      coarse_deflation_evals[k] = real(TensorRemove(innerProduct(coarse_deflation_vecs[k], Ac)));
      std::cout << GridLogMessage << "Coarse deflation eval[" << k << "] = "
                << coarse_deflation_evals[k] << std::endl;
    }
  }

  ///////////////////////////////////////////////////////////
  // Coarse CG solver with deflation guesser
  ///////////////////////////////////////////////////////////
  ConjugateGradient<CoarseVector> coarseCG(5.0e-2, 10000, false);
  DeflatedGuesser<CoarseVector>   coarseGuess(coarse_deflation_vecs, coarse_deflation_evals);
  HPDSolver<CoarseVector> CoarseSolve(LinOpCoarse, coarseCG, coarseGuess);

  ///////////////////////////////////////////////////////////
  // Fine-grid smoother: fixed-iteration CG on (M†M + lo*I)
  // lo=2.0 is the IR-shift high-pass parameter from mrhs-HDCG.
  ///////////////////////////////////////////////////////////
  RealD lo           = 2.0;
  int   smootherIters = 8;
  ShiftedHermOpLinearOperator<LatticeFermionD> ShiftedOp(HermFineOp, lo);
  CGSmoother<LatticeFermionD> Smoother(smootherIters, ShiftedOp);

  ///////////////////////////////////////////////////////////
  // Two-level V-cycle multigrid preconditioner
  ///////////////////////////////////////////////////////////
  MGPreconditioner<vSpinColourVector, vTComplex, nbasis>
    TwoLevelPrecon(Aggregates, HermFineOp, Smoother, Smoother, LinOpCoarse, CoarseSolve);

  ///////////////////////////////////////////////////////////
  // Outer solver: preconditioned GCR on M†M
  ///////////////////////////////////////////////////////////
  LatticeFermionD src(FGrid);    random(RNG5, src);
  LatticeFermionD result(FGrid); result = Zero();

  PrecGeneralisedConjugateResidual<LatticeFermionD>
    PGCR(1.0e-8, 10000, HermFineOp, TwoLevelPrecon, 16, 16);

  double t = -usecond();
  PGCR(src, result);
  t += usecond();
  std::cout << GridLogMessage << "Total solve time: " << t/1e6 << " s" << std::endl;

  Grid_finalize();
  return 0;
}
