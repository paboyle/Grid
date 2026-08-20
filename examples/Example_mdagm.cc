/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_mdagm.cc

    Copyright (C) 2023

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>

using namespace std;
using namespace Grid;

// Routes Op/AdjOp -> HermOp so that CoarsenOperator and CreateSubspace
// both see the HPD operator M†M rather than bare M.
template<class Field>
class HermOpAdaptor : public LinearOperatorBase<Field>
{
  LinearOperatorBase<Field> &wrapped;
public:
  HermOpAdaptor(LinearOperatorBase<Field> &wrapme) : wrapped(wrapme) {};
  void Op     (const Field &in, Field &out)   { wrapped.HermOp(in,out); }
  void HermOp (const Field &in, Field &out)   { wrapped.HermOp(in,out); }
  void AdjOp  (const Field &in, Field &out)   { wrapped.HermOp(in,out); }
  void OpDiag (const Field &in, Field &out)                  { GRID_ASSERT(0); }
  void OpDir  (const Field &in, Field &out,int dir,int disp) { GRID_ASSERT(0); }
  void OpDirAll(const Field &in, std::vector<Field> &out)    { GRID_ASSERT(0); }
  void HermOpAndNorm(const Field &in, Field &out, RealD &n1, RealD &n2) {
    wrapped.HermOp(in, out);
    ComplexD dot = innerProduct(in, out);
    n1 = real(dot);
    n2 = norm2(out);
  }
};

// Fixed-iteration CG smoother: runs exactly `iters` steps of CG on the
// shifted operator.  tolerance=0 so CG never exits early.
template<class Field>
class CGSmoother : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  typedef LinearOperatorBase<Field> FineOperator;
  FineOperator &_SmootherOperator;
  int iters;
  CGSmoother(int _iters, FineOperator &SmootherOperator)
    : _SmootherOperator(SmootherOperator), iters(_iters)
  {
    std::cout << GridLogMessage << " CGSmoother order " << iters << std::endl;
  }
  void operator()(const Field &in, Field &out)
  {
    ConjugateGradient<Field> CG(0.0, iters, false);
    out = Zero();
    CG(_SmootherOperator, in, out);
  }
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls    = 24;
  const int nbasis = 60;
  const int cb    = 0;
  RealD M5  = 1.8;
  RealD b   = 1.5;
  RealD c   = 0.5;

  RealD mass = 0.00078;
  { const char *e = getenv("MASS"); if (e && *e) mass = atof(e); }

  std::cout << GridLogMessage << "Mass: " << mass
            << "  Ls: " << Ls << "  b=" << b << " c=" << c << std::endl;
  std::cout << GridLogMessage << "nbasis: " << nbasis << std::endl;

  // ── Grids ──────────────────────────────────────────────────────────────
  std::vector<int> lat_size{48,48,48,96};

  GridCartesian         *UGrid   = SpaceTimeGrid::makeFourDimGrid(lat_size,
                                     GridDefaultSimd(Nd,vComplex::Nsimd()),
                                     GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         *FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian *FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  Coordinate Block({4,4,3,4});
  Coordinate clatt = lat_size;
  for (int d = 0; d < (int)clatt.size(); d++) clatt[d] /= Block[d];

  GridCartesian *Coarse4d = SpaceTimeGrid::makeFourDimGrid(clatt,
                              GridDefaultSimd(Nd,vComplex::Nsimd()),
                              GridDefaultMpi());
  GridCartesian *Coarse5d = SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);

  // ── RNGs ───────────────────────────────────────────────────────────────
  GridParallelRNG RNG5(FGrid);   RNG5.SeedFixedIntegers({5,6,7,8});
  GridParallelRNG RNG4(UGrid);   RNG4.SeedFixedIntegers({1,2,3,4});

  // ── Gauge field ────────────────────────────────────────────────────────
  LatticeGaugeField Umu(UGrid);
  FieldMetaData header;
  NerscIO::readConfiguration(Umu, header, std::string("/ccs/home/poare/ckpoint_lat.1000"));

  // ── Fermion operator ───────────────────────────────────────────────────
  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);

  MdagMLinearOperator<MobiusFermionD,LatticeFermionD> MdagMOp(Ddwf);
  HermOpAdaptor<LatticeFermionD> HermFineOp(MdagMOp);

  // ── Coarse geometry ────────────────────────────────────────────────────
  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector CoarseVector;
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis> Subspace;

  NextToNearestStencilGeometry5D geom(Coarse5d);

  // ── Power method: estimate upper end of M†M spectrum ───────────────────
  LatticeFermionD pm_src(FGrid); random(RNG5, pm_src);
  PowerMethod<LatticeFermionD> PM;
  RealD hi = PM(HermFineOp, pm_src);
  std::cout << GridLogMessage << "Power method: hi = " << hi << std::endl;

  // ── Smoother: fixed-iteration CG on (M†M + lo) ─────────────────────────
  // lo/hi ~ 2/95 matches the HDCG ratio; tune empirically.
  RealD lo  = hi / 40.0;
  int   ord = 12;
  std::cout << GridLogMessage << "Smoother shift lo = " << lo
            << "  order = " << ord << std::endl;
  ShiftedHermOpLinearOperator<LatticeFermionD> ShiftedFineOp(HermFineOp, lo);
  CGSmoother<LatticeFermionD> Smoother(ord, ShiftedFineOp);

  // ── Subspace via CG inverse iteration ──────────────────────────────────
  Subspace Aggregates(Coarse5d, FGrid, cb);
  Aggregates.CreateSubspace(RNG5, HermFineOp, nbasis);

  // ── Cheap coarse deflation: diagonalise W BEFORE block-GS ──────────────
  // Orthogonalise() applies block-local Gram-Schmidt which rotates subspace[i]
  // into orthonormal block-local combinations, destroying the near-null
  // property of individual vectors.  We must extract the coarse zero-mode
  // combinations χₖ = Σᵢ V[i,k] ψᵢ from the pre-GS near-null vectors first.
  std::vector<LatticeFermionD> chi(nbasis, FGrid);
  {
    std::vector<LatticeFermionD> &psi = Aggregates.subspace;
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

  // ── Coarse operator ────────────────────────────────────────────────────
  LittleDiracOperator LittleDiracOp(geom, FGrid, Coarse5d);
  LittleDiracOp.CoarsenOperator(HermFineOp, Aggregates);

  // ── Coarse linear operator ─────────────────────────────────────────────
  HermitianLinearOperator<LittleDiracOperator,CoarseVector> LinOpCoarse(LittleDiracOp);

  // ── Project χₖ to coarse grid; Rayleigh quotients give deflation evals ─
  // ProjectToSubspace uses the post-GS basis, correctly mapping the pre-GS
  // near-null combinations to coarse vectors via the block-local U†(x_c).
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

  // ── Coarse solve: CG + deflation guesser ──────────────────────────────
  ConjugateGradient<CoarseVector>  coarseCG(5.0e-2, 10000, false);
  DeflatedGuesser<CoarseVector>    coarseGuess(coarse_deflation_vecs, coarse_deflation_evals);
  HPDSolver<CoarseVector> CoarseSolve(LinOpCoarse, coarseCG, coarseGuess);

  // ── ADEF2 outer solve ──────────────────────────────────────────────────
  LatticeFermionD src(FGrid);    random(RNG5, src);
  LatticeFermionD result(FGrid); result = Zero();

  TwoLevelADEF2<LatticeFermionD, CoarseVector, Subspace>
    HDCG(1.0e-8, 1000,
         HermFineOp,
         Smoother,
         CoarseSolve,   // used in PcgM1
         CoarseSolve,   // used in Vstart
         Aggregates);

  HDCG(src, result);

  // ── Reference RBCG ─────────────────────────────────────────────────────
#if 0
  {
    SchurDiagMooeeOperator<MobiusFermionD,LatticeFermion> HermOpEO(Ddwf);
    LatticeFermionD rb_src(FrbGrid);    random(RNG5, rb_src);
    LatticeFermionD rb_res(FrbGrid);    rb_res = Zero();
    ConjugateGradient<LatticeFermionD> CG(1.0e-8, 30000, false);
    CG(HermOpEO, rb_src, rb_res);
  }
#endif

  Grid_finalize();
  return 0;
}
