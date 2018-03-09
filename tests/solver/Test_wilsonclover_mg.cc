/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/solver/Test_wilsonclover_mg.cc

    Copyright (C) 2017

    Author: Daniel Richtmann <daniel.richtmann@ur.de>

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
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidual.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

template<class Field, int nbasis> class TestVectorAnalyzer {
public:
  void operator()(LinearOperatorBase<Field> &Linop, std::vector<Field> const &vectors, int nn = nbasis) {

    auto positiveOnes = 0;

    std::vector<Field> tmp(4, vectors[0]._grid);
    Gamma              g5(Gamma::Algebra::Gamma5);

    std::cout << GridLogMessage << "Test vector analysis:" << std::endl;

    for(auto i = 0; i < nn; ++i) {

      Linop.Op(vectors[i], tmp[3]);

      tmp[0] = g5 * tmp[3];

      auto lambda = innerProduct(vectors[i], tmp[0]) / innerProduct(vectors[i], vectors[i]);

      tmp[1] = tmp[0] - lambda * vectors[i];

      auto mu = ::sqrt(norm2(tmp[1]) / norm2(vectors[i]));

      auto nrm = ::sqrt(norm2(vectors[i]));

      if(real(lambda) > 0)
        positiveOnes++;

      std::cout << GridLogMessage << std::scientific << std::setprecision(2) << std::setw(2) << std::showpos << "vector " << i << ": "
                << "singular value: " << lambda << ", singular vector precision: " << mu << ", norm: " << nrm << std::endl;
    }
    std::cout << GridLogMessage << std::scientific << std::setprecision(2) << std::setw(2) << std::showpos << positiveOnes << " out of "
              << nn << " vectors were positive" << std::endl;
  }
};

class myclass : Serializable {
public:
  // clang-format off
  GRID_SERIALIZABLE_CLASS_MEMBERS(myclass,
                                  int, domaindecompose,
                                  int, domainsize,
                                  int, coarsegrids,
                                  int, order,
                                  int, Ls,
                                  double, mq,
                                  double, lo,
                                  double, hi,
                                  int, steps);
  // clang-format on
  myclass(){};
};
myclass params;

template<int nbasis> struct CoarseGrids {
public:
  std::vector<std::vector<int>> LattSizes;
  std::vector<std::vector<int>> Seeds;
  std::vector<GridCartesian *>  Grids;
  std::vector<GridParallelRNG>  PRNGs;

  CoarseGrids(std::vector<std::vector<int>> const &blockSizes, int coarsegrids) {

    assert(blockSizes.size() == coarsegrids);

    std::cout << GridLogMessage << "Constructing " << coarsegrids << " CoarseGrids" << std::endl;

    for(int cl = 0; cl < coarsegrids; ++cl) { // may be a bit ugly and slow but not perf critical
      // need to differentiate between first and other coarse levels in size calculation
      LattSizes.push_back({cl == 0 ? GridDefaultLatt() : LattSizes[cl - 1]});
      Seeds.push_back(std::vector<int>(LattSizes[cl].size()));

      for(int d = 0; d < LattSizes[cl].size(); ++d) {
        LattSizes[cl][d] = LattSizes[cl][d] / blockSizes[cl][d];
        Seeds[cl][d]     = (cl + 1) * LattSizes[cl].size() + d + 1;
        // calculation unimportant, just to get. e.g., {5, 6, 7, 8} for first coarse level and so on
      }

      Grids.push_back(SpaceTimeGrid::makeFourDimGrid(LattSizes[cl], GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi()));
      PRNGs.push_back(GridParallelRNG(Grids[cl]));

      PRNGs[cl].SeedFixedIntegers(Seeds[cl]);

      std::cout << GridLogMessage << "cl = " << cl << ": LattSize = " << LattSizes[cl] << std::endl;
      std::cout << GridLogMessage << "cl = " << cl << ":    Seeds = " << Seeds[cl] << std::endl;
    }
  }
};

template<class Field> void testLinearOperator(LinearOperatorBase<Field> &LinOp, GridBase *Grid, std::string const &name = "") {

  std::vector<int> seeds({1, 2, 3, 4});
  GridParallelRNG  RNG(Grid);
  RNG.SeedFixedIntegers(seeds);

  {
    std::cout << GridLogMessage << "Testing that Mdiag + Σ_μ Mdir_μ == M for operator " << name << ":" << std::endl;

    // clang-format off
    Field src(Grid);    random(RNG, src);
    Field ref(Grid);    ref    = zero;
    Field result(Grid); result = zero;
    Field diag(Grid);   diag   = zero;
    Field sumDir(Grid); sumDir = zero;
    Field tmp(Grid);
    Field err(Grid);
    // clang-format on

    std::cout << setprecision(9);

    std::cout << GridLogMessage << " norm2(src)\t\t\t\t= " << norm2(src) << std::endl;

    LinOp.OpDiag(src, diag);
    std::cout << GridLogMessage << " norm2(Mdiag * src)\t\t\t= " << norm2(diag) << std::endl;

    for(int dir = 0; dir < 4; dir++) {
      for(auto disp : {+1, -1}) {
        LinOp.OpDir(src, tmp, dir, disp);
        std::cout << GridLogMessage << " norm2(Mdir_{" << dir << "," << disp << "} * src)\t\t= " << norm2(tmp) << std::endl;
        sumDir = sumDir + tmp;
      }
    }
    std::cout << GridLogMessage << " norm2(Σ_μ Mdir_μ * src)\t\t= " << norm2(sumDir) << std::endl;

    result = diag + sumDir;
    std::cout << GridLogMessage << " norm2((Mdiag + Σ_μ Mdir_μ) * src)\t= " << norm2(result) << std::endl;

    LinOp.Op(src, ref);
    std::cout << GridLogMessage << " norm2(M * src)\t\t\t= " << norm2(ref) << std::endl;

    err = ref - result;
    std::cout << GridLogMessage << " Absolute deviation\t\t\t= " << norm2(err) << std::endl;
    std::cout << GridLogMessage << " Relative deviation\t\t\t= " << norm2(err) / norm2(ref) << std::endl;
  }

  {
    std::cout << GridLogMessage << "Testing hermiticity stochastically for operator " << name << ":" << std::endl;

    // clang-format off
    Field phi(Grid); random(RNG, phi);
    Field chi(Grid); random(RNG, chi);
    Field MPhi(Grid);
    Field MdagChi(Grid);
    // clang-format on

    LinOp.Op(phi, MPhi);
    LinOp.AdjOp(chi, MdagChi);

    ComplexD chiMPhi    = innerProduct(chi, MPhi);
    ComplexD phiMdagChi = innerProduct(phi, MdagChi);

    ComplexD phiMPhi    = innerProduct(phi, MPhi);
    ComplexD chiMdagChi = innerProduct(chi, MdagChi);

    std::cout << GridLogMessage << " chiMPhi = " << chiMPhi << " phiMdagChi = " << phiMdagChi
              << " difference = " << chiMPhi - conjugate(phiMdagChi) << std::endl;

    std::cout << GridLogMessage << " phiMPhi = " << phiMPhi << " chiMdagChi = " << chiMdagChi << " <- should be real if hermitian"
              << std::endl;
  }

  {
    std::cout << GridLogMessage << "Testing linearity for operator " << name << ":" << std::endl;

    // clang-format off
    Field phi(Grid); random(RNG, phi);
    Field chi(Grid); random(RNG, chi);
    Field phiPlusChi(Grid);
    Field MPhi(Grid);
    Field MChi(Grid);
    Field MPhiPlusChi(Grid);
    Field linearityError(Grid);
    // clang-format on

    LinOp.Op(phi, MPhi);
    LinOp.Op(chi, MChi);

    phiPlusChi = phi + chi;

    LinOp.Op(phiPlusChi, MPhiPlusChi);

    linearityError = MPhiPlusChi - MPhi;
    linearityError = linearityError - MChi;

    std::cout << GridLogMessage << " norm2(linearityError) = " << norm2(linearityError) << std::endl;
  }
}

// template < class Fobj, class CComplex, int coarseSpins, int nbasis, class Matrix >
// class MultiGridPreconditioner : public LinearFunction< Lattice< Fobj > > {
template<class Fobj, class CComplex, int nbasis, class Matrix> class MultiGridPreconditioner : public LinearFunction<Lattice<Fobj>> {
public:
  typedef Aggregation<Fobj, CComplex, nbasis>     Aggregates;
  typedef CoarsenedMatrix<Fobj, CComplex, nbasis> CoarseOperator;

  typedef typename Aggregation<Fobj, CComplex, nbasis>::siteVector   siteVector;
  typedef typename Aggregation<Fobj, CComplex, nbasis>::CoarseScalar CoarseScalar;
  typedef typename Aggregation<Fobj, CComplex, nbasis>::CoarseVector CoarseVector;
  typedef typename Aggregation<Fobj, CComplex, nbasis>::CoarseMatrix CoarseMatrix;
  typedef typename Aggregation<Fobj, CComplex, nbasis>::FineField    FineField;
  typedef LinearOperatorBase<FineField>                              FineOperator;

  Aggregates &    _Aggregates;
  CoarseOperator &_CoarseOperator;
  Matrix &        _FineMatrix;
  FineOperator &  _FineOperator;
  Matrix &        _SmootherMatrix;
  FineOperator &  _SmootherOperator;

  // Constructor
  MultiGridPreconditioner(Aggregates &    Agg,
                          CoarseOperator &Coarse,
                          FineOperator &  Fine,
                          Matrix &        FineMatrix,
                          FineOperator &  Smooth,
                          Matrix &        SmootherMatrix)
    : _Aggregates(Agg)
    , _CoarseOperator(Coarse)
    , _FineOperator(Fine)
    , _FineMatrix(FineMatrix)
    , _SmootherOperator(Smooth)
    , _SmootherMatrix(SmootherMatrix) {}

  void operator()(const FineField &in, FineField &out) {

    CoarseVector coarseSrc(_CoarseOperator.Grid());
    CoarseVector coarseTmp(_CoarseOperator.Grid());
    CoarseVector coarseSol(_CoarseOperator.Grid());
    coarseSol = zero;

    GeneralisedMinimalResidual<CoarseVector> coarseGMRES(5.0e-2, 100, 25, false);
    GeneralisedMinimalResidual<FineField>    fineGMRES(5.0e-2, 100, 25, false);

    HermitianLinearOperator<CoarseOperator, CoarseVector> coarseHermOp(_CoarseOperator);
    MdagMLinearOperator<CoarseOperator, CoarseVector>     coarseMdagMOp(_CoarseOperator);
    MdagMLinearOperator<Matrix, FineField>                fineMdagMOp(_SmootherMatrix);

    FineField fineTmp1(in._grid);
    FineField fineTmp2(in._grid);

    RealD Ni = norm2(in);

    // no pre smoothing for now
    auto  preSmootherNorm     = 0;
    auto  preSmootherResidual = 0;
    RealD r;

    // Project to coarse grid, solve, project back to fine grid
    _Aggregates.ProjectToSubspace(coarseSrc, in);
    coarseGMRES(coarseMdagMOp, coarseSrc, coarseSol);
    _Aggregates.PromoteFromSubspace(coarseSol, out);

    // Recompute error
    _FineOperator.Op(out, fineTmp1);
    fineTmp1            = in - fineTmp1;
    r                   = norm2(fineTmp1);
    auto coarseResidual = std::sqrt(r / Ni);

    // Apply smoother, use GMRES for the moment
    fineGMRES(fineMdagMOp, in, out);

    // Recompute error
    _FineOperator.Op(out, fineTmp1);
    fineTmp1                  = in - fineTmp1;
    r                         = norm2(fineTmp1);
    auto postSmootherResidual = std::sqrt(r / Ni);

    std::cout << GridLogIterative << "Input norm = " << Ni << " Pre-Smoother norm " << preSmootherNorm
              << " Pre-Smoother residual = " << preSmootherResidual << " Coarse residual = " << coarseResidual
              << " Post-Smoother residual = " << postSmootherResidual << std::endl;
  }

  void runChecks(CoarseGrids<nbasis> &cGrids, int whichCoarseGrid) {

    /////////////////////////////////////////////
    // Some stuff we need for the checks below //
    /////////////////////////////////////////////
    auto tolerance = 1e-13; // TODO: this obviously depends on the precision we use, current value is for double

    std::vector<CoarseVector> cTmps(4, _CoarseOperator.Grid());
    std::vector<FineField>    fTmps(2, _Aggregates.subspace[0]._grid); // atm only for one coarser grid

    // need to construct an operator, since _CoarseOperator is not a LinearOperator but only a matrix (the name is a bit misleading)
    MdagMLinearOperator<CoarseOperator, CoarseVector> MdagMOp(_CoarseOperator);

    std::cout << GridLogMessage << "**************************************************" << std::endl;
    std::cout << GridLogMessage << "MG correctness check: 0 == (1 - P R) v" << std::endl;
    std::cout << GridLogMessage << "**************************************************" << std::endl;

    for(auto i = 0; i < _Aggregates.subspace.size(); ++i) {
      _Aggregates.ProjectToSubspace(cTmps[0], _Aggregates.subspace[i]); //   R v_i
      _Aggregates.PromoteFromSubspace(cTmps[0], fTmps[0]);              // P R v_i

      fTmps[1]       = _Aggregates.subspace[i] - fTmps[0]; // v_i - P R v_i
      auto deviation = std::sqrt(norm2(fTmps[1]) / norm2(_Aggregates.subspace[i]));

      std::cout << GridLogMessage << "Vector " << i << ": norm2(v_i) = " << norm2(_Aggregates.subspace[i])
                << " | norm2(R v_i) = " << norm2(cTmps[0]) << " | norm2(P R v_i) = " << norm2(fTmps[0])
                << " | relative deviation = " << deviation;

      if(deviation > tolerance) {
        std::cout << " > " << tolerance << " -> check failed" << std::endl;
        // abort();
      } else {
        std::cout << " < " << tolerance << " -> check passed" << std::endl;
      }
    }

    std::cout << GridLogMessage << "**************************************************" << std::endl;
    std::cout << GridLogMessage << "MG correctness check: 0 == (1 - R P) v_c" << std::endl;
    std::cout << GridLogMessage << "**************************************************" << std::endl;

    random(cGrids.PRNGs[whichCoarseGrid], cTmps[0]);

    _Aggregates.PromoteFromSubspace(cTmps[0], fTmps[0]); //   P v_c
    _Aggregates.ProjectToSubspace(cTmps[1], fTmps[0]);   // R P v_c

    cTmps[2]       = cTmps[0] - cTmps[1]; // v_c - R P v_c
    auto deviation = std::sqrt(norm2(cTmps[2]) / norm2(cTmps[0]));

    std::cout << GridLogMessage << "norm2(v_c) = " << norm2(cTmps[0]) << " | norm2(R P v_c) = " << norm2(cTmps[1])
              << " | norm2(P v_c) = " << norm2(fTmps[0]) << " | relative deviation = " << deviation;

    if(deviation > tolerance) {
      std::cout << " > " << tolerance << " -> check failed" << std::endl;
      // abort();
    } else {
      std::cout << " < " << tolerance << " -> check passed" << std::endl;
    }

    std::cout << GridLogMessage << "**************************************************" << std::endl;
    std::cout << GridLogMessage << "MG correctness check: 0 == (R D P - D_c) v_c" << std::endl;
    std::cout << GridLogMessage << "**************************************************" << std::endl;

    random(cGrids.PRNGs[whichCoarseGrid], cTmps[0]);

    _Aggregates.PromoteFromSubspace(cTmps[0], fTmps[0]); //     P v_c
    _FineOperator.Op(fTmps[0], fTmps[1]);                //   D P v_c
    _Aggregates.ProjectToSubspace(cTmps[1], fTmps[1]);   // R D P v_c

    MdagMOp.Op(cTmps[0], cTmps[2]); // D_c v_c

    cTmps[3]  = cTmps[1] - cTmps[2]; // R D P v_c - D_c v_c
    deviation = std::sqrt(norm2(cTmps[3]) / norm2(cTmps[1]));

    std::cout << GridLogMessage << "norm2(R D P v_c) = " << norm2(cTmps[1]) << " | norm2(D_c v_c) = " << norm2(cTmps[2])
              << " | relative deviation = " << deviation;

    if(deviation > tolerance) {
      std::cout << " > " << tolerance << " -> check failed" << std::endl;
      // abort();
    } else {
      std::cout << " < " << tolerance << " -> check passed" << std::endl;
    }

    std::cout << GridLogMessage << "**************************************************" << std::endl;
    std::cout << GridLogMessage << "MG correctness check: 0 == |(Im(v_c^dag D_c^dag D_c v_c)|" << std::endl;
    std::cout << GridLogMessage << "**************************************************" << std::endl;

    random(cGrids.PRNGs[whichCoarseGrid], cTmps[0]);

    MdagMOp.Op(cTmps[0], cTmps[1]);    //         D_c v_c
    MdagMOp.AdjOp(cTmps[1], cTmps[2]); // D_c^dag D_c v_c

    auto dot  = innerProduct(cTmps[0], cTmps[2]); //v_c^dag D_c^dag D_c v_c
    deviation = abs(imag(dot)) / abs(real(dot));

    std::cout << GridLogMessage << "Re(v_c^dag D_c^dag D_c v_c) = " << real(dot) << " | Im(v_c^dag D_c^dag D_c v_c) = " << imag(dot)
              << " | relative deviation = " << deviation;

    if(deviation > tolerance) {
      std::cout << " > " << tolerance << " -> check failed" << std::endl;
      // abort();
    } else {
      std::cout << " < " << tolerance << " -> check passed" << std::endl;
    }
  }
};

int main(int argc, char **argv) {

  Grid_init(&argc, &argv);

  params.domainsize      = 1;
  params.coarsegrids     = 1;
  params.domaindecompose = 0;
  params.order           = 30;
  params.Ls              = 1;
  params.mq              = -0.5;
  params.lo              = 0.5;
  params.hi              = 70.0;
  params.steps           = 1;

  typedef typename WilsonCloverFermionR::FermionField FermionField;
  typename WilsonCloverFermionR::ImplParams           wcImplparams;
  WilsonAnisotropyCoefficients                        wilsonAnisCoeff;

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Params: " << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  std::cout << params << std::endl;

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Set up some fine level stuff: " << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  GridCartesian *FGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *FrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(FGrid);

  std::vector<int> fSeeds({1, 2, 3, 4});
  GridParallelRNG  fPRNG(FGrid);
  fPRNG.SeedFixedIntegers(fSeeds);

  Gamma g5(Gamma::Algebra::Gamma5);

  // clang-format off
  FermionField      src(FGrid); gaussian(fPRNG, src);
  FermionField   result(FGrid); result = zero;
  LatticeGaugeField Umu(FGrid); SU3::HotConfiguration(fPRNG, Umu);
  // clang-format on

  RealD mass = params.mq;

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Set up some coarser levels stuff: " << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  const int nbasis = 20; // fix the number of test vector to the same
                         // number on every level for now

  //////////////////////////////////////////
  // toggle to run two/three level method
  //////////////////////////////////////////

  // two-level algorithm
  std::vector<std::vector<int>> blockSizes({{2, 2, 2, 2}});
  CoarseGrids<nbasis>           coarseGrids(blockSizes, 1);

  // // three-level algorithm
  // std::vector<std::vector<int>> blockSizes({{2, 2, 2, 2}, {2, 2, 1, 1}});
  // CoarseGrids<nbasis>           coarseGrids(blockSizes, 2);

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Some typedefs" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  // typedefs for transition from fine to first coarsened grid
  typedef vSpinColourVector                                                                       FineSiteVector;
  typedef vTComplex                                                                               CoarseSiteScalar;
  typedef Aggregation<FineSiteVector, CoarseSiteScalar, nbasis>                                   Subspace;
  typedef CoarsenedMatrix<FineSiteVector, CoarseSiteScalar, nbasis>                               CoarseOperator;
  typedef CoarseOperator::CoarseVector                                                            CoarseVector;
  typedef CoarseOperator::siteVector                                                              CoarseSiteVector;
  typedef TestVectorAnalyzer<FermionField, nbasis>                                                FineTVA;
  typedef MultiGridPreconditioner<FineSiteVector, CoarseSiteScalar, nbasis, WilsonCloverFermionR> FineMGPreconditioner;
  typedef TrivialPrecon<FermionField>                                                             FineTrivialPreconditioner;

  // typedefs for transition from a coarse to the next coarser grid (some defs remain the same)
  typedef Aggregation<CoarseSiteVector, CoarseSiteScalar, nbasis>                             SubSubSpace;
  typedef CoarsenedMatrix<CoarseSiteVector, CoarseSiteScalar, nbasis>                         CoarseCoarseOperator;
  typedef CoarseCoarseOperator::CoarseVector                                                  CoarseCoarseVector;
  typedef CoarseCoarseOperator::siteVector                                                    CoarseCoarseSiteVector;
  typedef TestVectorAnalyzer<CoarseVector, nbasis>                                            CoarseTVA;
  typedef MultiGridPreconditioner<CoarseSiteVector, CoarseSiteScalar, nbasis, CoarseOperator> CoarseMGPreconditioner;
  typedef TrivialPrecon<CoarseVector>                                                         CoarseTrivialPreconditioner;

  static_assert(std::is_same<CoarseVector, CoarseCoarseVector>::value, "CoarseVector and CoarseCoarseVector must be of the same type");

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Building the wilson clover operator on the fine grid" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  RealD                csw_r = 1.0;
  RealD                csw_t = 1.0;
  WilsonCloverFermionR Dwc(Umu, *FGrid, *FrbGrid, mass, csw_r, csw_t, wilsonAnisCoeff, wcImplparams);

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Setting up linear operators" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  MdagMLinearOperator<WilsonCloverFermionR, FermionField> FineMdagMOp(Dwc);

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Calling Aggregation class to build subspaces" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  Subspace FineAggregates(coarseGrids.Grids[0], FGrid, 0);

  assert((nbasis & 0x1) == 0);
  int nb = nbasis / 2;
  std::cout << GridLogMessage << " nbasis/2 = " << nb << std::endl;

  FineAggregates.CreateSubspace(fPRNG, FineMdagMOp /*, nb */); // Don't specify nb to see the orthogonalization check

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Test vector analysis after initial creation of subspace" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  FineTVA fineTVA;
  fineTVA(FineMdagMOp, FineAggregates.subspace);

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Projecting subspace to definite chirality" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  for(int n = 0; n < nb; n++) {
    FineAggregates.subspace[n + nb] = g5 * FineAggregates.subspace[n];
  }

  auto coarseSites = 1;
  for(auto const &elem : coarseGrids.LattSizes[0]) coarseSites *= elem;

  std::cout << GridLogMessage << "Norms of MG test vectors after chiral projection (coarse sites = " << coarseSites << ")" << std::endl;
  for(int n = 0; n < nbasis; n++) {
    std::cout << GridLogMessage << "vec[" << n << "] = " << norm2(FineAggregates.subspace[n]) << std::endl;
  }

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Building coarse representation of Dirac operator" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  CoarseOperator Dc(*coarseGrids.Grids[0]);

  Dc.CoarsenOperator(FGrid, FineMdagMOp, FineAggregates);

  MdagMLinearOperator<CoarseOperator, CoarseVector> CoarseMdagMOp(Dc);

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Test vector analysis after construction of coarse Dirac operator" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  fineTVA(FineMdagMOp, FineAggregates.subspace);

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Testing the linear operators" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  // clang-format off
  testLinearOperator(FineMdagMOp,   FGrid,                "FineMdagMOp");   std::cout << GridLogMessage << std::endl;
  testLinearOperator(CoarseMdagMOp, coarseGrids.Grids[0], "CoarseMdagMOp"); std::cout << GridLogMessage << std::endl;
  // clang-format on

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Building coarse vectors" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  CoarseVector coarseSource(coarseGrids.Grids[0]);
  CoarseVector coarseResult(coarseGrids.Grids[0]);
  gaussian(coarseGrids.PRNGs[0], coarseSource);
  coarseResult = zero;

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Building some coarse space solvers" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  std::vector<std::unique_ptr<OperatorFunction<CoarseVector>>> dummyCoarseSolvers;
  dummyCoarseSolvers.emplace_back(new GeneralisedMinimalResidual<CoarseVector>(5.0e-2, 100, 8, false));
  dummyCoarseSolvers.emplace_back(new MinimalResidual<CoarseVector>(5.0e-2, 100, 0.8, false));
  dummyCoarseSolvers.emplace_back(new ConjugateGradient<CoarseVector>(5.0e-2, 100, false));

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Testing some coarse space solvers" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  std::cout << GridLogMessage << "checking norm of coarse src " << norm2(coarseSource) << std::endl;

  for(auto const &solver : dummyCoarseSolvers) {
    coarseResult = zero;
    (*solver)(CoarseMdagMOp, coarseSource, coarseResult);
  }

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Building a multigrid preconditioner" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  FineMGPreconditioner      FineMGPrecon(FineAggregates, Dc, FineMdagMOp, Dwc, FineMdagMOp, Dwc);
  FineTrivialPreconditioner FineSimplePrecon;

  FineMGPrecon.runChecks(coarseGrids, 0);

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Building krylov subspace solvers w/ & w/o MG Preconditioner" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  std::vector<std::unique_ptr<OperatorFunction<FermionField>>> solvers;
  solvers.emplace_back(new FlexibleGeneralisedMinimalResidual<FermionField>(1.0e-12, 4000000, FineSimplePrecon, 25, false));
  solvers.emplace_back(new FlexibleGeneralisedMinimalResidual<FermionField>(1.0e-12, 100, FineMGPrecon, 25, false));
  solvers.emplace_back(new PrecGeneralisedConjugateResidual<FermionField>(1.0e-12, 4000000, FineSimplePrecon, 25, 25));

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Testing the (un)?preconditioned solvers" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  for(auto const &solver : solvers) {
    std::cout << GridLogMessage << "checking norm of fine src " << norm2(src) << std::endl;
    result = zero;
    (*solver)(FineMdagMOp, src, result);
    std::cout << std::endl;
  }

#if 0
  if(coarseGrids.LattSizes.size() == 2) {

    std::cout << GridLogMessage << "**************************************************" << std::endl;
    std::cout << GridLogMessage << "Some testing for construction of a second coarse level" << std::endl;
    std::cout << GridLogMessage << "**************************************************" << std::endl;

    // std::cout << GridLogMessage << "**************************************************" << std::endl;
    // std::cout << GridLogMessage << "Calling Aggregation class to build subspaces" << std::endl;
    // std::cout << GridLogMessage << "**************************************************" << std::endl;

    SubSubSpace CoarseAggregates(coarseGrids.Grids[1], coarseGrids.Grids[0], 0);
    CoarseAggregates.CreateSubspace(coarseGrids.PRNGs[0], CoarseMdagMOp);

    // std::cout << GridLogMessage << "**************************************************" << std::endl;
    // std::cout << GridLogMessage << "Test vector analysis after initial creation of subspace" << std::endl;
    // std::cout << GridLogMessage << "**************************************************" << std::endl;

    // // this doesn't work because this function applies g5 to a vector, which
    // // doesn't work for coarse vectors atm -> FIXME
    // CoarseTVA coarseTVA;
    // coarseTVA(CoarseMdagMOp, CoarseAggregates.subspace);

    // std::cout << GridLogMessage << "**************************************************" << std::endl;
    // std::cout << GridLogMessage << "Projecting subspace to definite chirality" << std::endl;
    // std::cout << GridLogMessage << "**************************************************" << std::endl;

    // // cannot apply g5 to coarse vectors atm -> FIXME
    // for(int n=0;n<nb;n++){
    //   CoarseAggregates.subspace[n+nb] = g5 * CoarseAggregates.subspace[n];
    //   std::cout<<GridLogMessage<<n<<" subspace "<<norm2(CoarseAggregates.subspace[n+nb])<<" "<<norm2(CoarseAggregates.subspace[n]) <<std::endl;
    // }

    auto coarseCoarseSites = 1;
    for(auto const &elem : coarseGrids.LattSizes[1]) coarseCoarseSites *= elem;

    std::cout << GridLogMessage << "Norms of MG test vectors after chiral projection (coarse coarse sites = " << coarseCoarseSites << ")"
              << std::endl;
    for(int n = 0; n < nbasis; n++) {
      std::cout << GridLogMessage << "vec[" << n << "] = " << norm2(CoarseAggregates.subspace[n]) << std::endl;
    }

    // std::cout << GridLogMessage << "**************************************************" << std::endl;
    // std::cout << GridLogMessage << "Building coarse coarse representation of Dirac operator" << std::endl;
    // std::cout << GridLogMessage << "**************************************************" << std::endl;

    CoarseCoarseOperator Dcc(*coarseGrids.Grids[1]);

    Dcc.CoarsenOperator(coarseGrids.Grids[0], CoarseMdagMOp, CoarseAggregates);

    MdagMLinearOperator<CoarseCoarseOperator, CoarseCoarseVector> CoarseCoarseMdagMOp(Dcc);

    // std::cout << GridLogMessage << "**************************************************" << std::endl;
    // std::cout << GridLogMessage << "Test vector analysis after construction of coarse Dirac operator" << std::endl;
    // std::cout << GridLogMessage << "**************************************************" << std::endl;

    // // this doesn't work because this function applies g5 to a vector, which
    // // doesn't work for coarse vectors atm -> FIXME
    // coarseTVA(CoarseMdagMOp, CoarseAggregates.subspace);

    // std::cout << GridLogMessage << "**************************************************" << std::endl;
    // std::cout << GridLogMessage << "Testing the linear operators" << std::endl;
    // std::cout << GridLogMessage << "**************************************************" << std::endl;

    // clang-format off
    testLinearOperator(CoarseMdagMOp,       coarseGrids.Grids[0], "CoarseMdagMOp");
    testLinearOperator(CoarseCoarseMdagMOp, coarseGrids.Grids[1], "CoarseCoarseMdagMOp");
    // clang-format on

    // std::cout << GridLogMessage << "**************************************************" << std::endl;
    // std::cout << GridLogMessage << "Building coarse coarse vectors" << std::endl;
    // std::cout << GridLogMessage << "**************************************************" << std::endl;

    CoarseCoarseVector coarseCoarseSource(coarseGrids.Grids[1]);
    CoarseCoarseVector coarseCoarseResult(coarseGrids.Grids[1]);
    gaussian(coarseGrids.PRNGs[1], coarseCoarseSource);
    coarseCoarseResult = zero;

    // std::cout << GridLogMessage << "**************************************************" << std::endl;
    // std::cout << GridLogMessage << "Building some coarse space solvers" << std::endl;
    // std::cout << GridLogMessage << "**************************************************" << std::endl;

    std::vector<std::unique_ptr<OperatorFunction<CoarseCoarseVector>>> dummyCoarseCoarseSolvers;
    dummyCoarseCoarseSolvers.emplace_back(new GeneralisedMinimalResidual<CoarseCoarseVector>(5.0e-2, 100, 8, false));
    dummyCoarseCoarseSolvers.emplace_back(new MinimalResidual<CoarseCoarseVector>(5.0e-2, 100, 0.8, false));
    dummyCoarseCoarseSolvers.emplace_back(new ConjugateGradient<CoarseCoarseVector>(5.0e-2, 100, false));

    // std::cout << GridLogMessage << "**************************************************" << std::endl;
    // std::cout << GridLogMessage << "Testing some coarse coarse space solvers" << std::endl;
    // std::cout << GridLogMessage << "**************************************************" << std::endl;

    std::cout << GridLogMessage << "checking norm of coarse coarse src " << norm2(coarseCoarseSource) << std::endl;

    for(auto const &solver : dummyCoarseCoarseSolvers) {
      coarseCoarseResult = zero;
      (*solver)(CoarseCoarseMdagMOp, coarseCoarseSource, coarseCoarseResult);
    }

    // std::cout << GridLogMessage << "**************************************************" << std::endl;
    // std::cout << GridLogMessage << "Building a multigrid preconditioner" << std::endl;
    // std::cout << GridLogMessage << "**************************************************" << std::endl;

    CoarseMGPreconditioner      CoarseMGPrecon(CoarseAggregates, Dcc, CoarseMdagMOp, Dc, CoarseMdagMOp, Dc);
    CoarseTrivialPreconditioner CoarseSimplePrecon;

    CoarseMGPrecon.runChecks(coarseGrids, 1);

    // std::cout << GridLogMessage << "**************************************************" << std::endl;
    // std::cout << GridLogMessage << "Building krylov subspace solvers w/ & w/o MG Preconditioner" << std::endl;
    // std::cout << GridLogMessage << "**************************************************" << std::endl;

    std::vector<std::unique_ptr<OperatorFunction<CoarseVector>>> solvers;
    solvers.emplace_back(new FlexibleGeneralisedMinimalResidual<CoarseVector>(1.0e-12, 4000000, CoarseSimplePrecon, 25, false));
    solvers.emplace_back(new FlexibleGeneralisedMinimalResidual<CoarseVector>(1.0e-12, 100, CoarseMGPrecon, 25, false));
    solvers.emplace_back(new PrecGeneralisedConjugateResidual<CoarseVector>(1.0e-12, 4000000, CoarseSimplePrecon, 25, 25));

    // std::cout << GridLogMessage << "**************************************************" << std::endl;
    // std::cout << GridLogMessage << "Testing the (un)?preconditioned solvers" << std::endl;
    // std::cout << GridLogMessage << "**************************************************" << std::endl;

    for(auto const &solver : solvers) {
      std::cout << GridLogMessage << "checking norm of fine src " << norm2(coarseSource) << std::endl;
      coarseResult = zero;
      (*solver)(CoarseMdagMOp, coarseSource, coarseResult);
      std::cout << std::endl;
    }

  }
#endif

  Grid_finalize();
}
