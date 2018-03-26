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

// clang-format off
struct MultiGridParams : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(MultiGridParams,
                                  int, nLevels,
                                  std::vector<std::vector<int>>, blockSizes,
                                  bool, kCycle);
  MultiGridParams(){};
};
MultiGridParams mgParams;
// clang-format on

struct LevelInfo {
public:
  std::vector<std::vector<int>> Seeds;
  std::vector<GridCartesian *>  Grids;
  std::vector<GridParallelRNG>  PRNGs;

  LevelInfo(GridCartesian *FineGrid, MultiGridParams const &mgParams) {

    auto nCoarseLevels = mgParams.blockSizes.size();

    assert(nCoarseLevels == mgParams.nLevels - 1);

    // set up values for finest grid
    Grids.push_back(FineGrid);
    Seeds.push_back({1, 2, 3, 4});
    PRNGs.push_back(GridParallelRNG(Grids.back()));
    PRNGs.back().SeedFixedIntegers(Seeds.back());

    // set up values for coarser grids
    for(int level = 1; level < mgParams.nLevels; ++level) {
      auto Nd  = Grids[level - 1]->_ndimension;
      auto tmp = Grids[level - 1]->_fdimensions;
      assert(tmp.size() == Nd);

      Seeds.push_back(std::vector<int>(Nd));

      for(int d = 0; d < Nd; ++d) {
        tmp[d] /= mgParams.blockSizes[level - 1][d];
        Seeds[level][d] = (level)*Nd + d + 1;
      }

      Grids.push_back(SpaceTimeGrid::makeFourDimGrid(tmp, GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi()));
      PRNGs.push_back(GridParallelRNG(Grids[level]));

      PRNGs[level].SeedFixedIntegers(Seeds[level]);
    }

    std::cout << GridLogMessage << "Constructed " << mgParams.nLevels << " levels" << std::endl;

    // The construction above corresponds to the finest level having level == 0
    // (simply because it's not as ugly to implement), but we need it the
    // other way round (i.e., the coarsest level to have level == 0) for the MG
    // Preconditioner -> reverse the vectors

    std::reverse(Seeds.begin(), Seeds.end());
    std::reverse(Grids.begin(), Grids.end());
    std::reverse(PRNGs.begin(), PRNGs.end());

    for(int level = 0; level < mgParams.nLevels; ++level) {
      std::cout << GridLogMessage << "level = " << level << ":" << std::endl;
      Grids[level]->show_decomposition();
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

template<class Fobj, class CoarseScalar, int nCoarseSpins, int nBasis, int level, class Matrix>
class MultiGridPreconditioner : public LinearFunction<Lattice<Fobj>> {
public:
  /////////////////////////////////////////////
  // Type Definitions
  /////////////////////////////////////////////

  typedef Aggregation<Fobj, CoarseScalar, nBasis>                                                                Aggregates;
  typedef CoarsenedMatrix<Fobj, CoarseScalar, nBasis>                                                            CoarseMatrix;
  typedef typename Aggregates::CoarseVector                                                                      CoarseVector;
  typedef typename Aggregates::siteVector                                                                        CoarseSiteVector;
  typedef Matrix                                                                                                 FineMatrix;
  typedef typename Aggregates::FineField                                                                         FineVector;
  typedef MultiGridPreconditioner<CoarseSiteVector, CoarseScalar, nCoarseSpins, nBasis, level - 1, CoarseMatrix> NextPreconditionerLevel;

  /////////////////////////////////////////////
  // Member Data
  /////////////////////////////////////////////

  MultiGridParams &                        _MultiGridParams;
  LevelInfo &                              _LevelInfo;
  FineMatrix &                             _FineMatrix;
  FineMatrix &                             _SmootherMatrix;
  Aggregates                               _Aggregates;
  CoarseMatrix                             _CoarseMatrix;
  std::unique_ptr<NextPreconditionerLevel> _NextPreconditionerLevel;

  /////////////////////////////////////////////
  // Member Functions
  /////////////////////////////////////////////

  MultiGridPreconditioner(MultiGridParams &mgParams, LevelInfo &LvlInfo, FineMatrix &FineMat, FineMatrix &SmootherMat)
    : _MultiGridParams(mgParams)
    , _LevelInfo(LvlInfo)
    , _FineMatrix(FineMat)
    , _SmootherMatrix(SmootherMat)
    , _Aggregates(_LevelInfo.Grids[level - 1], _LevelInfo.Grids[level], 0)
    , _CoarseMatrix(*_LevelInfo.Grids[level - 1]) {
    _NextPreconditionerLevel
      = std::unique_ptr<NextPreconditionerLevel>(new NextPreconditionerLevel(_MultiGridParams, _LevelInfo, _CoarseMatrix, _CoarseMatrix));
  }

  void setup() {

    Gamma                                       g5(Gamma::Algebra::Gamma5);
    MdagMLinearOperator<FineMatrix, FineVector> fineMdagMOp(_FineMatrix);

    _Aggregates.CreateSubspace(_LevelInfo.PRNGs[level], fineMdagMOp /*, nb */); // NOTE: Don't specify nb to see the orthogonalization check

    // TestVectorAnalyzer<FineVector, nbasis> fineTVA;
    // fineTVA(fineMdagMOp, _Aggregates.subspace);

    static_assert((nBasis & 0x1) == 0, "MG Preconditioner only supports an even number of basis vectors");
    int nb = nBasis / 2;

    // TODO: to get this to work for more than two levels, I would need to either implement coarse spins or have a template specialization of this class also for the finest level
    for(int n = 0; n < nb; n++) {
      _Aggregates.subspace[n + nb] = g5 * _Aggregates.subspace[n];
    }

    _CoarseMatrix.CoarsenOperator(_LevelInfo.Grids[level], fineMdagMOp, _Aggregates);

    _NextPreconditionerLevel->setup();
  }

  virtual void operator()(Lattice<Fobj> const &in, Lattice<Fobj> &out) {

    // TODO: implement a W-cycle
    if(_MultiGridParams.kCycle)
      kCycle(in, out);
    else
      vCycle(in, out);
  }

  void vCycle(Lattice<Fobj> const &in, Lattice<Fobj> &out) {

    RealD inputNorm = norm2(in);

    CoarseVector coarseSrc(_LevelInfo.Grids[level - 1]);
    CoarseVector coarseSol(_LevelInfo.Grids[level - 1]);
    coarseSol = zero;

    FineVector fineTmp(in._grid);

    TrivialPrecon<FineVector>                      fineTrivialPreconditioner;
    FlexibleGeneralisedMinimalResidual<FineVector> fineFGMRES(1.0e-14, 1, fineTrivialPreconditioner, 1, false);

    MdagMLinearOperator<FineMatrix, FineVector> fineMdagMOp(_FineMatrix);
    MdagMLinearOperator<FineMatrix, FineVector> fineSmootherMdagMOp(_SmootherMatrix);

    _Aggregates.ProjectToSubspace(coarseSrc, in);
    (*_NextPreconditionerLevel)(coarseSrc, coarseSol);
    _Aggregates.PromoteFromSubspace(coarseSol, out);

    fineMdagMOp.Op(out, fineTmp);
    fineTmp                                = in - fineTmp;
    auto r                                 = norm2(fineTmp);
    auto residualAfterCoarseGridCorrection = std::sqrt(r / inputNorm);

    fineFGMRES(fineSmootherMdagMOp, in, out);

    fineMdagMOp.Op(out, fineTmp);
    fineTmp                        = in - fineTmp;
    r                              = norm2(fineTmp);
    auto residualAfterPostSmoother = std::sqrt(r / inputNorm);

    std::cout << GridLogMG << " Level " << level << ": V-cycle: Input norm = " << std::sqrt(inputNorm)
              << " Coarse residual = " << residualAfterCoarseGridCorrection << " Post-Smoother residual = " << residualAfterPostSmoother
              << std::endl;
  }

  void kCycle(Lattice<Fobj> const &in, Lattice<Fobj> &out) {

    RealD inputNorm = norm2(in);

    CoarseVector coarseSrc(_LevelInfo.Grids[level - 1]);
    CoarseVector coarseSol(_LevelInfo.Grids[level - 1]);
    coarseSol = zero;

    FineVector fineTmp(in._grid);

    TrivialPrecon<FineVector>                        fineTrivialPreconditioner;
    FlexibleGeneralisedMinimalResidual<FineVector>   fineFGMRES(1.0e-14, 1, fineTrivialPreconditioner, 1, false);
    FlexibleGeneralisedMinimalResidual<CoarseVector> coarseFGMRES(1.0e-14, 1, *_NextPreconditionerLevel, 1, false);

    MdagMLinearOperator<FineMatrix, FineVector>     fineMdagMOp(_FineMatrix);
    MdagMLinearOperator<FineMatrix, FineVector>     fineSmootherMdagMOp(_SmootherMatrix);
    MdagMLinearOperator<CoarseMatrix, CoarseVector> coarseMdagMOp(_CoarseMatrix);

    _Aggregates.ProjectToSubspace(coarseSrc, in);
    coarseFGMRES(coarseMdagMOp, coarseSrc, coarseSol);
    _Aggregates.PromoteFromSubspace(coarseSol, out);

    fineMdagMOp.Op(out, fineTmp);
    fineTmp                                = in - fineTmp;
    auto r                                 = norm2(fineTmp);
    auto residualAfterCoarseGridCorrection = std::sqrt(r / inputNorm);

    fineFGMRES(fineSmootherMdagMOp, in, out);

    fineMdagMOp.Op(out, fineTmp);
    fineTmp                        = in - fineTmp;
    r                              = norm2(fineTmp);
    auto residualAfterPostSmoother = std::sqrt(r / inputNorm);

    std::cout << GridLogMG << " Level " << level << ": K-cycle: Input norm = " << std::sqrt(inputNorm)
              << " Coarse residual = " << residualAfterCoarseGridCorrection << " Post-Smoother residual = " << residualAfterPostSmoother
              << std::endl;
  }

  void runChecks() {

    auto tolerance   = 1e-13; // TODO: this obviously depends on the precision we use, current value is for double
    auto coarseLevel = level - 1;

    std::vector<FineVector>   fineTmps(2, _LevelInfo.Grids[level]);
    std::vector<CoarseVector> coarseTmps(4, _LevelInfo.Grids[level - 1]);

    MdagMLinearOperator<FineMatrix, FineVector>     fineMdagMOp(_FineMatrix);
    MdagMLinearOperator<CoarseMatrix, CoarseVector> coarseMdagMOp(_CoarseMatrix);

    std::cout << GridLogMG << " Level " << level << ": **************************************************" << std::endl;
    std::cout << GridLogMG << " Level " << level << ": MG correctness check: 0 == (1 - P R) v" << std::endl;
    std::cout << GridLogMG << " Level " << level << ": **************************************************" << std::endl;

    for(auto i = 0; i < _Aggregates.subspace.size(); ++i) {
      _Aggregates.ProjectToSubspace(coarseTmps[0], _Aggregates.subspace[i]); //   R v_i
      _Aggregates.PromoteFromSubspace(coarseTmps[0], fineTmps[0]);           // P R v_i

      fineTmps[1]    = _Aggregates.subspace[i] - fineTmps[0]; // v_i - P R v_i
      auto deviation = std::sqrt(norm2(fineTmps[1]) / norm2(_Aggregates.subspace[i]));

      std::cout << GridLogMG << " Level " << level << ": Vector " << i << ": norm2(v_i) = " << norm2(_Aggregates.subspace[i])
                << " | norm2(R v_i) = " << norm2(coarseTmps[0]) << " | norm2(P R v_i) = " << norm2(fineTmps[0])
                << " | relative deviation = " << deviation;

      if(deviation > tolerance) {
        std::cout << " > " << tolerance << " -> check failed" << std::endl;
        // abort();
      } else {
        std::cout << " < " << tolerance << " -> check passed" << std::endl;
      }
    }

    std::cout << GridLogMG << " Level " << level << ": **************************************************" << std::endl;
    std::cout << GridLogMG << " Level " << level << ": MG correctness check: 0 == (1 - R P) v_c" << std::endl;
    std::cout << GridLogMG << " Level " << level << ": **************************************************" << std::endl;

    random(_LevelInfo.PRNGs[coarseLevel], coarseTmps[0]);

    _Aggregates.PromoteFromSubspace(coarseTmps[0], fineTmps[0]); //   P v_c
    _Aggregates.ProjectToSubspace(coarseTmps[1], fineTmps[0]);   // R P v_c

    coarseTmps[2]  = coarseTmps[0] - coarseTmps[1]; // v_c - R P v_c
    auto deviation = std::sqrt(norm2(coarseTmps[2]) / norm2(coarseTmps[0]));

    std::cout << GridLogMG << " Level " << level << ": norm2(v_c) = " << norm2(coarseTmps[0])
              << " | norm2(R P v_c) = " << norm2(coarseTmps[1]) << " | norm2(P v_c) = " << norm2(fineTmps[0])
              << " | relative deviation = " << deviation;

    if(deviation > tolerance) {
      std::cout << " > " << tolerance << " -> check failed" << std::endl;
      // abort();
    } else {
      std::cout << " < " << tolerance << " -> check passed" << std::endl;
    }

    std::cout << GridLogMG << " Level " << level << ": **************************************************" << std::endl;
    std::cout << GridLogMG << " Level " << level << ": MG correctness check: 0 == (R D P - D_c) v_c" << std::endl;
    std::cout << GridLogMG << " Level " << level << ": **************************************************" << std::endl;

    random(_LevelInfo.PRNGs[coarseLevel], coarseTmps[0]);

    _Aggregates.PromoteFromSubspace(coarseTmps[0], fineTmps[0]); //     P v_c
    fineMdagMOp.Op(fineTmps[0], fineTmps[1]);                    //   D P v_c
    _Aggregates.ProjectToSubspace(coarseTmps[1], fineTmps[1]);   // R D P v_c

    coarseMdagMOp.Op(coarseTmps[0], coarseTmps[2]); // D_c v_c

    coarseTmps[3] = coarseTmps[1] - coarseTmps[2]; // R D P v_c - D_c v_c
    deviation     = std::sqrt(norm2(coarseTmps[3]) / norm2(coarseTmps[1]));

    std::cout << GridLogMG << " Level " << level << ": norm2(R D P v_c) = " << norm2(coarseTmps[1])
              << " | norm2(D_c v_c) = " << norm2(coarseTmps[2]) << " | relative deviation = " << deviation;

    if(deviation > tolerance) {
      std::cout << " > " << tolerance << " -> check failed" << std::endl;
      // abort();
    } else {
      std::cout << " < " << tolerance << " -> check passed" << std::endl;
    }

    std::cout << GridLogMG << " Level " << level << ": **************************************************" << std::endl;
    std::cout << GridLogMG << " Level " << level << ": MG correctness check: 0 == |(Im(v_c^dag D_c^dag D_c v_c)|" << std::endl;
    std::cout << GridLogMG << " Level " << level << ": **************************************************" << std::endl;

    random(_LevelInfo.PRNGs[coarseLevel], coarseTmps[0]);

    coarseMdagMOp.Op(coarseTmps[0], coarseTmps[1]);    //         D_c v_c
    coarseMdagMOp.AdjOp(coarseTmps[1], coarseTmps[2]); // D_c^dag D_c v_c

    auto dot  = innerProduct(coarseTmps[0], coarseTmps[2]); //v_c^dag D_c^dag D_c v_c
    deviation = abs(imag(dot)) / abs(real(dot));

    std::cout << GridLogMG << " Level " << level << ": Re(v_c^dag D_c^dag D_c v_c) = " << real(dot)
              << " | Im(v_c^dag D_c^dag D_c v_c) = " << imag(dot) << " | relative deviation = " << deviation;

    if(deviation > tolerance) {
      std::cout << " > " << tolerance << " -> check failed" << std::endl;
      // abort();
    } else {
      std::cout << " < " << tolerance << " -> check passed"
                << std::endl; // TODO: this check will work only when I got Mdag in CoarsenedMatrix to work
    }

    _NextPreconditionerLevel->runChecks();
  }
};

// Specialize the coarsest level, this corresponds to counting downwards with level: coarsest = 0, finest = N
template<class Fobj, class CoarseScalar, int nCoarseSpins, int nbasis, class Matrix>
class MultiGridPreconditioner<Fobj, CoarseScalar, nCoarseSpins, nbasis, 0, Matrix> : public LinearFunction<Lattice<Fobj>> {
public:
  /////////////////////////////////////////////
  // Type Definitions
  /////////////////////////////////////////////

  typedef Matrix        FineMatrix;
  typedef Lattice<Fobj> FineVector;

  /////////////////////////////////////////////
  // Member Data
  /////////////////////////////////////////////

  MultiGridParams &_MultiGridParams;
  LevelInfo &      _LevelInfo;
  FineMatrix &     _FineMatrix;
  FineMatrix &     _SmootherMatrix;

  /////////////////////////////////////////////
  // Member Functions
  /////////////////////////////////////////////

  MultiGridPreconditioner(MultiGridParams &mgParams, LevelInfo &LvlInfo, FineMatrix &FineMat, FineMatrix &SmootherMat)
    : _MultiGridParams(mgParams), _LevelInfo(LvlInfo), _FineMatrix(FineMat), _SmootherMatrix(SmootherMat) {}

  void setup() {}

  virtual void operator()(Lattice<Fobj> const &in, Lattice<Fobj> &out) {

    TrivialPrecon<FineVector>                      fineTrivialPreconditioner;
    FlexibleGeneralisedMinimalResidual<FineVector> fineFGMRES(1.0e-14, 1, fineTrivialPreconditioner, 1, false);

    MdagMLinearOperator<FineMatrix, FineVector> fineMdagMOp(_FineMatrix);

    fineFGMRES(fineMdagMOp, in, out);
  }

  void runChecks() {}
};

template<class Fobj, class CoarseScalar, int nCoarseSpins, int nbasis, class Matrix>
using FourLevelMGPreconditioner = MultiGridPreconditioner<Fobj, CoarseScalar, nCoarseSpins, nbasis, 4 - 1, Matrix>;

template<class Fobj, class CoarseScalar, int nCoarseSpins, int nbasis, class Matrix>
using ThreeLevelMGPreconditioner = MultiGridPreconditioner<Fobj, CoarseScalar, nCoarseSpins, nbasis, 3 - 1, Matrix>;

template<class Fobj, class CoarseScalar, int nCoarseSpins, int nbasis, class Matrix>
using TwoLevelMGPreconditioner = MultiGridPreconditioner<Fobj, CoarseScalar, nCoarseSpins, nbasis, 2 - 1, Matrix>;

template<class Fobj, class CoarseScalar, int nCoarseSpins, int nbasis, int nlevel, class Matrix>
using NLevelMGPreconditioner = MultiGridPreconditioner<Fobj, CoarseScalar, nCoarseSpins, nbasis, nlevel - 1, Matrix>;

int main(int argc, char **argv) {

  Grid_init(&argc, &argv);

  typename WilsonCloverFermionR::ImplParams wcImplparams;
  WilsonAnisotropyCoefficients              wilsonAnisCoeff;

  GridCartesian *FGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *FrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(FGrid);

  std::vector<int> fSeeds({1, 2, 3, 4});
  GridParallelRNG  fPRNG(FGrid);
  fPRNG.SeedFixedIntegers(fSeeds);

  Gamma g5(Gamma::Algebra::Gamma5);

  // clang-format off
  LatticeFermion    src(FGrid); gaussian(fPRNG, src);
  LatticeFermion result(FGrid); result = zero;
  LatticeGaugeField Umu(FGrid); SU3::HotConfiguration(fPRNG, Umu);
  // clang-format on

  RealD mass  = 0.5;
  RealD csw_r = 1.0;
  RealD csw_t = 1.0;

  const int nbasis = 20;

  WilsonFermionR       Dw(Umu, *FGrid, *FrbGrid, mass);
  WilsonCloverFermionR Dwc(Umu, *FGrid, *FrbGrid, mass, csw_r, csw_t, wilsonAnisCoeff, wcImplparams);

  // mgParams.blockSizes = {{2, 2, 2, 2}, {2, 2, 1, 1}, {1, 1, 2, 1}};
  // mgParams.blockSizes = {{2, 2, 2, 2}, {2, 2, 1, 1}};
  mgParams.blockSizes = {{2, 2, 2, 2}};
  mgParams.nLevels    = mgParams.blockSizes.size() + 1;
  mgParams.kCycle     = true;

  std::cout << mgParams << std::endl;

  LevelInfo levelInfo(FGrid, mgParams);

  static_assert(std::is_same<LatticeFermion, typename WilsonFermionR::FermionField>::value, "");
  static_assert(std::is_same<LatticeFermion, typename WilsonCloverFermionR::FermionField>::value, "");

  MdagMLinearOperator<WilsonFermionR, LatticeFermion>       MdagMOpDw(Dw);
  MdagMLinearOperator<WilsonCloverFermionR, LatticeFermion> MdagMOpDwc(Dwc);

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Testing Multigrid for Wilson" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  TrivialPrecon<LatticeFermion>                                                     TrivialPrecon;
  TwoLevelMGPreconditioner<vSpinColourVector, vTComplex, 1, nbasis, WilsonFermionR> TwoLevelMGPreconDw(mgParams, levelInfo, Dw, Dw);
  // ThreeLevelMGPreconditioner<vSpinColourVector, vTComplex, 1, nbasis, WilsonFermionR> ThreeLevelMGPreconDw(levelInfo, Dw, Dw);
  // FourLevelMGPreconditioner<vSpinColourVector, vTComplex, 1, nbasis, WilsonFermionR> FourLevelMGPreconDw(levelInfo, Dw, Dw);
  // NLevelMGPreconditioner<vSpinColourVector, vTComplex, 1, nbasis, 4, WilsonFermionR> NLevelMGPreconDw(levelInfo, Dw, Dw);

  TwoLevelMGPreconDw.setup();
  TwoLevelMGPreconDw.runChecks();

  // ThreeLevelMGPreconDw.setup();
  // ThreeLevelMGPreconDw.runChecks();

  // FourLevelMGPreconDw.setup();
  // FourLevelMGPreconDw.runChecks();

  // NLevelMGPreconDw.setup();
  // NLevelMGPreconDw.runChecks();

  std::vector<std::unique_ptr<OperatorFunction<LatticeFermion>>> solversDw;

  solversDw.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 50000, TrivialPrecon, 1000, false));
  solversDw.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 50000, TwoLevelMGPreconDw, 1000, false));
  // solversDw.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 50000, ThreeLevelMGPreconDw, 1000, false));
  // solversDw.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 50000, FourLevelMGPreconDw, 1000, false));
  // solversDw.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 50000, NLevelMGPreconDw, 1000, false));

  for(auto const &solver : solversDw) {
    std::cout << "Starting with a new solver" << std::endl;
    result = zero;
    (*solver)(MdagMOpDw, src, result);
    std::cout << std::endl;
  }

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Testing Multigrid for Wilson Clover" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  TwoLevelMGPreconditioner<vSpinColourVector, vTComplex, 1, nbasis, WilsonCloverFermionR> TwoLevelMGPreconDwc(
    mgParams, levelInfo, Dwc, Dwc);
  // ThreeLevelMGPreconditioner<vSpinColourVector, vTComplex, 1, nbasis, WilsonCloverFermionR> ThreeLevelMGPreconDwc(levelInfo, Dwc, Dwc);
  // FourLevelMGPreconditioner<vSpinColourVector, vTComplex, 1, nbasis, WilsonCloverFermionR> FourLevelMGPreconDwc(levelInfo, Dwc, Dwc);
  // NLevelMGPreconditioner<vSpinColourVector, vTComplex, 1, nbasis, 4, WilsonCloverFermionR> NLevelMGPreconDwc(levelInfo, Dwc, Dwc);

  TwoLevelMGPreconDwc.setup();
  TwoLevelMGPreconDwc.runChecks();

  // ThreeLevelMGPreconDwc.setup();
  // ThreeLevelMGPreconDwc.runChecks();

  // FourLevelMGPreconDwc.setup();
  // FourLevelMGPreconDwc.runChecks();

  // NLevelMGPreconDwc.setup();
  // NLevelMGPreconDwc.runChecks();

  std::vector<std::unique_ptr<OperatorFunction<LatticeFermion>>> solversDwc;

  solversDwc.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 50000, TrivialPrecon, 1000, false));
  solversDwc.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 50000, TwoLevelMGPreconDwc, 1000, false));
  // solversDwc.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 50000, ThreeLevelMGPreconDwc, 1000, false));
  // solversDwc.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 50000, FourLevelMGPreconDwc, 1000, false));
  // solversDwc.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 50000, NLevelMGPreconDwc, 1000, false));

  for(auto const &solver : solversDwc) {
    std::cout << "Starting with a new solver" << std::endl;
    result = zero;
    (*solver)(MdagMOpDwc, src, result);
    std::cout << std::endl;
  }

  Grid_finalize();
}
