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

// TODO: Can think about having one parameter struct per level and then a
// vector of these structs. How well would that work together with the
// serialization strategy of Grid?

// clang-format off
struct MultiGridParams : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(MultiGridParams,
                                  int,                           nLevels,
                                  std::vector<std::vector<int>>, blockSizes,           // size == nLevels - 1
                                  std::vector<double>,           smootherTol,          // size == nLevels - 1
                                  std::vector<int>,              smootherMaxOuterIter, // size == nLevels - 1
                                  std::vector<int>,              smootherMaxInnerIter, // size == nLevels - 1
                                  bool,                          kCycle,
                                  std::vector<double>,           kCycleTol,            // size == nLevels - 1
                                  std::vector<int>,              kCycleMaxOuterIter,   // size == nLevels - 1
                                  std::vector<int>,              kCycleMaxInnerIter,   // size == nLevels - 1
                                  double,                        coarseSolverTol,
                                  int,                           coarseSolverMaxOuterIter,
                                  int,                           coarseSolverMaxInnerIter);
  MultiGridParams(){};
};
// clang-format on

void checkParameterValidity(MultiGridParams const &params) {

  auto correctSize = params.nLevels - 1;

  assert(correctSize == params.blockSizes.size());
  assert(correctSize == params.smootherTol.size());
  assert(correctSize == params.smootherMaxOuterIter.size());
  assert(correctSize == params.smootherMaxInnerIter.size());
  assert(correctSize == params.kCycleTol.size());
  assert(correctSize == params.kCycleMaxOuterIter.size());
  assert(correctSize == params.kCycleMaxInnerIter.size());
}

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

    for(int level = 0; level < mgParams.nLevels; ++level) {
      std::cout << GridLogMessage << "level = " << level << ":" << std::endl;
      Grids[level]->show_decomposition();
    }
  }
};

template<class Field> class MultiGridPreconditionerBase : public LinearFunction<Field> {
public:
  virtual ~MultiGridPreconditionerBase()               = default;
  virtual void setup()                                 = 0;
  virtual void operator()(Field const &in, Field &out) = 0;
  virtual void runChecks()                             = 0;
};

template<class Fobj, class CoarseScalar, int nCoarseSpins, int nBasis, int nCoarserLevels, class Matrix>
class MultiGridPreconditioner : public MultiGridPreconditionerBase<Lattice<Fobj>> {
public:
  /////////////////////////////////////////////
  // Type Definitions
  /////////////////////////////////////////////

  // clang-format off
  typedef Aggregation<Fobj, CoarseScalar, nBasis>                                                                         Aggregates;
  typedef CoarsenedMatrix<Fobj, CoarseScalar, nBasis>                                                                     CoarseMatrix;
  typedef typename Aggregates::CoarseVector                                                                               CoarseVector;
  typedef typename Aggregates::siteVector                                                                                 CoarseSiteVector;
  typedef Matrix                                                                                                          FineMatrix;
  typedef typename Aggregates::FineField                                                                                  FineVector;
  typedef MultiGridPreconditioner<CoarseSiteVector, CoarseScalar, nCoarseSpins, nBasis, nCoarserLevels - 1, CoarseMatrix> NextPreconditionerLevel;
  // clang-format on

  /////////////////////////////////////////////
  // Member Data
  /////////////////////////////////////////////

  int                                      _CurrentLevel;
  int                                      _NextCoarserLevel;
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
    : _CurrentLevel(mgParams.nLevels - (nCoarserLevels + 1)) // _Level = 0 corresponds to finest
    , _NextCoarserLevel(_CurrentLevel + 1)                   // incremented for instances on coarser levels
    , _MultiGridParams(mgParams)
    , _LevelInfo(LvlInfo)
    , _FineMatrix(FineMat)
    , _SmootherMatrix(SmootherMat)
    , _Aggregates(_LevelInfo.Grids[_NextCoarserLevel], _LevelInfo.Grids[_CurrentLevel], 0)
    , _CoarseMatrix(*_LevelInfo.Grids[_NextCoarserLevel]) {
    _NextPreconditionerLevel
      = std::unique_ptr<NextPreconditionerLevel>(new NextPreconditionerLevel(_MultiGridParams, _LevelInfo, _CoarseMatrix, _CoarseMatrix));
  }

  void setup() {

    Gamma                                       g5(Gamma::Algebra::Gamma5);
    MdagMLinearOperator<FineMatrix, FineVector> fineMdagMOp(_FineMatrix);

    // NOTE: Don't specify nb here to see the orthogonalization check
    _Aggregates.CreateSubspace(_LevelInfo.PRNGs[_CurrentLevel], fineMdagMOp /*, nb */);

    // TestVectorAnalyzer<FineVector, nbasis> fineTVA;
    // fineTVA(fineMdagMOp, _Aggregates.subspace);

    static_assert((nBasis & 0x1) == 0, "MG Preconditioner only supports an even number of basis vectors");
    int nb = nBasis / 2;

    // // TODO: to get this to work for more than two levels, I would need to either implement coarse spins or have a template specialization of this class also for the finest level
    // for(int n = 0; n < nb; n++) {
    //   _Aggregates.subspace[n + nb] = g5 * _Aggregates.subspace[n];
    // }

    _CoarseMatrix.CoarsenOperator(_LevelInfo.Grids[_CurrentLevel], fineMdagMOp, _Aggregates);

    _NextPreconditionerLevel->setup();
  }

  virtual void operator()(Lattice<Fobj> const &in, Lattice<Fobj> &out) {

    conformable(_LevelInfo.Grids[_CurrentLevel], in._grid);
    conformable(in, out);

    // TODO: implement a W-cycle
    if(_MultiGridParams.kCycle)
      kCycle(in, out);
    else
      vCycle(in, out);
  }

  void vCycle(Lattice<Fobj> const &in, Lattice<Fobj> &out) {

    RealD inputNorm = norm2(in);

    CoarseVector coarseSrc(_LevelInfo.Grids[_NextCoarserLevel]);
    CoarseVector coarseSol(_LevelInfo.Grids[_NextCoarserLevel]);
    coarseSol = zero;

    FineVector fineTmp(in._grid);

    auto maxSmootherIter = _MultiGridParams.smootherMaxOuterIter[_CurrentLevel] * _MultiGridParams.smootherMaxInnerIter[_CurrentLevel];

    TrivialPrecon<FineVector>                      fineTrivialPreconditioner;
    FlexibleGeneralisedMinimalResidual<FineVector> fineFGMRES(_MultiGridParams.smootherTol[_CurrentLevel],
                                                              maxSmootherIter,
                                                              fineTrivialPreconditioner,
                                                              _MultiGridParams.smootherMaxInnerIter[_CurrentLevel],
                                                              false);

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

    std::cout << GridLogMG << " Level " << _CurrentLevel << ": V-cycle: Input norm = " << std::sqrt(inputNorm)
              << " Coarse residual = " << residualAfterCoarseGridCorrection << " Post-Smoother residual = " << residualAfterPostSmoother
              << std::endl;
  }

  void kCycle(Lattice<Fobj> const &in, Lattice<Fobj> &out) {

    RealD inputNorm = norm2(in);

    CoarseVector coarseSrc(_LevelInfo.Grids[_NextCoarserLevel]);
    CoarseVector coarseSol(_LevelInfo.Grids[_NextCoarserLevel]);
    coarseSol = zero;

    FineVector fineTmp(in._grid);

    auto smootherMaxIter = _MultiGridParams.smootherMaxOuterIter[_CurrentLevel] * _MultiGridParams.smootherMaxInnerIter[_CurrentLevel];
    auto kCycleMaxIter   = _MultiGridParams.kCycleMaxOuterIter[_CurrentLevel] * _MultiGridParams.kCycleMaxInnerIter[_CurrentLevel];

    TrivialPrecon<FineVector>                        fineTrivialPreconditioner;
    FlexibleGeneralisedMinimalResidual<FineVector>   fineFGMRES(_MultiGridParams.smootherTol[_CurrentLevel],
                                                              smootherMaxIter,
                                                              fineTrivialPreconditioner,
                                                              _MultiGridParams.smootherMaxInnerIter[_CurrentLevel],
                                                              false);
    FlexibleGeneralisedMinimalResidual<CoarseVector> coarseFGMRES(_MultiGridParams.kCycleTol[_CurrentLevel],
                                                                  kCycleMaxIter,
                                                                  *_NextPreconditionerLevel,
                                                                  _MultiGridParams.kCycleMaxInnerIter[_CurrentLevel],
                                                                  false);

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

    std::cout << GridLogMG << " Level " << _CurrentLevel << ": K-cycle: Input norm = " << std::sqrt(inputNorm)
              << " Coarse residual = " << residualAfterCoarseGridCorrection << " Post-Smoother residual = " << residualAfterPostSmoother
              << std::endl;
  }

  void runChecks() {

    auto tolerance = 1e-13; // TODO: this obviously depends on the precision we use, current value is for double

    std::vector<FineVector>   fineTmps(7, _LevelInfo.Grids[_CurrentLevel]);
    std::vector<CoarseVector> coarseTmps(4, _LevelInfo.Grids[_NextCoarserLevel]);

    MdagMLinearOperator<FineMatrix, FineVector>     fineMdagMOp(_FineMatrix);
    MdagMLinearOperator<CoarseMatrix, CoarseVector> coarseMdagMOp(_CoarseMatrix);

    std::cout << GridLogMG << " Level " << _CurrentLevel << ": **************************************************" << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": MG correctness check: 0 == (M - (Mdiag + Σ_μ Mdir_μ)) * v" << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": **************************************************" << std::endl;

    random(_LevelInfo.PRNGs[_CurrentLevel], fineTmps[0]);

    fineMdagMOp.Op(fineTmps[0], fineTmps[1]);     //     M * v
    fineMdagMOp.OpDiag(fineTmps[0], fineTmps[2]); // Mdiag * v

    fineTmps[4] = zero;
    for(int dir = 0; dir < 4; dir++) { //       Σ_μ Mdir_μ * v
      for(auto disp : {+1, -1}) {
        fineMdagMOp.OpDir(fineTmps[0], fineTmps[3], dir, disp);
        fineTmps[4] = fineTmps[4] + fineTmps[3];
      }
    }

    fineTmps[5] = fineTmps[2] + fineTmps[4]; // (Mdiag + Σ_μ Mdir_μ) * v

    fineTmps[6]    = fineTmps[1] - fineTmps[5];
    auto deviation = std::sqrt(norm2(fineTmps[6]) / norm2(fineTmps[1]));

    std::cout << GridLogMG << " Level " << _CurrentLevel << ": norm2(M * v)                    = " << norm2(fineTmps[1]) << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": norm2(Mdiag * v)                = " << norm2(fineTmps[2]) << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": norm2(Σ_μ Mdir_μ * v)           = " << norm2(fineTmps[4]) << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": norm2((Mdiag + Σ_μ Mdir_μ) * v) = " << norm2(fineTmps[5]) << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": relative deviation              = " << deviation;

    if(deviation > tolerance) {
      std::cout << " > " << tolerance << " -> check failed" << std::endl;
      // abort();
    } else {
      std::cout << " < " << tolerance << " -> check passed" << std::endl;
    }

    std::cout << GridLogMG << " Level " << _CurrentLevel << ": **************************************************" << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": MG correctness check: 0 == (1 - P R) v" << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": **************************************************" << std::endl;

    for(auto i = 0; i < _Aggregates.subspace.size(); ++i) {
      _Aggregates.ProjectToSubspace(coarseTmps[0], _Aggregates.subspace[i]); //   R v_i
      _Aggregates.PromoteFromSubspace(coarseTmps[0], fineTmps[0]);           // P R v_i

      fineTmps[1] = _Aggregates.subspace[i] - fineTmps[0]; // v_i - P R v_i
      deviation   = std::sqrt(norm2(fineTmps[1]) / norm2(_Aggregates.subspace[i]));

      std::cout << GridLogMG << " Level " << _CurrentLevel << ": Vector " << i << ": norm2(v_i) = " << norm2(_Aggregates.subspace[i])
                << " | norm2(R v_i) = " << norm2(coarseTmps[0]) << " | norm2(P R v_i) = " << norm2(fineTmps[0])
                << " | relative deviation = " << deviation;

      if(deviation > tolerance) {
        std::cout << " > " << tolerance << " -> check failed" << std::endl;
        // abort();
      } else {
        std::cout << " < " << tolerance << " -> check passed" << std::endl;
      }
    }

    std::cout << GridLogMG << " Level " << _CurrentLevel << ": **************************************************" << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": MG correctness check: 0 == (1 - R P) v_c" << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": **************************************************" << std::endl;

    random(_LevelInfo.PRNGs[_NextCoarserLevel], coarseTmps[0]);

    _Aggregates.PromoteFromSubspace(coarseTmps[0], fineTmps[0]); //   P v_c
    _Aggregates.ProjectToSubspace(coarseTmps[1], fineTmps[0]);   // R P v_c

    coarseTmps[2] = coarseTmps[0] - coarseTmps[1]; // v_c - R P v_c
    deviation     = std::sqrt(norm2(coarseTmps[2]) / norm2(coarseTmps[0]));

    std::cout << GridLogMG << " Level " << _CurrentLevel << ": norm2(v_c) = " << norm2(coarseTmps[0])
              << " | norm2(R P v_c) = " << norm2(coarseTmps[1]) << " | norm2(P v_c) = " << norm2(fineTmps[0])
              << " | relative deviation = " << deviation;

    if(deviation > tolerance) {
      std::cout << " > " << tolerance << " -> check failed" << std::endl;
      // abort();
    } else {
      std::cout << " < " << tolerance << " -> check passed" << std::endl;
    }

    std::cout << GridLogMG << " Level " << _CurrentLevel << ": **************************************************" << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": MG correctness check: 0 == (R D P - D_c) v_c" << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": **************************************************" << std::endl;

    random(_LevelInfo.PRNGs[_NextCoarserLevel], coarseTmps[0]);

    _Aggregates.PromoteFromSubspace(coarseTmps[0], fineTmps[0]); //     P v_c
    fineMdagMOp.Op(fineTmps[0], fineTmps[1]);                    //   D P v_c
    _Aggregates.ProjectToSubspace(coarseTmps[1], fineTmps[1]);   // R D P v_c

    coarseMdagMOp.Op(coarseTmps[0], coarseTmps[2]); // D_c v_c

    coarseTmps[3] = coarseTmps[1] - coarseTmps[2]; // R D P v_c - D_c v_c
    deviation     = std::sqrt(norm2(coarseTmps[3]) / norm2(coarseTmps[1]));

    std::cout << GridLogMG << " Level " << _CurrentLevel << ": norm2(R D P v_c) = " << norm2(coarseTmps[1])
              << " | norm2(D_c v_c) = " << norm2(coarseTmps[2]) << " | relative deviation = " << deviation;

    if(deviation > tolerance) {
      std::cout << " > " << tolerance << " -> check failed" << std::endl;
      // abort();
    } else {
      std::cout << " < " << tolerance << " -> check passed" << std::endl;
    }

    std::cout << GridLogMG << " Level " << _CurrentLevel << ": **************************************************" << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": MG correctness check: 0 == |(Im(v_c^dag D_c^dag D_c v_c)|" << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": **************************************************" << std::endl;

    random(_LevelInfo.PRNGs[_NextCoarserLevel], coarseTmps[0]);

    coarseMdagMOp.Op(coarseTmps[0], coarseTmps[1]);    //         D_c v_c
    coarseMdagMOp.AdjOp(coarseTmps[1], coarseTmps[2]); // D_c^dag D_c v_c

    auto dot  = innerProduct(coarseTmps[0], coarseTmps[2]); //v_c^dag D_c^dag D_c v_c
    deviation = abs(imag(dot)) / abs(real(dot));

    std::cout << GridLogMG << " Level " << _CurrentLevel << ": Re(v_c^dag D_c^dag D_c v_c) = " << real(dot)
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

// Specialization for the coarsest level
template<class Fobj, class CoarseScalar, int nCoarseSpins, int nbasis, class Matrix>
class MultiGridPreconditioner<Fobj, CoarseScalar, nCoarseSpins, nbasis, 0, Matrix> : public MultiGridPreconditionerBase<Lattice<Fobj>> {
public:
  /////////////////////////////////////////////
  // Type Definitions
  /////////////////////////////////////////////

  typedef Matrix        FineMatrix;
  typedef Lattice<Fobj> FineVector;

  /////////////////////////////////////////////
  // Member Data
  /////////////////////////////////////////////

  int              _CurrentLevel;
  MultiGridParams &_MultiGridParams;
  LevelInfo &      _LevelInfo;
  FineMatrix &     _FineMatrix;
  FineMatrix &     _SmootherMatrix;

  /////////////////////////////////////////////
  // Member Functions
  /////////////////////////////////////////////

  MultiGridPreconditioner(MultiGridParams &mgParams, LevelInfo &LvlInfo, FineMatrix &FineMat, FineMatrix &SmootherMat)
    : _CurrentLevel(mgParams.nLevels - (0 + 1))
    , _MultiGridParams(mgParams)
    , _LevelInfo(LvlInfo)
    , _FineMatrix(FineMat)
    , _SmootherMatrix(SmootherMat) {}

  void setup() {}

  virtual void operator()(Lattice<Fobj> const &in, Lattice<Fobj> &out) {

    conformable(_LevelInfo.Grids[_CurrentLevel], in._grid);
    conformable(in, out);

    auto coarseSolverMaxIter = _MultiGridParams.coarseSolverMaxOuterIter * _MultiGridParams.coarseSolverMaxInnerIter;

    // On the coarsest level we only have a fine what I above call the fine level, no coarse one
    TrivialPrecon<FineVector>                      fineTrivialPreconditioner;
    FlexibleGeneralisedMinimalResidual<FineVector> fineFGMRES(
      _MultiGridParams.coarseSolverTol, coarseSolverMaxIter, fineTrivialPreconditioner, _MultiGridParams.coarseSolverMaxInnerIter, false);

    MdagMLinearOperator<FineMatrix, FineVector> fineMdagMOp(_FineMatrix);

    fineFGMRES(fineMdagMOp, in, out);
  }

  void runChecks() {}
};

template<class Fobj, class CoarseScalar, int nCoarseSpins, int nbasis, int nLevels, class Matrix>
using NLevelMGPreconditioner = MultiGridPreconditioner<Fobj, CoarseScalar, nCoarseSpins, nbasis, nLevels - 1, Matrix>;

template<class Fobj, class CoarseScalar, int nCoarseSpins, int nbasis, class Matrix>
std::unique_ptr<MultiGridPreconditionerBase<Lattice<Fobj>>>
createMGInstance(MultiGridParams &mgParams, LevelInfo &levelInfo, Matrix &FineMat, Matrix &SmootherMat) {

  // clang-format off
  #define CASE_FOR_N_LEVELS(nLevels)                                                                                                       \
    case nLevels:                                                                                                                          \
      return std::unique_ptr<NLevelMGPreconditioner<Fobj, CoarseScalar, nCoarseSpins, nbasis, nLevels, Matrix>>(                           \
        new NLevelMGPreconditioner<Fobj, CoarseScalar, nCoarseSpins, nbasis, nLevels, Matrix>(mgParams, levelInfo, FineMat, SmootherMat)); \
      break;
  // clang-format on

  switch(mgParams.nLevels) {
    CASE_FOR_N_LEVELS(2);
    CASE_FOR_N_LEVELS(3);
    CASE_FOR_N_LEVELS(4);
    default:
      std::cout << GridLogError << "We currently only support nLevels ∈ {2, 3, 4}" << std::endl;
      exit(EXIT_FAILURE);
      break;
  }
}

int main(int argc, char **argv) {

  Grid_init(&argc, &argv);

  MultiGridParams mgParams;

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

  // Params for two-level MG preconditioner
  mgParams.nLevels                  = 2;
  mgParams.blockSizes               = {{2, 2, 2, 2}};
  mgParams.smootherTol              = {1e-14};
  mgParams.smootherMaxOuterIter     = {1};
  mgParams.smootherMaxInnerIter     = {1};
  mgParams.kCycle                   = true;
  mgParams.kCycleTol                = {1e-14};
  mgParams.kCycleMaxOuterIter       = {1};
  mgParams.kCycleMaxInnerIter       = {1};
  mgParams.coarseSolverTol          = 1e-14;
  mgParams.coarseSolverMaxOuterIter = 1;
  mgParams.coarseSolverMaxInnerIter = 1;

  // // Params for three-level MG preconditioner
  // mgParams.nLevels                  = 3;
  // mgParams.blockSizes               = {{2, 2, 2, 2}, {2, 2, 1, 1}};
  // mgParams.smootherTol              = {1e-14, 1e-14};
  // mgParams.smootherMaxOuterIter     = {1, 1};
  // mgParams.smootherMaxInnerIter     = {1, 1};
  // mgParams.kCycle                   = true;
  // mgParams.kCycleTol                = {1e-14, 1e-14};
  // mgParams.kCycleMaxOuterIter       = {1, 1};
  // mgParams.kCycleMaxInnerIter       = {1, 1};
  // mgParams.coarseSolverTol          = 1e-14;
  // mgParams.coarseSolverMaxOuterIter = 1;
  // mgParams.coarseSolverMaxInnerIter = 1;

  // // // Params for four-level MG preconditioner
  // mgParams.nLevels                  = 4;
  // mgParams.blockSizes               = {{2, 2, 2, 2}, {2, 2, 1, 1}, {1, 1, 2, 1}};
  // mgParams.smootherTol              = {1e-14, 1e-14, 1e-14};
  // mgParams.smootherMaxOuterIter     = {1, 1, 1};
  // mgParams.smootherMaxInnerIter     = {1, 1, 1};
  // mgParams.kCycle                   = true;
  // mgParams.kCycleTol                = {1e-14, 1e-14, 1e-14};
  // mgParams.kCycleMaxOuterIter       = {1, 1, 1};
  // mgParams.kCycleMaxInnerIter       = {1, 1, 1};
  // mgParams.coarseSolverTol          = 1e-14;
  // mgParams.coarseSolverMaxOuterIter = 1;
  // mgParams.coarseSolverMaxInnerIter = 1;

  checkParameterValidity(mgParams);

  std::cout << mgParams << std::endl;

  LevelInfo levelInfo(FGrid, mgParams);

  static_assert(std::is_same<LatticeFermion, typename WilsonFermionR::FermionField>::value, "");
  static_assert(std::is_same<LatticeFermion, typename WilsonCloverFermionR::FermionField>::value, "");

  MdagMLinearOperator<WilsonFermionR, LatticeFermion>       MdagMOpDw(Dw);
  MdagMLinearOperator<WilsonCloverFermionR, LatticeFermion> MdagMOpDwc(Dwc);

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Testing Multigrid for Wilson" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  TrivialPrecon<LatticeFermion> TrivialPrecon;
  auto MGPreconDw = createMGInstance<vSpinColourVector, vTComplex, 1, nbasis, WilsonFermionR>(mgParams, levelInfo, Dw, Dw);

  MGPreconDw->setup();
  MGPreconDw->runChecks();

  std::vector<std::unique_ptr<OperatorFunction<LatticeFermion>>> solversDw;

  solversDw.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 50000, TrivialPrecon, 100, false));
  solversDw.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 50000, *MGPreconDw, 100, false));

  for(auto const &solver : solversDw) {
    std::cout << "Starting with a new solver" << std::endl;
    result = zero;
    (*solver)(MdagMOpDw, src, result);
    std::cout << std::endl;
  }

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Testing Multigrid for Wilson Clover" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  auto MGPreconDwc = createMGInstance<vSpinColourVector, vTComplex, 1, nbasis, WilsonCloverFermionR>(mgParams, levelInfo, Dwc, Dwc);

  MGPreconDwc->setup();
  MGPreconDwc->runChecks();

  std::vector<std::unique_ptr<OperatorFunction<LatticeFermion>>> solversDwc;

  solversDwc.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 50000, TrivialPrecon, 100, false));
  solversDwc.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 50000, *MGPreconDwc, 100, false));

  for(auto const &solver : solversDwc) {
    std::cout << std::endl << "Starting with a new solver" << std::endl;
    result = zero;
    (*solver)(MdagMOpDwc, src, result);
    std::cout << std::endl;
  }

  Grid_finalize();
}
