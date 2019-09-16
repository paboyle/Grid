/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/solver/Test_multigrid_common.h

    Copyright (C) 2015-2018

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
#ifndef GRID_TEST_MULTIGRID_COMMON_H
#define GRID_TEST_MULTIGRID_COMMON_H

namespace Grid {

Zero zero;

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

  // constructor with default values
  MultiGridParams(int                           _nLevels                  = 2,
                  std::vector<std::vector<int>> _blockSizes               = {{4, 4, 4, 4}},
                  std::vector<double>           _smootherTol              = {1e-14},
                  std::vector<int>              _smootherMaxOuterIter     = {4},
                  std::vector<int>              _smootherMaxInnerIter     = {4},
                  bool                          _kCycle                   = true,
                  std::vector<double>           _kCycleTol                = {1e-1},
                  std::vector<int>              _kCycleMaxOuterIter       = {2},
                  std::vector<int>              _kCycleMaxInnerIter       = {5},
                  double                        _coarseSolverTol          = 5e-2,
                  int                           _coarseSolverMaxOuterIter = 10,
                  int                           _coarseSolverMaxInnerIter = 500)
  : nLevels(_nLevels)
  , blockSizes(_blockSizes)
  , smootherTol(_smootherTol)
  , smootherMaxOuterIter(_smootherMaxOuterIter)
  , smootherMaxInnerIter(_smootherMaxInnerIter)
  , kCycle(_kCycle)
  , kCycleTol(_kCycleTol)
  , kCycleMaxOuterIter(_kCycleMaxOuterIter)
  , kCycleMaxInnerIter(_kCycleMaxInnerIter)
  , coarseSolverTol(_coarseSolverTol)
  , coarseSolverMaxOuterIter(_coarseSolverMaxOuterIter)
  , coarseSolverMaxInnerIter(_coarseSolverMaxInnerIter)
  {}
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

      Grids.push_back(SpaceTimeGrid::makeFourDimGrid(tmp, Grids[level - 1]->_simd_layout, GridDefaultMpi()));
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
  virtual void runChecks(RealD tolerance)              = 0;
  virtual void reportTimings()                         = 0;
};

template<class Fobj, class CComplex, int nBasis, int nCoarserLevels, class Matrix>
class MultiGridPreconditioner : public MultiGridPreconditionerBase<Lattice<Fobj>> {
public:
  /////////////////////////////////////////////
  // Type Definitions
  /////////////////////////////////////////////

  // clang-format off
  typedef Aggregation<Fobj, CComplex, nBasis>                                                                         Aggregates;
  typedef CoarsenedMatrix<Fobj, CComplex, nBasis>                                                                     CoarseDiracMatrix;
  typedef typename Aggregates::CoarseVector                                                                           CoarseVector;
  typedef typename Aggregates::siteVector                                                                             CoarseSiteVector;
  typedef Matrix                                                                                                      FineDiracMatrix;
  typedef typename Aggregates::FineField                                                                              FineVector;
  typedef MultiGridPreconditioner<CoarseSiteVector, iScalar<CComplex>, nBasis, nCoarserLevels - 1, CoarseDiracMatrix> NextPreconditionerLevel;
  // clang-format on

  /////////////////////////////////////////////
  // Member Data
  /////////////////////////////////////////////

  int _CurrentLevel;
  int _NextCoarserLevel;

  MultiGridParams &_MultiGridParams;
  LevelInfo &      _LevelInfo;

  FineDiracMatrix & _FineMatrix;
  FineDiracMatrix & _SmootherMatrix;
  Aggregates        _Aggregates;
  CoarseDiracMatrix _CoarseMatrix;

  std::unique_ptr<NextPreconditionerLevel> _NextPreconditionerLevel;

  GridStopWatch _SetupTotalTimer;
  GridStopWatch _SetupCreateSubspaceTimer;
  GridStopWatch _SetupProjectToChiralitiesTimer;
  GridStopWatch _SetupCoarsenOperatorTimer;
  GridStopWatch _SetupNextLevelTimer;
  GridStopWatch _SolveTotalTimer;
  GridStopWatch _SolveRestrictionTimer;
  GridStopWatch _SolveProlongationTimer;
  GridStopWatch _SolveSmootherTimer;
  GridStopWatch _SolveNextLevelTimer;

  /////////////////////////////////////////////
  // Member Functions
  /////////////////////////////////////////////

  MultiGridPreconditioner(MultiGridParams &mgParams, LevelInfo &LvlInfo, FineDiracMatrix &FineMat, FineDiracMatrix &SmootherMat)
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

    resetTimers();
  }

  void setup() {

    _SetupTotalTimer.Start();

    static_assert((nBasis & 0x1) == 0, "MG Preconditioner only supports an even number of basis vectors");
    int nb = nBasis / 2;

    MdagMLinearOperator<FineDiracMatrix, FineVector> fineMdagMOp(_FineMatrix);

    _SetupCreateSubspaceTimer.Start();
    _Aggregates.CreateSubspace(_LevelInfo.PRNGs[_CurrentLevel], fineMdagMOp, nb);
    _SetupCreateSubspaceTimer.Stop();

    _SetupProjectToChiralitiesTimer.Start();
    FineVector tmp1(_Aggregates.subspace[0].Grid());
    FineVector tmp2(_Aggregates.subspace[0].Grid());
    for(int n = 0; n < nb; n++) {
      auto tmp1 = _Aggregates.subspace[n];
      G5C(tmp2, _Aggregates.subspace[n]);
      axpby(_Aggregates.subspace[n], 0.5, 0.5, tmp1, tmp2);
      axpby(_Aggregates.subspace[n + nb], 0.5, -0.5, tmp1, tmp2);
      std::cout << GridLogMG << " Level " << _CurrentLevel << ": Chirally doubled vector " << n << ". "
                << "norm2(vec[" << n << "]) = " << norm2(_Aggregates.subspace[n]) << ". "
                << "norm2(vec[" << n + nb << "]) = " << norm2(_Aggregates.subspace[n + nb]) << std::endl;
    }
    _SetupProjectToChiralitiesTimer.Stop();

    _SetupCoarsenOperatorTimer.Start();
    _CoarseMatrix.CoarsenOperator(_LevelInfo.Grids[_CurrentLevel], fineMdagMOp, _Aggregates);
    _SetupCoarsenOperatorTimer.Stop();

    _SetupNextLevelTimer.Start();
    _NextPreconditionerLevel->setup();
    _SetupNextLevelTimer.Stop();

    _SetupTotalTimer.Stop();
  }

  virtual void operator()(FineVector const &in, FineVector &out) {

    conformable(_LevelInfo.Grids[_CurrentLevel], in.Grid());
    conformable(in, out);

    // TODO: implement a W-cycle
    if(_MultiGridParams.kCycle)
      kCycle(in, out);
    else
      vCycle(in, out);
  }

  void vCycle(FineVector const &in, FineVector &out) {

    _SolveTotalTimer.Start();

    RealD inputNorm = norm2(in);

    CoarseVector coarseSrc(_LevelInfo.Grids[_NextCoarserLevel]);
    CoarseVector coarseSol(_LevelInfo.Grids[_NextCoarserLevel]);
    coarseSol = zero;

    FineVector fineTmp(in.Grid());

    auto maxSmootherIter = _MultiGridParams.smootherMaxOuterIter[_CurrentLevel] * _MultiGridParams.smootherMaxInnerIter[_CurrentLevel];

    TrivialPrecon<FineVector>                      fineTrivialPreconditioner;
    FlexibleGeneralisedMinimalResidual<FineVector> fineFGMRES(_MultiGridParams.smootherTol[_CurrentLevel],
                                                              maxSmootherIter,
                                                              fineTrivialPreconditioner,
                                                              _MultiGridParams.smootherMaxInnerIter[_CurrentLevel],
                                                              false);

    MdagMLinearOperator<FineDiracMatrix, FineVector> fineMdagMOp(_FineMatrix);
    MdagMLinearOperator<FineDiracMatrix, FineVector> fineSmootherMdagMOp(_SmootherMatrix);

    _SolveRestrictionTimer.Start();
    _Aggregates.ProjectToSubspace(coarseSrc, in);
    _SolveRestrictionTimer.Stop();

    _SolveNextLevelTimer.Start();
    (*_NextPreconditionerLevel)(coarseSrc, coarseSol);
    _SolveNextLevelTimer.Stop();

    _SolveProlongationTimer.Start();
    _Aggregates.PromoteFromSubspace(coarseSol, out);
    _SolveProlongationTimer.Stop();

    fineMdagMOp.Op(out, fineTmp);
    fineTmp                                = in - fineTmp;
    auto r                                 = norm2(fineTmp);
    auto residualAfterCoarseGridCorrection = std::sqrt(r / inputNorm);

    _SolveSmootherTimer.Start();
    fineFGMRES(fineSmootherMdagMOp, in, out);
    _SolveSmootherTimer.Stop();

    fineMdagMOp.Op(out, fineTmp);
    fineTmp                        = in - fineTmp;
    r                              = norm2(fineTmp);
    auto residualAfterPostSmoother = std::sqrt(r / inputNorm);

    std::cout << GridLogMG << " Level " << _CurrentLevel << ": V-cycle: Input norm = " << std::sqrt(inputNorm)
              << " Coarse residual = " << residualAfterCoarseGridCorrection << " Post-Smoother residual = " << residualAfterPostSmoother
              << std::endl;

    _SolveTotalTimer.Stop();
  }

  void kCycle(FineVector const &in, FineVector &out) {

    _SolveTotalTimer.Start();

    RealD inputNorm = norm2(in);

    CoarseVector coarseSrc(_LevelInfo.Grids[_NextCoarserLevel]);
    CoarseVector coarseSol(_LevelInfo.Grids[_NextCoarserLevel]);
    coarseSol = zero;

    FineVector fineTmp(in.Grid());

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

    MdagMLinearOperator<FineDiracMatrix, FineVector>     fineMdagMOp(_FineMatrix);
    MdagMLinearOperator<FineDiracMatrix, FineVector>     fineSmootherMdagMOp(_SmootherMatrix);
    MdagMLinearOperator<CoarseDiracMatrix, CoarseVector> coarseMdagMOp(_CoarseMatrix);

    _SolveRestrictionTimer.Start();
    _Aggregates.ProjectToSubspace(coarseSrc, in);
    _SolveRestrictionTimer.Stop();

    _SolveNextLevelTimer.Start();
    coarseFGMRES(coarseMdagMOp, coarseSrc, coarseSol);
    _SolveNextLevelTimer.Stop();

    _SolveProlongationTimer.Start();
    _Aggregates.PromoteFromSubspace(coarseSol, out);
    _SolveProlongationTimer.Stop();

    fineMdagMOp.Op(out, fineTmp);
    fineTmp                                = in - fineTmp;
    auto r                                 = norm2(fineTmp);
    auto residualAfterCoarseGridCorrection = std::sqrt(r / inputNorm);

    _SolveSmootherTimer.Start();
    fineFGMRES(fineSmootherMdagMOp, in, out);
    _SolveSmootherTimer.Stop();

    fineMdagMOp.Op(out, fineTmp);
    fineTmp                        = in - fineTmp;
    r                              = norm2(fineTmp);
    auto residualAfterPostSmoother = std::sqrt(r / inputNorm);

    std::cout << GridLogMG << " Level " << _CurrentLevel << ": K-cycle: Input norm = " << std::sqrt(inputNorm)
              << " Coarse residual = " << residualAfterCoarseGridCorrection << " Post-Smoother residual = " << residualAfterPostSmoother
              << std::endl;

    _SolveTotalTimer.Stop();
  }

  void runChecks(RealD tolerance) {

    std::vector<FineVector>   fineTmps(7, _LevelInfo.Grids[_CurrentLevel]);
    std::vector<CoarseVector> coarseTmps(4, _LevelInfo.Grids[_NextCoarserLevel]);

    MdagMLinearOperator<FineDiracMatrix, FineVector>     fineMdagMOp(_FineMatrix);
    MdagMLinearOperator<CoarseDiracMatrix, CoarseVector> coarseMdagMOp(_CoarseMatrix);

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
      abort();
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
        abort();
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
      abort();
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
      abort();
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
    deviation = std::abs(imag(dot)) / std::abs(real(dot));

    std::cout << GridLogMG << " Level " << _CurrentLevel << ": Re(v_c^dag D_c^dag D_c v_c) = " << real(dot)
              << " | Im(v_c^dag D_c^dag D_c v_c) = " << imag(dot) << " | relative deviation = " << deviation;

    if(deviation > tolerance) {
      std::cout << " > " << tolerance << " -> check failed" << std::endl;
      abort();
    } else {
      std::cout << " < " << tolerance << " -> check passed" << std::endl;
    }

    _NextPreconditionerLevel->runChecks(tolerance);
  }

  void reportTimings() {

    // clang-format off
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": Time elapsed: Sum   total            " <<                _SetupTotalTimer.Elapsed() + _SolveTotalTimer.Elapsed() << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": Time elapsed: Setup total            " <<                _SetupTotalTimer.Elapsed() << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": Time elapsed: Setup create subspace  " <<       _SetupCreateSubspaceTimer.Elapsed() << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": Time elapsed: Setup project chiral   " << _SetupProjectToChiralitiesTimer.Elapsed() << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": Time elapsed: Setup coarsen operator " <<      _SetupCoarsenOperatorTimer.Elapsed() << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": Time elapsed: Setup next level       " <<            _SetupNextLevelTimer.Elapsed() << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": Time elapsed: Solve total            " <<                _SolveTotalTimer.Elapsed() << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": Time elapsed: Solve restriction      " <<          _SolveRestrictionTimer.Elapsed() << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": Time elapsed: Solve prolongation     " <<         _SolveProlongationTimer.Elapsed() << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": Time elapsed: Solve smoother         " <<             _SolveSmootherTimer.Elapsed() << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": Time elapsed: Solve next level       " <<            _SolveNextLevelTimer.Elapsed() << std::endl;
    // clang-format on

    _NextPreconditionerLevel->reportTimings();
  }

  void resetTimers() {

    _SetupTotalTimer.Reset();
    _SetupCreateSubspaceTimer.Reset();
    _SetupProjectToChiralitiesTimer.Reset();
    _SetupCoarsenOperatorTimer.Reset();
    _SetupNextLevelTimer.Reset();
    _SolveTotalTimer.Reset();
    _SolveRestrictionTimer.Reset();
    _SolveProlongationTimer.Reset();
    _SolveSmootherTimer.Reset();
    _SolveNextLevelTimer.Reset();

    _NextPreconditionerLevel->resetTimers();
  }
};

// Specialization for the coarsest level
template<class Fobj, class CComplex, int nBasis, class Matrix>
class MultiGridPreconditioner<Fobj, CComplex, nBasis, 0, Matrix> : public MultiGridPreconditionerBase<Lattice<Fobj>> {
public:
  /////////////////////////////////////////////
  // Type Definitions
  /////////////////////////////////////////////

  typedef Matrix        FineDiracMatrix;
  typedef Lattice<Fobj> FineVector;

  /////////////////////////////////////////////
  // Member Data
  /////////////////////////////////////////////

  int _CurrentLevel;

  MultiGridParams &_MultiGridParams;
  LevelInfo &      _LevelInfo;

  FineDiracMatrix &_FineMatrix;
  FineDiracMatrix &_SmootherMatrix;

  GridStopWatch _SolveTotalTimer;
  GridStopWatch _SolveSmootherTimer;

  /////////////////////////////////////////////
  // Member Functions
  /////////////////////////////////////////////

  MultiGridPreconditioner(MultiGridParams &mgParams, LevelInfo &LvlInfo, FineDiracMatrix &FineMat, FineDiracMatrix &SmootherMat)
    : _CurrentLevel(mgParams.nLevels - (0 + 1))
    , _MultiGridParams(mgParams)
    , _LevelInfo(LvlInfo)
    , _FineMatrix(FineMat)
    , _SmootherMatrix(SmootherMat) {

    resetTimers();
  }

  void setup() {}

  virtual void operator()(FineVector const &in, FineVector &out) {

    _SolveTotalTimer.Start();

    conformable(_LevelInfo.Grids[_CurrentLevel], in.Grid());
    conformable(in, out);

    auto coarseSolverMaxIter = _MultiGridParams.coarseSolverMaxOuterIter * _MultiGridParams.coarseSolverMaxInnerIter;

    // On the coarsest level we only have what I above call the fine level, no coarse one
    TrivialPrecon<FineVector>                      fineTrivialPreconditioner;
    FlexibleGeneralisedMinimalResidual<FineVector> fineFGMRES(
      _MultiGridParams.coarseSolverTol, coarseSolverMaxIter, fineTrivialPreconditioner, _MultiGridParams.coarseSolverMaxInnerIter, false);

    MdagMLinearOperator<FineDiracMatrix, FineVector> fineMdagMOp(_FineMatrix);

    _SolveSmootherTimer.Start();
    fineFGMRES(fineMdagMOp, in, out);
    _SolveSmootherTimer.Stop();

    _SolveTotalTimer.Stop();
  }

  void runChecks(RealD tolerance) {}

  void reportTimings() {

    // clang-format off
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": Time elapsed: Solve total            " <<    _SolveTotalTimer.Elapsed() << std::endl;
    std::cout << GridLogMG << " Level " << _CurrentLevel << ": Time elapsed: Solve smoother         " << _SolveSmootherTimer.Elapsed() << std::endl;
    // clang-format on
  }

  void resetTimers() {

    _SolveTotalTimer.Reset();
    _SolveSmootherTimer.Reset();
  }
};

template<class Fobj, class CComplex, int nBasis, int nLevels, class Matrix>
using NLevelMGPreconditioner = MultiGridPreconditioner<Fobj, CComplex, nBasis, nLevels - 1, Matrix>;

template<class Fobj, class CComplex, int nBasis, class Matrix>
std::unique_ptr<MultiGridPreconditionerBase<Lattice<Fobj>>>
createMGInstance(MultiGridParams &mgParams, LevelInfo &levelInfo, Matrix &FineMat, Matrix &SmootherMat) {

#define CASE_FOR_N_LEVELS(nLevels)                                                                                     \
  case nLevels:                                                                                                        \
    return std::unique_ptr<NLevelMGPreconditioner<Fobj, CComplex, nBasis, nLevels, Matrix>>(                           \
      new NLevelMGPreconditioner<Fobj, CComplex, nBasis, nLevels, Matrix>(mgParams, levelInfo, FineMat, SmootherMat)); \
    break;

  switch(mgParams.nLevels) {
    CASE_FOR_N_LEVELS(2);
    CASE_FOR_N_LEVELS(3);
    CASE_FOR_N_LEVELS(4);
    default:
      std::cout << GridLogError << "We currently only support nLevels ∈ {2, 3, 4}" << std::endl;
      exit(EXIT_FAILURE);
      break;
  }
#undef CASE_FOR_N_LEVELS
}

}
#endif
