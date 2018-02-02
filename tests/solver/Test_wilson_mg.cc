/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/solver/Test_wilson_mg.cc

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

RealD InverseApproximation(RealD x) {
  return 1.0 / x;
}

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

    LinOp.Op(src, ref);
    std::cout << GridLogMessage << " norm2(M * src)            = " << norm2(ref) << std::endl;

    LinOp.OpDiag(src, diag);
    std::cout << GridLogMessage << " norm2(Mdiag * src)        = " << norm2(diag) << std::endl;

    for(int dir = 0; dir < 4; dir++) {
      for(auto disp : {+1, -1}) {
        LinOp.OpDir(src, tmp, dir, disp);
        std::cout << GridLogMessage << " norm2(Mdir_{" << dir << "," << disp << "} * src) = " << norm2(tmp) << std::endl;
        sumDir = sumDir + tmp;
      }
    }
    std::cout << GridLogMessage << " norm2(Σ_μ Mdir_μ * src)   = " << norm2(sumDir) << std::endl;

    result = diag + sumDir;
    err    = ref - result;

    std::cout << GridLogMessage << " Absolute deviation        = " << norm2(err) << std::endl;
    std::cout << GridLogMessage << " Relative deviation        = " << norm2(err) / norm2(ref) << std::endl;
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

  void PowerMethod(const FineField &in) {

    FineField p1(in._grid);
    FineField p2(in._grid);

    MdagMLinearOperator<Matrix, FineField> fMdagMOp(_FineMatrix);

    p1 = in;
    RealD absp2;
    for(int i = 0; i < 20; i++) {
      RealD absp1 = std::sqrt(norm2(p1));
      fMdagMOp.HermOp(p1, p2); // this is the G5 herm bit
      // _FineOperator.Op(p1,p2); // this is the G5 herm bit
      RealD absp2 = std::sqrt(norm2(p2));
      if(i % 10 == 9)
        std::cout << GridLogMessage << "Power method on mdagm " << i << " " << absp2 / absp1 << std::endl;
      p1 = p2 * (1.0 / std::sqrt(absp2));
    }
  }

  void operator()(const FineField &in, FineField &out) {
    if(params.domaindecompose) {
      operatorSAP(in, out);
    } else {
      operatorCheby(in, out);
    }
  }

    ////////////////////////////////////////////////////////////////////////
    // ADEF2: [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]
    // ADEF1: [MP+Q ] in = M [1 - A Q] in + Q in
    ////////////////////////////////////////////////////////////////////////
#if 1
  void operatorADEF2(const FineField &in, FineField &out) {

    CoarseVector Csrc(_CoarseOperator.Grid());
    CoarseVector Ctmp(_CoarseOperator.Grid());
    CoarseVector Csol(_CoarseOperator.Grid());

    ConjugateGradient<CoarseVector> CG(1.0e-10, 100000);
    ConjugateGradient<FineField>    fCG(3.0e-2, 1000);

    HermitianLinearOperator<CoarseOperator, CoarseVector> HermOp(_CoarseOperator);
    MdagMLinearOperator<CoarseOperator, CoarseVector>     MdagMOp(_CoarseOperator);
    MdagMLinearOperator<Matrix, FineField>                fMdagMOp(_FineMatrix);

    FineField tmp(in._grid);
    FineField res(in._grid);
    FineField Min(in._grid);

    // Monitor completeness of low mode space
    _Aggregates.ProjectToSubspace(Csrc, in);
    _Aggregates.PromoteFromSubspace(Csrc, out);
    std::cout << GridLogMessage << "Coarse Grid Preconditioner\nCompleteness in: " << std::sqrt(norm2(out) / norm2(in)) << std::endl;

    // [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]
    _FineOperator.Op(in, tmp); // this is the G5 herm bit
    fCG(fMdagMOp, tmp, Min);   // solves MdagM = g5 M g5M

    // Monitor completeness of low mode space
    _Aggregates.ProjectToSubspace(Csrc, Min);
    _Aggregates.PromoteFromSubspace(Csrc, out);
    std::cout << GridLogMessage << "Completeness Min: " << std::sqrt(norm2(out) / norm2(Min)) << std::endl;

    _FineOperator.Op(Min, tmp);
    tmp = in - tmp; // in - A Min

    Csol = zero;
    _Aggregates.ProjectToSubspace(Csrc, tmp);
    HermOp.AdjOp(Csrc, Ctmp); // Normal equations
    CG(MdagMOp, Ctmp, Csol);

    HermOp.Op(Csol, Ctmp);
    Ctmp = Ctmp - Csrc;
    std::cout << GridLogMessage << "coarse space true residual " << std::sqrt(norm2(Ctmp) / norm2(Csrc)) << std::endl;
    _Aggregates.PromoteFromSubspace(Csol, out);

    _FineOperator.Op(out, res);
    res = res - tmp;
    std::cout << GridLogMessage << "promoted sol residual " << std::sqrt(norm2(res) / norm2(tmp)) << std::endl;
    _Aggregates.ProjectToSubspace(Csrc, res);
    std::cout << GridLogMessage << "coarse space proj of residual " << norm2(Csrc) << std::endl;

    out = out + Min; // additive coarse space correction
    //    out = Min; // no additive coarse space correction

    _FineOperator.Op(out, tmp);
    tmp = tmp - in; // tmp is new residual

    std::cout << GridLogMessage << " Preconditioner in  " << norm2(in) << std::endl;
    std::cout << GridLogMessage << " Preconditioner out " << norm2(out) << std::endl;
    std::cout << GridLogMessage << "preconditioner thinks residual is " << std::sqrt(norm2(tmp) / norm2(in)) << std::endl;
  }
#endif
    // ADEF1: [MP+Q ] in = M [1 - A Q] in + Q in
#if 1
  void operatorADEF1(const FineField &in, FineField &out) {

    CoarseVector Csrc(_CoarseOperator.Grid());
    CoarseVector Ctmp(_CoarseOperator.Grid());
    CoarseVector Csol(_CoarseOperator.Grid());
    Csol = zero;

    ConjugateGradient<CoarseVector> CG(1.0e-10, 100000);
    ConjugateGradient<FineField>    fCG(3.0e-2, 1000);

    HermitianLinearOperator<CoarseOperator, CoarseVector> HermOp(_CoarseOperator);
    MdagMLinearOperator<CoarseOperator, CoarseVector>     MdagMOp(_CoarseOperator);
    ShiftedMdagMLinearOperator<Matrix, FineField>         fMdagMOp(_FineMatrix, 0.1);

    FineField tmp(in._grid);
    FineField res(in._grid);
    FineField Qin(in._grid);

    // Monitor completeness of low mode space
    //    _Aggregates.ProjectToSubspace  (Csrc,in);
    //    _Aggregates.PromoteFromSubspace(Csrc,out);
    //    std::cout<<GridLogMessage<<"Coarse Grid Preconditioner\nCompleteness in: "<<std::sqrt(norm2(out)/norm2(in))<<std::endl;

    _Aggregates.ProjectToSubspace(Csrc, in);
    HermOp.AdjOp(Csrc, Ctmp); // Normal equations
    CG(MdagMOp, Ctmp, Csol);
    _Aggregates.PromoteFromSubspace(Csol, Qin);

    //    Qin=0;
    _FineOperator.Op(Qin, tmp); // A Q in
    tmp = in - tmp;             // in - A Q in

    _FineOperator.Op(tmp, res); // this is the G5 herm bit
    fCG(fMdagMOp, res, out);    // solves  MdagM = g5 M g5M

    out = out + Qin;

    _FineOperator.Op(out, tmp);
    tmp = tmp - in; // tmp is new residual

    std::cout << GridLogMessage << "preconditioner thinks residual is " << std::sqrt(norm2(tmp) / norm2(in)) << std::endl;
  }
#endif

  void SAP(const FineField &src, FineField &psi) {

    Lattice<iScalar<vInteger>> coor(src._grid);
    Lattice<iScalar<vInteger>> subset(src._grid);

    FineField r(src._grid);
    FineField zz(src._grid);
    zz = zero;
    FineField vec1(src._grid);
    FineField vec2(src._grid);

    const Integer block = params.domainsize;

    subset = zero;
    for(int mu = 0; mu < Nd; mu++) {
      LatticeCoordinate(coor, mu + 1);
      coor   = div(coor, block);
      subset = subset + coor;
    }
    subset = mod(subset, (Integer)2);

    ShiftedMdagMLinearOperator<Matrix, FineField> fMdagMOp(_SmootherMatrix, 0.0);
    Chebyshev<FineField>                          Cheby(params.lo, params.hi, params.order, InverseApproximation);

    RealD resid;
    for(int i = 0; i < params.steps; i++) {

      // Even domain residual
      _FineOperator.Op(psi, vec1); // this is the G5 herm bit
      r     = src - vec1;
      resid = norm2(r) / norm2(src);
      std::cout << "SAP " << i << " resid " << resid << std::endl;

      // Even domain solve
      r = where(subset == (Integer)0, r, zz);
      _SmootherOperator.AdjOp(r, vec1);
      Cheby(fMdagMOp, vec1, vec2); // solves  MdagM = g5 M g5M
      psi = psi + vec2;

      // Odd domain residual
      _FineOperator.Op(psi, vec1); // this is the G5 herm bit
      r = src - vec1;
      r = where(subset == (Integer)1, r, zz);

      resid = norm2(r) / norm2(src);
      std::cout << "SAP " << i << " resid " << resid << std::endl;

      // Odd domain solve
      _SmootherOperator.AdjOp(r, vec1);
      Cheby(fMdagMOp, vec1, vec2); // solves  MdagM = g5 M g5M
      psi = psi + vec2;

      _FineOperator.Op(psi, vec1); // this is the G5 herm bit
      r     = src - vec1;
      resid = norm2(r) / norm2(src);
      std::cout << "SAP " << i << " resid " << resid << std::endl;
    }
  };

  void SmootherTest(const FineField &in) {

    FineField vec1(in._grid);
    FineField vec2(in._grid);

    RealD lo[3] = {0.5, 1.0, 2.0};

    //    MdagMLinearOperator<Matrix,FineField>        fMdagMOp(_FineMatrix);
    ShiftedMdagMLinearOperator<Matrix, FineField> fMdagMOp(_SmootherMatrix, 0.0);

    RealD Ni, r;

    Ni = norm2(in);

    for(int ilo = 0; ilo < 3; ilo++) {
      for(int ord = 5; ord < 50; ord *= 2) {

        _SmootherOperator.AdjOp(in, vec1);

        Chebyshev<FineField> Cheby(lo[ilo], 70.0, ord, InverseApproximation);
        Cheby(fMdagMOp, vec1, vec2); // solves  MdagM = g5 M g5M

        _FineOperator.Op(vec2, vec1); // this is the G5 herm bit
        vec1 = in - vec1;             // tmp  = in - A Min
        r    = norm2(vec1);
        std::cout << GridLogMessage << "Smoother resid " << std::sqrt(r / Ni) << std::endl;
      }
    }
  }

  void operatorCheby(const FineField &in, FineField &out) {

    CoarseVector Csrc(_CoarseOperator.Grid());
    CoarseVector Ctmp(_CoarseOperator.Grid());
    CoarseVector Csol(_CoarseOperator.Grid());
    Csol = zero;

    ConjugateGradient<CoarseVector> CG(3.0e-3, 100000);
    //    ConjugateGradient<FineField>    fCG(3.0e-2,1000);

    HermitianLinearOperator<CoarseOperator, CoarseVector> HermOp(_CoarseOperator);
    MdagMLinearOperator<CoarseOperator, CoarseVector>     MdagMOp(_CoarseOperator);
    //    MdagMLinearOperator<Matrix,FineField>        fMdagMOp(_FineMatrix);
    ShiftedMdagMLinearOperator<Matrix, FineField> fMdagMOp(_SmootherMatrix, 0.0);

    FineField vec1(in._grid);
    FineField vec2(in._grid);

    //    Chebyshev<FineField> Cheby    (0.5,70.0,30,InverseApproximation);
    //    Chebyshev<FineField> ChebyAccu(0.5,70.0,30,InverseApproximation);
    Chebyshev<FineField> Cheby(params.lo, params.hi, params.order, InverseApproximation);
    Chebyshev<FineField> ChebyAccu(params.lo, params.hi, params.order, InverseApproximation);
    //    Cheby.JacksonSmooth();
    //    ChebyAccu.JacksonSmooth();

    //    _Aggregates.ProjectToSubspace  (Csrc,in);
    //    _Aggregates.PromoteFromSubspace(Csrc,out);
    //    std::cout<<GridLogMessage<<"Completeness: "<<std::sqrt(norm2(out)/norm2(in))<<std::endl;

    //    ofstream fout("smoother");
    //    Cheby.csv(fout);

    // V11 multigrid.
    // Use a fixed chebyshev and hope hermiticity helps.

    // To make a working smoother for indefinite operator
    // must multiply by "Mdag" (ouch loses all low mode content)
    // and apply to poly approx of (mdagm)^-1.
    // so that we end up with an odd polynomial.

    RealD Ni = norm2(in);

    _SmootherOperator.AdjOp(in, vec1); // this is the G5 herm bit
    ChebyAccu(fMdagMOp, vec1, out);    // solves  MdagM = g5 M g5M

    std::cout << GridLogMessage << "Smoother norm " << norm2(out) << std::endl;

    // Update with residual for out
    _FineOperator.Op(out, vec1); // this is the G5 herm bit
    vec1 = in - vec1;            // tmp  = in - A Min

    RealD r = norm2(vec1);

    std::cout << GridLogMessage << "Smoother resid " << std::sqrt(r / Ni) << " " << r << " " << Ni << std::endl;

    _Aggregates.ProjectToSubspace(Csrc, vec1);
    HermOp.AdjOp(Csrc, Ctmp); // Normal equations
    CG(MdagMOp, Ctmp, Csol);
    _Aggregates.PromoteFromSubspace(Csol, vec1); // Ass^{-1} [in - A Min]_s
                                                 // Q = Q[in - A Min]
    out = out + vec1;

    // Three preconditioner smoothing -- hermitian if C3 = C1
    // Recompute error
    _FineOperator.Op(out, vec1); // this is the G5 herm bit
    vec1 = in - vec1;            // tmp  = in - A Min
    r    = norm2(vec1);

    std::cout << GridLogMessage << "Coarse resid " << std::sqrt(r / Ni) << std::endl;

    // Reapply smoother
    _SmootherOperator.Op(vec1, vec2); // this is the G5 herm bit
    ChebyAccu(fMdagMOp, vec2, vec1);  // solves  MdagM = g5 M g5M

    out  = out + vec1;
    vec1 = in - vec1; // tmp  = in - A Min
    r    = norm2(vec1);
    std::cout << GridLogMessage << "Smoother resid " << std::sqrt(r / Ni) << std::endl;
  }

  void operatorSAP(const FineField &in, FineField &out) {

    CoarseVector Csrc(_CoarseOperator.Grid());
    CoarseVector Ctmp(_CoarseOperator.Grid());
    CoarseVector Csol(_CoarseOperator.Grid());
    Csol = zero;

    ConjugateGradient<CoarseVector> CG(1.0e-3, 100000);

    HermitianLinearOperator<CoarseOperator, CoarseVector> HermOp(_CoarseOperator);
    MdagMLinearOperator<CoarseOperator, CoarseVector>     MdagMOp(_CoarseOperator);

    FineField vec1(in._grid);
    FineField vec2(in._grid);

    _Aggregates.ProjectToSubspace(Csrc, in);
    _Aggregates.PromoteFromSubspace(Csrc, out);
    std::cout << GridLogMessage << "Completeness: " << std::sqrt(norm2(out) / norm2(in)) << std::endl;

    // To make a working smoother for indefinite operator
    // must multiply by "Mdag" (ouch loses all low mode content)
    // and apply to poly approx of (mdagm)^-1.
    // so that we end up with an odd polynomial.
    SAP(in, out);

    // Update with residual for out
    _FineOperator.Op(out, vec1); // this is the G5 herm bit
    vec1 = in - vec1;            // tmp  = in - A Min

    RealD r  = norm2(vec1);
    RealD Ni = norm2(in);
    std::cout << GridLogMessage << "SAP resid " << std::sqrt(r / Ni) << " " << r << " " << Ni << std::endl;

    _Aggregates.ProjectToSubspace(Csrc, vec1);
    HermOp.AdjOp(Csrc, Ctmp); // Normal equations
    CG(MdagMOp, Ctmp, Csol);
    _Aggregates.PromoteFromSubspace(Csol, vec1); // Ass^{-1} [in - A Min]_s
                                                 // Q = Q[in - A Min]
    out = out + vec1;

    // Three preconditioner smoothing -- hermitian if C3 = C1
    // Recompute error
    _FineOperator.Op(out, vec1); // this is the G5 herm bit
    vec1 = in - vec1;            // tmp  = in - A Min
    r    = norm2(vec1);

    std::cout << GridLogMessage << "Coarse resid " << std::sqrt(r / Ni) << std::endl;

    // Reapply smoother
    SAP(vec1, vec2);
    out = out + vec2;

    // Update with residual for out
    _FineOperator.Op(out, vec1); // this is the G5 herm bit
    vec1 = in - vec1;            // tmp  = in - A Min

    r  = norm2(vec1);
    Ni = norm2(in);
    std::cout << GridLogMessage << "SAP resid(post) " << std::sqrt(r / Ni) << " " << r << " " << Ni << std::endl;
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
                << " | relative deviation = " << deviation << std::endl;

      if(deviation > tolerance) {
        std::cout << GridLogError << "Vector " << i << ": relative deviation check failed " << deviation << " > " << tolerance << std::endl;
        abort();
      }
    }
    std::cout << GridLogMessage << "Check passed!" << std::endl;

    std::cout << GridLogMessage << "**************************************************" << std::endl;
    std::cout << GridLogMessage << "MG correctness check: 0 == (1 - R P) v_c" << std::endl;
    std::cout << GridLogMessage << "**************************************************" << std::endl;

    random(cGrids.PRNGs[whichCoarseGrid], cTmps[0]);

    _Aggregates.PromoteFromSubspace(cTmps[0], fTmps[0]); //   P v_c
    _Aggregates.ProjectToSubspace(cTmps[1], fTmps[0]);   // R P v_c

    cTmps[2]       = cTmps[0] - cTmps[1]; // v_c - R P v_c
    auto deviation = std::sqrt(norm2(cTmps[2]) / norm2(cTmps[0]));

    std::cout << GridLogMessage << "norm2(v_c) = " << norm2(cTmps[0]) << " | norm2(R P v_c) = " << norm2(cTmps[1])
              << " | norm2(P v_c) = " << norm2(fTmps[0]) << " | relative deviation = " << deviation << std::endl;

    if(deviation > tolerance) {
      std::cout << GridLogError << "relative deviation check failed " << deviation << " > " << tolerance << std::endl;
      abort();
    }
    std::cout << GridLogMessage << "Check passed!" << std::endl;

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
              << " | relative deviation = " << deviation << std::endl;

    if(deviation > tolerance) {
      std::cout << GridLogError << "relative deviation check failed " << deviation << " > " << tolerance << std::endl;
      abort();
    }
    std::cout << GridLogMessage << "Check passed!" << std::endl;

    std::cout << GridLogMessage << "**************************************************" << std::endl;
    std::cout << GridLogMessage << "MG correctness check: 0 == |(Im(v_c^dag D_c^dag D_c v_c)|" << std::endl;
    std::cout << GridLogMessage << "**************************************************" << std::endl;

    random(cGrids.PRNGs[whichCoarseGrid], cTmps[0]);

    MdagMOp.Op(cTmps[0], cTmps[1]);    //         D_c v_c
    MdagMOp.AdjOp(cTmps[1], cTmps[2]); // D_c^dag D_c v_c

    // // alternative impl, which is better?
    // MdagMOp.HermOp(cTmps[0], cTmps[2]); // D_c^dag D_c v_c

    auto dot  = innerProduct(cTmps[0], cTmps[2]); //v_c^dag D_c^dag D_c v_c
    deviation = abs(imag(dot)) / abs(real(dot));

    std::cout << GridLogMessage << "Re(v_c^dag D_c^dag D_c v_c) = " << real(dot) << " | Im(v_c^dag D_c^dag D_c v_c) = " << imag(dot)
              << " | relative deviation = " << deviation << std::endl;

    if(deviation > tolerance) {
      std::cout << GridLogError << "relative deviation check failed " << deviation << " > " << tolerance << std::endl;
      abort();
    }
    std::cout << GridLogMessage << "Check passed!" << std::endl;
  }
};

int main(int argc, char **argv) {

  Grid_init(&argc, &argv);

  params.domainsize      = 1;
  params.coarsegrids     = 1;
  params.domaindecompose = 0;
  params.order           = 30;
  params.Ls              = 1;
  // params.mq = .13;
  params.mq    = .5;
  params.lo    = 0.5;
  params.hi    = 70.0;
  params.steps = 1;

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
  LatticeFermion        src(FGrid); gaussian(fPRNG, src); // src=src + g5 * src;
  LatticeFermion     result(FGrid); result = zero;
  LatticeFermion        ref(FGrid); ref = zero;
  LatticeFermion        tmp(FGrid);
  LatticeFermion        err(FGrid);
  LatticeGaugeField     Umu(FGrid); SU3::HotConfiguration(fPRNG, Umu);
  LatticeGaugeField   UmuDD(FGrid);
  LatticeColourMatrix     U(FGrid);
  LatticeColourMatrix     zz(FGrid);
  // clang-format on

  if(params.domaindecompose) {
    Lattice<iScalar<vInteger>> coor(FGrid);
    zz = zero;
    for(int mu = 0; mu < Nd; mu++) {
      LatticeCoordinate(coor, mu);
      U = PeekIndex<LorentzIndex>(Umu, mu);
      U = where(mod(coor, params.domainsize) == (Integer)0, zz, U);
      PokeIndex<LorentzIndex>(UmuDD, U, mu);
    }
  } else {
    UmuDD = Umu;
  }

  RealD mass = params.mq;

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Set up some coarser levels stuff: " << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  const int nbasis = 20; // we fix the number of test vector to the same
                         // number on every level for now

  //////////////////////////////////////////
  // toggle to run two/three level method
  //////////////////////////////////////////

  // // two-level algorithm
  // std::vector<std::vector<int>> blockSizes({{2, 2, 2, 2}});
  // CoarseGrids<nbasis>           coarseGrids(blockSizes, 1);

  // three-level algorithm
  std::vector<std::vector<int>> blockSizes({{2, 2, 2, 2}, {2, 2, 1, 1}});
  CoarseGrids<nbasis>           coarseGrids(blockSizes, 2);

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Building the wilson operator on the fine grid" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  WilsonFermionR Dw(Umu, *FGrid, *FrbGrid, mass);
  WilsonFermionR DwDD(UmuDD, *FGrid, *FrbGrid, mass);

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Some typedefs" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  // typedefs for transition from fine to first coarsened grid
  typedef vSpinColourVector                                                                 FineSiteVector;
  typedef vTComplex                                                                         CoarseSiteScalar;
  typedef Aggregation<FineSiteVector, CoarseSiteScalar, nbasis>                             Subspace;
  typedef CoarsenedMatrix<FineSiteVector, CoarseSiteScalar, nbasis>                         CoarseOperator;
  typedef CoarseOperator::CoarseVector                                                      CoarseVector;
  typedef CoarseOperator::siteVector                                                        CoarseSiteVector;
  typedef TestVectorAnalyzer<LatticeFermion, nbasis>                                        FineTVA;
  typedef MultiGridPreconditioner<FineSiteVector, CoarseSiteScalar, nbasis, WilsonFermionR> FineMGPreconditioner;
  typedef TrivialPrecon<LatticeFermion>                                                     FineTrivialPreconditioner;

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
  std::cout << GridLogMessage << "Calling Aggregation class to build subspaces" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  MdagMLinearOperator<WilsonFermionR, LatticeFermion> FineHermPosdefOp(Dw);
  Subspace                                            FineAggregates(coarseGrids.Grids[0], FGrid, 0);

  assert((nbasis & 0x1) == 0);
  int nb = nbasis / 2;
  std::cout << GridLogMessage << " nbasis/2 = " << nb << std::endl;

  FineAggregates.CreateSubspace(fPRNG, FineHermPosdefOp /*, nb */); // Don't specify nb to see the orthogonalization check

  std::cout << GridLogMessage << "Test vector analysis after initial creation of MG test vectors" << std::endl;
  FineTVA fineTVA;
  fineTVA(FineHermPosdefOp, FineAggregates.subspace, nb);

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

  // using Gamma5HermitianLinearOperator corresponds to working with H = g5 * D
  Gamma5HermitianLinearOperator<WilsonFermionR, LatticeFermion> FineHermIndefOp(Dw);
  Gamma5HermitianLinearOperator<WilsonFermionR, LatticeFermion> FineHermIndefOpDD(DwDD);
  CoarseOperator                                                Dc(*coarseGrids.Grids[0]);
  Dc.CoarsenOperator(FGrid, FineHermIndefOp, FineAggregates); // uses only linop.OpDiag & linop.OpDir

  std::cout << GridLogMessage << "Test vector analysis after construction of D_c" << std::endl;
  fineTVA(FineHermPosdefOp, FineAggregates.subspace, nb);

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Building coarse vectors" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  CoarseVector coarseSource(coarseGrids.Grids[0]);
  CoarseVector coarseResult(coarseGrids.Grids[0]);
  gaussian(coarseGrids.PRNGs[0], coarseSource);
  coarseResult = zero;

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Testing some coarse space solvers" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  MdagMLinearOperator<CoarseOperator, CoarseVector> CoarseHermPosdefOp(Dc);

  std::vector<std::unique_ptr<OperatorFunction<CoarseVector>>> coarseSolvers;
  coarseSolvers.emplace_back(new GeneralisedMinimalResidual<CoarseVector>(5.0e-2, 100, 8, false));
  coarseSolvers.emplace_back(new MinimalResidual<CoarseVector>(5.0e-2, 100, false));
  coarseSolvers.emplace_back(new ConjugateGradient<CoarseVector>(5.0e-2, 100, false));

  for(auto const &solver : coarseSolvers) {
    coarseResult = zero;
    (*solver)(CoarseHermPosdefOp, coarseSource, coarseResult);
  }

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Testing the operators" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  std::cout << GridLogMessage << "MdagMLinearOperator<WilsonFermionR, LatticeFermion> FineHermPosdefOp(Dw);" << std::endl;
  testOperator(FineHermPosdefOp, FGrid);
  std::cout << GridLogMessage << "Gamma5HermitianLinearOperator<WilsonFermionR, LatticeFermion> FineHermIndefOp(Dw);" << std::endl;
  testOperator(FineHermIndefOp, FGrid);
  std::cout << GridLogMessage << "Gamma5HermitianLinearOperator<WilsonFermionR, LatticeFermion> FineHermIndefOpDD(DwDD);" << std::endl;
  testOperator(FineHermIndefOpDD, FGrid);
  std::cout << GridLogMessage << "MdagMLinearOperator<CoarseOperator, CoarseVector> CoarseHermPosdefOp(Dc);" << std::endl;
  testOperator(CoarseHermPosdefOp, coarseGrids.Grids[0]);

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Building deflation preconditioner " << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  FineMGPreconditioner FineMGPrecon(FineAggregates, Dc, FineHermIndefOp, Dw, FineHermIndefOp, Dw);

  FineMGPreconditioner FineMGPreconDD(FineAggregates, Dc, FineHermIndefOp, Dw, FineHermIndefOpDD, DwDD);

  FineTrivialPreconditioner FineSimplePrecon;

  FineMGPrecon.runChecks(coarseGrids, 0);

  if(coarseGrids.LattSizes.size() == 2) {

    std::cout << GridLogMessage << "**************************************************" << std::endl;
    std::cout << GridLogMessage << "Dummy testing for building a second coarse level" << std::endl;
    std::cout << GridLogMessage << "**************************************************" << std::endl;

    SubSubSpace CoarseAggregates(coarseGrids.Grids[1], coarseGrids.Grids[0], 0);
    CoarseAggregates.CreateSubspace(coarseGrids.PRNGs[0], CoarseHermPosdefOp);

    // // this doesn't work because this function applies g5 to a vector, which
    // // doesn't work for coarse vectors atm -> FIXME
    // CoarseTVA coarseTVA;
    // coarseTVA(CoarseHermPosdefOp, CoarseAggregates.subspace, nb);

    // // cannot apply g5 to coarse vectors atm -> FIXME
    // for(int n=0;n<nb;n++){
    //   CoarseAggregates.subspace[n+nb] = g5 * CoarseAggregates.subspace[n]; // multiply with g5 normally instead of G5R5 since this specific to DWF
    //   std::cout<<GridLogMessage<<n<<" subspace "<<norm2(CoarseAggregates.subspace[n+nb])<<" "<<norm2(CoarseAggregates.subspace[n]) <<std::endl;
    // }

    auto coarseCoarseSites = 1;
    for(auto const &elem : coarseGrids.LattSizes[1]) coarseCoarseSites *= elem;

    std::cout << GridLogMessage << "Norms of MG test vectors after chiral projection (coarse coarse sites = " << coarseCoarseSites << ")"
              << std::endl;
    for(int n = 0; n < nbasis; n++) {
      std::cout << GridLogMessage << "vec[" << n << "] = " << norm2(CoarseAggregates.subspace[n]) << std::endl;
    }

    CoarseCoarseOperator Dcc(*coarseGrids.Grids[1]);
    Dcc.CoarsenOperator(coarseGrids.Grids[0], CoarseHermPosdefOp, CoarseAggregates); // uses only linop.OpDiag & linop.OpDir

    // // this doesn't work because this function applies g5 to a vector, which
    // // doesn't work for coarse vectors atm -> FIXME
    // std::cout << GridLogMessage << "Test vector analysis after construction of D_c_c" << std::endl;
    // coarseTVA(CoarseHermPosdefOp, CoarseAggregates.subspace, nb);

    CoarseCoarseVector coarseCoarseSource(coarseGrids.Grids[1]);
    CoarseCoarseVector coarseCoarseResult(coarseGrids.Grids[1]);
    gaussian(coarseGrids.PRNGs[1], coarseCoarseSource);
    coarseCoarseResult = zero;

    MdagMLinearOperator<CoarseCoarseOperator, CoarseCoarseVector> CoarseCoarseHermPosdefOp(Dcc);

    std::vector<std::unique_ptr<OperatorFunction<CoarseCoarseVector>>> coarseCoarseSolvers;
    coarseSolvers.emplace_back(new GeneralisedMinimalResidual<CoarseCoarseVector>(5.0e-2, 100, 8, false));
    coarseSolvers.emplace_back(new MinimalResidual<CoarseCoarseVector>(5.0e-2, 100, false));
    coarseSolvers.emplace_back(new ConjugateGradient<CoarseCoarseVector>(5.0e-2, 100, false));

    for(auto const &solver : coarseCoarseSolvers) {
      coarseCoarseResult = zero;
      (*solver)(CoarseCoarseHermPosdefOp, coarseCoarseSource, coarseCoarseResult);
    }

    CoarseMGPreconditioner      CoarseMGPrecon(CoarseAggregates, Dcc, CoarseHermPosdefOp, Dc, CoarseHermPosdefOp, Dc);

    CoarseMGPrecon.runChecks(coarseGrids, 1);

    std::cout << GridLogMessage << "ARTIFICIAL ABORT" << std::endl;
    abort();
  }

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Building VPGCR and FGMRES solvers w/ & w/o MG Preconditioner" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  std::vector<std::unique_ptr<OperatorFunction<LatticeFermion>>> solvers;
  solvers.emplace_back(new PrecGeneralisedConjugateResidual<LatticeFermion>(1.0e-12, 100, FineMGPrecon, 8, 8));
  solvers.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 100, FineMGPrecon, 8));
  solvers.emplace_back(new PrecGeneralisedConjugateResidual<LatticeFermion>(1.0e-12, 4000000, FineSimplePrecon, 8, 8));
  solvers.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(1.0e-12, 4000000, FineSimplePrecon, 8));

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Testing the solvers" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

  std::cout << GridLogMessage << "checking norm src " << norm2(src) << std::endl;

  for(auto const &solver : solvers) {
    result = zero;
    (*solver)(FineHermIndefOp, src, result);
  }

  Grid_finalize();
}
