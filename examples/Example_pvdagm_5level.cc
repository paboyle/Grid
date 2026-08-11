/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_pvdagm_5level.cc

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

#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidual.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>
#include <Grid/algorithms/iterative/BiCGSTAB.h>

using namespace std;
using namespace Grid;

template <class T> void readFile(T& out, std::string const fname){
  #ifdef HAVE_LIME
    std::cout << Grid::GridLogMessage << "Reading: " << fname << std::endl;
    Grid::emptyUserRecord record;
    Grid::ScidacReader SR;
    SR.open(fname);
    SR.readScidacFieldRecord(out, record);
    SR.close();
  #endif
}

template <class T> void writeFile(T& in, std::string const fname){
  #ifdef HAVE_LIME
    std::cout << Grid::GridLogMessage << "Writing: " << fname << std::endl;
    Grid::emptyUserRecord record;
    Grid::ScidacWriter SW(in.Grid()->IsBoss());
    SW.open(fname);
    SW.writeScidacFieldRecord(in, record);
    SW.close();
  #endif
}

template <class Field>
void saveSubspace(std::vector<Field> &subspace, std::string const fname){
  #ifdef HAVE_LIME
    std::cout << Grid::GridLogMessage << "Saving subspace (" << subspace.size() << " vectors) to: " << fname << std::endl;
    Grid::emptyUserRecord record;
    Grid::ScidacWriter SW(subspace[0].Grid()->IsBoss());
    SW.open(fname);
    for (int k = 0; k < (int)subspace.size(); k++)
      SW.writeScidacFieldRecord(subspace[k], record);
    SW.close();
  #endif
}

template <class Field>
void loadSubspace(std::vector<Field> &subspace, std::string const fname){
  #ifdef HAVE_LIME
    std::cout << Grid::GridLogMessage << "Loading subspace (" << subspace.size() << " vectors) from: " << fname << std::endl;
    Grid::emptyUserRecord record;
    Grid::ScidacReader SR;
    SR.open(fname);
    for (int k = 0; k < (int)subspace.size(); k++)
      SR.readScidacFieldRecord(subspace[k], record);
    SR.close();
  #endif
}

template<class Matrix,class Field>
class PVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
  int nApp;
  int nAppDag;
public:
  PVdagMLinearOperator(Matrix &Mat,Matrix &PV): _Mat(Mat),_PV(PV), nApp(0), nAppDag(0) {};

  void OpDiag (const Field &in, Field &out) {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
  void Op     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _Mat.M(in,tmp);
    _PV.Mdag(tmp,out);
    nApp++;
  }
  void AdjOp     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _PV.M(in,tmp);
    _Mat.Mdag(tmp,out);
    nAppDag++;
  }
  void clear() { nApp = 0; nAppDag = 0; }
  void getApplications() {
    std::cout << GridLogMessage << "# applications of PVdagM: " << nApp << std::endl;
    std::cout << GridLogMessage << "# applications of PVdagM^dag: " << nAppDag << std::endl;
    std::cout << GridLogMessage << "# applications total: " << nApp + nAppDag << std::endl;
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    HermOp(in,out);
    ComplexD dot = innerProduct(in,out);
    n1=real(dot);
    n2=norm2(out);
  }
  void HermOp(const Field &in, Field &out){
    Field tmp(in.Grid());
    Op(in,tmp);
    AdjOp(tmp,out);
  }
};

template<class Matrix,class Field>
class MdagPVLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
public:
  MdagPVLinearOperator(Matrix &Mat,Matrix &PV): _Mat(Mat),_PV(PV){};

  void OpDiag (const Field &in, Field &out) {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
  void Op     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _PV.M(in,tmp);
    _Mat.Mdag(tmp,out);
  }
  void AdjOp     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _Mat.M(in,tmp);
    _PV.Mdag(tmp,out);
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    ComplexD dot = innerProduct(in,out);
    n1=real(dot);
    n2=norm2(out);
  }
  void HermOp(const Field &in, Field &out){
    Field tmp(in.Grid());
    Op(in,tmp);
    AdjOp(tmp,out);
  }
};

template<class Matrix,class Field>
class ShiftedPVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
  RealD shift;
public:
  ShiftedPVdagMLinearOperator(RealD _shift,Matrix &Mat,Matrix &PV): shift(_shift),_Mat(Mat),_PV(PV){};

  void OpDiag (const Field &in, Field &out) {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
  void Op     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _Mat.M(in,tmp);
    _PV.Mdag(tmp,out);
    out = out + shift * in;
  }
  void AdjOp     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _PV.M(tmp,out);
    _Mat.Mdag(in,tmp);
    out = out + shift * in;
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ assert(0); }
  void HermOp(const Field &in, Field &out){
    Field tmp(in.Grid());
    Op(in,tmp);
    AdjOp(tmp,out);
  }
};

// Lüscher deflated guesser (arXiv:0706.2298 Sec A.3) for a non-Hermitian solve.
// C_{st} = <psi[s] | LinOp | psi[t]>;  guess = sum_s c_s psi[s]  where c = C^{-1} psi† src.
template<class Field>
class LuscherGuesser : public LinearFunction<Field> {
  const std::vector<Field> &psi;
  Eigen::MatrixXcd          C_inv;
public:
  using LinearFunction<Field>::operator();
  LuscherGuesser(const std::vector<Field> &psi_, const Eigen::MatrixXcd &Cinv_)
    : psi(psi_), C_inv(Cinv_) {}
  virtual void operator()(const Field &src, Field &guess) {
    int N = psi.size();
    Eigen::VectorXcd b(N);
    for (int t = 0; t < N; t++)
      b(t) = TensorRemove(innerProduct(psi[t], src));
    Eigen::VectorXcd c = C_inv * b;
    guess = Zero();
    for (int s = 0; s < N; s++)
      guess += ComplexD(c(s)) * psi[s];
  }
};

template<class Fobj,class CComplex,int nbasis>
class MGPreconditioner : public LinearFunction< Lattice<Fobj> > {
public:
  using LinearFunction<Lattice<Fobj> >::operator();

  typedef Aggregation<Fobj,CComplex,nbasis> Aggregates;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::FineField    FineField;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseVector CoarseVector;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseMatrix CoarseMatrix;
  typedef LinearOperatorBase<FineField>                            FineOperator;
  typedef LinearFunction    <FineField>                            FineSmoother;
  typedef LinearOperatorBase<CoarseVector>                         CoarseOperator;
  typedef LinearFunction    <CoarseVector>                         CoarseSolver;

  Aggregates     & _Aggregates;
  FineOperator   & _FineOperator;
  FineSmoother   & _PreSmoother;
  FineSmoother   & _PostSmoother;
  CoarseOperator & _CoarseOperator;
  CoarseSolver   & _CoarseSolve;
  CoarseSolver   & _CoarseGuesser;

  int    level;  void Level(int lv) {level = lv; };

  MGPreconditioner(Aggregates &Agg,
		   FineOperator &Fine,
		   FineSmoother &PreSmoother,
		   FineSmoother &PostSmoother,
		   CoarseOperator &CoarseOperator_,
		   CoarseSolver &CoarseSolve_,
		   CoarseSolver &CoarseGuesser_)
    : _Aggregates(Agg),
      _FineOperator(Fine),
      _PreSmoother(PreSmoother),
      _PostSmoother(PostSmoother),
      _CoarseOperator(CoarseOperator_),
      _CoarseSolve(CoarseSolve_),
      _CoarseGuesser(CoarseGuesser_),
      level(1)  {  }

  virtual void operator()(const FineField &in, FineField & out)
  {
    GridBase *CoarseGrid = _Aggregates.CoarseGrid;
    CoarseVector Csrc(CoarseGrid);
    CoarseVector Csol(CoarseGrid);
    FineField vec1(in.Grid());
    FineField vec2(in.Grid());

    double t;
    out = Zero();
    t=-usecond();
    _PreSmoother(in,out);
    t+=usecond();
    std::cout<<GridLogMessage << "PreSmoother took "<< t/1000.0<< "ms" <<std::endl;

    _FineOperator.Op(out,vec1);  sub(vec1, in ,vec1);

    t=-usecond();
    _Aggregates.ProjectToSubspace(Csrc,vec1);
    t+=usecond();
    std::cout<<GridLogMessage << "Project to coarse took "<< t/1000.0<< "ms" <<std::endl;

    t=-usecond();
    _CoarseGuesser(Csrc,Csol);
    _CoarseSolve(Csrc,Csol);
    t+=usecond();
    std::cout<<GridLogMessage << "Coarse solve took "<< t/1000.0<< "ms" <<std::endl;

    t=-usecond();
    _Aggregates.PromoteFromSubspace(Csol,vec1);
    add(out,out,vec1);
    t+=usecond();
    std::cout<<GridLogMessage << "Promote to this level took "<< t/1000.0<< "ms" <<std::endl;

    _FineOperator.Op(out,vec1);  sub(vec1 ,in , vec1);

    t=-usecond();
    vec2=Zero();
    _PostSmoother(vec1,vec2);
    t+=usecond();
    std::cout<<GridLogMessage << "PostSmoother took "<< t/1000.0<< "ms" <<std::endl;

    add(out,out,vec2);
  }
};

// Generic shifted linear operator: wraps any LinearOperatorBase and adds shift*I.
// Used to condition the coarse-level GCR smoother, analogous to ShiftedPVdagMLinearOperator
// at the fine level.
template<class Field>
class ShiftedLinearOperator : public LinearOperatorBase<Field> {
  LinearOperatorBase<Field> &_Op;
  RealD shift;
public:
  ShiftedLinearOperator(RealD _shift, LinearOperatorBase<Field> &Op) : shift(_shift), _Op(Op) {}
  void OpDiag   (const Field &in, Field &out)                         { assert(0); }
  void OpDir    (const Field &in, Field &out, int dir, int disp)      { assert(0); }
  void OpDirAll (const Field &in, std::vector<Field> &out)            { assert(0); }
  void Op       (const Field &in, Field &out) { _Op.Op(in, out);    out = out + shift * in; }
  void AdjOp    (const Field &in, Field &out) { _Op.AdjOp(in, out); out = out + shift * in; }
  void HermOpAndNorm(const Field &in, Field &out, RealD &n1, RealD &n2) { assert(0); }
  void HermOp   (const Field &in, Field &out) { Field tmp(in.Grid()); Op(in,tmp); AdjOp(tmp,out); }
};

template<int NB, class PVdagM_t, class ShiftedPVdagM_t, class Subspace, class LittleDiracOperator, class CoarseVector, class TwoLevelMG>
void runMG(
  GridCartesian *FGrid,
  GridCartesian *Coarse5d,
  GridCartesian *CoarseCoarse5d,
  GridCartesian *CoarseCoarseCoarse5d,
  GridCartesian *CoarseCoarseCoarseCoarse5d,
  NextToNearestStencilGeometry5D geom,
  PVdagM_t &PVdagM,
  ShiftedPVdagM_t &ShiftedPVdagM,
  Subspace &AggregatesPD
) {
  std::vector<LatticeFermion> subspace = AggregatesPD.subspace;
  assert((int)subspace.size() == NB);
  const int nbasis = NB;
  const int cb = 0;

  CoarseVector c_src(Coarse5d);
  CoarseVector c_res(Coarse5d);
  Complex one(1.0);

  LatticeFermionD f_src(FGrid);
  LatticeFermionD f_res(FGrid);

  TrivialPrecon<CoarseVector>      simpleC;
  TrivialPrecon<LatticeFermionD>   simple_fine;

  //////////////////////////////////////////////////////////////////////
  // Level 0→1: coarsen PVdagM, build LinOpCoarse
  //////////////////////////////////////////////////////////////////////
  LittleDiracOperator LittleDiracOpPV(geom, FGrid, Coarse5d);
  LittleDiracOpPV.CoarsenOperator(PVdagM, AggregatesPD);

  NonHermitianLinearOperator<LittleDiracOperator,CoarseVector> LinOpCoarse(LittleDiracOpPV);

  //////////////////////////////////////////////////////////////////////
  // Baseline: plain PGCR on LinOpCoarse (reference for comparison)
  //////////////////////////////////////////////////////////////////////
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<" Level 1 solve: plain PGCR baseline"<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;

  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector> L2PGCR_baseline(3.0e-2,1100,LinOpCoarse,simpleC,10,10);
  L2PGCR_baseline.Level(2);
  L2PGCR_baseline.Name("Cbaseline");
  c_src = one;
  c_res = Zero();
  L2PGCR_baseline(c_src,c_res);

  //////////////////////////////////////////////////////////////////////
  // psi_coarse: coarse projections of pre-GS fine null vectors.
  // These are the Level 1 near-null vectors, promoted from Level 0.
  // Used as the aggregation basis for Level 1→2 coarsening.
  //////////////////////////////////////////////////////////////////////
  std::vector<CoarseVector> psi_coarse(nbasis, Coarse5d);
  for (int k = 0; k < nbasis; k++)
    AggregatesPD.ProjectToSubspace(psi_coarse[k], subspace[k]);

  //////////////////////////////////////////////////////////////////////
  // Diagnostics: W (fine projected matrix) and C (Galerkin check)
  //////////////////////////////////////////////////////////////////////
  {
    Eigen::MatrixXcd W = Eigen::MatrixXcd::Zero(nbasis, nbasis);
    LatticeFermion ftmp(FGrid);
    for (int j = 0; j < nbasis; j++) {
      PVdagM.Op(subspace[j], ftmp);
      for (int i = 0; i < nbasis; i++)
        W(i,j) = TensorRemove(innerProduct(subspace[i], ftmp));
    }
    RealD normW = W.norm();
    std::cout << GridLogMessage << "Fine projected matrix ||W|| = " << normW << std::endl;

    Eigen::MatrixXcd C = Eigen::MatrixXcd::Zero(nbasis, nbasis);
    CoarseVector Ac(Coarse5d);
    for (int l = 0; l < nbasis; l++) {
      LinOpCoarse.Op(psi_coarse[l], Ac);
      for (int k = 0; k < nbasis; k++)
        C(k,l) = TensorRemove(innerProduct(psi_coarse[k], Ac));
    }
    RealD normC      = C.norm();
    RealD normCmCdag = (C - C.adjoint()).norm();
    std::cout << GridLogMessage << "Coarse null matrix ||C||            = " << normC << std::endl;
    std::cout << GridLogMessage << "Coarse null matrix ||C - C†||/||C|| = " << normCmCdag/normC << std::endl;
    std::cout << GridLogMessage << "Galerkin check ||C||/||W||          = " << normC/normW << std::endl;
  }

  //////////////////////////////////////////////////////////////////////
  // Level 1→2: set up aggregation using psi_coarse as subspace.
  // Block factor 2,2,3,2 (removes odd local sublattice in z given MPI
  // geometry 3×6×4×4 where z-local at Level 1 is 6).
  // psi_coarse are assigned directly; CoarsenOperator performs
  // block-GS orthogonalisation before building LinOpCoarseCoarse.
  //////////////////////////////////////////////////////////////////////
  // innerProduct(CoarseSiteObj, CoarseSiteObj) returns iScalar<vTComplex>, so CComplex
  // for the L1→L2 level must be iScalar<vTComplex>, not vTComplex.
  typedef typename CoarseVector::vector_object                            CoarseSiteObj;
  typedef iScalar<vTComplex>                                              vTTComplex;
  typedef GeneralCoarsenedMatrix<CoarseSiteObj,vTTComplex,NB>            LittleDiracOperatorL2;
  typedef typename LittleDiracOperatorL2::CoarseVector                   CoarseCoarseVector;
  typedef Aggregation<CoarseSiteObj,vTTComplex,NB>                       SubspaceL2;
  typedef MGPreconditioner<CoarseSiteObj,vTTComplex,NB>                  L1to2MG;

  SubspaceL2 AggregatesL2(CoarseCoarse5d, Coarse5d, cb);
  for (int k = 0; k < nbasis; k++)
    AggregatesL2.subspace[k] = psi_coarse[k];

  NextToNearestStencilGeometry5D geom2(CoarseCoarse5d);
  LittleDiracOperatorL2 LittleDiracOpL2(geom2, Coarse5d, CoarseCoarse5d);
  LittleDiracOpL2.CoarsenOperator(LinOpCoarse, AggregatesL2);

  NonHermitianLinearOperator<LittleDiracOperatorL2,CoarseCoarseVector> LinOpCC(LittleDiracOpL2);

  TrivialPrecon<CoarseCoarseVector> simpleCC;

  //////////////////////////////////////////////////////////////////////
  // Lüscher deflation guesser for L3PGCR.
  // Step 1: project psi_coarse[k] (promoted fine null vectors) to
  //         CoarseCoarseVector space — these cover the zero-momentum
  //         component of the near-null space of LinOpCC.
  // Step 2: breed Nextra additional null vectors directly on LinOpCC
  //         using GCR with random sources — these pick up near-null
  //         modes at all spatial frequencies not spanned by step 1.
  // Step 3: build C_{st} = <psi_cc[s]|LinOpCC|psi_cc[t]> over the
  //         full augmented basis and invert directly via Eigen LU.
  //////////////////////////////////////////////////////////////////////
  std::vector<CoarseCoarseVector> psi_cc(nbasis, CoarseCoarse5d);
  for (int k = 0; k < nbasis; k++)
    AggregatesL2.ProjectToSubspace(psi_cc[k], psi_coarse[k]);

  {
    int Nextra = nbasis;  // breed as many extra as we have promoted ones
    if ( getenv("CC_NEXTRA") ) Nextra = atoi(getenv("CC_NEXTRA"));
    GridParallelRNG RNG_CC(CoarseCoarse5d);
    RNG_CC.SeedFixedIntegers({11,13,17,19});
    PrecGeneralisedConjugateResidualNonHermitian<CoarseCoarseVector>
      nullGCR(1e-2, 200, LinOpCC, simpleCC, 32, 32);
    CoarseCoarseVector tmp(CoarseCoarse5d);
    for (int k = 0; k < Nextra; k++) {
      CoarseCoarseVector src(CoarseCoarse5d);
      gaussian(RNG_CC, src);
      tmp = Zero();
      nullGCR(src, tmp);
      psi_cc.push_back(tmp);
    }
    std::cout << GridLogMessage << "LinOpCC deflation basis: " << nbasis
              << " promoted + " << Nextra << " bred = " << psi_cc.size() << " total" << std::endl;
  }

  const int Naug = psi_cc.size();
  Eigen::MatrixXcd Ccc = Eigen::MatrixXcd::Zero(Naug, Naug);
  {
    CoarseCoarseVector Acc(CoarseCoarse5d);
    for (int l = 0; l < Naug; l++) {
      LinOpCC.Op(psi_cc[l], Acc);
      for (int k = 0; k < Naug; k++)
        Ccc(k,l) = TensorRemove(innerProduct(psi_cc[k], Acc));
    }
  }
  {
    RealD normCcc      = Ccc.norm();
    RealD normCccmCdag = (Ccc - Ccc.adjoint()).norm();
    std::cout << GridLogMessage << "Coarse-coarse deflation matrix ||Ccc||              = " << normCcc << std::endl;
    std::cout << GridLogMessage << "Coarse-coarse deflation matrix ||Ccc-Ccc†||/||Ccc|| = " << normCccmCdag/normCcc << std::endl;
  }
  Eigen::MatrixXcd Ccc_inv = Ccc.inverse();
  LuscherGuesser<CoarseCoarseVector> CCDeflGuesser(psi_cc, Ccc_inv);

  //////////////////////////////////////////////////////////////////////
  // Level 2→3: coarsen LinOpCC using the RAW promoted psi_cc as aggregation
  // to build the Level 4 (coarse-coarse-coarse) operator.
  //   psi_cc[0..nbasis-1] are the coarse-coarse near-null vectors, projected
  //   from the RAW psi_coarse (themselves projected from the RAW fine null
  //   vectors) -- the pre-block-GS chain the whole construction depends on.
  //   CoarsenOperator block-GS orthogonalises AggregatesL3.subspace IN PLACE,
  //   so assign COPIES of psi_cc and keep psi_cc itself raw.
  //
  // Tensor depth deepens once more: innerProduct(CoarseCoarseSiteObj,...) returns
  // iScalar<vTTComplex>, so CComplex for the L2→L3 level is iScalar<iScalar<vTComplex>>.
  //////////////////////////////////////////////////////////////////////
  typedef typename CoarseCoarseVector::vector_object                     CoarseCoarseSiteObj;
  typedef iScalar<vTTComplex>                                            vTTTComplex;
  typedef GeneralCoarsenedMatrix<CoarseCoarseSiteObj,vTTTComplex,NB>     LittleDiracOperatorL3;
  typedef typename LittleDiracOperatorL3::CoarseVector                   CoarseCoarseCoarseVector;
  typedef Aggregation<CoarseCoarseSiteObj,vTTTComplex,NB>               SubspaceL3;
  typedef MGPreconditioner<CoarseCoarseSiteObj,vTTTComplex,NB>          L2to3MG;

  SubspaceL3 AggregatesL3(CoarseCoarseCoarse5d, CoarseCoarse5d, cb);
  for (int k = 0; k < nbasis; k++)
    AggregatesL3.subspace[k] = psi_cc[k];          // raw promoted; COPY, keeps psi_cc raw

  NextToNearestStencilGeometry5D geom3(CoarseCoarseCoarse5d);
  LittleDiracOperatorL3 LittleDiracOpL3(geom3, CoarseCoarse5d, CoarseCoarseCoarse5d);
  LittleDiracOpL3.CoarsenOperator(LinOpCC, AggregatesL3);   // block-GS's AggregatesL3.subspace in place

  NonHermitianLinearOperator<LittleDiracOperatorL3,CoarseCoarseCoarseVector> LinOpCCC(LittleDiracOpL3);
  TrivialPrecon<CoarseCoarseCoarseVector> simpleCCC;

  //////////////////////////////////////////////////////////////////////
  // Level 3→4: coarsen LinOpCCC to build the Level 5 operator, using a
  // TRUNCATED basis of only the first NB5 (< nbasis) raw promoted null vectors.
  //   psi_ccc[k] = raw psi_cc projected through the (block-GS'd) L3 aggregation
  //   -- the pre-block-GS chain continued one level deeper.  We keep only the
  //   leading NB5: after the global orthogonalisation of the original fine null
  //   vectors the early indices retain the most-null content (shared low-mode
  //   components are peeled in first), so the leading NB5 are the crudely-most-
  //   null slice.  This is the cheap "first 30" truncation test; a principled
  //   sigma-ordered rotation of psi_ccc would replace the slice, not the idea.
  //   NB: a positive result is conservative (sigma-ordering can only help); a
  //   negative one is inconclusive until the sigma-ordered NB5 is tried.
  //
  // Tensor depth deepens once more: CComplex for the L3→L4 level is
  // iScalar<vTTTComplex>.  NB5 (the coarse dimension) is independent of the
  // depth -- it just makes the coarsest site vector NB5-dimensional.
  //////////////////////////////////////////////////////////////////////
  const int NB5 = 30;   // compile-time: changing it re-instantiates the L4/L5 tensors
  std::cout << GridLogMessage << "PARAM NB5 (truncated coarsest basis) = " << NB5 << std::endl;
  assert(NB5 <= nbasis);

  std::vector<CoarseCoarseCoarseVector> psi_ccc(nbasis, CoarseCoarseCoarse5d);
  for (int k = 0; k < nbasis; k++)
    AggregatesL3.ProjectToSubspace(psi_ccc[k], psi_cc[k]);   // raw psi_cc -> L4 null vectors

  //////////////////////////////////////////////////////////////////////
  // Optional sigma-ordering of psi_ccc (SVD_REORDER set): replace the crude
  // first-NB5 slice with the NB5 genuinely-most-null directions of span(psi_ccc)
  // under LinOpCCC.  For a NON-NORMAL operator the nullness measure is the
  // singular value of A restricted to the span -- eig of Q†A†AQ -- NOT the
  // numerical range Q†AQ (which non-normality contaminates).  Robust route:
  // whiten by the Gram (drop near-dependent directions), Hermitian-eig the
  // whitened A†A, rotate.  The printed singular spectrum IS the SVD study: where
  // it falls off tells you the natural NB5, and the same numbers illuminate why
  // the earlier singular-subspace deflation re-entered.  Safe here because we
  // ORDER vectors that then feed a Galerkin projection, not REMOVE a subspace.
  // Default (unset) leaves psi_ccc in raw order == the "first 30" test.
  //////////////////////////////////////////////////////////////////////
  if ( getenv("SVD_REORDER") ) {
    std::cout << GridLogMessage << "SVD_REORDER: sigma-ordering psi_ccc under LinOpCCC" << std::endl;

    Eigen::MatrixXcd G(nbasis,nbasis);                    // Gram  = Psi^dag Psi
    for (int i=0;i<nbasis;i++)
    for (int j=0;j<nbasis;j++)
      G(i,j) = TensorRemove(innerProduct(psi_ccc[i],psi_ccc[j]));

    std::vector<CoarseCoarseCoarseVector> Apsi(nbasis, CoarseCoarseCoarse5d);
    for (int j=0;j<nbasis;j++) LinOpCCC.Op(psi_ccc[j], Apsi[j]);

    Eigen::MatrixXcd M(nbasis,nbasis);                    // A^dagA = Psi^dag A^dag A Psi
    for (int i=0;i<nbasis;i++)
    for (int j=0;j<nbasis;j++)
      M(i,j) = TensorRemove(innerProduct(Apsi[i],Apsi[j]));

    // Whiten by the Gram: G = Ug diag(g) Ug^dag; keep g > tol*max; T = Ug diag(1/sqrt g).
    // Q = Psi T is then orthonormal (Q^dag Q = T^dag G T = I).
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> esG(G);
    Eigen::VectorXd g = esG.eigenvalues();               // ascending, real
    RealD gmax = g(nbasis-1);
    RealD gtol = 1.0e-9 * gmax;
    int keep = 0; for (int i=0;i<nbasis;i++) if (g(i) > gtol) keep++;
    std::cout << GridLogMessage << "  Gram spectrum: min=" << g(0) << " max=" << gmax
              << " cond=" << gmax/std::max(g(0),1.0e-300) << " keep=" << keep << "/" << nbasis << std::endl;
    assert(keep >= NB5);

    Eigen::MatrixXcd T(nbasis, keep);                    // whitening (largest-g first)
    { int c=0;
      for (int i=nbasis-1;i>=0;i--) if (g(i) > gtol) { T.col(c) = esG.eigenvectors().col(i)/std::sqrt(g(i)); c++; }
    }

    Eigen::MatrixXcd Mw = T.adjoint() * M * T;           // whitened A^dagA (keep x keep, Hermitian)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> esM(Mw);
    Eigen::VectorXd s2 = esM.eigenvalues();              // ascending sigma^2 (most-null first)
    std::cout << GridLogMessage << "  Singular spectrum sigma_k (most-null first):" << std::endl;
    for (int k=0;k<keep;k++)
      std::cout << GridLogMessage << "    sigma[" << k << "] = " << std::sqrt(std::max(s2(k),0.0)) << std::endl;

    Eigen::MatrixXcd R = T * esM.eigenvectors();         // coeffs over Psi, sigma-ordered orthonormal dirs
    std::vector<CoarseCoarseCoarseVector> phi(keep, CoarseCoarseCoarse5d);
    for (int k=0;k<keep;k++) {
      phi[k] = Zero();
      for (int j=0;j<nbasis;j++)
        phi[k] = phi[k] + ComplexD(R(j,k)) * psi_ccc[j];
    }
    for (int k=0;k<keep;k++) psi_ccc[k] = phi[k];        // psi_ccc[0..NB5-1] now = most-null dirs
    std::cout << GridLogMessage << "SVD_REORDER: psi_ccc replaced by sigma-ordered directions" << std::endl;
  }

  typedef typename CoarseCoarseCoarseVector::vector_object              CoarseCoarseCoarseSiteObj;
  typedef iScalar<vTTTComplex>                                          vTTTTComplex;
  typedef GeneralCoarsenedMatrix<CoarseCoarseCoarseSiteObj,vTTTTComplex,NB5> LittleDiracOperatorL4;
  typedef typename LittleDiracOperatorL4::CoarseVector                  CoarseCoarseCoarseCoarseVector;
  typedef Aggregation<CoarseCoarseCoarseSiteObj,vTTTTComplex,NB5>       SubspaceL4;
  typedef MGPreconditioner<CoarseCoarseCoarseSiteObj,vTTTTComplex,NB5>  L3to4MG;

  SubspaceL4 AggregatesL4(CoarseCoarseCoarseCoarse5d, CoarseCoarseCoarse5d, cb);
  for (int k = 0; k < NB5; k++)
    AggregatesL4.subspace[k] = psi_ccc[k];         // FIRST NB5 raw promoted vectors (truncation)

  NextToNearestStencilGeometry5D geom4(CoarseCoarseCoarseCoarse5d);
  LittleDiracOperatorL4 LittleDiracOpL4(geom4, CoarseCoarseCoarse5d, CoarseCoarseCoarseCoarse5d);
  LittleDiracOpL4.CoarsenOperator(LinOpCCC, AggregatesL4);   // block-GS's AggregatesL4.subspace in place

  NonHermitianLinearOperator<LittleDiracOperatorL4,CoarseCoarseCoarseCoarseVector> LinOpCCCC(LittleDiracOpL4);
  TrivialPrecon<CoarseCoarseCoarseCoarseVector> simpleCCCC;

  //////////////////////////////////////////////////////////////////////
  // Level 5 bottom solve: GCR on a SHIFTED LinOpCCCC (the coarsest, most
  // non-normal operator).  l5_shift slides its field of values off the origin;
  // defaults to 0.0 (bare LinOpCCCC) until opted in.  This is the level a dense
  // direct inverse would eventually replace: rank = NB5 * sites(clatt4).
  //////////////////////////////////////////////////////////////////////
  RealD l5_shift = 0.0;
  if(getenv("l5_shift")) l5_shift = atof(getenv("l5_shift"));
  std::cout << GridLogMessage << "PARAM l5_shift = " << l5_shift << std::endl;

  ShiftedLinearOperator<CoarseCoarseCoarseCoarseVector> ShiftedLinOpCCCC(l5_shift, LinOpCCCC);
  PrecGeneralisedConjugateResidualNonHermitian<CoarseCoarseCoarseCoarseVector> L5PGCR(1.0e-1,200,ShiftedLinOpCCCC,simpleCCCC,16,16);
  L5PGCR.Level(5);
  L5PGCR.Name("CCCCouter");

  //////////////////////////////////////////////////////////////////////
  // Level 3→4 V-cycle: depth-2 SHIFTED smoother on LinOpCCC + Level 5 bottom.
  // Level 4 is no longer the bottom -- it is smoothed shallowly and recursed to
  // Level 5, mirroring how Level 3 recurses to Level 4.
  //////////////////////////////////////////////////////////////////////
  RealD ccc_smoother_shift = 0.05;
  int   ccc_smoother_nstep = 2;
  if(getenv("ccc_smoother_shift")) ccc_smoother_shift = atof(getenv("ccc_smoother_shift"));
  if(getenv("ccc_smoother_nstep")) ccc_smoother_nstep = atoi(getenv("ccc_smoother_nstep"));

  ShiftedLinearOperator<CoarseCoarseCoarseVector> ShiftedLinOpCCC(ccc_smoother_shift, LinOpCCC);
  PrecGeneralisedConjugateResidualNonHermitian<CoarseCoarseCoarseVector>
    CoarseCoarseCoarseSmootherGCR(0.01,1,ShiftedLinOpCCC,simpleCCC,ccc_smoother_nstep,ccc_smoother_nstep);
  CoarseCoarseCoarseSmootherGCR.SetZeroGuess(1);  // smoother slot: caller zeroes guess
  CoarseCoarseCoarseSmootherGCR.Level(4);
  CoarseCoarseCoarseSmootherGCR.Name("CCCsmoother");

  L3to4MG L3to4Precon(AggregatesL4,
                      LinOpCCC,
                      simpleCCC,                     // no pre-smoother
                      CoarseCoarseCoarseSmootherGCR, // post-smoother: depth-2 shifted GCR
                      LinOpCCCC,
                      L5PGCR,
                      simpleCCCC);                   // trivial guesser at the bottom

  //////////////////////////////////////////////////////////////////////
  // Level 4 (coarse-coarse-coarse) solve: GCR preconditioned by the L3→L4 V-cycle.
  //////////////////////////////////////////////////////////////////////
  PrecGeneralisedConjugateResidualNonHermitian<CoarseCoarseCoarseVector> L4MGsolver(1.0e-1,200,LinOpCCC,L3to4Precon,16,16);
  L4MGsolver.Level(4);
  L4MGsolver.Name("CCCouter");

  //////////////////////////////////////////////////////////////////////
  // Level 2→3 V-cycle: depth-2 SHIFTED smoother on LinOpCC + Level 4 solve.
  // The shift slides the coarse-coarse field of values off the origin so a
  // 2-step smoother has something to bite on a non-normal operator (IRS idea).
  //////////////////////////////////////////////////////////////////////
  RealD cc_smoother_shift = 0.01;
  int   cc_smoother_nstep = 2;
  if(getenv("cc_smoother_shift")) cc_smoother_shift = atof(getenv("cc_smoother_shift"));
  if(getenv("cc_smoother_nstep")) cc_smoother_nstep = atoi(getenv("cc_smoother_nstep"));

  ShiftedLinearOperator<CoarseCoarseVector> ShiftedLinOpCC(cc_smoother_shift, LinOpCC);
  PrecGeneralisedConjugateResidualNonHermitian<CoarseCoarseVector>
    CoarseCoarseSmootherGCR(0.01,1,ShiftedLinOpCC,simpleCC,cc_smoother_nstep,cc_smoother_nstep);
  CoarseCoarseSmootherGCR.SetZeroGuess(1);  // smoother slot: caller zeroes guess
  CoarseCoarseSmootherGCR.Level(3);
  CoarseCoarseSmootherGCR.Name("CCsmoother");

  L2to3MG L2to3Precon(AggregatesL3,
                      LinOpCC,
                      simpleCC,                  // no pre-smoother
                      CoarseCoarseSmootherGCR,   // post-smoother: depth-2 shifted GCR
                      LinOpCCC,
                      L4MGsolver,                // coarse solve is now the L3→L4 V-cycle
                      simpleCCC);                // trivial guesser

  //////////////////////////////////////////////////////////////////////
  // Level 3 (coarse-coarse) solve: GCR preconditioned by the L2→L3 V-cycle.
  // Replaces the plain L3PGCR of the 3-level build -- the coarse-coarse level
  // is now smoothed shallowly and recursed rather than solved deeply.
  //////////////////////////////////////////////////////////////////////
  PrecGeneralisedConjugateResidualNonHermitian<CoarseCoarseVector> L3MGsolver(1.0e-1,200,LinOpCC,L2to3Precon,16,16);
  L3MGsolver.Level(3);
  L3MGsolver.Name("CCouter");

  //////////////////////////////////////////////////////////////////////
  // Coarse-level GCR smoother for Level 1→2 V-cycle.
  // Mirrors fine-grid SmootherGCR: shifted operator + fixed step count.
  // coarse_smoother_shift and coarse_smoother_nstep are the tuning knobs.
  //////////////////////////////////////////////////////////////////////
  RealD coarse_smoother_shift = 0.01;
  int   coarse_smoother_nstep = 2;   // depth-2 smoother on the coarse level
  if(getenv("coarse_smoother_shift")) coarse_smoother_shift = atof(getenv("coarse_smoother_shift"));
  if(getenv("coarse_smoother_nstep")) coarse_smoother_nstep = atoi(getenv("coarse_smoother_nstep"));
     
  ShiftedLinearOperator<CoarseVector> ShiftedLinOpCoarse(coarse_smoother_shift, LinOpCoarse);
  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector> CoarseSmootherGCR(0.01,1,ShiftedLinOpCoarse,simpleC,coarse_smoother_nstep,coarse_smoother_nstep);
  CoarseSmootherGCR.SetZeroGuess(1);   // smoother slot: caller zeroes guess
  CoarseSmootherGCR.Level(2);
  CoarseSmootherGCR.Name("Csmoother");

  //////////////////////////////////////////////////////////////////////
  // Level 1→2 V-cycle preconditioner.
  //////////////////////////////////////////////////////////////////////
  L1to2MG L1to2Precon(AggregatesL2,
                       LinOpCoarse,
                       simpleC,            // no pre-smoother (matches fine-grid setup)
                       CoarseSmootherGCR,  // post-smoother: depth-2 shifted GCR
                       LinOpCC,
                       L3MGsolver,         // coarse-coarse solve is now the L2→L3 V-cycle
                       CCDeflGuesser);     // Lüscher guesser: psi_cc C^{-1} psi_cc†

  //////////////////////////////////////////////////////////////////////
  // Standalone Level 1 two-level solve test.
  // Compare against plain PGCR baseline above.
  //////////////////////////////////////////////////////////////////////
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<" Level 1 solve: two-level MG preconditioned PGCR"<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;

  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector> L2MGsolver(3.0e-2,200,LinOpCoarse,L1to2Precon,16,16);
  L2MGsolver.Level(2);
  L2MGsolver.Name("Couter");
  c_res = Zero();
  L2MGsolver(c_src,c_res);

  std::cout << GridLogMessage << "Level 1 two-level test: PVdagM operator uses:" << std::endl;
  PVdagM.getApplications();
  PVdagM.clear();

  //////////////////////////////////////////////////////////////////////
  // Full five-level outer solve
  //////////////////////////////////////////////////////////////////////
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<" Five-level outer solve"<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;

  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD> SmootherGCR(0.01,1,ShiftedPVdagM,simple_fine,16,16);
  SmootherGCR.SetZeroGuess(1);         // pre/post smoother slots zero their guess
  SmootherGCR.Level(1);
  SmootherGCR.Name("Fsmoother");

  f_src = one;

  // Pre-smoother: none (TrivialPrecon); post-smoother: shifted PGCR.
  // Coarse solver: L2MGsolver (PGCR preconditioned by Level 1→2 V-cycle).
  TwoLevelMG ThreeLevelPrecon(AggregatesPD,
                               PVdagM,
                               simple_fine,
                               SmootherGCR,
                               LinOpCoarse,
                               L2MGsolver,
                               simpleC);

  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermion> L1PGCR(1.0e-8,1000,PVdagM,ThreeLevelPrecon,16,16);
  L1PGCR.Level(1);
  L1PGCR.Name("Fouter");

  f_res = Zero();
  L1PGCR(f_src,f_res);

  std::cout << GridLogMessage << "Five-level outer solve: PVdagM operator uses:" << std::endl;
  PVdagM.getApplications();
  PVdagM.clear();
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls   = 24;
  RealD     M5   = 1.8;
  RealD     b    = 1.5;
  RealD     c    = 0.5;
  RealD     mass = 0.00078;
  if ( getenv("MASS") ) mass = atof(getenv("MASS"));

  const int nbasis = 60;

  std::cout << GridLogMessage << "Mass: " << mass << ", Ls: " << Ls << ", b=" << b << ", c=" << c << std::endl;
  std::cout << GridLogMessage << "nbasis: " << nbasis << std::endl;

  std::vector<int> lat_size {48, 48, 48, 96};

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(lat_size, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Level 1 coarse grid: block 2^4 from fine (48×48×48×96 → 24×24×24×48, Ls=1)
  Coordinate clatt = lat_size;
  for (int d = 0; d < 4; d++) clatt[d] /= 2;
  std::cout << GridLogMessage << "Level 1 coarse lattice: " << clatt << std::endl;

  GridCartesian *Coarse4d  = SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian *Coarse5d  = SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);

  // Level 2 coarse-coarse grid: block 2,2,3,3 from Level 1 (24×24×24×48 → 12×12×8×16, Ls=1).
  // MPI geometry 3.6.4.4 (288 ranks): fine local {16,8,12,24}.
  // Level 1 local {8,4,6,12}; Level 2 local {4,2,2,4}.
  // z blocked by 3: z-Level1-local=6; 6/3=2 (even), 6/2=3 (odd) → must use 3.
  // t blocked by 3: t-Level1-local=12; 12/3=4 divisible by Nsimd=4 (gen-simd-width=64).
  //   t-block=2 gives t2-local=6, 6 mod 4 ≠ 0, fails Grid SIMD assertion. ✓
  // With {4,2,2,4}: Nsimd=4 goes into x or t (both =4). ✓
  Coordinate clatt2 = clatt;
  clatt2[0] /= 2;
  clatt2[1] /= 2;
  clatt2[2] /= 3;
  clatt2[3] /= 3;
  std::cout << GridLogMessage << "Level 2 coarse-coarse lattice: " << clatt2 << std::endl;

  GridCartesian *CoarseCoarse4d = SpaceTimeGrid::makeFourDimGrid(clatt2, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian *CoarseCoarse5d = SpaceTimeGrid::makeFiveDimGrid(1,CoarseCoarse4d);

  // Level 3 coarse-coarse-coarse grid: block clatt2 = {12,12,8,16} -> {6,12,8,8}.
  // GEOMETRY (mpi 3.6.4.4, Nsimd=4 => SIMD layout {1,1,2,2}, factor 2 on z and t):
  //   every grid needs z-local and t-local EVEN.  clatt2-local is {4,2,2,4}, so
  //   z-local=2 is already at its minimum even value and CANNOT be blocked (2->1
  //   is odd and trips the SIMD assertion); y-local=2 would go to 1 (degenerate).
  //   Only x and t have room, so block {2,1,1,2}: clatt3 {6,12,8,8}, L4-local
  //   {2,2,2,2} -- all dims even and >=2.  z stays unblocked by construction.
  Coordinate clatt3 = clatt2;
  clatt3[0] /= 2;   // x: 12 -> 6   (x-local 4 -> 2)
  // clatt3[1] (y) unblocked: y-local 2 -> blocking gives 1 (degenerate)
  // clatt3[2] (z) unblocked: z-local 2 is SIMD-pinned even, cannot halve
  clatt3[3] /= 2;   // t: 16 -> 8   (t-local 4 -> 2)
  std::cout << GridLogMessage << "Level 3 coarse-coarse-coarse lattice: " << clatt3 << std::endl;

  GridCartesian *CoarseCoarseCoarse4d = SpaceTimeGrid::makeFourDimGrid(clatt3, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian *CoarseCoarseCoarse5d = SpaceTimeGrid::makeFiveDimGrid(1,CoarseCoarseCoarse4d);

  // Level 4 coarse^4 grid: block clatt3 = {6,12,8,8} -> {3,6,8,8} via {2,2,1,1}.
  //   mpi 3.6.4.4 => clatt4-local {1,1,2,2}: z-local=2, t-local=2 stay EVEN (SIMD
  //   factor 2 pins them), so z,t are unblocked; x,y (SIMD factor 1) halve to
  //   local 1 -- fully distributed but legal for the halo-depth-1 NextToNearest
  //   stencil.  1152 sites; with NB5=30 that is the 34,560-rank coarsest operator
  //   a dense direct inverse would target.
  Coordinate clatt4 = clatt3;
  clatt4[0] /= 2;   // x: 6  -> 3   (x-local 2 -> 1)
  clatt4[1] /= 2;   // y: 12 -> 6   (y-local 2 -> 1)
  // clatt4[2] (z) unblocked: z-local 2 is SIMD-pinned even
  // clatt4[3] (t) unblocked: t-local 2 is SIMD-pinned even
  std::cout << GridLogMessage << "Level 4 coarse^4 lattice: " << clatt4 << std::endl;

  GridCartesian *CoarseCoarseCoarseCoarse4d = SpaceTimeGrid::makeFourDimGrid(clatt4, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian *CoarseCoarseCoarseCoarse5d = SpaceTimeGrid::makeFiveDimGrid(1,CoarseCoarseCoarseCoarse4d);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG  RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG  RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid);
  std::cout << GridLogMessage << "Reading gauge field" << std::endl;
  FieldMetaData header;
  std::string file("/ccs/home/poare/ckpoint_lat.1000");
  NerscIO::readConfiguration(Umu,header,file);

  RealD b_ = 1.5;
  RealD c_ = 0.5;
  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b_,c_);
  MobiusFermionD Dpv (Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,1.0, M5,b_,c_);

  typedef PVdagMLinearOperator<MobiusFermionD,LatticeFermionD>        PVdagM_t;
  typedef ShiftedPVdagMLinearOperator<MobiusFermionD,LatticeFermionD> ShiftedPVdagM_t;
  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>  LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector                           CoarseVector;
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis>             Subspace;
  typedef MGPreconditioner<vSpinColourVector,vTComplex,nbasis>        TwoLevelMG;

  PVdagM_t        PVdagM(Ddwf,Dpv);
  ShiftedPVdagM_t ShiftedPVdagM(0.01,Ddwf,Dpv);

  NextToNearestStencilGeometry5D geom(Coarse5d);

  // Subspace cache: save after generation, reload on subsequent runs to skip expensive setup.
  // Set SUBSPACE_FILE to override the default path.
  std::string subspace_file = "/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/subspace_nb"
                              + std::to_string(nbasis) + ".scidac";
  if ( getenv("SUBSPACE_FILE") ) subspace_file = std::string(getenv("SUBSPACE_FILE"));

  // Check if subspace file exists (boss rank checks, result broadcast via GlobalSum).
  uint64_t file_exists = 0;
  if ( UGrid->IsBoss() ) {
    std::ifstream f(subspace_file);
    file_exists = f.good() ? 1 : 0;
  }
  UGrid->GlobalSum(file_exists);

  const int cb = 0;
  Subspace AggregatesGCR(Coarse5d,FGrid,cb);

  if ( file_exists ) {
    std::cout << GridLogMessage << "*** Loading subspace from disk ***" << std::endl;
    loadSubspace(AggregatesGCR.subspace, subspace_file);
    // Insurance: GLOBAL (whole-lattice) orthonormalise, in case the cached file
    // predates the GlobalOrthonormalise() that CreateSubspaceGCR now applies
    // (Aggregates.h:196).  It is span-preserving and makes the vectors globally
    // orthonormal -- it is NOT the block Orthogonalise() below, so it does NOT
    // cause the psi_coarse->e_k trap.  It also (re)establishes the weak nullness
    // gradient (shared most-null components peeled into the early indices) that
    // the "first NB5" truncation relies on.  Idempotent if the file was already
    // globally orthonormal.  The RAW subspace copy in runMG happens AFTER this
    // call, so the raw-null (pre-block-GS) discipline is preserved.
    AggregatesGCR.GlobalOrthonormalise();
    // DO NOT block-orthogonalise here: runMG copies subspace[] as the RAW
    // (pre-block-GS) basis and CoarsenOperator block-GS's it in place later.
    // Orthogonalising now defeats the raw-null discipline (psi_coarse -> e_k)
    // and poisons L2/L3/L4.  See project_block_orthogonalise_leak.
    // AggregatesGCR.Orthogonalise();
    std::cout << GridLogMessage << "Subspace loaded, globally orthonormalised (raw block basis preserved)." << std::endl;
  } else {
    std::cout << GridLogMessage << "*** GCR subspace generation ***" << std::endl;
    AggregatesGCR.CreateSubspaceGCR(RNG5,PVdagM,nbasis);
    std::cout << GridLogMessage << "Subspace generation: PVdagM operator uses:" << std::endl;
    PVdagM.getApplications();
    PVdagM.clear();
    saveSubspace(AggregatesGCR.subspace, subspace_file);
    std::cout << GridLogMessage << "Subspace saved to: " << subspace_file << std::endl;
  }

  runMG<nbasis,PVdagM_t,ShiftedPVdagM_t,Subspace,LittleDiracOperator,CoarseVector,TwoLevelMG>(
    FGrid,
    Coarse5d,
    CoarseCoarse5d,
    CoarseCoarseCoarse5d,
    CoarseCoarseCoarseCoarse5d,
    geom,
    PVdagM,
    ShiftedPVdagM,
    AggregatesGCR
  );

  std::cout << GridLogMessage << "Done" << std::endl;
  Grid_finalize();
  return 0;
}
