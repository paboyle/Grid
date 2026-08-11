/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_pvdagm_3level.cc

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
  std::string      _name;
  
  int    level;  void Level(int lv) {level = lv; };

  MGPreconditioner(Aggregates &Agg,
		   FineOperator &Fine,
		   FineSmoother &PreSmoother,
		   FineSmoother &PostSmoother,
		   CoarseOperator &CoarseOperator_,
		   CoarseSolver &CoarseSolve_,
		   CoarseSolver &CoarseGuesser_,
		   std::string name = std::string("unnamed"))
    : _Aggregates(Agg),
      _FineOperator(Fine),
      _PreSmoother(PreSmoother),
      _PostSmoother(PostSmoother),
      _CoarseOperator(CoarseOperator_),
      _CoarseSolve(CoarseSolve_),
      _CoarseGuesser(CoarseGuesser_),
      _name(name),
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
    std::cout<<GridLogMessage << _name << "PreSmoother took "<< t/1000.0<< "ms" <<std::endl;

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
    std::cout<<GridLogMessage << _name << "Promote to this level took "<< t/1000.0<< "ms" <<std::endl;

    _FineOperator.Op(out,vec1);  sub(vec1 ,in , vec1);

    t=-usecond();
    vec2=Zero();
    _PostSmoother(vec1,vec2);
    t+=usecond();
    std::cout<<GridLogMessage << _name <<"PostSmoother took "<< t/1000.0<< "ms" <<std::endl;

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
  c_src = one;
  c_res = Zero();
  L2PGCR_baseline(c_src,c_res);
  PVdagM.getApplications();
  PVdagM.clear();

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

  //////////////////////////////////////////////////////////////////////
  // Level 2 solver: plain GCR, no further coarsening
  //////////////////////////////////////////////////////////////////////
  TrivialPrecon<CoarseCoarseVector> simpleCC;
  // L3PGCR is an inner solver inside the L1→2 V-cycle; does not need to converge
  // to fine-grid precision. Loose tolerance (3e-2) and large restart (64) to allow
  // the Krylov space to span enough of the near-null spectrum of LinOpCC per cycle.
  PrecGeneralisedConjugateResidualNonHermitian<CoarseCoarseVector> L3PGCR(1.0e-4,5,LinOpCC,simpleCC,64,64);
  L3PGCR.Level(3);

  //////////////////////////////////////////////////////////////////////
  // Coarse-level GCR smoother for Level 1→2 V-cycle.
  // Mirrors fine-grid SmootherGCR: shifted operator + fixed step count.
  // coarse_smoother_shift and coarse_smoother_nstep are the tuning knobs.
  //////////////////////////////////////////////////////////////////////
  RealD coarse_smoother_shift = 0.0;
  int   coarse_smoother_nstep = 8;
  if(getenv("coarse_smoother_shift")) coarse_smoother_shift = atof(getenv("coarse_smoother_shift"));
  if(getenv("coarse_smoother_nstep")) coarse_smoother_nstep = atoi(getenv("coarse_smoother_nstep"));
     
  ShiftedLinearOperator<CoarseVector> ShiftedLinOpCoarse(coarse_smoother_shift, LinOpCoarse);
  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector> CoarseSmootherGCR(0.0,1,ShiftedLinOpCoarse,simpleC,coarse_smoother_nstep,coarse_smoother_nstep);
  CoarseSmootherGCR.Level(2);

  //////////////////////////////////////////////////////////////////////
  // Level 1→2 V-cycle preconditioner.
  //////////////////////////////////////////////////////////////////////
  L1to2MG L1to2Precon(AggregatesL2,
		      LinOpCoarse,
		      simpleC,            // no pre-smoother (matches fine-grid setup)
		      CoarseSmootherGCR,  // post-smoother: 12 GCR steps
		      LinOpCC,
		      L3PGCR,
		      simpleCC,
		      std::string("LinOpC"));

  //////////////////////////////////////////////////////////////////////
  // Standalone Level 1 two-level solve test.
  // Compare against plain PGCR baseline above.
  //////////////////////////////////////////////////////////////////////
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<" Level 1 solve: two-level MG preconditioned PGCR"<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;

  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector> L2MGsolver(3.0e-2,200,LinOpCoarse,L1to2Precon,16,16);
  L2MGsolver.Level(2);
  c_res = Zero();
  L2MGsolver(c_src,c_res);

  std::cout << GridLogMessage << "Level 1 two-level test: PVdagM operator uses:" << std::endl;
  PVdagM.getApplications();
  PVdagM.clear();

  //////////////////////////////////////////////////////////////////////
  // Full three-level outer solve
  //////////////////////////////////////////////////////////////////////
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<" Three-level outer solve"<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;

  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD> SmootherGCR(0.00,1,ShiftedPVdagM,simple_fine,16,16);
  SmootherGCR.Level(1);

  f_src = one;

  // Pre-smoother: none (TrivialPrecon); post-smoother: shifted PGCR.
  // Coarse solver: L2MGsolver (PGCR preconditioned by Level 1→2 V-cycle).
  TwoLevelMG ThreeLevelPrecon(AggregatesPD,
			      PVdagM,
			      simple_fine,
			      SmootherGCR,
			      LinOpCoarse,
			      L2MGsolver,
			      simpleC,
			      std::string("PVdagM"));

  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermion> L1PGCR(1.0e-8,1000,PVdagM,ThreeLevelPrecon,16,16);
  L1PGCR.Level(1);

  f_res = Zero();
  L1PGCR(f_src,f_res);

  std::cout << GridLogMessage << "Three-level outer solve: PVdagM operator uses:" << std::endl;
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

  RealD madj = 1.0;
  if ( getenv("MADJ") ) madj=atof(getenv("MADJ"));
  std::cout << "PV mass set to "<<madj<<std::endl;

  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b_,c_);
  MobiusFermionD Dpv (Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,madj, M5,b_,c_);

  typedef PVdagMLinearOperator<MobiusFermionD,LatticeFermionD>        PVdagM_t;
  typedef ShiftedPVdagMLinearOperator<MobiusFermionD,LatticeFermionD> ShiftedPVdagM_t;
  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>  LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector                           CoarseVector;
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis>             Subspace;
  typedef MGPreconditioner<vSpinColourVector,vTComplex,nbasis>        TwoLevelMG;

  PVdagM_t        PVdagM(Ddwf,Dpv);
  ShiftedPVdagM_t ShiftedPVdagM(0.00,Ddwf,Dpv);

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
    // Re-orthogonalise after loading to ensure block-GS condition holds.
    //    AggregatesGCR.Orthogonalise();
    std::cout << GridLogMessage << "Subspace loaded and re-orthogonalised." << std::endl;
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
    geom,
    PVdagM,
    ShiftedPVdagM,
    AggregatesGCR
  );

  std::cout << GridLogMessage << "Done" << std::endl;
  Grid_finalize();
  return 0;
}
