/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_pvdagm_mrhs_3level_DenseCoarseMatrix.cc

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */

// MultiRHS (valence) THREE-level multigrid for PVdagM with a DENSE, EXACT,
// non-iterative coarse-coarse bottom -- the LIBRARY-CLASS successor of
// Example_pvdagm_mrhs_3level_dense.cc, which is FROZEN as the regression
// baseline / champion-provenance artifact (21.7 s/RHS at BLOCK=2.2.3.3,
// BLOCK2=8.4.2.4, nb60, CSO3/FSO6/CST0.04 on 36 Frontier nodes).
//
// The dense bottom is now Grid/algorithms/multigrid/DenseCoarseMatrix.h:
//  - stencil -> dense DIRECT import (no probe assembly: rows are local data)
//    + IMPORT CERTIFICATE  (DENSE_IMPORT_SIGN=-1 flips convention, no rebuild)
//  - split-K apply via GridBLAS.gemmBatched with explicit leading dimensions
//    (DENSE_SPLITK chunks, default 32) -- the fig-11 software split-K
//  - deviceVector / GridBLAS throughout the apply: platform-agnostic
//
// INTERCHANGE: same SLAB_FILE per-rank format as the frozen example (stem MUST
// encode cfg/mass/blocking/nbasis; the header guards only N/nrows/nbasis) and
// the same env-var set, so existing sbatch scripts drive either binary.
//
// A/B acceptance (old binary = control):
//   slab-cached  : outer counts match EXACTLY (identical apply data; split-K
//                  changes only fp32 reduction order); wall delta = split-K gain.
//   fresh setup  : outer equal-or-+-1 (import vs probe = rounding); VERIFY
//                  ~7e-4 both; setup delta = import gain (~93 s probe retired).
//
// Level structure, solvers, and tuning knobs are UNCHANGED from the frozen
// example.  Env: MASS SUBSPACE_FILE NRHS BLOCK BLOCK2 FineSmootherShift/Order
// CoarseSmootherShift/Nstep CoarseSolverTol/Order DENSE_CC DENSE_CC_CHECK
// DENSE_SPLITK DENSE_DEVICE_SUM DENSE_IMPORT_SIGN DENSE_APPLY_PROFILE
// L3_TOL L3_MAXIT L3_NSTEP OuterMmax OuterNstep OuterTol

#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>
#include <Grid/algorithms/multigrid/DenseCoarseMatrix.h>

#include <memory>

using namespace std;
using namespace Grid;

RealD FineSmootherShift    = 0.1;
int   FineSmootherOrder    = 16;
RealD CoarseSmootherShift  = 0.1;
int   CoarseSmootherNstep  = 4;
RealD CoarseSolverTol      = 0.03;
int   CoarseSolverOrder    = 200;
RealD L3Tol                = 2.5e-1;
int   L3MaxIt              = 50;
int   L3Nstep              = 50;
RealD OuterTol             = 1.0e-8;
int   OuterMmax            = 8;
int   OuterNstep           = 8;
int   Nrhs                 = 12;
int   UseDenseCC           = 1;
RealD mass                 = 0.00078;

void ParseEnvironment(void)
{
  if(getenv("MASS"))               mass               = atof(getenv("MASS"));
  if(getenv("FineSmootherShift"))  FineSmootherShift  = atof(getenv("FineSmootherShift"));
  if(getenv("FineSmootherOrder"))  FineSmootherOrder  = atoi(getenv("FineSmootherOrder"));
  if(getenv("CoarseSmootherShift"))CoarseSmootherShift= atof(getenv("CoarseSmootherShift"));
  if(getenv("CoarseSmootherNstep"))CoarseSmootherNstep= atoi(getenv("CoarseSmootherNstep"));
  if(getenv("CoarseSolverTol"))    CoarseSolverTol    = atof(getenv("CoarseSolverTol"));
  if(getenv("CoarseSolverOrder"))  CoarseSolverOrder  = atoi(getenv("CoarseSolverOrder"));
  if(getenv("L3_TOL"))             L3Tol              = atof(getenv("L3_TOL"));
  if(getenv("L3_MAXIT"))           L3MaxIt            = atoi(getenv("L3_MAXIT"));
  if(getenv("L3_NSTEP"))           L3Nstep            = atoi(getenv("L3_NSTEP"));
  if(getenv("OuterTol"))           OuterTol           = atof(getenv("OuterTol"));
  if(getenv("OuterMmax"))          OuterMmax          = atoi(getenv("OuterMmax"));
  if(getenv("OuterNstep"))         OuterNstep         = atoi(getenv("OuterNstep"));
  if(getenv("NRHS"))               Nrhs               = atoi(getenv("NRHS"));
  if(getenv("DENSE_CC"))           UseDenseCC         = atoi(getenv("DENSE_CC"));

  std::cout << GridLogMessage << "PARAM: MASS               " << mass               << std::endl;
  std::cout << GridLogMessage << "PARAM: NRHS               " << Nrhs               << std::endl;
  std::cout << GridLogMessage << "PARAM: DENSE_CC           " << UseDenseCC         << std::endl;
  std::cout << GridLogMessage << "PARAM: FineSmootherShift  " << FineSmootherShift  << std::endl;
  std::cout << GridLogMessage << "PARAM: FineSmootherOrder  " << FineSmootherOrder  << std::endl;
  std::cout << GridLogMessage << "PARAM: CoarseSmootherShift" << CoarseSmootherShift<< std::endl;
  std::cout << GridLogMessage << "PARAM: CoarseSmootherNstep" << CoarseSmootherNstep<< std::endl;
  std::cout << GridLogMessage << "PARAM: CoarseSolverTol    " << CoarseSolverTol    << std::endl;
  std::cout << GridLogMessage << "PARAM: CoarseSolverOrder  " << CoarseSolverOrder  << std::endl;
  std::cout << GridLogMessage << "PARAM: L3_TOL             " << L3Tol              << std::endl;
  std::cout << GridLogMessage << "PARAM: OuterMmax          " << OuterMmax          << std::endl;
  std::cout << GridLogMessage << "PARAM: OuterNstep         " << OuterNstep         << std::endl;
}

template <class Field>
void saveSubspace(std::vector<Field> &subspace, std::string const fname){
#ifdef HAVE_LIME
  Grid::emptyUserRecord record;
  Grid::ScidacWriter SW(subspace[0].Grid()->IsBoss());
  SW.open(fname);
  for (int k = 0; k < (int)subspace.size(); k++) SW.writeScidacFieldRecord(subspace[k], record);
  SW.close();
#endif
}
template <class Field>
void loadSubspace(std::vector<Field> &subspace, std::string const fname){
#ifdef HAVE_LIME
  Grid::emptyUserRecord record;
  Grid::ScidacReader SR;
  SR.open(fname);
  for (int k = 0; k < (int)subspace.size(); k++) SR.readScidacFieldRecord(subspace[k], record);
  SR.close();
#endif
}

//////////////////////////////////////////////////////////////////////
// A = PV^dag M (non-Hermitian), and shifted variant for smoothers.
//////////////////////////////////////////////////////////////////////
template<class Matrix,class Field>
class PVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat; Matrix &_PV;
public:
  PVdagMLinearOperator(Matrix &Mat,Matrix &PV): _Mat(Mat),_PV(PV) {};
  void OpDiag (const Field &in, Field &out) { assert(0); }
  void OpDir  (const Field &in, Field &out,int dir,int disp) { assert(0); }
  void OpDirAll  (const Field &in, std::vector<Field> &out){ assert(0); };
  void Op     (const Field &in, Field &out){ Field tmp(in.Grid()); _Mat.M(in,tmp); _PV.Mdag(tmp,out); }
  void AdjOp  (const Field &in, Field &out){ Field tmp(in.Grid()); _PV.M(in,tmp); _Mat.Mdag(tmp,out); }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ HermOp(in,out); ComplexD d=innerProduct(in,out); n1=real(d); n2=norm2(out); }
  void HermOp(const Field &in, Field &out){ Field tmp(in.Grid()); Op(in,tmp); AdjOp(tmp,out); }
};

template<class Matrix,class Field>
class ShiftedPVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat; Matrix &_PV;
public:
  RealD shift;
  ShiftedPVdagMLinearOperator(RealD _shift,Matrix &Mat,Matrix &PV): shift(_shift),_Mat(Mat),_PV(PV){};
  void OpDiag (const Field &in, Field &out) { assert(0); }
  void OpDir  (const Field &in, Field &out,int dir,int disp) { assert(0); }
  void OpDirAll  (const Field &in, std::vector<Field> &out){ assert(0); };
  void Op     (const Field &in, Field &out){ Field tmp(in.Grid()); _Mat.M(in,tmp); _PV.Mdag(tmp,out); out = out + shift*in; }
  void AdjOp  (const Field &in, Field &out){ Field tmp(in.Grid()); _PV.M(tmp,out); _Mat.Mdag(in,tmp); out = out + shift*in; }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ assert(0); }
  void HermOp(const Field &in, Field &out){ Field tmp(in.Grid()); Op(in,tmp); AdjOp(tmp,out); }
};

// Generic shift wrapper (for the coarse-level smoother on the 6D mrhs coarse operator).
template<class Field>
class ShiftedLinearOperator : public LinearOperatorBase<Field> {
  LinearOperatorBase<Field> &_Op; RealD shift;
public:
  ShiftedLinearOperator(RealD _shift, LinearOperatorBase<Field> &Op) : _Op(Op), shift(_shift) {}
  void OpDiag  (const Field &in, Field &out) { assert(0); }
  void OpDir   (const Field &in, Field &out,int dir,int disp) { assert(0); }
  void OpDirAll (const Field &in, std::vector<Field> &out) { assert(0); }
  void Op      (const Field &in, Field &out) { _Op.Op(in,out);    out = out + shift*in; }
  void AdjOp   (const Field &in, Field &out) { _Op.AdjOp(in,out); out = out + shift*in; }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){ assert(0); }
  void HermOp  (const Field &in, Field &out) { Field tmp(in.Grid()); Op(in,tmp); AdjOp(tmp,out); }
};

//////////////////////////////////////////////////////////////////////
// Dense CC solve on the PACKED 6D mrhs coarse-coarse field: drop-in
// for the L3 PGCR, delegating to the library DenseCoarseMatrix.
//////////////////////////////////////////////////////////////////////
template<class DenseType, class CoarseCoarseField>
class MrhsDenseCCSolve : public LinearFunction<CoarseCoarseField> {
public:
  DenseType &_Dense;
  GridBase *_CoarseCoarse5d;
  int _nrhs;
  MrhsDenseCCSolve(DenseType &D, GridBase *cc5d, int nrhs)
    : _Dense(D), _CoarseCoarse5d(cc5d), _nrhs(nrhs) {}
  using LinearFunction<CoarseCoarseField>::operator();
  virtual void operator()(const CoarseCoarseField &in, CoarseCoarseField &out){
    if ( getenv("DENSE_CC_CHECK") ) {
      // Audit path: per-rhs 5D unpack so ApplyBatch can run the _Op defect
      // check per rhs.  ~50ms/call of slice/split overhead -- audit only.
      CoarseCoarseField tmp(in.Grid());
      tmp = in;
      std::vector<CoarseCoarseField> split_in (_nrhs,_CoarseCoarse5d);
      std::vector<CoarseCoarseField> split_out(_nrhs,_CoarseCoarse5d);
      for(int r=0;r<_nrhs;r++) ExtractSliceFast(split_in[r], tmp, r, 0);
      _Dense.ApplyBatch(split_in, split_out);
      for(int r=0;r<_nrhs;r++) InsertSliceFast(split_out[r], out, r, 0);
    } else {
      GRID_TRACE("MrhsDensCCSolve::ApplyBatch6D");
      _Dense.ApplyBatch6D(in, out, _nrhs);
    }
  }
};

//////////////////////////////////////////////////////////////////////
// mrhs interfaces + single-polynomial mrhs PGCR (verbatim from the
// frozen Example_pvdagm_mrhs_3level_dense.cc)
//////////////////////////////////////////////////////////////////////
template<class Field>
class MrhsLinearFunction {
public:
  virtual void operator()(std::vector<Field> &in, std::vector<Field> &out) = 0;
};

template<class Field>
class MrhsPGCRNonHermitian {
public:
  RealD Tolerance; Integer MaxIterations; int mmax,nstep,steps,level;
  int ZeroGuess = 0; int FirstCycle = 0;   // caller contract: zero guess => first-cycle r0 = src
  std::string name = "Level 1";
  LinearOperatorBase<Field> &Linop;
  MrhsLinearFunction<Field> &Preconditioner;
  void Level(int lv){ name = "Level " + std::to_string(lv); level=lv; }
  void Name(std::string n){ name = n; }
  void SetZeroGuess(int z){ ZeroGuess=z; }
  MrhsPGCRNonHermitian(RealD tol,Integer maxit,LinearOperatorBase<Field> &_Linop,MrhsLinearFunction<Field> &Prec,int _mmax,int _nstep)
    : Tolerance(tol),MaxIterations(maxit),Linop(_Linop),Preconditioner(Prec),mmax(_mmax),nstep(_nstep){ level=1; }
  static RealD vnorm2(std::vector<Field> &x){ GRID_TRACE("GCR-vnorm2"); RealD s=0; for(auto &f:x) s+=norm2(f); return s; }
  static ComplexD vinnerProduct(std::vector<Field> &x,std::vector<Field> &y){ GRID_TRACE("GCR-vinnerProduct"); ComplexD s(0); for(int r=0;r<(int)x.size();r++) s+=innerProduct(x[r],y[r]); return s; }
  static void vaxpy(std::vector<Field> &z,ComplexD a,std::vector<Field> &x,std::vector<Field> &y){ GRID_TRACE("GCR-vaxpy"); for(int r=0;r<(int)z.size();r++) axpy(z[r],a,x[r],y[r]); }
  void vOp(std::vector<Field> &in,std::vector<Field> &out){ GRID_TRACE("GCR-vOp"); for(int r=0;r<(int)in.size();r++) Linop.Op(in[r],out[r]); }  
  void operator()(std::vector<Field> &src,std::vector<Field> &psi){
    GRID_TRACE((name+"MrhsPGCRNonHermitian").c_str());
    RealD cp,ssq,rsq; int nrhs=src.size(); GridBase *grid=src[0].Grid();
    ssq=vnorm2(src); rsq=Tolerance*Tolerance*ssq;
    std::vector<Field> r(nrhs,grid);
    GridStopWatch T; T.Start(); steps=0; FirstCycle=1;
    for(int k=0;k<MaxIterations;k++){
      cp=GCRnStep(src,psi,rsq);
      std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name<<" MrhsPGCR("<<mmax<<","<<nstep<<") "<<steps<<" steps cp = "<<cp<<" target "<<rsq<<std::endl;
      if(cp<rsq){
        T.Stop(); vOp(psi,r); for(int rr=0;rr<nrhs;rr++) axpy(r[rr],-1.0,src[rr],r[rr]);
        RealD tr=vnorm2(r);
        std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name<<" MrhsPGCR: Converged on iteration "<<steps
                 <<" computed residual "<<std::sqrt(cp/ssq)<<" true residual "<<std::sqrt(tr/ssq)<<" target "<<Tolerance<<std::endl;
        std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name<<" MrhsPGCR Time elapsed: Total "<<T.Elapsed()<<std::endl;
        for(int rr=0;rr<nrhs;rr++){ RealD rn=std::sqrt(norm2(r[rr])/norm2(src[rr])); std::cout<<GridLogMessage<<"MrhsPGCR per-rhs true residual["<<rr<<"] = "<<rn<<std::endl; }
        return;
      }
    }
    std::cout<<GridLogMessage<<"MrhsPGCR: did not converge"<<std::endl;
  }
  RealD GCRnStep(std::vector<Field> &src,std::vector<Field> &psi,RealD rsq){
    RealD cp; ComplexD a,b,rq; RealD zAAz; int nrhs=src.size(); GridBase *grid=src[0].Grid();
    std::vector<Field> r(nrhs,grid),z(nrhs,grid),Az(nrhs,grid);
    std::vector< std::vector<Field> > q(mmax,std::vector<Field>(nrhs,grid));
    std::vector< std::vector<Field> > p(mmax,std::vector<Field>(nrhs,grid));
    std::vector<RealD> qq(mmax);
    std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name<<" MrhsPGCR nStep("<<nstep<<")"<<std::endl;
    if (ZeroGuess && FirstCycle) { for(int rr=0;rr<nrhs;rr++){ psi[rr]=Zero(); r[rr]=src[rr]; } }
    else                         { vOp(psi,Az); for(int rr=0;rr<nrhs;rr++) r[rr]=src[rr]-Az[rr]; }
    FirstCycle=0;
    Preconditioner(r,z); vOp(z,Az); zAAz=vnorm2(Az);
    p[0]=z; q[0]=Az; qq[0]=zAAz; cp=vnorm2(r);
    for(int k=0;k<nstep;k++){
      steps++; int kp=k+1, peri_k=k%mmax, peri_kp=kp%mmax;
      {
        GRID_TRACE("GCR-project");
        rq=vinnerProduct(q[peri_k],r); a=rq/qq[peri_k];
        vaxpy(psi,a,p[peri_k],psi); vaxpy(r,-a,q[peri_k],r); cp=vnorm2(r);
      }      
      //      rq=vinnerProduct(q[peri_k],r); a=rq/qq[peri_k];
      //      vaxpy(psi,a,p[peri_k],psi); vaxpy(r,-a,q[peri_k],r); cp=vnorm2(r);
      std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name<<" MrhsPGCR step["<<steps<<"]  resid "<<cp<<" target "<<rsq<<std::endl;
      if((k==nstep-1)||(cp<rsq)) return cp;
      Preconditioner(r,z); vOp(z,Az); zAAz=vnorm2(Az);

      {
        GRID_TRACE("GCR-orthog");
        q[peri_kp]=Az; p[peri_kp]=z;
        int northog=((kp)>(mmax-1))?(mmax-1):(kp);
        for(int back=0;back<northog;back++){ int peri_back=(k-back)%mmax; GRID_ASSERT((k-back)>=0);
          b=-real(vinnerProduct(q[peri_back],Az))/qq[peri_back];
          vaxpy(p[peri_kp],b,p[peri_back],p[peri_kp]); vaxpy(q[peri_kp],b,q[peri_back],q[peri_kp]); }
        qq[peri_kp]=vnorm2(q[peri_kp]);
      }      
      //      q[peri_kp]=Az; p[peri_kp]=z;
      //      int northog=((kp)>(mmax-1))?(mmax-1):(kp);
      //      for(int back=0;back<northog;back++){ int peri_back=(k-back)%mmax; GRID_ASSERT((k-back)>=0);
      //        b=-real(vinnerProduct(q[peri_back],Az))/qq[peri_back];
      //        vaxpy(p[peri_kp],b,p[peri_back],p[peri_kp]); vaxpy(q[peri_kp],b,q[peri_back],q[peri_kp]); }
      //      qq[peri_kp]=vnorm2(q[peri_kp]);
    }
    GRID_ASSERT(0); return cp;
  }
};

//////////////////////////////////////////////////////////////////////
// L2->L3 mrhs V-cycle: LinearFunction on the 6D mrhs COARSE field.
// The coarse-coarse solve slot takes EITHER the dense mrhs solve
// (DENSE_CC=1) or the L3 PGCR (DENSE_CC=0).
//////////////////////////////////////////////////////////////////////
template<class CoarseField, class CoarseCoarseField>
class MrhsCoarseThreeLevelPrec : public LinearFunction<CoarseField> {
public:
  LinearOperatorBase<CoarseField>          &_CoarseOp;          // mrhs coarse op (6D)
  LinearFunction<CoarseField>              &_CoarseSmoother;    // shifted 6D coarse smoother
  MultiRHSBlockProject<CoarseField>        &_Projector;         // L2->L3 (vector-based)
  LinearFunction<CoarseCoarseField>        &_CoarseCoarseSolve; // L3 solve (6D cc)
  GridBase *_Coarse5d, *_CoarseCoarse5d, *_CoarseCoarseMrhs;
  int _nrhs;

  MrhsCoarseThreeLevelPrec(LinearOperatorBase<CoarseField> &CoarseOp,
                           LinearFunction<CoarseField> &CoarseSmoother,
                           MultiRHSBlockProject<CoarseField> &Projector,
                           LinearFunction<CoarseCoarseField> &CoarseCoarseSolve,
                           GridBase *Coarse5d, GridBase *CoarseCoarse5d, GridBase *CoarseCoarseMrhs, int nrhs)
    : _CoarseOp(CoarseOp), _CoarseSmoother(CoarseSmoother), _Projector(Projector),
      _CoarseCoarseSolve(CoarseCoarseSolve),
      _Coarse5d(Coarse5d), _CoarseCoarse5d(CoarseCoarse5d), _CoarseCoarseMrhs(CoarseCoarseMrhs), _nrhs(nrhs) {}

  using LinearFunction<CoarseField>::operator();
  virtual void operator()(const CoarseField &in, CoarseField &out) {
    int nrhs=_nrhs; double t;
    CoarseField vec1(in.Grid());
    CoarseField vec2(in.Grid());

    // trivial pre-smoother
    out = in;

    // residual (6D coarse)
    _CoarseOp.Op(out,vec1);  sub(vec1,in,vec1);

    // restrict: unpack 6D coarse -> vector<Coarse> -> blockProject -> vector<CoarseCoarse> -> pack 6D cc
    std::vector<CoarseField>       csplit(nrhs,_Coarse5d);
    std::vector<CoarseCoarseField> ccsplit(nrhs,_CoarseCoarse5d);
    CoarseCoarseField CCsrc(_CoarseCoarseMrhs);
    CoarseCoarseField CCsol(_CoarseCoarseMrhs);

    t=-usecond();
    {
      GRID_TRACE("L2L3-Vcycle - Extract/blockProject/Insert");
      for(int r=0;r<nrhs;r++) ExtractSliceFast(csplit[r],vec1,r,0);
      _Projector.blockProject(csplit,ccsplit);
      for(int r=0;r<nrhs;r++) InsertSliceFast(ccsplit[r],CCsrc,r,0);
    }
    t+=usecond();
    std::cout<<GridLogMessage<<"L2->L3 restrict took "<<t/1000.0<<"ms"<<std::endl;

    // L3 solve (dense mrhs GEMM, or PGCR)
    t=-usecond();
    {
      GRID_TRACE("L2L3-Vcycle - CoarseCoarseSolve");
      CCsol=Zero();
      _CoarseCoarseSolve(CCsrc,CCsol);
    }
    t+=usecond();
    std::cout<<GridLogMessage<<"L3 coarse-coarse solve took "<<t/1000.0<<"ms  ("<<t/1000.0/nrhs<<"ms/rhs)"<<std::endl;

    // prolong: unpack 6D cc -> blockPromote -> pack 6D coarse; add correction
    t=-usecond();
    {
      GRID_TRACE("L2L3-Vcycle - Ext/blockPromote/Ins");
      for(int r=0;r<nrhs;r++) ExtractSliceFast(ccsplit[r],CCsol,r,0);
      _Projector.blockPromote(csplit,ccsplit);
      for(int r=0;r<nrhs;r++) InsertSliceFast(csplit[r],vec1,r,0);
      add(out,out,vec1);
    }
    t+=usecond();
    std::cout<<GridLogMessage<<"L2->L3 prolong took "<<t/1000.0<<"ms"<<std::endl;

    // residual + coarse smoother (6D coarse)
    {
      GRID_TRACE("L2L3-Vcycle - Resid+CoarseSmooth");
      _CoarseOp.Op(out,vec1);  sub(vec1,in,vec1);
      vec2=Zero();
      _CoarseSmoother(vec1,vec2);
      add(out,out,vec2);
    }
  }
};

//////////////////////////////////////////////////////////////////////
// L1->L2 mrhs V-cycle (verbatim from the frozen example)
//////////////////////////////////////////////////////////////////////
template<class FineField, class MrhsCoarseVector, class FineSmoother>
class MrhsTwoLevelMG : public MrhsLinearFunction<FineField> {
public:
  typedef MrhsCoarseVector CoarseVector;
  LinearOperatorBase<FineField>   &_FineOperator;
  FineSmoother                    &_PostSmoother;
  MultiRHSBlockProject<FineField> &_Projector;
  LinearFunction<CoarseVector>    &_CoarseSolve;
  GridBase *_CoarseGrid, *_CoarseGridMrhs;
  MrhsTwoLevelMG(LinearOperatorBase<FineField> &FineOp, FineSmoother &Post,
                 MultiRHSBlockProject<FineField> &Projector, LinearFunction<CoarseVector> &CoarseSolve,
                 GridBase *CoarseGrid, GridBase *CoarseGridMrhs)
    : _FineOperator(FineOp),_PostSmoother(Post),_Projector(Projector),_CoarseSolve(CoarseSolve),
      _CoarseGrid(CoarseGrid),_CoarseGridMrhs(CoarseGridMrhs){}
  virtual void operator()(std::vector<FineField> &in, std::vector<FineField> &out){
    int nrhs=in.size(); GridBase *fgrid=in[0].Grid(); double t;

    std::vector<FineField> vec1(nrhs,fgrid),vec2(nrhs,fgrid);
    {
      GRID_TRACE("L1L2-Vcycle - Resid");
      for(int r=0;r<nrhs;r++) out[r]=in[r];
      for(int r=0;r<nrhs;r++){
	_FineOperator.Op(out[r],vec1[r]);
	sub(vec1[r],in[r],vec1[r]);
      }
    }
    std::vector<CoarseVector> Csrc_split(nrhs,_CoarseGrid), Csol_split(nrhs,_CoarseGrid);
    CoarseVector CsrcMrhs(_CoarseGridMrhs), CsolMrhs(_CoarseGridMrhs);
    t=-usecond();
    {
      GRID_TRACE("L1L2-Vcycle - blockProject");
      _Projector.blockProject(vec1,Csrc_split);
      for(int r=0;r<nrhs;r++) InsertSliceFast(Csrc_split[r],CsrcMrhs,r,0);
    }
    t+=usecond(); std::cout<<GridLogMessage<<"Mrhs project+pack took "<<t/1000.0<<"ms"<<std::endl;
    t=-usecond();
    {
      GRID_TRACE("L1L2-Vcycle - CoarseSolve");
      CsolMrhs=Zero();
      _CoarseSolve(CsrcMrhs,CsolMrhs);
    }
    t+=usecond();
    std::cout<<GridLogMessage<<"Mrhs coarse solve took "<<t/1000.0<<"ms  ("<<t/1000.0/nrhs<<"ms/rhs)"<<std::endl;
    t=-usecond();
    {
      GRID_TRACE("L1L2-Vcycle - ExtractSlice/blockPromote/add delta");
      for(int r=0;r<nrhs;r++) ExtractSliceFast(Csol_split[r],CsolMrhs,r,0);
      _Projector.blockPromote(vec1,Csol_split);
      for(int r=0;r<nrhs;r++) add(out[r],out[r],vec1[r]);
    }
    t+=usecond(); std::cout<<GridLogMessage<<"Mrhs unpack+promote took "<<t/1000.0<<"ms"<<std::endl;
    {
      GRID_TRACE("L1L2-Vcycle - Residuals");
      for(int r=0;r<nrhs;r++){ _FineOperator.Op(out[r],vec1[r]); sub(vec1[r],in[r],vec1[r]); }
    }
    t=-usecond();
    {
      GRID_TRACE("L1L2-Vcycle - Smoothers");
      for(int r=0;r<nrhs;r++){ vec2[r]=Zero(); _PostSmoother(vec1[r],vec2[r]); add(out[r],out[r],vec2[r]); }
    }
    t+=usecond(); std::cout<<GridLogMessage<<"Mrhs post-smooth took "<<t/1000.0<<"ms  ("<<t/1000.0/nrhs<<"ms/rhs)"<<std::endl;
  }
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);
  ParseEnvironment();

  const int Ls=24; RealD M5=1.8, b=1.5, c=0.5;
  const int nbasis=60; const int nrhs=Nrhs;
  GRID_ASSERT(nrhs % vComplex::Nsimd() == 0);

  std::vector<int> lat_size {48,48,48,96};
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(lat_size, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Level 1 blocking (default 2^4)
  Coordinate clatt = lat_size;
  Coordinate Block({2,2,2,2});
  if ( getenv("BLOCK") ){ GridCmdOptionIntVector(std::string(getenv("BLOCK")),Block); GRID_ASSERT(Block.size()==4); }
  for(int d=0;d<4;d++){ GRID_ASSERT(lat_size[d]%Block[d]==0); clatt[d]=lat_size[d]/Block[d]; }
  std::cout << GridLogMessage << "Block  " << Block  << "  coarse lattice       " << clatt << std::endl;

  // Level 2 blocking: SUPERCOARSE default 8,4,3,6 -> CC [3,6,8,8], the dense floor.
  Coordinate cclatt = clatt;
  Coordinate Block2({4,4,3,6});
  if ( getenv("BLOCK2") ){ GridCmdOptionIntVector(std::string(getenv("BLOCK2")),Block2); GRID_ASSERT(Block2.size()==4); }
  for(int d=0;d<4;d++){ GRID_ASSERT(clatt[d]%Block2[d]==0); cclatt[d]=clatt[d]/Block2[d]; }
  std::cout << GridLogMessage << "Block2 " << Block2 << "  coarse-coarse lattice " << cclatt << std::endl;

  GridCartesian *Coarse4d = SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian *Coarse5d = SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);
  GridCartesian *CoarseCoarse4d = SpaceTimeGrid::makeFourDimGrid(cclatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian *CoarseCoarse5d = SpaceTimeGrid::makeFiveDimGrid(1,CoarseCoarse4d);

  // 6D mrhs grids: rhs is dim 0, SIMD across rhs
  Coordinate mpi=GridDefaultMpi();
  Coordinate rhMpi ({1,1,mpi[0],mpi[1],mpi[2],mpi[3]});
  Coordinate rhSimd({vComplex::Nsimd(),1,1,1,1,1});
  Coordinate rhLatt ({nrhs,1,clatt[0], clatt[1], clatt[2], clatt[3]});
  Coordinate rhLatt2({nrhs,1,cclatt[0],cclatt[1],cclatt[2],cclatt[3]});
  GridCartesian *CoarseMrhs       = new GridCartesian(rhLatt, rhSimd,rhMpi);
  GridCartesian *CoarseCoarseMrhs = new GridCartesian(rhLatt2,rhSimd,rhMpi);

  GridParallelRNG RNG5(FGrid); RNG5.SeedFixedIntegers({5,6,7,8});

  LatticeGaugeField Umu(UGrid);
  std::cout << GridLogMessage << "Reading gauge field" << std::endl;
  FieldMetaData header;
  std::string file("/ccs/home/poare/ckpoint_lat.1000");
  NerscIO::readConfiguration(Umu,header,file);

  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);
  MobiusFermionD Dpv (Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,1.0, M5,b,c);

  typedef PVdagMLinearOperator<MobiusFermionD,LatticeFermionD>        PVdagM_t;
  typedef ShiftedPVdagMLinearOperator<MobiusFermionD,LatticeFermionD> ShiftedPVdagM_t;

  // Level 1 tensor types
  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>       LittleDiracOperator;
  typedef MultiGeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>  MrhsLittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector                                CoarseVector;
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis>                  Subspace;

  // Level 2 tensor types (coarsening deepens the nest by one iScalar)
  typedef CoarseVector::vector_object                                      CoarseSiteObj;
  typedef iScalar<vTComplex>                                              vTTComplex;
  typedef GeneralCoarsenedMatrix<CoarseSiteObj,vTTComplex,nbasis>          LittleDiracOperatorL2;
  typedef MultiGeneralCoarsenedMatrix<CoarseSiteObj,vTTComplex,nbasis>     MrhsLittleDiracOperatorL2;
  typedef LittleDiracOperatorL2::CoarseVector                              CoarseCoarseVector;
  typedef Aggregation<CoarseSiteObj,vTTComplex,nbasis>                     SubspaceL2;

  // The library dense bottom over the L2 coarse operator
  typedef DenseCoarseMatrix<vTTComplex,nbasis>                             DenseCC_t;

  PVdagM_t        PVdagM(Ddwf,Dpv);
  ShiftedPVdagM_t ShiftedPVdagM(FineSmootherShift,Ddwf,Dpv);

  NextToNearestStencilGeometry5D geom (Coarse5d);
  NextToNearestStencilGeometry5D geom2(CoarseCoarse5d);

  //////////////////////////////////////////////////////////////////////
  // Subspace: load RAW (no Orthogonalise!), or generate.
  //////////////////////////////////////////////////////////////////////
  std::string subspace_file = "/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/subspace_nb"
                            + std::to_string(nbasis) + ".scidac";
  if ( getenv("SUBSPACE_FILE") ) subspace_file = std::string(getenv("SUBSPACE_FILE"));
  uint64_t file_exists=0;
  if ( UGrid->IsBoss() ){ std::ifstream f(subspace_file); file_exists=f.good()?1:0; }
  UGrid->GlobalSum(file_exists);

  const int cb=0;
  Subspace AggregatesGCR(Coarse5d,FGrid,cb);
  if ( file_exists ){
    std::cout << GridLogMessage << "*** Loading subspace from disk (kept RAW) ***" << std::endl;
    loadSubspace(AggregatesGCR.subspace, subspace_file);
  } else {
    std::cout << GridLogMessage << "*** GCR subspace generation ***" << std::endl;
    AggregatesGCR.CreateSubspaceGCR(RNG5,PVdagM,nbasis);
    saveSubspace(AggregatesGCR.subspace, subspace_file);
  }

  // RAW copy of the fine null vectors BEFORE CoarsenOperator block-orthonormalises in place.
  std::vector<LatticeFermionD> rawNull(nbasis,FGrid);
  for(int k=0;k<nbasis;k++) rawNull[k]=AggregatesGCR.subspace[k];

  //////////////////////////////////////////////////////////////////////
  // Coarsen L1->L2 and L2->L3 with SINGLE-RHS machinery; import to mrhs via
  // CopyMatrix.  The L2 (coarse-coarse) single-RHS operator is HOISTED to
  // main scope: DenseCoarseMatrix imports its stencil and uses its M for
  // certificates, so it must stay alive for the whole run.
  //////////////////////////////////////////////////////////////////////
  MrhsLittleDiracOperator   mrhsLittleDiracOpPV(geom,  CoarseMrhs);
  MrhsLittleDiracOperatorL2 mrhsLittleDiracOpL2(geom2, CoarseCoarseMrhs);
  MultiRHSBlockProject<LatticeFermionD> MrhsProjector;
  MultiRHSBlockProject<CoarseVector>    MrhsProjectorL2;

  LittleDiracOperatorL2 LittleDiracOpL2(geom2,Coarse5d,CoarseCoarse5d);
  NonHermitianLinearOperator<LittleDiracOperatorL2,CoarseCoarseVector> LinOpCC5d(LittleDiracOpL2);

  {
    // --- L1->L2 single-RHS coarse operator (scoped: its padded _A is the memory peak) ---
    LittleDiracOperator LittleDiracOpPV(geom,FGrid,Coarse5d);
    LittleDiracOpPV.CoarsenOperator(PVdagM, AggregatesGCR);   // orthonormalises AggregatesGCR.subspace in place
    mrhsLittleDiracOpPV.CopyMatrix(LittleDiracOpPV);
    MrhsProjector.Allocate(nbasis,FGrid,Coarse5d);
    MrhsProjector.ImportBasis(AggregatesGCR.subspace);
    NonHermitianLinearOperator<LittleDiracOperator,CoarseVector> LinOpCoarse(LittleDiracOpPV);

    // --- psi_coarse = P^dag (RAW fine null) -> Galerkin images, NOT e_k ---
    std::vector<CoarseVector> psi_coarse(nbasis,Coarse5d);
    for(int k=0;k<nbasis;k++) AggregatesGCR.ProjectToSubspace(psi_coarse[k], rawNull[k]);
    rawNull.clear(); rawNull.shrink_to_fit();
    {
      RealD s2=0.0;
      for(int i=0;i<nbasis;i++) for(int j=0;j<nbasis;j++){
        ComplexD sij=TensorRemove(innerProduct(psi_coarse[i],psi_coarse[j]));
        ComplexD d=sij-(i==j?ComplexD(1.0):ComplexD(0.0)); s2+=real(d)*real(d)+imag(d)*imag(d);
      }
      std::cout<<GridLogMessage<<"GUARD: ||<psi_coarse|psi_coarse> - I||_F = "<<std::sqrt(s2)
               <<"   (~0.23 good; ~sqrt(N_coarse)="<<std::sqrt((double)Coarse5d->gSites())<<" = e_k leak)"<<std::endl;
    }

    // --- L2->L3 single-RHS coarsening ---
    SubspaceL2 AggregatesL2(CoarseCoarse5d,Coarse5d,cb);
    for(int k=0;k<nbasis;k++) AggregatesL2.subspace[k]=psi_coarse[k];
    LittleDiracOpL2.CoarsenOperator(LinOpCoarse, AggregatesL2);
    mrhsLittleDiracOpL2.CopyMatrix(LittleDiracOpL2);
    MrhsProjectorL2.Allocate(nbasis,Coarse5d,CoarseCoarse5d);
    MrhsProjectorL2.ImportBasis(AggregatesL2.subspace);

    // --- guard psi_cc ---
    {
      std::vector<CoarseCoarseVector> psi_cc(nbasis,CoarseCoarse5d);
      for(int k=0;k<nbasis;k++) AggregatesL2.ProjectToSubspace(psi_cc[k], psi_coarse[k]);
      RealD s2=0.0;
      for(int i=0;i<nbasis;i++) for(int j=0;j<nbasis;j++){
        ComplexD sij=TensorRemove(innerProduct(psi_cc[i],psi_cc[j]));
        ComplexD d=sij-(i==j?ComplexD(1.0):ComplexD(0.0)); s2+=real(d)*real(d)+imag(d)*imag(d);
      }
      std::cout<<GridLogMessage<<"GUARD: ||<psi_cc|psi_cc> - I||_F = "<<std::sqrt(s2)
               <<"   (~0.23 good; ~sqrt(N_cc)="<<std::sqrt((double)CoarseCoarse5d->gSites())<<" = e_k leak)"<<std::endl;
    }
  } // single-RHS FINE op + padded _A + AggregatesL2 + psi_coarse freed here

  NonHermitianLinearOperator<MrhsLittleDiracOperator,CoarseVector>         mrhsLinOpCoarse(mrhsLittleDiracOpPV);
  NonHermitianLinearOperator<MrhsLittleDiracOperatorL2,CoarseCoarseVector> mrhsLinOpCC(mrhsLittleDiracOpL2);

  //////////////////////////////////////////////////////////////////////
  // DENSE coarse-coarse bottom: the LIBRARY class, constructed AFTER the
  // fine coarsening frees its memory peak.  Imports the stencil of the
  // hoisted single-RHS LittleDiracOpL2 directly (no probing).
  //////////////////////////////////////////////////////////////////////
  std::unique_ptr<DenseCC_t> DenseCC;
  std::unique_ptr<MrhsDenseCCSolve<DenseCC_t,CoarseCoarseVector>> MrhsDenseCC;
  if (UseDenseCC) {
    std::cout << GridLogMessage << "**********************************************" << std::endl;
    std::cout << GridLogMessage << " Dense CC inverse setup (library DenseCoarseMatrix)" << std::endl;
    std::cout << GridLogMessage << "**********************************************" << std::endl;
    DenseCC.reset(new DenseCC_t(CoarseCoarse5d));
    DenseCC->Import(LittleDiracOpL2);
    MrhsDenseCC.reset(new MrhsDenseCCSolve<DenseCC_t,CoarseCoarseVector>(*DenseCC, CoarseCoarse5d, nrhs));
  }

  //////////////////////////////////////////////////////////////////////
  // Solvers, innermost first.
  //////////////////////////////////////////////////////////////////////
  TrivialPrecon<CoarseVector>       simpleC;
  TrivialPrecon<CoarseCoarseVector> simpleCC;
  TrivialPrecon<LatticeFermionD>    simple_fine;

  // L3 (coarse-coarse) iterative solve: PGCR on the 6D cc operator (DENSE_CC=0 branch)
  PrecGeneralisedConjugateResidualNonHermitian<CoarseCoarseVector>
    L3PGCR(L3Tol,L3MaxIt,mrhsLinOpCC,simpleCC,L3Nstep,L3Nstep);
  L3PGCR.Level(3);
  L3PGCR.Name("CCouter");

  LinearFunction<CoarseCoarseVector> *ccSolve;
  if (UseDenseCC) ccSolve = MrhsDenseCC.get();
  else            ccSolve = &L3PGCR;

  // L2 coarse smoother: shifted 6D coarse op, fixed nstep
  ShiftedLinearOperator<CoarseVector> ShiftedMrhsCoarse(CoarseSmootherShift, mrhsLinOpCoarse);
  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>
    CoarseSmootherGCR(0.01,1,ShiftedMrhsCoarse,simpleC,CoarseSmootherNstep,CoarseSmootherNstep);
  CoarseSmootherGCR.Level(2);
  CoarseSmootherGCR.Name("Csmoother");
  CoarseSmootherGCR.SetZeroGuess(1);   // caller zeroes vec2: skip r0 apply every L2 iteration

  // L2->L3 V-cycle preconditioner (operates on 6D coarse field)
  MrhsCoarseThreeLevelPrec<CoarseVector,CoarseCoarseVector>
    L2to3Precon(mrhsLinOpCoarse, CoarseSmootherGCR, MrhsProjectorL2, *ccSolve,
                Coarse5d, CoarseCoarse5d, CoarseCoarseMrhs, nrhs);

  // L2 coarse solve: PGCR on 6D coarse op, preconditioned by the L2->L3 V-cycle
  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>
    L2PGCR(CoarseSolverTol, CoarseSolverOrder/16, mrhsLinOpCoarse, L2to3Precon, 16, 16);
  L2PGCR.Level(2);
  L2PGCR.Name("Couter");
  L2PGCR.SetZeroGuess(1);              // caller zeroes CsolMrhs; restarts still recompute r

  // Fine smoother (per-rhs, looped in the L1->L2 V-cycle)
  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD>
    SmootherGCR(0.0,1,ShiftedPVdagM,simple_fine,FineSmootherOrder,FineSmootherOrder);
  SmootherGCR.Level(1);
  SmootherGCR.Name("Fsmoother");
  SmootherGCR.SetZeroGuess(1);         // caller zeroes vec2[r]: saves 12 fine mults/outer

  // L1->L2 V-cycle (fine); its coarse solve is the three-level L2PGCR
  typedef PrecGeneralisedConjugateResidualNonHermitian<LatticeFermionD> FineSmoother_t;
  MrhsTwoLevelMG<LatticeFermionD,CoarseVector,FineSmoother_t>
    ThreeLevelPrecon(PVdagM, SmootherGCR, MrhsProjector, L2PGCR, Coarse5d, CoarseMrhs);

  // Outer mrhs solve
  MrhsPGCRNonHermitian<LatticeFermionD>
    L1PGCR(OuterTol,1000,PVdagM,ThreeLevelPrecon,OuterMmax,OuterNstep);
  L1PGCR.Level(1);
  L1PGCR.Name("Fouter");
  L1PGCR.SetZeroGuess(1);              // sol[r]=Zero() below; restarts recompute r as always

  //////////////////////////////////////////////////////////////////////
  // Sources and solve
  //////////////////////////////////////////////////////////////////////
  std::vector<LatticeFermionD> src(nrhs,FGrid), sol(nrhs,FGrid);
  for(int r=0;r<nrhs;r++){ gaussian(RNG5,src[r]); sol[r]=Zero(); }

  std::cout << GridLogMessage << "**********************************************" << std::endl;
  std::cout << GridLogMessage << " MultiRHS THREE-level solve (DenseCoarseMatrix bottom): " << nrhs << " RHS " << std::endl;
  std::cout << GridLogMessage << "**********************************************" << std::endl;

  GridStopWatch w; w.Start();
  L1PGCR(src,sol);
  w.Stop();
  std::cout << GridLogMessage << "MultiRHS 3-level dense solve total " << w.Elapsed()
            << "  (per RHS: " << w.useconds()/1.0e6/nrhs << " s)" << std::endl;

  { LatticeFermionD Ax(FGrid); RealD worst=0.0;
    for(int r=0;r<nrhs;r++){ PVdagM.Op(sol[r],Ax); Ax=Ax-src[r];
      RealD rn=std::sqrt(norm2(Ax)/norm2(src[r]));
      std::cout << GridLogMessage << "FINAL: rhs["<<r<<"] true residual = " << rn << std::endl;
      worst=std::max(worst,rn); }
    std::cout << GridLogMessage << "FINAL: worst-case residual = " << worst << std::endl;
  }

  std::cout << GridLogMessage << "Done" << std::endl;
  Grid_finalize();
  return 0;
}
