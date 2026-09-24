/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/multigrid/MrhsMultiGrid.h

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See the full license in the file "LICENSE" in the top level distribution
    directory
*************************************************************************************/
/*  END LEGAL */
#pragma once

#include <Grid/algorithms/deflation/MultiRHSBlockProject.h>
#include <Grid/algorithms/multigrid/MrhsPreconditioner.h>

NAMESPACE_BEGIN(Grid);


//////////////////////////////////////////////////////////////////////
// Single-polynomial mrhs PGCR: one step length and one set of
// orthogonalisation coefficients SHARED across the right-hand sides
// (vnorm2/vinnerProduct sum over rhs).  For DWF the dense spectrum
// makes sharing near-free and the summed reductions amortise; see the
// mrhs-HDCG paper (arXiv:2409.03904).  Orthogonalisation is classical
// Gram-Schmidt, batched per rhs (rankInnerProductMulti/axpyMulti) with
// ONE GlobalSumVector for all coefficients.
//
// OnStep is a programmatic per-step hook for algorithmic studies
// (e.g. driving the GCR coefficient recorder); it is not a consumer
// parameter.
//////////////////////////////////////////////////////////////////////
template<class Field>
class MrhsPGCRNonHermitian {
public:
  RealD Tolerance; Integer MaxIterations; int mmax,nstep,steps,level;
  int ZeroGuess = 0; int FirstCycle = 0;
  std::string name = "Level 1";
  // Trace range names carry the instance name, so a profile separates the
  // solvers that otherwise nest as one shared range.
  std::string trace_op = "MrhsPGCR::vOp", trace_orthog = "MrhsPGCR orthog";
  LinearOperatorBase<Field> &Linop;
  MrhsLinearFunction<Field> &Preconditioner;
  std::function<void(int)> OnStep;      // called with the outer step count after every step
  void Level(int lv){ name = "Level " + std::to_string(lv); level=lv; }
  void Name(std::string n){ name = n; trace_op = name+" MrhsPGCR::vOp"; trace_orthog = name+" MrhsPGCR orthog"; }
  void SetZeroGuess(int z){ ZeroGuess=z; }
  MrhsPGCRNonHermitian(RealD tol,Integer maxit,LinearOperatorBase<Field> &_Linop,MrhsLinearFunction<Field> &Prec,int _mmax,int _nstep)
    : Tolerance(tol),MaxIterations(maxit),Linop(_Linop),Preconditioner(Prec),mmax(_mmax),nstep(_nstep){ level=1; }
  static RealD vnorm2(std::vector<Field> &x){ RealD s=0; for(auto &f:x) s+=norm2(f); return s; }
  static ComplexD vinnerProduct(std::vector<Field> &x,std::vector<Field> &y){ ComplexD s(0); for(int r=0;r<(int)x.size();r++) s+=innerProduct(x[r],y[r]); return s; }
  static void vaxpy(std::vector<Field> &z,ComplexD a,std::vector<Field> &x,std::vector<Field> &y){ for(int r=0;r<(int)z.size();r++) axpy(z[r],a,x[r],y[r]); }
  void vOp(std::vector<Field> &in,std::vector<Field> &out){ GRID_TRACE(trace_op.c_str()); for(int r=0;r<(int)in.size();r++) Linop.Op(in[r],out[r]); }
  void operator()(std::vector<Field> &src,std::vector<Field> &psi){
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
        return;
      }
    }
    std::cout<<GridLogMessage<<"MrhsPGCR: did not converge"<<std::endl;
  }
  RealD GCRnStep(std::vector<Field> &src,std::vector<Field> &psi,RealD rsq){
    RealD cp; ComplexD a,rq; int nrhs=src.size(); GridBase *grid=src[0].Grid();
    std::vector<Field> r(nrhs,grid),Az(nrhs,grid);   // Az: restart residual scratch only
    std::vector< std::vector<Field> > q(mmax,std::vector<Field>(nrhs,grid));
    std::vector< std::vector<Field> > p(mmax,std::vector<Field>(nrhs,grid));
    std::vector<RealD> qq(mmax);
    if (ZeroGuess && FirstCycle) { for(int rr=0;rr<nrhs;rr++){ psi[rr]=Zero(); r[rr]=src[rr]; } }
    else                         { vOp(psi,Az); for(int rr=0;rr<nrhs;rr++) r[rr]=src[rr]-Az[rr]; }
    FirstCycle=0;
    // p[0]=Prec(r), q[0]=A p[0], produced directly in the history slots (no copies)
    Preconditioner(r,p[0]); vOp(p[0],q[0]); qq[0]=vnorm2(q[0]); cp=vnorm2(r);
    for(int k=0;k<nstep;k++){
      steps++; int kp=k+1, peri_k=k%mmax, peri_kp=kp%mmax;
      if ( OnStep ) OnStep(steps);
      rq=vinnerProduct(q[peri_k],r); a=rq/qq[peri_k];
      vaxpy(psi,a,p[peri_k],psi); vaxpy(r,-a,q[peri_k],r); cp=vnorm2(r);
      std::cout<<GridLogMessage<<std::string(level,'\t')<<" "<<name<<" MrhsPGCR step["<<steps<<"]  resid "<<cp<<" target "<<rsq<<std::endl;
      if((k==nstep-1)||(cp<rsq)) return cp;
      // New direction straight into its history slot: p=Prec(r), q=A p.
      Preconditioner(r,p[peri_kp]);
      vOp(p[peri_kp],q[peri_kp]);
      int northog=((kp)>(mmax-1))?(mmax-1):(kp);
      {
        GRID_TRACE(trace_orthog.c_str());
        // Classical Gram-Schmidt: all coefficients against the UN-updated new q
        // (independent, batchable), then apply.  Complex coefficient: the
        // operator is non-Hermitian, real(<q_j,Aq>) alone left q's non-orthogonal.
        // Batched per rhs (one fused kernel + one reduction each), the shared
        // coefficient summed over rhs on the host, ONE GlobalSumVector.
        std::vector<ComplexD> bcoef(northog,ComplexD(0.0)), part;
        for(int rr=0;rr<nrhs;rr++){
          std::vector<const Field*> qwin(northog);
          for(int back=0;back<northog;back++){ int peri_back=(k-back)%mmax; GRID_ASSERT((k-back)>=0); qwin[back]=&q[peri_back][rr]; }
          rankInnerProductMulti(part,qwin,q[peri_kp][rr]);
          for(int back=0;back<northog;back++) bcoef[back]+=part[back];
        }
        if(northog) grid->GlobalSumVector(&bcoef[0],northog);
        for(int back=0;back<northog;back++){ int peri_back=(k-back)%mmax; bcoef[back]=-bcoef[back]/qq[peri_back]; }
        for(int rr=0;rr<nrhs;rr++){
          std::vector<const Field*> qwin(northog), pwin(northog);
          for(int back=0;back<northog;back++){ int peri_back=(k-back)%mmax; qwin[back]=&q[peri_back][rr]; pwin[back]=&p[peri_back][rr]; }
          axpyMulti(p[peri_kp][rr],bcoef,pwin);
          axpyMulti(q[peri_kp][rr],bcoef,qwin);
        }
      }
      qq[peri_kp]=vnorm2(q[peri_kp]);
    }
    GRID_ASSERT(0); return cp;
  }
};

//////////////////////////////////////////////////////////////////////
// Dense L3 solve on the packed D+1 coarse-coarse field.  Adapts any
// class exposing ApplyBatch6D (DenseCoarseMatrix) to a LinearFunction.
//////////////////////////////////////////////////////////////////////
template<class DenseType, class CoarseCoarseField>
class MrhsDenseCCSolve : public LinearFunction<CoarseCoarseField> {
public:
  DenseType &_Dense;
  int _nrhs;
  MrhsDenseCCSolve(DenseType &D, int nrhs) : _Dense(D), _nrhs(nrhs) {}
  using LinearFunction<CoarseCoarseField>::operator();
  virtual void operator()(const CoarseCoarseField &in, CoarseCoarseField &out){
    _Dense.ApplyBatch6D(in, out, _nrhs);
  }
};

//////////////////////////////////////////////////////////////////////
// L2->L3 half V-cycle on the D+1 coarse field: coarse-coarse correction
// through the mixed blockProject, then post-smooth.
//////////////////////////////////////////////////////////////////////
template<class CoarseField, class CoarseCoarseField>
class MrhsCoarseThreeLevelPrec : public LinearFunction<CoarseField> {
public:
  LinearOperatorBase<CoarseField>          &_CoarseOp;
  LinearFunction<CoarseField>              &_CoarseSmoother;
  MultiRHSBlockProject<CoarseField>        &_Projector;
  LinearFunction<CoarseCoarseField>        &_CoarseCoarseSolve;
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
    CoarseField vec1(in.Grid());
    CoarseField vec2(in.Grid());
    out = in;
    _CoarseOp.Op(out,vec1);  sub(vec1,in,vec1);

    // restrict, through the mixed blockProject: D+1 coarse in, D+1 cc out
    CoarseCoarseField CCsrc(_CoarseCoarseMrhs);
    CoarseCoarseField CCsol(_CoarseCoarseMrhs);
    _Projector.blockProject(vec1,CCsrc);

    _CoarseCoarseSolve(CCsrc,CCsol);

    _Projector.blockPromote(vec1,CCsol);
    add(out,out,vec1);

    _CoarseOp.Op(out,vec1);  sub(vec1,in,vec1);
    _CoarseSmoother(vec1,vec2);
    add(out,out,vec2);
  }
};

//////////////////////////////////////////////////////////////////////
// L1->L2 mrhs V-cycle.  The whole V-cycle is preconditioner: its fine
// residuals and the smoother may run with sloppy halos; the caller
// (the outer Krylov) gets the exact operator back on exit.  SetSloppy
// is wired by the composer to PVdagMLinearOperator::SloppyComms (a
// no-op by default), replacing the file-scope global the example used.
//////////////////////////////////////////////////////////////////////
// Projector_t defaults to a transfer operator whose STORE matches FineField,
// which is the case whenever the coarse sector and the fine level share a
// precision; pass it explicitly when they differ.
template<class FineField, class MrhsCoarseVector, class FineSmoother,
         class Projector_t = MultiRHSBlockProject<FineField> >
class MrhsTwoLevelMG : public MrhsPreconditioner<FineField> {
public:
  typedef MrhsCoarseVector CoarseVector;
  LinearOperatorBase<FineField>   &_FineOperator;
  FineSmoother                    &_PostSmoother;
  Projector_t                     &_Projector;   // store precision need not match FineField
  LinearFunction<CoarseVector>    &_CoarseSolve;
  GridBase *_CoarseGrid, *_CoarseGridMrhs;
  std::function<void(int)> SetSloppy = [](int){};
  int SloppyComms = 0;                 // value passed to SetSloppy on entry
  MrhsTwoLevelMG(LinearOperatorBase<FineField> &FineOp, FineSmoother &Post,
                 Projector_t &Projector, LinearFunction<CoarseVector> &CoarseSolve,
                 GridBase *CoarseGrid, GridBase *CoarseGridMrhs)
    : _FineOperator(FineOp),_PostSmoother(Post),_Projector(Projector),_CoarseSolve(CoarseSolve),
      _CoarseGrid(CoarseGrid),_CoarseGridMrhs(CoarseGridMrhs){}
  virtual void operator()(std::vector<FineField> &in, std::vector<FineField> &out){
    GRID_TRACE("MGVcycle");
    SetSloppy(SloppyComms);
    int nrhs=in.size(); GridBase *fgrid=in[0].Grid();
    std::vector<FineField> vec1(nrhs,fgrid),vec2(nrhs,fgrid);
    for(int r=0;r<nrhs;r++) out[r]=in[r];
    { GRID_TRACE("MGFineResidual");
      for(int r=0;r<nrhs;r++){ _FineOperator.Op(out[r],vec1[r]); sub(vec1[r],in[r],vec1[r]); }
    }
    // fine vector -> D+1 coarse, via the mixed blockProject
    CoarseVector CsrcMrhs(_CoarseGridMrhs), CsolMrhs(_CoarseGridMrhs);
    { GRID_TRACE("MGProject");
      _Projector.blockProject(vec1,CsrcMrhs);
    }
    CsolMrhs=Zero();
    { GRID_TRACE("MGCoarseSolve");
      _CoarseSolve(CsrcMrhs,CsolMrhs);
    }
    { GRID_TRACE("MGPromote");
      _Projector.blockPromote(vec1,CsolMrhs);
      for(int r=0;r<nrhs;r++) add(out[r],out[r],vec1[r]);
    }
    { GRID_TRACE("MGFineResidual2");
      for(int r=0;r<nrhs;r++){ _FineOperator.Op(out[r],vec1[r]); sub(vec1[r],in[r],vec1[r]); }
    }
    { GRID_TRACE("MGPostSmooth");
      for(int r=0;r<nrhs;r++){
        _PostSmoother(vec1[r],vec2[r]); add(out[r],out[r],vec2[r]);
      }
    }
    SetSloppy(0);
  }
};

//////////////////////////////////////////////////////////////////////
// The fp64/fp32 seam of the solve chain.  The outer Krylov hands fp64
// residuals to its preconditioner; this adapter converts them to fp32,
// runs an fp32 preconditioner (the whole V-cycle), and converts the
// correction back.  Two precisionChange per outer step per rhs, on
// workspaces built once.  The outer operator never sees fp32.
//
// The fp32 scratch is allocated ON FIRST USE, not in the constructor: the
// solver builds this seam whichever precision is selected, and at Nrhs 12
// on a 48^3x96 Ls=24 rank the scratch is 2 GB of device memory.  A seam
// that is never called must cost nothing.  Same idiom as the GCR history
// vectors: re-made only if the grid or the rhs count changes.
//////////////////////////////////////////////////////////////////////
template<class FieldD, class FieldF>
class MrhsMixedPrecPreconditioner : public MrhsPreconditioner<FieldD> {
public:
  MrhsPreconditioner<FieldF> &_Inner;
  precisionChangeWorkspace    _ws_d2f;   // out fp32, in fp64
  precisionChangeWorkspace    _ws_f2d;   // out fp64, in fp32
  GridBase                   *_gridF;
  std::vector<FieldF>         _in_f, _out_f;
  MrhsMixedPrecPreconditioner(MrhsPreconditioner<FieldF> &Inner, GridBase *gridD, GridBase *gridF, int nrhs)
    : _Inner(Inner), _ws_d2f(gridF,gridD), _ws_f2d(gridD,gridF), _gridF(gridF) {}
  void Scratch(int nrhs){
    if ( (int)_in_f.size() == nrhs ) return;
    _in_f.clear();  _in_f.reserve(nrhs);
    _out_f.clear(); _out_f.reserve(nrhs);
    for(int r=0;r<nrhs;r++){ _in_f.emplace_back(_gridF); _out_f.emplace_back(_gridF); }
  }
  virtual void operator()(std::vector<FieldD> &in, std::vector<FieldD> &out){
    GRID_TRACE("MGPrecisionSeam");
    int nrhs=in.size(); Scratch(nrhs);
    for(int r=0;r<nrhs;r++) precisionChange(_in_f[r],in[r],_ws_d2f);
    _Inner(_in_f,_out_f);
    for(int r=0;r<nrhs;r++) precisionChange(out[r],_out_f[r],_ws_f2d);
  }
  // The start and the timers are the inner preconditioner's
  virtual void Vstart(std::vector<FieldD> &x, std::vector<FieldD> &src){
    GRID_TRACE("MGPrecisionSeamVstart");
    int nrhs=src.size(); Scratch(nrhs);
    for(int r=0;r<nrhs;r++) precisionChange(_in_f[r],src[r],_ws_d2f);
    _Inner.Vstart(_out_f,_in_f);
    for(int r=0;r<nrhs;r++) precisionChange(x[r],_out_f[r],_ws_f2d);
  }
  virtual void ResetTimers(void)                       { _Inner.ResetTimers(); }
  virtual void ReportTimers(const std::string &prefix) { _Inner.ReportTimers(prefix); }
};

NAMESPACE_END(Grid);
