/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/multigrid/Smoothers.h

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
#include <Grid/algorithms/iterative/GCRCoefficients.h>

NAMESPACE_BEGIN(Grid);

//////////////////////////////////////////////////////////////////////////////
// Multigrid smoothers as LinearFunction<Field> objects.
//
// Gathered here from ~15 copies in tests/solver and tests/debug (HDCR /
// HDCG era), plus the fixed-polynomial smoothers of 2026:
//
//   ChebyshevSmoother          Chebyshev approx to 1/x on [lo,hi], applied
//                              through HermOp (the original HDCG smoother on
//                              a Hermitian shifted operator).
//   ChebyshevNonHermitianSmoother  same polynomial applied through Op(), for
//                              a non-Hermitian (near-normal, real-spectrum)
//                              smoother operator such as shifted PVdagM.
//   ChebyshevInverter          one Chebyshev-corrected step with residual print.
//   MirsSmoother               shifted-MdagM CG, HDCG arXiv:1402.2585.
//   GCRReplaySmoother          replays a GCR's recorded step lengths a_k and
//                              orthogonalisation coefficients b_kj with NO
//                              inner products: one matvec per step, zero
//                              reductions.  The "PreconditionerMirsPoly"
//                              idea of 1402.2585 p.13 applied to GCR.
//
// Recording: PrecGeneralisedConjugateResidualNonHermitian::SetCoefficientRecorder
// accumulates per-step means over calls into a GCRCoefficients; construct a
// GCRReplaySmoother from it.
//////////////////////////////////////////////////////////////////////////////

inline RealD InverseApproximation(RealD x){ return 1.0/x; }

//////////////////////////////////////////////////////////////////////////////
// HermOp-based Chebyshev smoother.  The second template parameter and the
// 5-argument constructor exist only so the historical call sites
//   ChebyshevSmoother<LatticeFermion,DomainWallFermionD> S(lo,hi,ord,HermOp,Ddwf);
// compile unchanged; the Matrix argument was never used.
//////////////////////////////////////////////////////////////////////////////
template<class Field,class Matrix=void> class ChebyshevSmoother : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  typedef LinearOperatorBase<Field> FineOperator;
  FineOperator   & _SmootherOperator;
  Chebyshev<Field> Cheby;
  ChebyshevSmoother(RealD _lo,RealD _hi,int _ord, FineOperator &SmootherOperator) :
    _SmootherOperator(SmootherOperator),
    Cheby(_lo,_hi,_ord,InverseApproximation)
  {
    std::cout << GridLogMessage<<" Chebyshev smoother order "<<_ord<<" ["<<_lo<<","<<_hi<<"]"<<std::endl;
  };
  template<class M>
  ChebyshevSmoother(RealD _lo,RealD _hi,int _ord, FineOperator &SmootherOperator, M &) :
    ChebyshevSmoother(_lo,_hi,_ord,SmootherOperator) {};
  void operator() (const Field &in, Field &out)
  {
    Cheby(_SmootherOperator,in,out);
  }
};

//////////////////////////////////////////////////////////////////////////////
// Op()-based Chebyshev smoother: x = S(A) r with S the Chebyshev fit to 1/x
// on [lo,hi].  Same three-term recurrence as Chebyshev<Field>::operator()
// but through Op, for the non-Hermitian smoother operators of the PVdagM
// multigrid (real coefficients / near-normal, as the recorded GCR
// coefficients show).  `order` matvecs, no reductions.
//////////////////////////////////////////////////////////////////////////////
template<class Field> class ChebyshevNonHermitianSmoother : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  LinearOperatorBase<Field> &Linop;
  RealD lo, hi; int order;
  std::vector<RealD> Coeffs;
  int Verbose = 0;
  std::string name = "cheb";
  ChebyshevNonHermitianSmoother(RealD _lo,RealD _hi,int _order,LinearOperatorBase<Field> &Op)
    : Linop(Op), lo(_lo), hi(_hi), order(_order)
  {
    GRID_ASSERT(order>=2);
    Coeffs.resize(order);
    for(int j=0;j<order;j++){
      RealD s=0;
      for(int k=0;k<order;k++){
        RealD y=std::cos(M_PI*(k+0.5)/order);
        RealD x=0.5*(y*(hi-lo)+(hi+lo));
        s=s+InverseApproximation(x)*std::cos( j*M_PI*(k+0.5)/order );
      }
      Coeffs[j] = s * 2.0/order;
    }
    std::cout << GridLogMessage<<" ChebyshevNonHermitian smoother order "<<order<<" ["<<lo<<","<<hi<<"]"<<std::endl;
  }
  void operator() (const Field &in, Field &out)
  {
    GRID_TRACE("ChebyshevSmoother");
    GridBase *grid=in.Grid();
    Field T0(grid); T0 = in;
    Field T1(grid), T2(grid), y(grid);
    Field *Tnm=&T0, *Tn=&T1, *Tnp=&T2;
    RealD xscale = 2.0/(hi-lo);
    RealD mscale = -(hi+lo)/(hi-lo);
    Linop.Op(T0,y);
    axpby(T1,xscale,mscale,y,in);
    axpby(out,0.5*Coeffs[0],Coeffs[1],T0,T1);
    for(int n=2;n<order;n++){
      Linop.Op(*Tn,y);
      axpby(y,xscale,mscale,y,(*Tn));
      axpby(*Tnp,2.0,-1.0,y,(*Tnm));
      if ( Coeffs[n] != 0.0 ) axpy(out,Coeffs[n],*Tnp,out);
      Field *swizzle=Tnm; Tnm=Tn; Tn=Tnp; Tnp=swizzle;
    }
    if ( Verbose ) { Linop.Op(out,y); y = y - in; std::cout << GridLogMessage << " " << name << " cheb |r|/|r0| = " << std::sqrt(norm2(y)/norm2(in)) << std::endl; }
  }
};

template<class Field> class ChebyshevInverter : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  typedef LinearOperatorBase<Field> FineOperator;
  FineOperator   & _Operator;
  Chebyshev<Field> Cheby;
  ChebyshevInverter(RealD _lo,RealD _hi,int _ord, FineOperator &Operator) :
    _Operator(Operator),
    Cheby(_lo,_hi,_ord,InverseApproximation)
  {
    std::cout << GridLogMessage<<" Chebyshev Inverter order "<<_ord<<" ["<<_lo<<","<<_hi<<"]"<<std::endl;
  };
  void operator() (const Field &in, Field &out)
  {
    Field r(in.Grid());
    Field AinvR(in.Grid());
    _Operator.HermOp(out,r);
    r = in - r; // b - A x
    Cheby(_Operator,r,AinvR); // A^{-1} ( b - A x ) ~ A^{-1} b - x
    out = out + AinvR;
    _Operator.HermOp(out,r);
    r = in - r; // b - A x
    RealD rr = norm2(r);
    RealD ss = norm2(in);
    std::cout << GridLogMessage << "ChebshevInverse resid " <<::sqrt(rr/ss)<<std::endl;
  }
};

//////////////////////////////////////////////////////////////////////////////
// MIRS: CG on the infra-red shifted MdagM (HDCG, arXiv:1402.2585).
//////////////////////////////////////////////////////////////////////////////
template<class Field,class Matrix> class MirsSmoother : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  typedef LinearOperatorBase<Field>                            FineOperator;
  Matrix         & SmootherMatrix;
  FineOperator   & SmootherOperator;
  RealD tol;
  RealD shift;
  int   maxit;
  MirsSmoother(RealD _shift,RealD _tol,int _maxit,FineOperator &_SmootherOperator,Matrix &_SmootherMatrix) :
    shift(_shift),tol(_tol),maxit(_maxit),
    SmootherOperator(_SmootherOperator),
    SmootherMatrix(_SmootherMatrix)
  {};
  void operator() (const Field &in, Field &out)
  {
    ZeroGuesser<Field> Guess;
    ConjugateGradient<Field>  CG(tol,maxit,false);
    Field src(in.Grid());
    ShiftedMdagMLinearOperator<SparseMatrixBase<Field>,Field> MdagMOp(SmootherMatrix,shift);
    SmootherOperator.AdjOp(in,src);
    Guess(src,out);
    CG(MdagMOp,src,out);
  }
};

//////////////////////////////////////////////////////////////////////////////
// Replay of a recorded GCR with a trivial preconditioner:
//   p_0 = r_0 ;  x_{k+1} = x_k + a_k p_k ;  r_{k+1} = r_k - a_k A p_k ;
//   p_{k+1} = r_{k+1} + sum_{j<northog(k)} b_kj p_{k-j}
// One matvec per step, no reductions, history of mmax p vectors held
// persistently (allocated on first use, per grid).
//////////////////////////////////////////////////////////////////////////////
template<class Field> class GCRReplaySmoother : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  LinearOperatorBase<Field> &Linop;
  int mmax, nstep;
  std::vector<ComplexD>               a;
  std::vector<std::vector<ComplexD> > b;
  GridBase *hist_grid = nullptr;
  std::vector<Field> p;
  int Verbose = 0;                       // 1: print |r_m|/|r_0| per call (one extra reduction)
  std::string name = "replay";
  GCRReplaySmoother(LinearOperatorBase<Field> &Op, const GCRCoefficients &c) : Linop(Op)
  {
    mmax  = c.mmax;  GRID_ASSERT(mmax>=1);
    nstep = c.Steps(); GRID_ASSERT(nstep>=1);
    a.resize(nstep); b.resize(nstep);
    // No std::abs / std::isfinite on ComplexD: it is thrust::complex on HIP
    // builds and neither overload exists.  Work on the real and imaginary
    // parts explicitly.
    auto finite = [](ComplexD z){ RealD x=real(z), y=imag(z); return (x==x) && (y==y) && (x-x==0.0) && (y-y==0.0); };
    auto cabs   = [](ComplexD z){ RealD x=real(z), y=imag(z); return std::sqrt(x*x+y*y); };
    RealD amax=0.0, bmax=0.0;
    for(int k=0;k<nstep;k++){
      a[k] = c.A(k);
      GRID_ASSERT( finite(a[k]) );
      amax = std::max(amax, cabs(a[k]));
      b[k].resize(c.NB(k));
      GRID_ASSERT( c.NB(k) <= mmax-1 );
      for(int j=0;j<c.NB(k);j++){
        b[k][j] = c.B(k,j);
        GRID_ASSERT( finite(b[k][j]) );
        bmax = std::max(bmax, cabs(b[k][j]));
      }
    }
    std::cout << GridLogMessage << " GCRReplaySmoother: max|a| " << amax << " max|b| " << bmax << std::endl;
    std::cout << GridLogMessage << " GCRReplaySmoother: " << nstep << " steps, mmax " << mmax
              << ", from " << c.Calls() << " recorded calls" << std::endl;
  }
  void operator() (const Field &src, Field &psi)
  {
    GRID_TRACE("GCRReplaySmoother");
    GridBase *grid = src.Grid();
    if ( hist_grid != grid ) {
      p.clear(); p.reserve(mmax);
      for(int i=0;i<mmax;i++) p.emplace_back(grid);
      hist_grid = grid;
    }
    Field r(grid), q(grid);
    r = src;
    psi = Zero();
    p[0] = r;
    RealD r0 = Verbose ? norm2(src) : 0.0;
    for(int k=0;k<nstep;k++){
      int kp=k+1, peri_k=k%mmax, peri_kp=kp%mmax;
      Linop.Op(p[peri_k],q);                 // q_k = A p_k
      axpy(psi,  a[k], p[peri_k], psi);
      if ( k==nstep-1 ) {
        if ( Verbose ) { axpy(r,-a[k],q,r); std::cout << GridLogMessage << " " << name << " replay |r|/|r0| = " << std::sqrt(norm2(r)/r0) << std::endl; }
        break;
      }
      axpy(r,   -a[k], q,         r);
      p[peri_kp] = r;
      for(int j=0;j<(int)b[k].size();j++){
        int peri_back=(k-j)%mmax;
        axpy(p[peri_kp], b[k][j], p[peri_back], p[peri_kp]);
      }
    }
  }
};

//////////////////////////////////////////////////////////////////////////////
// A LinearFunction that forwards to a replaceable target: lets a V-cycle be
// built once around a smoother slot whose implementation is swapped at run
// time (record with the adaptive GCR, then replay the polynomial).
//////////////////////////////////////////////////////////////////////////////
template<class Field> class SwitchableSmoother : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  LinearFunction<Field> *current;
  std::string label;
  SwitchableSmoother(LinearFunction<Field> &initial, const std::string &l="initial") : current(&initial), label(l) {}
  void Set(LinearFunction<Field> &f, const std::string &l)
  {
    current = &f; label = l;
    std::cout << GridLogMessage << " SwitchableSmoother -> " << l << std::endl;
  }
  void operator() (const Field &in, Field &out) { (*current)(in,out); }
};

NAMESPACE_END(Grid);
