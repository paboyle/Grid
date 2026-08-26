/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/PrecGeneralisedConjugateResidual.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
#ifndef GRID_PREC_GCR_NON_HERM_H
#define GRID_PREC_GCR_NON_HERM_H
#include <Grid/algorithms/iterative/GCRCoefficients.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////
//VPGCR Abe and Zhang, 2005.
//INTERNATIONAL JOURNAL OF NUMERICAL ANALYSIS AND MODELING
//Computing and Information Volume 2, Number 2, Pages 147-161
//NB. Likely not original reference since they are focussing on a preconditioner variant.
//    but VPGCR was nicely written up in their paper
///////////////////////////////////////////////////////////////////////////////////////////////////////
NAMESPACE_BEGIN(Grid);

#define GCRLogLevel std::cout << GridLogMessage <<std::string(level,'\t')<< name<<" "

template<class Field>
class PrecGeneralisedConjugateResidualNonHermitian : public LinearFunction<Field> {
public:                                                
  using LinearFunction<Field>::operator();
  RealD   Tolerance;
  RealD   SSQ;
  Integer MaxIterations;
  int verbose;
  int mmax;
  int nstep;
  int steps;
  int level;
  GridStopWatch PrecTimer;
  GridStopWatch MatTimer;
  GridStopWatch LinalgTimer;
  std::string name;
  int ZeroGuess  = 0;  // caller contract: guess is always zero => first-cycle r0 = src, skip the apply
  // persistent GCR history (see GCRnStep)
  GridBase           *hist_grid = nullptr;
  std::vector<Field>  q;
  std::vector<Field>  p;
  std::vector<RealD>  qq;
  int FirstCycle = 0;

  LinearFunction<Field>     &Preconditioner;
  LinearOperatorBase<Field> &Linop;

  void Name(std::string _name) { name = _name; };

  void Level(int n) { Name("Level " + std::to_string(n)); level = n; }

  void SetZeroGuess(int z) { ZeroGuess = z; };
  // Coefficient logging: one line per step with the step length a_k and the
  // orthogonalisation coefficients b_{k,j}.  These are the data from which a
  // FIXED polynomial smoother can be harvested: if they are stable from call
  // to call, the adaptive GCR can be replaced by a stationary p(A) with the
  // same applies and no reductions.  Off by default; boss rank prints.
  int  LogCoeffs = 0;
  void LogCoefficients(int l) { LogCoeffs = l; };
  // Optional recorder of the per-step coefficients (means over calls), for
  // replay by GCRReplaySmoother (Smoothers.h).  Records only; no effect on
  // the iteration.
  GCRCoefficients *Recorder = nullptr;
  void SetCoefficientRecorder(GCRCoefficients *r) { Recorder = r; if(r) r->mmax = mmax; };
  // Free the persistent history (e.g. when this solver is replaced by a
  // replayed polynomial): re-made on the next call if ever needed again.
  void ReleaseHistory(void) { q.clear(); p.clear(); qq.clear(); hist_grid = nullptr; };

  PrecGeneralisedConjugateResidualNonHermitian(RealD tol,Integer maxit,LinearOperatorBase<Field> &_Linop,LinearFunction<Field> &Prec,int _mmax,int _nstep) : 
    Tolerance(tol), 
    MaxIterations(maxit),
    Linop(_Linop),
    Preconditioner(Prec),
    mmax(_mmax),
    nstep(_nstep)
  {
    Level(1);
    verbose=1;
  };

  void operator() (const Field &src, Field &psi){

    //    psi=Zero();
    RealD cp, ssq,rsq;
    ssq=norm2(src);
    SSQ=ssq;
    rsq=Tolerance*Tolerance*ssq;
      
    Field r(src.Grid());

    PrecTimer.Reset();
    MatTimer.Reset();
    LinalgTimer.Reset();

    GridStopWatch SolverTimer;
    SolverTimer.Start();

    steps=0;
    FirstCycle=1;
    for(int k=0;k<MaxIterations;k++){

      cp=GCRnStep(src,psi,rsq);

      GCRLogLevel <<"PGCR("<<mmax<<","<<nstep<<") "<< steps <<" steps cp = "<<sqrt(cp/ssq)<<" target "<<sqrt(rsq/ssq) <<std::endl;

      if(cp<rsq) {

	SolverTimer.Stop();

	Linop.Op(psi,r);
	axpy(r,-1.0,src,r);
	RealD tr = norm2(r);
	GCRLogLevel<<"PGCR: Converged on iteration " <<steps
		 << " computed residual "<<sqrt(cp/ssq)
		 << " true residual "    <<sqrt(tr/ssq)
		 << " target "           <<Tolerance <<std::endl;

	GCRLogLevel<<"PGCR Time elapsed: Total  "<< SolverTimer.Elapsed() <<std::endl;
	return;
      }

    }
    GCRLogLevel<<"Variable Preconditioned GCR did not converge"<<std::endl;
    //    GRID_ASSERT(0);
  }

  RealD GCRnStep(const Field &src, Field &psi,RealD rsq){

    RealD cp;
    ComplexD a, b;
    ComplexD rq;

    GridBase *grid = src.Grid();

    // Only r and one scratch for the restart residual; the new p/q directions
    // are produced directly in their persistent history slots.
    Field r(grid);
    Field Az(grid);

    ////////////////////////////////
    // history for flexible orthog
    ////////////////////////////////
    // History arrays are PERSISTENT across calls (allocated once per grid,
    // re-made only if the grid or mmax changes).  The per-call form
    // std::vector<Field>(mmax,grid) built a temporary and copy-constructed it
    // mmax times on every restart cycle -- measured ~5 ms per fine-smoother
    // call.  Safe: every entry is written before it is read within a cycle
    // (q[kp],p[kp] assigned before the northog loop can reach them), so no
    // stale content is ever consumed.
    if ( hist_grid != grid || (int)q.size() != mmax ) {
      q.clear(); p.clear();
      q.reserve(mmax); p.reserve(mmax);
      for(int i=0;i<mmax;i++){ q.emplace_back(grid); p.emplace_back(grid); }
      qq.assign(mmax, 0.0);
      hist_grid = grid;
    }
      
    GCRLogLevel<< "PGCR nStep("<<nstep<<")"<<std::endl;

    //////////////////////////////////
    // r0 = src - A x0.  ZeroGuess: on the first cycle x0==0 by caller
    // contract (enforced here), so r0 = src exactly; skip the apply.
    // Restart cycles (psi!=0) always do the full computation.
    //////////////////////////////////
    if (ZeroGuess && FirstCycle) {
      psi = Zero();
      LinalgTimer.Start();
      r = src;
      LinalgTimer.Stop();
    } else {
      MatTimer.Start();
      Linop.Op(psi,Az);
      MatTimer.Stop();
      LinalgTimer.Start();
      r=src-Az;
      LinalgTimer.Stop();
    }
    FirstCycle=0;

    /////////////////////
    // p = Prec(r)
    /////////////////////

    // p[0] = Prec(r), q[0] = A p[0], written straight into the history slots
    PrecTimer.Start();
    Preconditioner(r,p[0]);
    PrecTimer.Stop();

    MatTimer.Start();
    Linop.Op(p[0],q[0]);
    MatTimer.Stop();

    LinalgTimer.Start();

    qq[0]= norm2(q[0]);

    cp =norm2(r);
    LinalgTimer.Stop();
    GCRLogLevel<< "PGCR true residual "<< sqrt(cp/SSQ)     <<std::endl;

    for(int k=0;k<nstep;k++){
      GRID_TRACE("PGCR_step");
      steps++;

      int kp     = k+1;
      int peri_k = k %mmax;
      int peri_kp= kp%mmax;

      LinalgTimer.Start();
      rq= innerProduct(q[peri_k],r); // what if rAr not real?
      a = rq/qq[peri_k];
      if ( Recorder ) Recorder->RecordA(k,a);

      axpy(psi,a,p[peri_k],psi);         

      cp = axpy_norm(r,-a,q[peri_k],r);
      LinalgTimer.Stop();
      if ( LogCoeffs ) {
        GCRLogLevel<<"coeff["<<k<<"] a=("<<real(a)<<","<<imag(a)<<")"
                   <<" |r|/|r0|="<<sqrt(cp/SSQ)<<std::endl;
      }

      GCRLogLevel<< "PGCR step["<<steps<<"]  resid " << sqrt(cp/SSQ)<<std::endl;

      if((k==nstep-1)||(cp<rsq)){
	return cp;
      }

      // New direction written straight into its history slot: p = Prec(r), q = A p.
      PrecTimer.Start();
      Preconditioner(r,p[peri_kp]);
      PrecTimer.Stop();

      MatTimer.Start();
      Linop.Op(p[peri_kp],q[peri_kp]);
      MatTimer.Stop();

      LinalgTimer.Start();

      int northog = ((kp)>(mmax-1))?(mmax-1):(kp);  // if more than mmax done, we orthog all mmax history.
      std::ostringstream bs;

      // Classical Gram-Schmidt: every coefficient is taken against the
      // UN-updated new q (so all northog inner products are independent and
      // batchable), then the window is applied.  The coefficient is complex:
      // for a non-Hermitian operator <q_j,Aq> is complex and keeping only the
      // real part left the q's non-orthogonal.
      // Batched: one fused kernel + one reduction for all coefficients,
      // one fused pass per update (independent of mmax).
      std::vector<const Field*> qwin(northog), pwin(northog);
      for(int back=0;back<northog;back++){
	int peri_back=(k-back)%mmax;   	  GRID_ASSERT((k-back)>=0);
	GRID_ASSERT(peri_back!=peri_kp);
	qwin[back] = &q[peri_back];
	pwin[back] = &p[peri_back];
      }
      std::vector<ComplexD> bcoef;
      innerProductMulti(bcoef,qwin,q[peri_kp]);
      for(int back=0;back<northog;back++){
	int peri_back=(k-back)%mmax;
	bcoef[back] = -bcoef[back]/qq[peri_back];
        if ( LogCoeffs ) bs<<" b["<<back<<"]="<<bcoef[back];
      }
      if ( Recorder && northog ) Recorder->RecordB(k,bcoef);
      if ( northog ) {
	axpyMulti(p[peri_kp],bcoef,pwin);
	qq[peri_kp]=axpyMultiNorm(q[peri_kp],bcoef,qwin);
      } else {
	qq[peri_kp]=norm2(q[peri_kp]);
      }
      if ( LogCoeffs && northog ) {
        GCRLogLevel<<"coeff["<<k<<"]"<<bs.str()<<std::endl;
      }
      qq[peri_kp]=norm2(q[peri_kp]); // could use axpy_norm
      LinalgTimer.Stop();
    }
    GRID_ASSERT(0); // never reached
    return cp;
  }
};
NAMESPACE_END(Grid);

#undef GCRLogLevel
#endif
