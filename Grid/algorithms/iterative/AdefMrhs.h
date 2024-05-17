    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/AdefGeneric.h

    Copyright (C) 2015

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
#pragma once


  /*
   * Compared to Tang-2009:  P=Pleft. P^T = PRight Q=MssInv. 
   * Script A = SolverMatrix 
   * Script P = Preconditioner
   *
   * Implement ADEF-2
   *
   * Vstart = P^Tx + Qb
   * M1 = P^TM + Q
   * M2=M3=1
   */
NAMESPACE_BEGIN(Grid);


template<class Field>
class TwoLevelCGmrhs
{
 public:
  RealD   Tolerance;
  Integer MaxIterations;
  GridBase *grid;

  // Fine operator, Smoother, CoarseSolver
  LinearOperatorBase<Field>   &_FineLinop;
  LinearFunction<Field>   &_Smoother;

  GridStopWatch ProjectTimer;
  GridStopWatch PromoteTimer;
  GridStopWatch DeflateTimer;
  GridStopWatch CoarseTimer;
  GridStopWatch FineTimer;
  GridStopWatch SmoothTimer;
  GridStopWatch InsertTimer;

  
  // more most opertor functions
  TwoLevelCGmrhs(RealD tol,
		 Integer maxit,
		 LinearOperatorBase<Field>   &FineLinop,
		 LinearFunction<Field>       &Smoother,
		 GridBase *fine) : 
    Tolerance(tol), 
    MaxIterations(maxit),
    _FineLinop(FineLinop),
    _Smoother(Smoother)
  {
    grid       = fine;
  };
  
  // Vector case
  virtual void operator() (std::vector<Field> &src, std::vector<Field> &x)
  {
    std::cout << GridLogMessage<<"HDCG: mrhs fPcg starting"<<std::endl;
    src[0].Grid()->Barrier();
    int nrhs = src.size();
    std::vector<RealD> f(nrhs);
    std::vector<RealD> rtzp(nrhs);
    std::vector<RealD> rtz(nrhs);
    std::vector<RealD> a(nrhs);
    std::vector<RealD> d(nrhs);
    std::vector<RealD> b(nrhs);
    std::vector<RealD> rptzp(nrhs);
    /////////////////////////////
    // Set up history vectors
    /////////////////////////////
    int mmax = 3;

    std::vector<std::vector<Field> > p(nrhs);   for(int r=0;r<nrhs;r++)  p[r].resize(mmax,grid);
    std::vector<std::vector<Field> > mmp(nrhs); for(int r=0;r<nrhs;r++) mmp[r].resize(mmax,grid);
    std::vector<std::vector<RealD> > pAp(nrhs); for(int r=0;r<nrhs;r++) pAp[r].resize(mmax);

    std::vector<Field> z(nrhs,grid);
    std::vector<Field>  mp (nrhs,grid);
    std::vector<Field>  r  (nrhs,grid);
    std::vector<Field>  mu (nrhs,grid);

    //Initial residual computation & set up
    std::vector<RealD> src_nrm(nrhs);
    for(int rhs=0;rhs<nrhs;rhs++) {
      src_nrm[rhs]=norm2(src[rhs]);
      assert(src_nrm[rhs]!=0.0);
    }
    std::vector<RealD> tn(nrhs);

    GridStopWatch HDCGTimer;
    //////////////////////////
    // x0 = Vstart -- possibly modify guess
    //////////////////////////
    Vstart(x,src);

    for(int rhs=0;rhs<nrhs;rhs++){
      // r0 = b -A x0
      _FineLinop.HermOp(x[rhs],mmp[rhs][0]);
      axpy (r[rhs], -1.0,mmp[rhs][0], src[rhs]);    // Recomputes r=src-Ax0
    }

    //////////////////////////////////
    // Compute z = M1 x
    //////////////////////////////////
    // This needs a multiRHS version for acceleration
    PcgM1(r,z);

    std::vector<RealD> ssq(nrhs);
    std::vector<RealD> rsq(nrhs);
    std::vector<Field> pp(nrhs,grid);

    for(int rhs=0;rhs<nrhs;rhs++){
      rtzp[rhs] =real(innerProduct(r[rhs],z[rhs]));
      p[rhs][0]=z[rhs];
      ssq[rhs]=norm2(src[rhs]);
      rsq[rhs]=  ssq[rhs]*Tolerance*Tolerance;
      //      std::cout << GridLogMessage<<"mrhs HDCG: "<<rhs<<" k=0 residual "<<rtzp[rhs]<<" rsq "<<rsq[rhs]<<"\n";
    }

    ProjectTimer.Reset();
    PromoteTimer.Reset();
    DeflateTimer.Reset();
    CoarseTimer.Reset();
    SmoothTimer.Reset();
    FineTimer.Reset();
    InsertTimer.Reset();

    GridStopWatch M1Timer;
    GridStopWatch M2Timer;
    GridStopWatch M3Timer;
    GridStopWatch LinalgTimer;

    HDCGTimer.Start();

    std::vector<RealD> rn(nrhs);
    for (int k=0;k<=MaxIterations;k++){
    
      int peri_k  = k % mmax;
      int peri_kp = (k+1) % mmax;

      for(int rhs=0;rhs<nrhs;rhs++){
	rtz[rhs]=rtzp[rhs];
	M3Timer.Start();
	d[rhs]= PcgM3(p[rhs][peri_k],mmp[rhs][peri_k]);
	M3Timer.Stop();
	a[rhs] = rtz[rhs]/d[rhs];

	LinalgTimer.Start();
	// Memorise this
	pAp[rhs][peri_k] = d[rhs];

	axpy(x[rhs],a[rhs],p[rhs][peri_k],x[rhs]);
	rn[rhs] = axpy_norm(r[rhs],-a[rhs],mmp[rhs][peri_k],r[rhs]);
	LinalgTimer.Stop();
      }

      // Compute z = M x (for *all* RHS)
      M1Timer.Start();
      PcgM1(r,z);
      M1Timer.Stop();
      
      RealD max_rn=0.0;
      LinalgTimer.Start();
      for(int rhs=0;rhs<nrhs;rhs++){

	rtzp[rhs] =real(innerProduct(r[rhs],z[rhs]));

	//	std::cout << GridLogMessage<<"HDCG::fPcg rhs"<<rhs<<" iteration "<<k<<" : inner rtzp "<<rtzp[rhs]<<"\n";
	mu[rhs]=z[rhs];

	p[rhs][peri_kp]=mu[rhs];

	// Standard search direction p == z + b p 
	b[rhs] = (rtzp[rhs])/rtz[rhs];

	int northog = (k>mmax-1)?(mmax-1):k;        // This is the fCG-Tr(mmax-1) algorithm
	for(int back=0; back < northog; back++){
	  int peri_back = (k-back)%mmax;
	  RealD pbApk= real(innerProduct(mmp[rhs][peri_back],p[rhs][peri_kp]));
	  RealD beta = -pbApk/pAp[rhs][peri_back];
	  axpy(p[rhs][peri_kp],beta,p[rhs][peri_back],p[rhs][peri_kp]);
	}

	RealD rrn=sqrt(rn[rhs]/ssq[rhs]);
	RealD rtn=sqrt(rtz[rhs]/ssq[rhs]);
	RealD rtnp=sqrt(rtzp[rhs]/ssq[rhs]);
	
	std::cout<<GridLogMessage<<"HDCG:fPcg rhs "<<rhs<<" k= "<<k<<" residual = "<<rrn<<"\n";
	if ( rrn > max_rn ) max_rn = rrn;
      }
      LinalgTimer.Stop();

      // Stopping condition based on worst case
      if ( max_rn <= Tolerance ) { 

	HDCGTimer.Stop();
	std::cout<<GridLogMessage<<"HDCG: mrhs fPcg converged in "<<k<<" iterations and "<<HDCGTimer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs fPcg : Linalg  "<<LinalgTimer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs fPcg : fine M3 "<<M3Timer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs fPcg : prec M1 "<<M1Timer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"**** M1 breakdown:"<<std::endl;
	std::cout<<GridLogMessage<<"HDCG: mrhs fPcg : Project "<<ProjectTimer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs fPcg : Promote "<<PromoteTimer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs fPcg : Deflate "<<DeflateTimer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs fPcg : Coarse  "<<CoarseTimer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs fPcg : Fine    "<<FineTimer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs fPcg : Smooth  "<<SmoothTimer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs fPcg : Insert  "<<InsertTimer.Elapsed()<<std::endl;;

	for(int rhs=0;rhs<nrhs;rhs++){
	  _FineLinop.HermOp(x[rhs],mmp[rhs][0]);			  
	  Field tmp(grid);
	  axpy(tmp,-1.0,src[rhs],mmp[rhs][0]);
      
	  RealD  mmpnorm = sqrt(norm2(mmp[rhs][0]));
	  RealD  xnorm   = sqrt(norm2(x[rhs]));
	  RealD  srcnorm = sqrt(norm2(src[rhs]));
	  RealD  tmpnorm = sqrt(norm2(tmp));
	  RealD  true_residual = tmpnorm/srcnorm;
	  std::cout<<GridLogMessage
		   <<"HDCG: true residual ["<<rhs<<"] is "<<true_residual
		   <<" solution "<<xnorm
		   <<" source "<<srcnorm
		   <<" mmp "<<mmpnorm	  
		   <<std::endl;
	}
	return;
      }
      
    }
    HDCGTimer.Stop();
    std::cout<<GridLogMessage<<"HDCG: not converged "<<HDCGTimer.Elapsed()<<std::endl;
    for(int rhs=0;rhs<nrhs;rhs++){
      RealD  xnorm   = sqrt(norm2(x[rhs]));
      RealD  srcnorm = sqrt(norm2(src[rhs]));
      std::cout<<GridLogMessage<<"HDCG: non-converged solution "<<xnorm<<" source "<<srcnorm<<std::endl;
    }
  }
  

 public:

  virtual void PcgM1(std::vector<Field> & in,std::vector<Field> & out) = 0;
  virtual void Vstart(std::vector<Field> & x,std::vector<Field> & src) = 0;
  virtual void PcgM2(const Field & in, Field & out) {
    out=in;
  }

  virtual RealD PcgM3(const Field & p, Field & mmp){
    RealD dd;
    _FineLinop.HermOp(p,mmp);
    ComplexD dot = innerProduct(p,mmp);
    dd=real(dot);
    return dd;
  }

};

template<class Field, class CoarseField>
class TwoLevelADEF2mrhs : public TwoLevelCGmrhs<Field>
{
public:
  GridBase *coarsegrid;
  GridBase *coarsegridmrhs;
  LinearFunction<CoarseField> &_CoarseSolverMrhs;
  LinearFunction<CoarseField> &_CoarseSolverPreciseMrhs;
  MultiRHSBlockProject<Field>    &_Projector;
  MultiRHSDeflation<CoarseField> &_Deflator;

  
  TwoLevelADEF2mrhs(RealD tol,
		    Integer maxit,
		    LinearOperatorBase<Field>    &FineLinop,
		    LinearFunction<Field>        &Smoother,
		    LinearFunction<CoarseField>  &CoarseSolverMrhs,
		    LinearFunction<CoarseField>  &CoarseSolverPreciseMrhs,
		    MultiRHSBlockProject<Field>    &Projector,
		    MultiRHSDeflation<CoarseField> &Deflator,
		    GridBase *_coarsemrhsgrid) :
    TwoLevelCGmrhs<Field>(tol, maxit,FineLinop,Smoother,Projector.fine_grid),
    _CoarseSolverMrhs(CoarseSolverMrhs),
    _CoarseSolverPreciseMrhs(CoarseSolverPreciseMrhs),
    _Projector(Projector),
    _Deflator(Deflator)
  {
    coarsegrid = Projector.coarse_grid;
    coarsegridmrhs = _coarsemrhsgrid;// Thi could be in projector
  };

  // Override Vstart
  virtual void Vstart(std::vector<Field> & x,std::vector<Field> & src)
  {
    int nrhs=x.size();
    ///////////////////////////////////
    // Choose x_0 such that 
    // x_0 = guess +  (A_ss^inv) r_s = guess + Ass_inv [src -Aguess]
    //                               = [1 - Ass_inv A] Guess + Assinv src
    //                               = P^T guess + Assinv src 
    //                               = Vstart  [Tang notation]
    // This gives:
    // W^T (src - A x_0) = src_s - A guess_s - r_s
    //                   = src_s - (A guess)_s - src_s  + (A guess)_s 
    //                   = 0 
    ///////////////////////////////////
    std::vector<CoarseField> PleftProj(nrhs,this->coarsegrid);
    std::vector<CoarseField> PleftMss_proj(nrhs,this->coarsegrid);
    CoarseField PleftProjMrhs(this->coarsegridmrhs);
    CoarseField PleftMss_projMrhs(this->coarsegridmrhs);

    this->_Projector.blockProject(src,PleftProj);
    this->_Deflator.DeflateSources(PleftProj,PleftMss_proj);
    for(int rhs=0;rhs<nrhs;rhs++) {
      InsertSliceFast(PleftProj[rhs],PleftProjMrhs,rhs,0);
      InsertSliceFast(PleftMss_proj[rhs],PleftMss_projMrhs,rhs,0); // the guess
    }
    
    this->_CoarseSolverPreciseMrhs(PleftProjMrhs,PleftMss_projMrhs); // Ass^{-1} r_s

    for(int rhs=0;rhs<nrhs;rhs++) {
      ExtractSliceFast(PleftMss_proj[rhs],PleftMss_projMrhs,rhs,0);
    }
    this->_Projector.blockPromote(x,PleftMss_proj);
  }

  virtual void PcgM1(std::vector<Field> & in,std::vector<Field> & out){

    int nrhs=in.size();

    // [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]
    std::vector<Field> tmp(nrhs,this->grid);
    std::vector<Field> Min(nrhs,this->grid);

    std::vector<CoarseField> PleftProj(nrhs,this->coarsegrid);
    std::vector<CoarseField> PleftMss_proj(nrhs,this->coarsegrid);

    CoarseField PleftProjMrhs(this->coarsegridmrhs);
    CoarseField PleftMss_projMrhs(this->coarsegridmrhs);

    for(int rhs=0;rhs<nrhs;rhs++) {

      this->SmoothTimer.Start();
      this->_Smoother(in[rhs],Min[rhs]);
      this->SmoothTimer.Stop();

      this->FineTimer.Start();
      this->_FineLinop.HermOp(Min[rhs],out[rhs]);

      axpy(tmp[rhs],-1.0,out[rhs],in[rhs]);          // resid  = in - A Min
      this->FineTimer.Stop();

    }

    this->ProjectTimer.Start();
    this->_Projector.blockProject(tmp,PleftProj);
    this->ProjectTimer.Stop();
    this->DeflateTimer.Start();
    this->_Deflator.DeflateSources(PleftProj,PleftMss_proj);
    this->DeflateTimer.Stop();
    this->InsertTimer.Start();
    for(int rhs=0;rhs<nrhs;rhs++) {
      InsertSliceFast(PleftProj[rhs],PleftProjMrhs,rhs,0);
      InsertSliceFast(PleftMss_proj[rhs],PleftMss_projMrhs,rhs,0); // the guess
    }
    this->InsertTimer.Stop();

    this->CoarseTimer.Start();
    this->_CoarseSolverMrhs(PleftProjMrhs,PleftMss_projMrhs); // Ass^{-1} [in - A Min]_s
    this->CoarseTimer.Stop();

    this->InsertTimer.Start();
    for(int rhs=0;rhs<nrhs;rhs++) {
      ExtractSliceFast(PleftMss_proj[rhs],PleftMss_projMrhs,rhs,0);
    }
    this->InsertTimer.Stop();
    this->PromoteTimer.Start();
    this->_Projector.blockPromote(tmp,PleftMss_proj);// tmp= Q[in - A Min]  
    this->PromoteTimer.Stop();
    this->FineTimer.Start();
    for(int rhs=0;rhs<nrhs;rhs++) {
      axpy(out[rhs],1.0,Min[rhs],tmp[rhs]); // Min+tmp
    }
    this->FineTimer.Stop();
  }
};
  

NAMESPACE_END(Grid);


