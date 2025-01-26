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
  MultiRHSBlockCGLinalg<Field> _BlockCGLinalg;

  GridStopWatch ProjectTimer;
  GridStopWatch PromoteTimer;
  GridStopWatch DeflateTimer;
  GridStopWatch CoarseTimer;
  GridStopWatch FineTimer;
  GridStopWatch SmoothTimer;
  GridStopWatch InsertTimer;

  /*
    Field rrr;
  Field sss;
  Field qqq;
  Field zzz;
  */  
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
    /*
    rrr(fine),
    sss(fine),
    qqq(fine),
    zzz(fine)
*/
  {
    grid       = fine;
  };
  
  // Vector case
  virtual void operator() (std::vector<Field> &src, std::vector<Field> &x)
  {
    //    SolveSingleSystem(src,x);
    SolvePrecBlockCG(src,x);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////
// Thin QR factorisation (google it)
////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //Dimensions
  // R_{ferm x Nblock} =  Q_{ferm x Nblock} x  C_{Nblock x Nblock} -> ferm x Nblock
  //
  // Rdag R = m_rr = Herm = L L^dag        <-- Cholesky decomposition (LLT routine in Eigen)
  //
  //   Q  C = R => Q = R C^{-1}
  //
  // Want  Ident = Q^dag Q = C^{-dag} R^dag R C^{-1} = C^{-dag} L L^dag C^{-1} = 1_{Nblock x Nblock} 
  //
  // Set C = L^{dag}, and then Q^dag Q = ident 
  //
  // Checks:
  // Cdag C = Rdag R ; passes.
  // QdagQ  = 1      ; passes
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  void ThinQRfact (Eigen::MatrixXcd &m_zz,
		   Eigen::MatrixXcd &C,
		   Eigen::MatrixXcd &Cinv,
		   std::vector<Field> &  Q,
		   std::vector<Field> & MQ,
		   const std::vector<Field> & Z,
		   const std::vector<Field> & MZ)
  {
    RealD t0=usecond();
    _BlockCGLinalg.InnerProductMatrix(m_zz,MZ,Z);
    RealD t1=usecond();

    m_zz = 0.5*(m_zz+m_zz.adjoint());
    
    Eigen::MatrixXcd L    = m_zz.llt().matrixL(); 
    
    C    = L.adjoint();
    Cinv = C.inverse();
    
    RealD t3=usecond();
    _BlockCGLinalg.MulMatrix( Q,Cinv,Z);
    _BlockCGLinalg.MulMatrix(MQ,Cinv,MZ);
    RealD t4=usecond();
    std::cout << " ThinQRfact IP    :"<< t1-t0<<" us"<<std::endl;
    std::cout << " ThinQRfact Eigen :"<< t3-t1<<" us"<<std::endl;
    std::cout << " ThinQRfact MulMat:"<< t4-t3<<" us"<<std::endl;
  }

  virtual void SolvePrecBlockCG (std::vector<Field> &src, std::vector<Field> &X)
  {
    std::cout << GridLogMessage<<"HDCG: mrhs fPrecBlockcg starting"<<std::endl;
    src[0].Grid()->Barrier();
    int nrhs = src.size();
    //    std::vector<RealD> f(nrhs);
    //    std::vector<RealD> rtzp(nrhs);
    //    std::vector<RealD> rtz(nrhs);
    //    std::vector<RealD> a(nrhs);
    //    std::vector<RealD> d(nrhs);
    //    std::vector<RealD> b(nrhs);
    //    std::vector<RealD> rptzp(nrhs);

    ////////////////////////////////////////////
    //Initial residual computation & set up
    ////////////////////////////////////////////
    std::vector<RealD> ssq(nrhs);
    for(int rhs=0;rhs<nrhs;rhs++){
      ssq[rhs]=norm2(src[rhs]); assert(ssq[rhs]!=0.0);
    }      

    ///////////////////////////
    // Fields -- eliminate duplicates between fPcg and block cg
    ///////////////////////////
    std::vector<Field> Mtmp(nrhs,grid);
    std::vector<Field> tmp(nrhs,grid);
    std::vector<Field>   Z(nrhs,grid); // Rename Z to R
    std::vector<Field>  MZ(nrhs,grid); // Rename MZ to Z
    std::vector<Field>   Q(nrhs,grid); // 
    std::vector<Field>  MQ(nrhs,grid); // Rename to P
    std::vector<Field>   D(nrhs,grid);
    std::vector<Field>  AD(nrhs,grid);
    
    /************************************************************************
     * Preconditioned Block conjugate gradient rQ
     * Generalise Sebastien Birk Thesis, after Dubrulle 2001.
     * Introduce preconditioning following Saad Ch9
     ************************************************************************
     * Dimensions:
     *
     *   X,B etc... ==(Nferm x nrhs)
     *  Matrix A==(Nferm x Nferm)
     *  
     * Nferm = Nspin x Ncolour x Ncomplex x Nlattice_site
     * QC => Thin QR factorisation (google it)
     *
     * R = B-AX
     * Z = Mi R
     * QC = Z
     * D = Q 
     * for k: 
     *   R  = AD
     *   Z  = Mi R
     *   M  = [D^dag R]^{-1}
     *   X  = X + D M C
     *   QS = Q - Z.M
     *   D  = Q + D S^dag
     *   C  = S C
     */
    Eigen::MatrixXcd m_DZ     = Eigen::MatrixXcd::Identity(nrhs,nrhs);
    Eigen::MatrixXcd m_M      = Eigen::MatrixXcd::Identity(nrhs,nrhs);
    Eigen::MatrixXcd m_zz     = Eigen::MatrixXcd::Zero(nrhs,nrhs);
    Eigen::MatrixXcd m_rr     = Eigen::MatrixXcd::Zero(nrhs,nrhs);
    
    Eigen::MatrixXcd m_C      = Eigen::MatrixXcd::Zero(nrhs,nrhs);
    Eigen::MatrixXcd m_Cinv   = Eigen::MatrixXcd::Zero(nrhs,nrhs);
    Eigen::MatrixXcd m_S      = Eigen::MatrixXcd::Zero(nrhs,nrhs);
    Eigen::MatrixXcd m_Sinv   = Eigen::MatrixXcd::Zero(nrhs,nrhs);
    
    Eigen::MatrixXcd m_tmp    = Eigen::MatrixXcd::Identity(nrhs,nrhs);
    Eigen::MatrixXcd m_tmp1   = Eigen::MatrixXcd::Identity(nrhs,nrhs);

    GridStopWatch HDCGTimer;

    //////////////////////////
    // x0 = Vstart -- possibly modify guess
    //////////////////////////
    Vstart(X,src);

    //////////////////////////
    // R = B-AX
    //////////////////////////
    for(int rhs=0;rhs<nrhs;rhs++){
      // r0 = b -A x0
      _FineLinop.HermOp(X[rhs],tmp[rhs]);
      axpy (Z[rhs], -1.0,tmp[rhs], src[rhs]);    // Computes R=Z=src - A X0
    }

    //////////////////////////////////
    // Compute MZ = M1 Z = M1 B - M1 A x0
    //////////////////////////////////
    PcgM1(Z,MZ);  

    //////////////////////////////////
    // QC = Z
    //////////////////////////////////
    ThinQRfact (m_zz, m_C, m_Cinv, Q, MQ, Z, MZ);

    //////////////////////////////////
    // D=MQ
    //////////////////////////////////
    for(int b=0;b<nrhs;b++) D[b]=MQ[b]; // LLT rotation of the MZ basis of search dirs

    std::cout << GridLogMessage<<"PrecBlockCGrQ vec computed initial residual and QR fact " <<std::endl;

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
    GridStopWatch InnerProdTimer;

    HDCGTimer.Start();

    std::vector<RealD> rn(nrhs);
    for (int k=0;k<=MaxIterations;k++){

      ////////////////////
      // Z  = AD
      ////////////////////
      M3Timer.Start();
      for(int b=0;b<nrhs;b++) _FineLinop.HermOp(D[b], Z[b]);      
      M3Timer.Stop();

      ////////////////////
      // MZ  = M1 Z <==== the Multigrid preconditioner
      ////////////////////
      M1Timer.Start();
      PcgM1(Z,MZ);
      M1Timer.Stop();

      FineTimer.Start();
      ////////////////////
      // M  = [D^dag Z]^{-1} = (<Ddag MZ>_M)^{-1} inner prod, generalising Saad derivation of Precon CG
      ////////////////////
      InnerProdTimer.Start();
      _BlockCGLinalg.InnerProductMatrix(m_DZ,D,Z);
      InnerProdTimer.Stop();
      m_M       = m_DZ.inverse();

      ///////////////////////////
      // X  = X + D MC
      ///////////////////////////
      m_tmp     = m_M * m_C;
      LinalgTimer.Start();
      _BlockCGLinalg.MaddMatrix(X,m_tmp, D,X);     // D are the search directions and X takes the updates 
      LinalgTimer.Stop();

      ///////////////////////////
      // QS = Q - M Z
      // (MQ) S = MQ - M (M1Z)
      ///////////////////////////
      LinalgTimer.Start();
      _BlockCGLinalg.MaddMatrix(tmp ,m_M, Z, Q,-1.0);
      _BlockCGLinalg.MaddMatrix(Mtmp,m_M,MZ,MQ,-1.0);
      ThinQRfact (m_zz, m_S, m_Sinv, Q, MQ, tmp, Mtmp);
      LinalgTimer.Stop();

      ////////////////////////////
      // D  = MQ + D S^dag
      ////////////////////////////
      m_tmp = m_S.adjoint();
      LinalgTimer.Start();
      _BlockCGLinalg.MaddMatrix(D,m_tmp,D,MQ);
      LinalgTimer.Stop();

      ////////////////////////////
      // C  = S C
      ////////////////////////////
      m_C = m_S*m_C;
      
      ////////////////////////////
      // convergence monitor
      ////////////////////////////
      m_rr = m_C.adjoint() * m_C;
      
      FineTimer.Stop();

      RealD max_resid=0;
      RealD rrsum=0;
      RealD sssum=0;
      RealD rr;

      for(int b=0;b<nrhs;b++) {
	rrsum+=real(m_rr(b,b));
	sssum+=ssq[b];
	rr = real(m_rr(b,b))/ssq[b];
	if ( rr > max_resid ) max_resid = rr;
      }
      std::cout << GridLogMessage <<
	  "\t Prec BlockCGrQ Iteration "<<k<<" ave resid "<< std::sqrt(rrsum/sssum) << " max "<< std::sqrt(max_resid) <<std::endl;


      if ( max_resid < Tolerance*Tolerance ) { 

	HDCGTimer.Stop();
	std::cout<<GridLogMessage<<"HDCG: mrhs PrecBlockCGrQ converged in "<<k<<" iterations and "<<HDCGTimer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs PrecBlockCGrQ : Linalg  "<<LinalgTimer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs PrecBlockCGrQ : fine H  "<<M3Timer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs PrecBlockCGrQ : prec M1 "<<M1Timer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"**** M1 breakdown:"<<std::endl;
	std::cout<<GridLogMessage<<"HDCG: mrhs PrecBlockCGrQ : Project "<<ProjectTimer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs PrecBlockCGrQ : Promote "<<PromoteTimer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs PrecBlockCGrQ : Deflate "<<DeflateTimer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs PrecBlockCGrQ : Coarse  "<<CoarseTimer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs PrecBlockCGrQ : Fine    "<<FineTimer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs PrecBlockCGrQ : Smooth  "<<SmoothTimer.Elapsed()<<std::endl;;
	std::cout<<GridLogMessage<<"HDCG: mrhs PrecBlockCGrQ : Insert  "<<InsertTimer.Elapsed()<<std::endl;;

	for(int rhs=0;rhs<nrhs;rhs++){

	  _FineLinop.HermOp(X[rhs],tmp[rhs]);			  

	  Field mytmp(grid);
	  axpy(mytmp,-1.0,src[rhs],tmp[rhs]);
      
	  RealD  xnorm   = sqrt(norm2(X[rhs]));
	  RealD  srcnorm = sqrt(norm2(src[rhs]));
	  RealD  tmpnorm = sqrt(norm2(mytmp));
	  RealD  true_residual = tmpnorm/srcnorm;
	  std::cout<<GridLogMessage
		   <<"HDCG: true residual ["<<rhs<<"] is "<<true_residual
		   <<" solution "<<xnorm
		   <<" source "<<srcnorm
		   <<std::endl;
	}
	return;
      }
      
    }
    HDCGTimer.Stop();
    std::cout<<GridLogMessage<<"HDCG: PrecBlockCGrQ not converged "<<HDCGTimer.Elapsed()<<std::endl;
    assert(0);
  }

  virtual void SolveSingleSystem (std::vector<Field> &src, std::vector<Field> &x)
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

    //    this->rrr=in[0];

#undef SMOOTHER_BLOCK_SOLVE
#if SMOOTHER_BLOCK_SOLVE
    this->SmoothTimer.Start();
    this->_Smoother(in,Min);
    this->SmoothTimer.Stop();
#else
    for(int rhs=0;rhs<nrhs;rhs++) {
      this->SmoothTimer.Start();
      this->_Smoother(in[rhs],Min[rhs]);
      this->SmoothTimer.Stop();
    }
#endif
    //    this->sss=Min[0];
    
    for(int rhs=0;rhs<nrhs;rhs++) {
      
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
    //    this->qqq=tmp[0];
    for(int rhs=0;rhs<nrhs;rhs++) {
      axpy(out[rhs],1.0,Min[rhs],tmp[rhs]); // Min+tmp
    }
    //    this->zzz=out[0];
    this->FineTimer.Stop();
  }
};


NAMESPACE_END(Grid);


