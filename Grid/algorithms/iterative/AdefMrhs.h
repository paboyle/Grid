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

#include <Grid/algorithms/multigrid/MrhsPreconditioner.h>

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


//////////////////////////////////////////////////////////////////////
// Two-level CG family on a vector of right-hand sides, with the
// preconditioner (M1 and Vstart) as an MrhsPreconditioner object:
//   fPcg           flexible PCG, one Krylov per rhs, shared M1 call
//   PrecBlockCGrQ  preconditioned BlockCGrQ (arXiv:2409.03904 s2.2);
//                  assumes a stationary preconditioner; not the
//                  production choice (its linalg grows as nrhs^2)
// Also a LinearFunction<Field>: one rhs is the nrhs=1 case.
//////////////////////////////////////////////////////////////////////
enum class MrhsCGAlgorithm { fPcg, PrecBlockCGrQ };

template<class Field>
class TwoLevelCGmrhs : public LinearFunction<Field>
{
 public:
  using LinearFunction<Field>::operator();
  RealD   Tolerance;
  Integer MaxIterations;
  GridBase *grid;
  MrhsCGAlgorithm Algorithm;

  LinearOperatorBase<Field>    &_FineLinop;
  MrhsPreconditioner<Field>    &_Precon;
  MultiRHSBlockCGLinalg<Field>  _BlockCGLinalg;

  TwoLevelCGmrhs(RealD tol,
		 Integer maxit,
		 LinearOperatorBase<Field> &FineLinop,
		 MrhsPreconditioner<Field> &Precon,
		 GridBase *fine,
		 MrhsCGAlgorithm alg = MrhsCGAlgorithm::fPcg) :
    Tolerance(tol),
    MaxIterations(maxit),
    Algorithm(alg),
    _FineLinop(FineLinop),
    _Precon(Precon)
  {
    grid = fine;
  };

  // mrhs entry point
  virtual void operator() (std::vector<Field> &src, std::vector<Field> &x)
  {
    if ( Algorithm == MrhsCGAlgorithm::PrecBlockCGrQ ) SolvePrecBlockCG(src,x);
    else                                              SolveSingleSystem(src,x);
  }
  // LinearFunction entry points
  virtual void operator() (const Field &src, Field &x)
  {
    std::vector<Field> S(1,src);
    std::vector<Field> X(1,x);
    (*this)(S,X);
    x = X[0];
  }
  virtual void operator() (const std::vector<Field> &src, std::vector<Field> &x)
  {
    std::vector<Field> S(src);
    (*this)(S,x);
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
      ssq[rhs]=norm2(src[rhs]); GRID_ASSERT(ssq[rhs]!=0.0);
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
    _Precon.Vstart(X,src);

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
    _Precon(Z,MZ);  

    //////////////////////////////////
    // QC = Z
    //////////////////////////////////
    ThinQRfact (m_zz, m_C, m_Cinv, Q, MQ, Z, MZ);

    //////////////////////////////////
    // D=MQ
    //////////////////////////////////
    for(int b=0;b<nrhs;b++) D[b]=MQ[b]; // LLT rotation of the MZ basis of search dirs

    std::cout << GridLogMessage<<"PrecBlockCGrQ vec computed initial residual and QR fact " <<std::endl;

    _Precon.ResetTimers();

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
      _Precon(Z,MZ);
      M1Timer.Stop();

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
	_Precon.ReportTimers("HDCG: mrhs PrecBlockCGrQ : ");

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
    GRID_ASSERT(0);
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
      GRID_ASSERT(src_nrm[rhs]!=0.0);
    }
    std::vector<RealD> tn(nrhs);

    GridStopWatch HDCGTimer;
    //////////////////////////
    // x0 = Vstart -- possibly modify guess
    //////////////////////////
    _Precon.Vstart(x,src);

    for(int rhs=0;rhs<nrhs;rhs++){
      // r0 = b -A x0
      _FineLinop.HermOp(x[rhs],mmp[rhs][0]);
      axpy (r[rhs], -1.0,mmp[rhs][0], src[rhs]);    // Recomputes r=src-Ax0
    }

    //////////////////////////////////
    // Compute z = M1 x
    //////////////////////////////////
    // This needs a multiRHS version for acceleration
    _Precon(r,z);

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

    _Precon.ResetTimers();

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
      _Precon(r,z);
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
	_Precon.ReportTimers("HDCG: mrhs fPcg : ");

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

  RealD PcgM3(const Field & p, Field & mmp){
    RealD dd;
    _FineLinop.HermOp(p,mmp);
    ComplexD dot = innerProduct(p,mmp);
    dd=real(dot);
    return dd;
  }
};

//////////////////////////////////////////////////////////////////////
// ADEF-2 as an mrhs preconditioner object (Tang, Nabben, Vuik, Erlangga):
//   M1     = [1 - Q A] M + Q   with  Q = P A_c^{-1} P^dag  (deflated coarse solve)
//   Vstart = Q b               the coarse-corrected start (precise coarse solve)
// The coarse space is the D+1 multiRHS field throughout: the projector's
// mixed overloads take the vector of fine fields straight to it, and the
// deflator and coarse solver act on it.  No per-rhs slices anywhere.
// The D+1 grid follows the solver's Nrhs through SetCoarseGridMrhs.
//////////////////////////////////////////////////////////////////////
// Projector_t defaults to a transfer operator whose STORE matches Field.
template<class Field, class CoarseField, class Projector_t = MultiRHSBlockProject<Field> >
class MrhsADEF2Preconditioner : public MrhsPreconditioner<Field>
{
public:
  GridBase                       *grid;             // fine
  GridBase                       *coarsegridmrhs;   // D+1, rhs innermost
  LinearOperatorBase<Field>      &_FineLinop;
  LinearFunction<Field>          &_Smoother;
  LinearFunction<CoarseField>    &_CoarseSolver;        // in M1
  LinearFunction<CoarseField>    &_CoarseSolverPrecise; // in Vstart
  Projector_t                    &_Projector;   // store precision need not match Field
  MultiRHSDeflation<CoarseField> &_Deflator;            // nev==0: no deflation, zero coarse guess

  MrhsADEF2Preconditioner(LinearOperatorBase<Field>      &FineLinop,
			  LinearFunction<Field>          &Smoother,
			  LinearFunction<CoarseField>    &CoarseSolver,
			  LinearFunction<CoarseField>    &CoarseSolverPrecise,
			  Projector_t                    &Projector,
			  MultiRHSDeflation<CoarseField> &Deflator,
			  GridBase *CoarseGridMrhs,
			  GridBase *FineGrid) :
    _FineLinop(FineLinop), _Smoother(Smoother),
    _CoarseSolver(CoarseSolver), _CoarseSolverPrecise(CoarseSolverPrecise),
    _Projector(Projector), _Deflator(Deflator)
  {
    // The fine grid is given, NOT taken from the projector: the projector's
    // grid carries its STORE's layout, which need not be this Field's (an
    // fp32 preconditioner may project through an fp64 store, or the reverse).
    grid           = FineGrid;
    coarsegridmrhs = CoarseGridMrhs;
  };
  // Store matches Field: the projector's grid IS this Field's grid.
  MrhsADEF2Preconditioner(LinearOperatorBase<Field>      &FineLinop,
			  LinearFunction<Field>          &Smoother,
			  LinearFunction<CoarseField>    &CoarseSolver,
			  LinearFunction<CoarseField>    &CoarseSolverPrecise,
			  Projector_t                    &Projector,
			  MultiRHSDeflation<CoarseField> &Deflator,
			  GridBase *CoarseGridMrhs) :
    MrhsADEF2Preconditioner(FineLinop,Smoother,CoarseSolver,CoarseSolverPrecise,
			    Projector,Deflator,CoarseGridMrhs,Projector.fine_grid) {};

  void SetCoarseGridMrhs(GridBase *CoarseGridMrhs) { coarsegridmrhs = CoarseGridMrhs; }

  // The coarse solver's starting guess: deflated, or zero
  void CoarseGuess(CoarseField &src,CoarseField &guess)
  {
    this->DeflateTimer.Start();
    if ( _Deflator.nev > 0 ) _Deflator.DeflateSources(src,guess);
    else                     guess = Zero();
    this->DeflateTimer.Stop();
  }
  // x_0 = Q b
  virtual void Vstart(std::vector<Field> &x,std::vector<Field> &src)
  {
    CoarseField Csrc(coarsegridmrhs), Csol(coarsegridmrhs);
    this->ProjectTimer.Start();  _Projector.blockProject(src,Csrc);  this->ProjectTimer.Stop();
    CoarseGuess(Csrc,Csol);
    this->CoarseTimer.Start();   _CoarseSolverPrecise(Csrc,Csol);    this->CoarseTimer.Stop();
    this->PromoteTimer.Start();  _Projector.blockPromote(x,Csol);    this->PromoteTimer.Stop();
  }
  // [1 - Q A] M in + Q in = Min + Q [in - A Min]
  virtual void operator()(std::vector<Field> &in,std::vector<Field> &out)
  {
    int nrhs=in.size();
    std::vector<Field> tmp(nrhs,grid), Min(nrhs,grid);
    CoarseField Csrc(coarsegridmrhs), Csol(coarsegridmrhs);

    this->SmoothTimer.Start();
    for(int r=0;r<nrhs;r++) _Smoother(in[r],Min[r]);
    this->SmoothTimer.Stop();

    this->FineTimer.Start();
    for(int r=0;r<nrhs;r++){
      _FineLinop.HermOp(Min[r],out[r]);
      axpy(tmp[r],-1.0,out[r],in[r]);          // in - A Min
    }
    this->FineTimer.Stop();

    this->ProjectTimer.Start();  _Projector.blockProject(tmp,Csrc);  this->ProjectTimer.Stop();
    CoarseGuess(Csrc,Csol);
    this->CoarseTimer.Start();   _CoarseSolver(Csrc,Csol);           this->CoarseTimer.Stop();
    this->PromoteTimer.Start();  _Projector.blockPromote(tmp,Csol);  this->PromoteTimer.Stop();

    this->FineTimer.Start();
    for(int r=0;r<nrhs;r++) axpy(out[r],1.0,Min[r],tmp[r]);       // Min + Q[in - A Min]
    this->FineTimer.Stop();
  }
};

NAMESPACE_END(Grid);
