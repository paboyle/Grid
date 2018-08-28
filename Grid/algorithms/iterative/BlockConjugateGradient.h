/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/BlockConjugateGradient.h

Copyright (C) 2017

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef GRID_BLOCK_CONJUGATE_GRADIENT_H
#define GRID_BLOCK_CONJUGATE_GRADIENT_H


namespace Grid {

enum BlockCGtype { BlockCG, BlockCGrQ, CGmultiRHS };

//////////////////////////////////////////////////////////////////////////
// Block conjugate gradient. Dimension zero should be the block direction
//////////////////////////////////////////////////////////////////////////
template <class Field>
class BlockConjugateGradient : public OperatorFunction<Field> {
 public:


  typedef typename Field::scalar_type scomplex;

  int blockDim ;
  int Nblock;

  BlockCGtype CGtype;
  bool ErrorOnNoConverge;  // throw an assert when the CG fails to converge.
                           // Defaults true.
  RealD Tolerance;
  Integer MaxIterations;
  Integer IterationsToComplete; //Number of iterations the CG took to finish. Filled in upon completion
  
  BlockConjugateGradient(BlockCGtype cgtype,int _Orthog,RealD tol, Integer maxit, bool err_on_no_conv = true)
    : Tolerance(tol), CGtype(cgtype),   blockDim(_Orthog),  MaxIterations(maxit), ErrorOnNoConverge(err_on_no_conv)
  {};

////////////////////////////////////////////////////////////////////////////////////////////////////
// Thin QR factorisation (google it)
////////////////////////////////////////////////////////////////////////////////////////////////////
void ThinQRfact (Eigen::MatrixXcd &m_rr,
		 Eigen::MatrixXcd &C,
		 Eigen::MatrixXcd &Cinv,
		 Field & Q,
		 const Field & R)
{
  int Orthog = blockDim; // First dimension is block dim; this is an assumption
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
  sliceInnerProductMatrix(m_rr,R,R,Orthog);

  // Force manifest hermitian to avoid rounding related
  m_rr = 0.5*(m_rr+m_rr.adjoint());

#if 0
  std::cout << " Calling Cholesky  ldlt on m_rr "  << m_rr <<std::endl;
  Eigen::MatrixXcd L_ldlt = m_rr.ldlt().matrixL(); 
  std::cout << " Called Cholesky  ldlt on m_rr "  << L_ldlt <<std::endl;
  auto  D_ldlt = m_rr.ldlt().vectorD(); 
  std::cout << " Called Cholesky  ldlt on m_rr "  << D_ldlt <<std::endl;
#endif

  //  std::cout << " Calling Cholesky  llt on m_rr "  <<std::endl;
  Eigen::MatrixXcd L    = m_rr.llt().matrixL(); 
  //  std::cout << " Called Cholesky  llt on m_rr "  << L <<std::endl;
  C    = L.adjoint();
  Cinv = C.inverse();
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Q = R C^{-1}
  //
  // Q_j  = R_i Cinv(i,j) 
  //
  // NB maddMatrix conventions are Right multiplication X[j] a[j,i] already
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  sliceMulMatrix(Q,Cinv,R,Orthog);
}
////////////////////////////////////////////////////////////////////////////////////////////////////
// Call one of several implementations
////////////////////////////////////////////////////////////////////////////////////////////////////
void operator()(LinearOperatorBase<Field> &Linop, const Field &Src, Field &Psi) 
{
  if ( CGtype == BlockCGrQ ) {
    BlockCGrQsolve(Linop,Src,Psi);
  } else if (CGtype == BlockCG ) {
    BlockCGsolve(Linop,Src,Psi);
  } else if (CGtype == CGmultiRHS ) {
    CGmultiRHSsolve(Linop,Src,Psi);
  } else {
    assert(0);
  }
}

////////////////////////////////////////////////////////////////////////////
// BlockCGrQ implementation:
//--------------------------
// X is guess/Solution
// B is RHS
// Solve A X_i = B_i    ;        i refers to Nblock index
////////////////////////////////////////////////////////////////////////////
void BlockCGrQsolve(LinearOperatorBase<Field> &Linop, const Field &B, Field &X) 
{
  int Orthog = blockDim; // First dimension is block dim; this is an assumption
  Nblock = B._grid->_fdimensions[Orthog];

  std::cout<<GridLogMessage<<" Block Conjugate Gradient : Orthog "<<Orthog<<" Nblock "<<Nblock<<std::endl;

  X.checkerboard = B.checkerboard;
  conformable(X, B);

  Field tmp(B);
  Field Q(B);
  Field D(B);
  Field Z(B);
  Field AD(B);

  Eigen::MatrixXcd m_DZ     = Eigen::MatrixXcd::Identity(Nblock,Nblock);
  Eigen::MatrixXcd m_M      = Eigen::MatrixXcd::Identity(Nblock,Nblock);
  Eigen::MatrixXcd m_rr     = Eigen::MatrixXcd::Zero(Nblock,Nblock);

  Eigen::MatrixXcd m_C      = Eigen::MatrixXcd::Zero(Nblock,Nblock);
  Eigen::MatrixXcd m_Cinv   = Eigen::MatrixXcd::Zero(Nblock,Nblock);
  Eigen::MatrixXcd m_S      = Eigen::MatrixXcd::Zero(Nblock,Nblock);
  Eigen::MatrixXcd m_Sinv   = Eigen::MatrixXcd::Zero(Nblock,Nblock);

  Eigen::MatrixXcd m_tmp    = Eigen::MatrixXcd::Identity(Nblock,Nblock);
  Eigen::MatrixXcd m_tmp1   = Eigen::MatrixXcd::Identity(Nblock,Nblock);

  // Initial residual computation & set up
  std::vector<RealD> residuals(Nblock);
  std::vector<RealD> ssq(Nblock);

  sliceNorm(ssq,B,Orthog);
  RealD sssum=0;
  for(int b=0;b<Nblock;b++) sssum+=ssq[b];

  sliceNorm(residuals,B,Orthog);
  for(int b=0;b<Nblock;b++){ assert(std::isnan(residuals[b])==0); }

  sliceNorm(residuals,X,Orthog);
  for(int b=0;b<Nblock;b++){ assert(std::isnan(residuals[b])==0); }

  /************************************************************************
   * Block conjugate gradient rQ (Sebastien Birk Thesis, after Dubrulle 2001)
   ************************************************************************
   * Dimensions:
   *
   *   X,B==(Nferm x Nblock)
   *   A==(Nferm x Nferm)
   *  
   * Nferm = Nspin x Ncolour x Ncomplex x Nlattice_site
   * 
   * QC = R = B-AX, D = Q     ; QC => Thin QR factorisation (google it)
   * for k: 
   *   Z  = AD
   *   M  = [D^dag Z]^{-1}
   *   X  = X + D MC
   *   QS = Q - ZM
   *   D  = Q + D S^dag
   *   C  = S C
   */
  ///////////////////////////////////////
  // Initial block: initial search dir is guess
  ///////////////////////////////////////
  std::cout << GridLogMessage<<"BlockCGrQ algorithm initialisation " <<std::endl;

  //1.  QC = R = B-AX, D = Q     ; QC => Thin QR factorisation (google it)

  Linop.HermOp(X, AD);
  tmp = B - AD;  
  //std::cout << GridLogMessage << " initial tmp " << norm2(tmp)<< std::endl;
  ThinQRfact (m_rr, m_C, m_Cinv, Q, tmp);
  //std::cout << GridLogMessage << " initial Q " << norm2(Q)<< std::endl;
  //std::cout << GridLogMessage << " m_rr " << m_rr<<std::endl;
  //std::cout << GridLogMessage << " m_C " << m_C<<std::endl;
  //std::cout << GridLogMessage << " m_Cinv " << m_Cinv<<std::endl;
  D=Q;

  std::cout << GridLogMessage<<"BlockCGrQ computed initial residual and QR fact " <<std::endl;

  ///////////////////////////////////////
  // Timers
  ///////////////////////////////////////
  GridStopWatch sliceInnerTimer;
  GridStopWatch sliceMaddTimer;
  GridStopWatch QRTimer;
  GridStopWatch MatrixTimer;
  GridStopWatch SolverTimer;
  SolverTimer.Start();

  int k;
  for (k = 1; k <= MaxIterations; k++) {

    //3. Z  = AD
    MatrixTimer.Start();
    Linop.HermOp(D, Z);      
    MatrixTimer.Stop();
    //std::cout << GridLogMessage << " norm2 Z " <<norm2(Z)<<std::endl;

    //4. M  = [D^dag Z]^{-1}
    sliceInnerTimer.Start();
    sliceInnerProductMatrix(m_DZ,D,Z,Orthog);
    sliceInnerTimer.Stop();
    m_M       = m_DZ.inverse();
    //std::cout << GridLogMessage << " m_DZ " <<m_DZ<<std::endl;
    
    //5. X  = X + D MC
    m_tmp     = m_M * m_C;
    sliceMaddTimer.Start();
    sliceMaddMatrix(X,m_tmp, D,X,Orthog);     
    sliceMaddTimer.Stop();

    //6. QS = Q - ZM
    sliceMaddTimer.Start();
    sliceMaddMatrix(tmp,m_M,Z,Q,Orthog,-1.0);
    sliceMaddTimer.Stop();
    QRTimer.Start();
    ThinQRfact (m_rr, m_S, m_Sinv, Q, tmp);
    QRTimer.Stop();
    
    //7. D  = Q + D S^dag
    m_tmp = m_S.adjoint();
    sliceMaddTimer.Start();
    sliceMaddMatrix(D,m_tmp,D,Q,Orthog);
    sliceMaddTimer.Stop();

    //8. C  = S C
    m_C = m_S*m_C;
    
    /*********************
     * convergence monitor
     *********************
     */
    m_rr = m_C.adjoint() * m_C;

    RealD max_resid=0;
    RealD rrsum=0;
    RealD rr;

    for(int b=0;b<Nblock;b++) {
      rrsum+=real(m_rr(b,b));
      rr = real(m_rr(b,b))/ssq[b];
      if ( rr > max_resid ) max_resid = rr;
    }

    std::cout << GridLogIterative << "\titeration "<<k<<" rr_sum "<<rrsum<<" ssq_sum "<< sssum
	      <<" ave "<<std::sqrt(rrsum/sssum) << " max "<< max_resid <<std::endl;

    if ( max_resid < Tolerance*Tolerance ) { 

      SolverTimer.Stop();

      std::cout << GridLogMessage<<"BlockCGrQ converged in "<<k<<" iterations"<<std::endl;

      for(int b=0;b<Nblock;b++){
	std::cout << GridLogMessage<< "\t\tblock "<<b<<" computed resid "
		  << std::sqrt(real(m_rr(b,b))/ssq[b])<<std::endl;
      }
      std::cout << GridLogMessage<<"\tMax residual is "<<std::sqrt(max_resid)<<std::endl;

      Linop.HermOp(X, AD);
      AD = AD-B;
      std::cout << GridLogMessage <<"\t True residual is " << std::sqrt(norm2(AD)/norm2(B)) <<std::endl;

      std::cout << GridLogMessage << "Time Breakdown "<<std::endl;
      std::cout << GridLogMessage << "\tElapsed    " << SolverTimer.Elapsed()     <<std::endl;
      std::cout << GridLogMessage << "\tMatrix     " << MatrixTimer.Elapsed()     <<std::endl;
      std::cout << GridLogMessage << "\tInnerProd  " << sliceInnerTimer.Elapsed() <<std::endl;
      std::cout << GridLogMessage << "\tMaddMatrix " << sliceMaddTimer.Elapsed()  <<std::endl;
      std::cout << GridLogMessage << "\tThinQRfact " << QRTimer.Elapsed()  <<std::endl;
	    
      IterationsToComplete = k;
      return;
    }

  }
  std::cout << GridLogMessage << "BlockConjugateGradient(rQ) did NOT converge" << std::endl;

  if (ErrorOnNoConverge) assert(0);
  IterationsToComplete = k;
}
//////////////////////////////////////////////////////////////////////////
// Block conjugate gradient; Original O'Leary Dimension zero should be the block direction
//////////////////////////////////////////////////////////////////////////
void BlockCGsolve(LinearOperatorBase<Field> &Linop, const Field &Src, Field &Psi) 
{
  int Orthog = blockDim; // First dimension is block dim; this is an assumption
  Nblock = Src._grid->_fdimensions[Orthog];

  std::cout<<GridLogMessage<<" Block Conjugate Gradient : Orthog "<<Orthog<<" Nblock "<<Nblock<<std::endl;

  Psi.checkerboard = Src.checkerboard;
  conformable(Psi, Src);

  Field P(Src);
  Field AP(Src);
  Field R(Src);
  
  Eigen::MatrixXcd m_pAp    = Eigen::MatrixXcd::Identity(Nblock,Nblock);
  Eigen::MatrixXcd m_pAp_inv= Eigen::MatrixXcd::Identity(Nblock,Nblock);
  Eigen::MatrixXcd m_rr     = Eigen::MatrixXcd::Zero(Nblock,Nblock);
  Eigen::MatrixXcd m_rr_inv = Eigen::MatrixXcd::Zero(Nblock,Nblock);

  Eigen::MatrixXcd m_alpha      = Eigen::MatrixXcd::Zero(Nblock,Nblock);
  Eigen::MatrixXcd m_beta   = Eigen::MatrixXcd::Zero(Nblock,Nblock);

  // Initial residual computation & set up
  std::vector<RealD> residuals(Nblock);
  std::vector<RealD> ssq(Nblock);

  sliceNorm(ssq,Src,Orthog);
  RealD sssum=0;
  for(int b=0;b<Nblock;b++) sssum+=ssq[b];

  sliceNorm(residuals,Src,Orthog);
  for(int b=0;b<Nblock;b++){ assert(std::isnan(residuals[b])==0); }

  sliceNorm(residuals,Psi,Orthog);
  for(int b=0;b<Nblock;b++){ assert(std::isnan(residuals[b])==0); }

  // Initial search dir is guess
  Linop.HermOp(Psi, AP);
  

  /************************************************************************
   * Block conjugate gradient (Stephen Pickles, thesis 1995, pp 71, O Leary 1980)
   ************************************************************************
   * O'Leary : R = B - A X
   * O'Leary : P = M R ; preconditioner M = 1
   * O'Leary : alpha = PAP^{-1} RMR
   * O'Leary : beta  = RMR^{-1}_old RMR_new
   * O'Leary : X=X+Palpha
   * O'Leary : R_new=R_old-AP alpha
   * O'Leary : P=MR_new+P beta
   */

  R = Src - AP;  
  P = R;
  sliceInnerProductMatrix(m_rr,R,R,Orthog);

  GridStopWatch sliceInnerTimer;
  GridStopWatch sliceMaddTimer;
  GridStopWatch MatrixTimer;
  GridStopWatch SolverTimer;
  SolverTimer.Start();

  int k;
  for (k = 1; k <= MaxIterations; k++) {

    RealD rrsum=0;
    for(int b=0;b<Nblock;b++) rrsum+=real(m_rr(b,b));

    std::cout << GridLogIterative << "\titeration "<<k<<" rr_sum "<<rrsum<<" ssq_sum "<< sssum
	      <<" / "<<std::sqrt(rrsum/sssum) <<std::endl;

    MatrixTimer.Start();
    Linop.HermOp(P, AP);
    MatrixTimer.Stop();

    // Alpha
    sliceInnerTimer.Start();
    sliceInnerProductMatrix(m_pAp,P,AP,Orthog);
    sliceInnerTimer.Stop();
    m_pAp_inv = m_pAp.inverse();
    m_alpha   = m_pAp_inv * m_rr ;

    // Psi, R update
    sliceMaddTimer.Start();
    sliceMaddMatrix(Psi,m_alpha, P,Psi,Orthog);     // add alpha *  P to psi
    sliceMaddMatrix(R  ,m_alpha,AP,  R,Orthog,-1.0);// sub alpha * AP to resid
    sliceMaddTimer.Stop();

    // Beta
    m_rr_inv = m_rr.inverse();
    sliceInnerTimer.Start();
    sliceInnerProductMatrix(m_rr,R,R,Orthog);
    sliceInnerTimer.Stop();
    m_beta = m_rr_inv *m_rr;

    // Search update
    sliceMaddTimer.Start();
    sliceMaddMatrix(AP,m_beta,P,R,Orthog);
    sliceMaddTimer.Stop();
    P= AP;

    /*********************
     * convergence monitor
     *********************
     */
    RealD max_resid=0;
    RealD rr;
    for(int b=0;b<Nblock;b++){
      rr = real(m_rr(b,b))/ssq[b];
      if ( rr > max_resid ) max_resid = rr;
    }
    
    if ( max_resid < Tolerance*Tolerance ) { 

      SolverTimer.Stop();

      std::cout << GridLogMessage<<"BlockCG converged in "<<k<<" iterations"<<std::endl;
      for(int b=0;b<Nblock;b++){
	std::cout << GridLogMessage<< "\t\tblock "<<b<<" computed resid "
		  << std::sqrt(real(m_rr(b,b))/ssq[b])<<std::endl;
      }
      std::cout << GridLogMessage<<"\tMax residual is "<<std::sqrt(max_resid)<<std::endl;

      Linop.HermOp(Psi, AP);
      AP = AP-Src;
      std::cout << GridLogMessage <<"\t True residual is " << std::sqrt(norm2(AP)/norm2(Src)) <<std::endl;

      std::cout << GridLogMessage << "Time Breakdown "<<std::endl;
      std::cout << GridLogMessage << "\tElapsed    " << SolverTimer.Elapsed()     <<std::endl;
      std::cout << GridLogMessage << "\tMatrix     " << MatrixTimer.Elapsed()     <<std::endl;
      std::cout << GridLogMessage << "\tInnerProd  " << sliceInnerTimer.Elapsed() <<std::endl;
      std::cout << GridLogMessage << "\tMaddMatrix " << sliceMaddTimer.Elapsed()  <<std::endl;
	    
      IterationsToComplete = k;
      return;
    }

  }
  std::cout << GridLogMessage << "BlockConjugateGradient did NOT converge" << std::endl;

  if (ErrorOnNoConverge) assert(0);
  IterationsToComplete = k;
}
//////////////////////////////////////////////////////////////////////////
// multiRHS conjugate gradient. Dimension zero should be the block direction
// Use this for spread out across nodes
//////////////////////////////////////////////////////////////////////////
void CGmultiRHSsolve(LinearOperatorBase<Field> &Linop, const Field &Src, Field &Psi) 
{
  int Orthog = blockDim; // First dimension is block dim
  Nblock = Src._grid->_fdimensions[Orthog];

  std::cout<<GridLogMessage<<"MultiRHS Conjugate Gradient : Orthog "<<Orthog<<" Nblock "<<Nblock<<std::endl;

  Psi.checkerboard = Src.checkerboard;
  conformable(Psi, Src);

  Field P(Src);
  Field AP(Src);
  Field R(Src);
  
  std::vector<ComplexD> v_pAp(Nblock);
  std::vector<RealD> v_rr (Nblock);
  std::vector<RealD> v_rr_inv(Nblock);
  std::vector<RealD> v_alpha(Nblock);
  std::vector<RealD> v_beta(Nblock);

  // Initial residual computation & set up
  std::vector<RealD> residuals(Nblock);
  std::vector<RealD> ssq(Nblock);

  sliceNorm(ssq,Src,Orthog);
  RealD sssum=0;
  for(int b=0;b<Nblock;b++) sssum+=ssq[b];

  sliceNorm(residuals,Src,Orthog);
  for(int b=0;b<Nblock;b++){ assert(std::isnan(residuals[b])==0); }

  sliceNorm(residuals,Psi,Orthog);
  for(int b=0;b<Nblock;b++){ assert(std::isnan(residuals[b])==0); }

  // Initial search dir is guess
  Linop.HermOp(Psi, AP);

  R = Src - AP;  
  P = R;
  sliceNorm(v_rr,R,Orthog);

  GridStopWatch sliceInnerTimer;
  GridStopWatch sliceMaddTimer;
  GridStopWatch sliceNormTimer;
  GridStopWatch MatrixTimer;
  GridStopWatch SolverTimer;

  SolverTimer.Start();
  int k;
  for (k = 1; k <= MaxIterations; k++) {

    RealD rrsum=0;
    for(int b=0;b<Nblock;b++) rrsum+=real(v_rr[b]);

    std::cout << GridLogIterative << "\titeration "<<k<<" rr_sum "<<rrsum<<" ssq_sum "<< sssum
	      <<" / "<<std::sqrt(rrsum/sssum) <<std::endl;

    MatrixTimer.Start();
    Linop.HermOp(P, AP);
    MatrixTimer.Stop();

    // Alpha
    sliceInnerTimer.Start();
    sliceInnerProductVector(v_pAp,P,AP,Orthog);
    sliceInnerTimer.Stop();
    for(int b=0;b<Nblock;b++){
      v_alpha[b] = v_rr[b]/real(v_pAp[b]);
    }

    // Psi, R update
    sliceMaddTimer.Start();
    sliceMaddVector(Psi,v_alpha, P,Psi,Orthog);     // add alpha *  P to psi
    sliceMaddVector(R  ,v_alpha,AP,  R,Orthog,-1.0);// sub alpha * AP to resid
    sliceMaddTimer.Stop();

    // Beta
    for(int b=0;b<Nblock;b++){
      v_rr_inv[b] = 1.0/v_rr[b];
    }
    sliceNormTimer.Start();
    sliceNorm(v_rr,R,Orthog);
    sliceNormTimer.Stop();
    for(int b=0;b<Nblock;b++){
      v_beta[b] = v_rr_inv[b] *v_rr[b];
    }

    // Search update
    sliceMaddTimer.Start();
    sliceMaddVector(P,v_beta,P,R,Orthog);
    sliceMaddTimer.Stop();

    /*********************
     * convergence monitor
     *********************
     */
    RealD max_resid=0;
    for(int b=0;b<Nblock;b++){
      RealD rr = v_rr[b]/ssq[b];
      if ( rr > max_resid ) max_resid = rr;
    }
    
    if ( max_resid < Tolerance*Tolerance ) { 

      SolverTimer.Stop();

      std::cout << GridLogMessage<<"MultiRHS solver converged in " <<k<<" iterations"<<std::endl;
      for(int b=0;b<Nblock;b++){
	std::cout << GridLogMessage<< "\t\tBlock "<<b<<" computed resid "<< std::sqrt(v_rr[b]/ssq[b])<<std::endl;
      }
      std::cout << GridLogMessage<<"\tMax residual is "<<std::sqrt(max_resid)<<std::endl;

      Linop.HermOp(Psi, AP);
      AP = AP-Src;
      std::cout <<GridLogMessage << "\tTrue residual is " << std::sqrt(norm2(AP)/norm2(Src)) <<std::endl;

      std::cout << GridLogMessage << "Time Breakdown "<<std::endl;
      std::cout << GridLogMessage << "\tElapsed    " << SolverTimer.Elapsed()     <<std::endl;
      std::cout << GridLogMessage << "\tMatrix     " << MatrixTimer.Elapsed()     <<std::endl;
      std::cout << GridLogMessage << "\tInnerProd  " << sliceInnerTimer.Elapsed() <<std::endl;
      std::cout << GridLogMessage << "\tNorm       " << sliceNormTimer.Elapsed() <<std::endl;
      std::cout << GridLogMessage << "\tMaddMatrix " << sliceMaddTimer.Elapsed()  <<std::endl;


      IterationsToComplete = k;
      return;
    }

  }
  std::cout << GridLogMessage << "MultiRHSConjugateGradient did NOT converge" << std::endl;

  if (ErrorOnNoConverge) assert(0);
  IterationsToComplete = k;
}

};

}
#endif
