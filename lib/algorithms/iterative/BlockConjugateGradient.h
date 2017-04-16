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

//////////////////////////////////////////////////////////////////////////
// Block conjugate gradient. Dimension zero should be the block direction
//////////////////////////////////////////////////////////////////////////
template <class Field>
class BlockConjugateGradient : public OperatorFunction<Field> {
 public:

  typedef typename Field::scalar_type scomplex;

  const int blockDim = 0;

  int Nblock;
  bool ErrorOnNoConverge;  // throw an assert when the CG fails to converge.
                           // Defaults true.
  RealD Tolerance;
  Integer MaxIterations;
  Integer IterationsToComplete; //Number of iterations the CG took to finish. Filled in upon completion
  
  BlockConjugateGradient(RealD tol, Integer maxit, bool err_on_no_conv = true)
    : Tolerance(tol),
    MaxIterations(maxit),
    ErrorOnNoConverge(err_on_no_conv){};

void operator()(LinearOperatorBase<Field> &Linop, const Field &Src, Field &Psi) 
{
  int Orthog = 0; // First dimension is block dim
  Nblock = Src._grid->_fdimensions[Orthog];
  std::cout<<GridLogMessage<<" Block Conjugate Gradient : Orthog "<<Orthog<<std::endl;
  std::cout<<GridLogMessage<<" Block Conjugate Gradient : Nblock "<<Nblock<<std::endl;

  Psi.checkerboard = Src.checkerboard;
  conformable(Psi, Src);

  Field P(Src);
  Field AP(Src);
  Field R(Src);
  
  GridStopWatch LinalgTimer;
  GridStopWatch MatrixTimer;
  GridStopWatch SolverTimer;

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

  int k;
  for (k = 1; k <= MaxIterations; k++) {

    RealD rrsum=0;
    for(int b=0;b<Nblock;b++) rrsum+=real(m_rr(b,b));

    std::cout << GridLogIterative << " iteration "<<k<<" rr_sum "<<rrsum<<" ssq_sum "<< sssum
	      <<" / "<<std::sqrt(rrsum/sssum) <<std::endl;

    Linop.HermOp(P, AP);

    // Alpha
    sliceInnerProductMatrix(m_pAp,P,AP,Orthog);
    m_pAp_inv = m_pAp.inverse();
    m_alpha   = m_pAp_inv * m_rr ;

    // Psi, R update
    sliceMaddMatrix(Psi,m_alpha, P,Psi,Orthog);     // add alpha *  P to psi
    sliceMaddMatrix(R  ,m_alpha,AP,  R,Orthog,-1.0);// sub alpha * AP to resid

    // Beta
    m_rr_inv = m_rr.inverse();
    sliceInnerProductMatrix(m_rr,R,R,Orthog);
    m_beta = m_rr_inv *m_rr;

    // Search update
    sliceMaddMatrix(AP,m_beta,P,R,Orthog);
    P= AP;

    /*********************
     * convergence monitor
     *********************
     */
    RealD max_resid=0;
    for(int b=0;b<Nblock;b++){
      RealD rr = real(m_rr(b,b))/ssq[b];
      if ( rr > max_resid ) max_resid = rr;
    }
    
    if ( max_resid < Tolerance*Tolerance ) { 

      std::cout << GridLogMessage<<" Block solver has converged in "
		<<k<<" iterations; max residual is "<<std::sqrt(max_resid)<<std::endl;

      for(int b=0;b<Nblock;b++){
	std::cout << GridLogMessage<< " block "<<b<<" resid "<< std::sqrt(real(m_rr(b,b))/ssq[b])<<std::endl;
      }

      Linop.HermOp(Psi, AP);
      AP = AP-Src;
      std::cout << " Block solver true residual is " << std::sqrt(norm2(AP)/norm2(Src)) <<std::endl;
      IterationsToComplete = k;
      return;
    }

  }
  std::cout << GridLogMessage << "BlockConjugateGradient did NOT converge" << std::endl;

  if (ErrorOnNoConverge) assert(0);
  IterationsToComplete = k;
}
};


//////////////////////////////////////////////////////////////////////////
// multiRHS conjugate gradient. Dimension zero should be the block direction
//////////////////////////////////////////////////////////////////////////
template <class Field>
class MultiRHSConjugateGradient : public OperatorFunction<Field> {
 public:

  typedef typename Field::scalar_type scomplex;

  const int blockDim = 0;

  int Nblock;
  bool ErrorOnNoConverge;  // throw an assert when the CG fails to converge.
                           // Defaults true.
  RealD Tolerance;
  Integer MaxIterations;
  Integer IterationsToComplete; //Number of iterations the CG took to finish. Filled in upon completion
  
   MultiRHSConjugateGradient(RealD tol, Integer maxit, bool err_on_no_conv = true)
    : Tolerance(tol),
    MaxIterations(maxit),
    ErrorOnNoConverge(err_on_no_conv){};

void operator()(LinearOperatorBase<Field> &Linop, const Field &Src, Field &Psi) 
{
  int Orthog = 0; // First dimension is block dim
  Nblock = Src._grid->_fdimensions[Orthog];
  std::cout<<GridLogMessage<<" MultiRHS Conjugate Gradient : Orthog "<<Orthog<<std::endl;
  std::cout<<GridLogMessage<<" MultiRHS Conjugate Gradient : Nblock "<<Nblock<<std::endl;

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

  int k;
  for (k = 1; k <= MaxIterations; k++) {

    RealD rrsum=0;
    for(int b=0;b<Nblock;b++) rrsum+=real(v_rr[b]);

    std::cout << GridLogIterative << " iteration "<<k<<" rr_sum "<<rrsum<<" ssq_sum "<< sssum
	      <<" / "<<std::sqrt(rrsum/sssum) <<std::endl;

    Linop.HermOp(P, AP);

    // Alpha
    sliceInnerProductVector(v_pAp,P,AP,Orthog);
    for(int b=0;b<Nblock;b++){
      v_alpha[b] = v_rr[b]/real(v_pAp[b]);
    }

    // Psi, R update
    sliceMaddVector(Psi,v_alpha, P,Psi,Orthog);     // add alpha *  P to psi
    sliceMaddVector(R  ,v_alpha,AP,  R,Orthog,-1.0);// sub alpha * AP to resid

    // Beta
    for(int b=0;b<Nblock;b++){
      v_rr_inv[b] = 1.0/v_rr[b];
    }
    sliceNorm(v_rr,R,Orthog);
    for(int b=0;b<Nblock;b++){
      v_beta[b] = v_rr_inv[b] *v_rr[b];
    }

    // Search update
    sliceMaddVector(P,v_beta,P,R,Orthog);

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
      std::cout << GridLogMessage<<" MultiRHS solver has converged in "
		<<k<<" iterations; max residual is "<<std::sqrt(max_resid)<<std::endl;
      for(int b=0;b<Nblock;b++){
	std::cout << GridLogMessage<< " block "<<b<<" resid "<< std::sqrt(v_rr[b]/ssq[b])<<std::endl;
      }

      Linop.HermOp(Psi, AP);
      AP = AP-Src;
      std::cout << " MultiRHS solver true residual is " << std::sqrt(norm2(AP)/norm2(Src)) <<std::endl;
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
