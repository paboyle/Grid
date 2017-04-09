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

#include <Grid/Eigen/Dense>

namespace Grid {

GridBase         *makeSubSliceGrid(const GridBase *BlockSolverGrid,int Orthog)
{
  int NN    = BlockSolverGrid->_ndimension;
  int nsimd = BlockSolverGrid->Nsimd();

  std::vector<int> latt_phys(0);
  std::vector<int> simd_phys(0);
  std::vector<int>  mpi_phys(0);
  
  for(int d=0;d<NN;d++){
    if( d!=Orthog ) { 
    latt_phys.push_back(BlockSolverGrid->_fdimensions[d]);
    simd_phys.push_back(BlockSolverGrid->_simd_layout[d]);
     mpi_phys.push_back(BlockSolverGrid->_processors[d]);
    }
  }
  return (GridBase *)new GridCartesian(latt_phys,simd_phys,mpi_phys); 
}
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Need to move sliceInnerProduct, sliceAxpy, sliceNorm etc... into lattice sector along with sliceSum
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class vobj>
static void sliceMaddMatrix (Lattice<vobj> &R,Eigen::MatrixXcd &aa,const Lattice<vobj> &X,const Lattice<vobj> &Y,int Orthog,RealD scale=1.0) 
{    
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  int Nblock = X._grid->GlobalDimensions()[Orthog];
    
  GridBase *FullGrid  = X._grid;
  GridBase *SliceGrid = makeSubSliceGrid(FullGrid,Orthog);
  
  Lattice<vobj> Xslice(SliceGrid);
  Lattice<vobj> Rslice(SliceGrid);
  // If we based this on Cshift it would work for spread out
  // but it would be even slower
  for(int i=0;i<Nblock;i++){
    ExtractSlice(Rslice,Y,i,Orthog);
    for(int j=0;j<Nblock;j++){
      ExtractSlice(Xslice,X,j,Orthog);
      Rslice = Rslice + Xslice*(scale*aa(j,i));
    }
    InsertSlice(Rslice,R,i,Orthog);
  }
};
template<class vobj>
static void sliceInnerProductMatrix(  Eigen::MatrixXcd &mat, const Lattice<vobj> &lhs,const Lattice<vobj> &rhs,int Orthog) 
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  GridBase *FullGrid  = lhs._grid;
  GridBase *SliceGrid = makeSubSliceGrid(FullGrid,Orthog);

  int Nblock = FullGrid->GlobalDimensions()[Orthog];
  
  Lattice<vobj> Lslice(SliceGrid);
  Lattice<vobj> Rslice(SliceGrid);

  mat = Eigen::MatrixXcd::Zero(Nblock,Nblock);

  for(int i=0;i<Nblock;i++){
    ExtractSlice(Lslice,lhs,i,Orthog);
    for(int j=0;j<Nblock;j++){
      ExtractSlice(Rslice,rhs,j,Orthog);
      mat(i,j) = innerProduct(Lslice,Rslice);
    }
  }
#undef FORCE_DIAG
#ifdef FORCE_DIAG
  for(int i=0;i<Nblock;i++){
    for(int j=0;j<Nblock;j++){
      if ( i != j ) mat(i,j)=0.0;
    }
  }
#endif
  return;
}
template<class vobj>
static void sliceInnerProductVector( std::vector<ComplexD> & vec, const Lattice<vobj> &lhs,const Lattice<vobj> &rhs,int Orthog) 
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::tensor_reduced scalar;
  typedef typename scalar::scalar_object  scomplex;
  
  int Nblock = lhs._grid->GlobalDimensions()[Orthog];

  vec.resize(Nblock);
  std::vector<scomplex> sip(Nblock);
  Lattice<scalar> IP(lhs._grid); 

  IP=localInnerProduct(lhs,rhs);
  sliceSum(IP,sip,Orthog);
  
  for(int ss=0;ss<Nblock;ss++){
    vec[ss] = TensorRemove(sip[ss]);
  }
}
template<class vobj>
static void sliceNorm (std::vector<RealD> &sn,const Lattice<vobj> &rhs,int Orthog) {

  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  
  int Nblock = rhs._grid->GlobalDimensions()[Orthog];
  std::vector<ComplexD> ip(Nblock);
  sn.resize(Nblock);
  
  sliceInnerProductVector(ip,rhs,rhs,Orthog);
  for(int ss=0;ss<Nblock;ss++){
    sn[ss] = real(ip[ss]);
  }
};
/*
template<class vobj>
static void sliceInnerProductMatrixOld(  Eigen::MatrixXcd &mat, const Lattice<vobj> &lhs,const Lattice<vobj> &rhs,int Orthog) 
{
  typedef typename vobj::scalar_object  sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::tensor_reduced scalar;
  typedef typename scalar::scalar_object  scomplex;

  int Nblock = lhs._grid->GlobalDimensions()[Orthog];

  std::cout << " sliceInnerProductMatrix Dim "<<Orthog<<" Nblock " << Nblock<<std::endl;

  Lattice<scalar> IP(lhs._grid); 
  std::vector<scomplex> sip(Nblock);
    
  mat = Eigen::MatrixXcd::Zero(Nblock,Nblock);

  Lattice<vobj> tmp = rhs;
  
  for(int s1=0;s1<Nblock;s1++){
    
    IP=localInnerProduct(lhs,tmp);
    sliceSum(IP,sip,Orthog);

    std::cout << "InnerProductMatrix ["<<s1<<"] = ";
    for(int ss=0;ss<Nblock;ss++){
      std::cout << TensorRemove(sip[ss])<<" ";
    }
    std::cout << std::endl;

    for(int ss=0;ss<Nblock;ss++){
      mat(ss,(s1+ss)%Nblock) = TensorRemove(sip[ss]);
    }
    if ( s1!=(Nblock-1) ) { 
      tmp = Cshift(tmp,Orthog,1);
    }
  }
}
*/

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
}
#endif
