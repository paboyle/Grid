    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/ImplicitlyRestartedLanczos.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Chulwoo Jung <chulwoo@bnl.gov>
Author: Christoph Lehner <clehner@bnl.gov>

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
#ifndef GRID_BIRL_H
#define GRID_BIRL_H

#include <string.h> //memset
//#include <zlib.h>
#include <sys/stat.h>

#include <Grid/algorithms/iterative/BlockImplicitlyRestartedLanczos/BlockedGrid.h>
#include <Grid/algorithms/iterative/BlockImplicitlyRestartedLanczos/FieldBasisVector.h>
#include <Grid/algorithms/iterative/BlockImplicitlyRestartedLanczos/BlockProjector.h>

namespace Grid { 

template<class Field>
void basisOrthogonalize(std::vector<Field> &basis,Field &w,int k) 
{
  for(int j=0; j<k; ++j){
    auto ip = innerProduct(basis[j],w);
    w = w - ip*basis[j];
  }
}

template<class Field>
void basisRotate(std::vector<Field> &basis,Eigen::MatrixXd& Qt,int j0, int j1, int k0,int k1,int Nm) 
{
  typedef typename Field::vector_object vobj;
  GridBase* grid = basis[0]._grid;
      
  parallel_region
  {
    std::vector < vobj > B(Nm); // Thread private
        
    parallel_for_internal(int ss=0;ss < grid->oSites();ss++){
      for(int j=j0; j<j1; ++j) B[j]=0.;
      
      for(int j=j0; j<j1; ++j){
	for(int k=k0; k<k1; ++k){
	  B[j] +=Qt(j,k) * basis[k]._odata[ss];
	}
      }
      for(int j=j0; j<j1; ++j){
	  basis[j]._odata[ss] = B[j];
      }
    }
  }
}

template<class Field>
void basisReorderInPlace(std::vector<Field> &_v,std::vector<RealD>& sort_vals, std::vector<int>& idx) 
{
  int vlen = idx.size();

  assert(vlen>=1);
  assert(vlen<=sort_vals.size());
  assert(vlen<=_v.size());

  for (size_t i=0;i<vlen;i++) {

    if (idx[i] != i) {

      assert(idx[i] > i);
      //////////////////////////////////////
      // idx[i] is a table of desired sources giving a permutation.
      //
      // Swap v[i] with v[idx[i]].
      //
      // Find  j>i for which _vnew[j] = _vold[i],
      // track the move idx[j] => idx[i]
      // track the move idx[i] => i
      //////////////////////////////////////
      size_t j;
      for (j=i;j<idx.size();j++)
	if (idx[j]==i)
	  break;

      assert(j!=idx.size());
      assert(idx[j]==i);

      std::swap(_v[i]._odata,_v[idx[i]]._odata); // should use vector move constructor, no data copy
      std::swap(sort_vals[i],sort_vals[idx[i]]);

      idx[j] = idx[i];
      idx[i] = i;
    }
  }
}

std::vector<int> basisSortGetIndex(std::vector<RealD>& sort_vals) 
{
  std::vector<int> idx(sort_vals.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(), [&sort_vals](int i1, int i2) {
    return ::fabs(sort_vals[i1]) < ::fabs(sort_vals[i2]);
  });
  return idx;
}

template<class Field>
void basisSortInPlace(std::vector<Field> & _v,std::vector<RealD>& sort_vals, bool reverse) 
{
  std::vector<int> idx = basisSortGetIndex(sort_vals);
  if (reverse)
    std::reverse(idx.begin(), idx.end());
  
  basisReorderInPlace(_v,sort_vals,idx);
}

// PAB: faster to compute the inner products first then fuse loops.
// If performance critical can improve.
template<class Field>
void basisDeflate(const std::vector<Field> &_v,const std::vector<RealD>& eval,const Field& src_orig,Field& result) {
  result = zero;
  assert(_v.size()==eval.size());
  int N = (int)_v.size();
  for (int i=0;i<N;i++) {
    Field& tmp = _v[i];
    axpy(result,TensorRemove(innerProduct(tmp,src_orig)) / eval[i],tmp,result);
  }
}

  /*  enum IRLdiagonalisation { 
    IRLdiagonaliseWithDSTEGR,
    IRLdiagonaliseWithQR,
    IRLdiagonaliseWithEigen
    };*/

/////////////////////////////////////////////////////////////
// Implicitly restarted lanczos
/////////////////////////////////////////////////////////////

template<class Field> 
class BlockImplicitlyRestartedLanczos {
 private:
  const RealD small = 1.0e-8;
  int MaxIter;
  int MinRestart; // Minimum number of restarts; only check for convergence after
  int Nstop;   // Number of evecs checked for convergence
  int Nk;      // Number of converged sought
  //  int Np;      // Np -- Number of spare vecs in krylov space //  == Nm - Nk
  int Nm;      // Nm -- total number of vectors
  IRLdiagonalisation diagonalisation;
  int orth_period;
    
  RealD OrthoTime;
  RealD eresid, betastp;
  ////////////////////////////////
  // Embedded objects
  ////////////////////////////////
  SortEigen<Field> _sort;
  LinearFunction<Field> &_HermOp;
  LinearFunction<Field> &_HermOpTest;
  /////////////////////////
  // Constructor
  /////////////////////////
public:       
 BlockImplicitlyRestartedLanczos(LinearFunction<Field> & HermOp,
				 LinearFunction<Field> & HermOpTest,
				 int _Nstop, // sought vecs
				 int _Nk, // sought vecs
				 int _Nm, // spare vecs
				 RealD _eresid, // resid in lmdue deficit 
				 RealD _betastp, // if beta(k) < betastp: converged
				 int _MaxIter, // Max iterations
				 int _MinRestart, int _orth_period = 1,
				 IRLdiagonalisation _diagonalisation= IRLdiagonaliseWithEigen) :
  _HermOp(HermOp),      _HermOpTest(HermOpTest),
    Nstop(_Nstop)  ,      Nk(_Nk),      Nm(_Nm),
    eresid(_eresid),      betastp(_betastp),
    MaxIter(_MaxIter)  ,      MinRestart(_MinRestart),
    orth_period(_orth_period), diagonalisation(_diagonalisation)  { };

  ////////////////////////////////
  // Helpers
  ////////////////////////////////
  template<typename T>  static RealD normalise(T& v) 
  {
    RealD nn = norm2(v);
    nn = sqrt(nn);
    v = v * (1.0/nn);
    return nn;
  }

  void orthogonalize(Field& w, BasisFieldVector<Field>& evec,int k)
  {
    OrthoTime-=usecond()/1e6;
    //evec.orthogonalize(w,k);
    basisOrthogonalize(evec._v,w,k);
    normalise(w);
    OrthoTime+=usecond()/1e6;
  }

/* Rudy Arthur's thesis pp.137
------------------------
Require: M > K P = M − K †
Compute the factorization AVM = VM HM + fM eM 
repeat
  Q=I
  for i = 1,...,P do
    QiRi =HM −θiI Q = QQi
    H M = Q †i H M Q i
  end for
  βK =HM(K+1,K) σK =Q(M,K)
  r=vK+1βK +rσK
  VK =VM(1:M)Q(1:M,1:K)
  HK =HM(1:K,1:K)
  →AVK =VKHK +fKe†K † Extend to an M = K + P step factorization AVM = VMHM + fMeM
until convergence
*/
  void calc(std::vector<RealD>& eval, BasisFieldVector<Field>& evec,  const Field& src, int& Nconv, bool reverse, int SkipTest)
  {
    GridBase *grid = src._grid;
    assert(grid == evec[0]._grid);
    
    GridLogIRL.TimingMode(1);
    std::cout << GridLogIRL <<"**************************************************************************"<< std::endl;
    std::cout << GridLogIRL <<" ImplicitlyRestartedLanczos::calc() starting iteration 0 /  "<< MaxIter<< std::endl;
    std::cout << GridLogIRL <<"**************************************************************************"<< std::endl;
    std::cout << GridLogIRL <<" -- seek   Nk    = " << Nk    <<" vectors"<< std::endl;
    std::cout << GridLogIRL <<" -- accept Nstop = " << Nstop <<" vectors"<< std::endl;
    std::cout << GridLogIRL <<" -- total  Nm    = " << Nm    <<" vectors"<< std::endl;
    std::cout << GridLogIRL <<" -- size of eval = " << eval.size() << std::endl;
    std::cout << GridLogIRL <<" -- size of evec = " << evec.size() << std::endl;
    if ( diagonalisation == IRLdiagonaliseWithDSTEGR ) {
      std::cout << GridLogIRL << "Diagonalisation is DSTEGR "<<std::endl;
    } else if ( diagonalisation == IRLdiagonaliseWithQR ) { 
      std::cout << GridLogIRL << "Diagonalisation is QR "<<std::endl;
    }  else if ( diagonalisation == IRLdiagonaliseWithEigen ) { 
      std::cout << GridLogIRL << "Diagonalisation is Eigen "<<std::endl;
    }
    std::cout << GridLogIRL <<"**************************************************************************"<< std::endl;
	
    assert(Nm <= evec.size() && Nm <= eval.size());
    
    // quickly get an idea of the largest eigenvalue to more properly normalize the residuum
    RealD evalMaxApprox = 0.0;
    {
      auto src_n = src;
      auto tmp = src;
      const int _MAX_ITER_IRL_MEVAPP_ = 50;
      for (int i=0;i<_MAX_ITER_IRL_MEVAPP_;i++) {
	_HermOpTest(src_n,tmp);
	RealD vnum = real(innerProduct(src_n,tmp)); // HermOp.
	RealD vden = norm2(src_n);
	RealD na = vnum/vden;
	if (fabs(evalMaxApprox/na - 1.0) < 0.05)
	  i=_MAX_ITER_IRL_MEVAPP_;
	evalMaxApprox = na;
	std::cout << GridLogIRL << " Approximation of largest eigenvalue: " << evalMaxApprox << std::endl;
	src_n = tmp;
      }
    }
	
    std::vector<RealD> lme(Nm);  
    std::vector<RealD> lme2(Nm);
    std::vector<RealD> eval2(Nm);
    std::vector<RealD> eval2_copy(Nm);
    Eigen::MatrixXd Qt = Eigen::MatrixXd::Zero(Nm,Nm);

    Field f(grid);
    Field v(grid);
    int k1 = 1;
    int k2 = Nk;
    RealD beta_k;

    Nconv = 0;
  
    // Set initial vector
    evec[0] = src;
    normalise(evec[0]);
	
    // Initial Nk steps
    OrthoTime=0.;
    for(int k=0; k<Nk; ++k) step(eval,lme,evec,f,Nm,k);
    std::cout<<GridLogIRL <<"Initial "<< Nk <<"steps done "<<std::endl;
    std::cout<<GridLogIRL <<"Initial steps:OrthoTime "<<OrthoTime<< "seconds"<<std::endl;

    //////////////////////////////////
    // Restarting loop begins
    //////////////////////////////////
    int iter;
    for(iter = 0; iter<MaxIter; ++iter){
      
      OrthoTime=0.;

      std::cout<< GridLogMessage <<" **********************"<< std::endl;
      std::cout<< GridLogMessage <<" Restart iteration = "<< iter << std::endl;
      std::cout<< GridLogMessage <<" **********************"<< std::endl;

      std::cout<<GridLogIRL <<" running "<<Nm-Nk <<" steps: "<<std::endl;
      for(int k=Nk; k<Nm; ++k) step(eval,lme,evec,f,Nm,k);
      f *= lme[Nm-1];

      std::cout<<GridLogIRL <<" "<<Nm-Nk <<" steps done "<<std::endl;
      std::cout<<GridLogIRL <<"Initial steps:OrthoTime "<<OrthoTime<< "seconds"<<std::endl;
	  
      //////////////////////////////////
      // getting eigenvalues
      //////////////////////////////////
      for(int k=0; k<Nm; ++k){
	eval2[k] = eval[k+k1-1];
	lme2[k] = lme[k+k1-1];
      }
      Qt = Eigen::MatrixXd::Identity(Nm,Nm);
      diagonalize(eval2,lme2,Nm,Nm,Qt,grid);
      std::cout<<GridLogIRL <<" diagonalized "<<std::endl;

      //////////////////////////////////
      // sorting
      //////////////////////////////////
      eval2_copy = eval2;

      _sort.push(eval2,Nm);

      std::cout<<GridLogIRL <<" evals sorted "<<std::endl;
      for(int ip=0; ip<k2; ++ip) std::cout<<GridLogIRL << "eval "<< ip << " "<< eval2[ip] << std::endl;

      //////////////////////////////////
      // Implicitly shifted QR transformations
      //////////////////////////////////
      Qt = Eigen::MatrixXd::Identity(Nm,Nm);
      std::cout<<GridLogIRL << "QR decompose " << std::endl;
      for(int ip=k2; ip<Nm; ++ip){ 
	QR_decomp(eval,lme,Nm,Nm,Qt,eval2[ip],k1,Nm);
      }
      std::cout<<GridLogIRL <<"QR decompose done "<<std::endl;

      assert(k2<Nm);
      assert(k2<Nm);
      assert(k1>0);
      //      evec.rotate(Qt,k1-1,k2+1,0,Nm,Nm); /// big constraint on the basis
      basisRotate(evec._v,Qt,k1-1,k2+1,0,Nm,Nm); /// big constraint on the basis

      std::cout<<GridLogIRL <<"QR rotation done "<<std::endl;
      
      ////////////////////////////////////////////////////
      // Compressed vector f and beta(k2)
      ////////////////////////////////////////////////////
      f *= Qt(k2-1,Nm-1);
      f += lme[k2-1] * evec[k2];
      beta_k = norm2(f);
      beta_k = sqrt(beta_k);
      std::cout<<GridLogIRL<<" beta(k) = "<<beta_k<<std::endl;
	  
      RealD betar = 1.0/beta_k;
      evec[k2] = betar * f;
      lme[k2-1] = beta_k;
	  
      ////////////////////////////////////////////////////
      // Convergence test
      ////////////////////////////////////////////////////
      for(int k=0; k<Nm; ++k){    
	eval2[k] = eval[k];
	lme2[k] = lme[k];
	//	std::cout<<GridLogIRL << "eval2[" << k << "] = " << eval2[k] << std::endl;
      }
      Qt = Eigen::MatrixXd::Identity(Nm,Nm);
      diagonalize(eval2,lme2,Nk,Nm,Qt,grid);
      std::cout<<GridLogIRL <<" Diagonalized "<<std::endl;
	  
      Nconv = 0;
      if (iter >= MinRestart) {
	std::cout << GridLogIRL << "Rotation to test convergence " << std::endl;
	
	Field ev0_orig(grid);
	ev0_orig = evec[0];
	    
	//	evec.rotate(Qt,0,Nk,0,Nk,Nm);
	basisRotate(evec._v,Qt,0,Nk,0,Nk,Nm);	    

	{
	  std::cout << GridLogIRL << "Test convergence" << std::endl;
	  Field B(grid);
	      
	  for(int j = 0; j<Nk; j+=SkipTest){
	    B=evec[j];

	    //std::cout << "Checkerboard: " << evec[j].checkerboard << std::endl; 
	    B.checkerboard = evec[0].checkerboard;

	    _HermOpTest(B,v);
		
	    RealD vnum = real(innerProduct(B,v)); // HermOp.
	    RealD vden = norm2(B);
	    RealD vv0 = norm2(v);
	    eval2[j] = vnum/vden;
	    v -= eval2[j]*B;
	    RealD vv = norm2(v) / ::pow(evalMaxApprox,2.0);
	    std::cout.precision(13);
	    std::cout<<GridLogIRL << "[" << std::setw(3)<< std::setiosflags(std::ios_base::right) <<j<<"] "
		     <<"eval = "<<std::setw(25)<< std::setiosflags(std::ios_base::left)<< eval2[j] << " (" << eval2_copy[j] << ")"
		     <<" |H B[i] - eval[i]B[i]|^2 / evalMaxApprox^2 " << std::setw(25)<< std::setiosflags(std::ios_base::right)<< vv
		     <<" "<< vnum/(sqrt(vden)*sqrt(vv0))
		     <<std::endl;
		
	    // change the criteria as evals are supposed to be sorted, all evals smaller(larger) than Nstop should have converged
	    if((vv<eresid*eresid) && (j == Nconv) ){
	      Nconv+=SkipTest;
	    }
	  }
	      
	  // test if we converged, if so, terminate
	  std::cout<<GridLogIRL<<" #modes converged: "<<Nconv<<std::endl;
	  if( Nconv>=Nstop || beta_k < betastp){
	    goto converged;
	  }
	      
	  std::cout << GridLogIRL << "Convergence testing: Rotating back" << std::endl;
	  //B[j] +=Qt[k+_Nm*j] * _v[k]._odata[ss];
	  {
	    Eigen::MatrixXd qm = Eigen::MatrixXd::Zero(Nk,Nk); // Restrict Qt to Nk x Nk
	    for (int k=0;k<Nk;k++)
	      for (int j=0;j<Nk;j++)
		qm(j,k) = Qt(j,k);

	    Eigen::MatrixXd qmI = qm.inverse();
	    std::cout << GridLogIRL << "Inverted ("<<Nk<<"x"<<Nk<<") matrix " << std::endl;

	    
	    RealD res_check_rotate_inverse = (qm*qmI - Eigen::MatrixXd::Identity(Nk,Nk)).norm(); // sqrt( |X|^2 )
	    assert(res_check_rotate_inverse < 1e-7);
	    //evec.rotate(qmI,0,Nk,0,Nk,Nm);
	    basisRotate(evec._v,qmI,0,Nk,0,Nk,Nm);
		
	    axpy(ev0_orig,-1.0,evec[0],ev0_orig);
	    std::cout << GridLogIRL << "Rotation done ; error = " << res_check_rotate_inverse << ");"<<std::endl;
	    std::cout << GridLogIRL << " | evec[0] - evec[0]_orig | = "    << ::sqrt(norm2(ev0_orig)) << std::endl;
	  }
	}
      } else {
	std::cout << GridLogIRL << "iter < MinRestart: do not yet test for convergence\n";
      } // end of iter loop
    }

    std::cout<<GridLogError<<"\n NOT converged.\n";
    abort();
	
  converged:

    if (SkipTest == 1) {
      eval = eval2;
    } else {
      // test quickly
      for (int j=0;j<Nstop;j+=SkipTest) {
	std::cout<<GridLogIRL << "Eigenvalue[" << j << "] = " << eval2[j] << " (" << eval2_copy[j] << ")" << std::endl;
      }
      eval2_copy.resize(eval2.size());
      eval = eval2_copy;
    }
    //    evec.sortInPlace(eval,reverse);
    basisSortInPlace(evec._v,eval,reverse);
    // test // PAB -- what does this test ?
    for (int j=0;j<Nstop;j++) {
      std::cout<<GridLogIRL << " |e[" << j << "]|^2 = " << norm2(evec[j]) << std::endl;
    }
       
    std::cout << GridLogIRL <<"**************************************************************************"<< std::endl;
    std::cout << GridLogIRL << "ImplicitlyRestartedLanczos CONVERGED ; Summary :\n";
    std::cout << GridLogIRL <<"**************************************************************************"<< std::endl;
    std::cout << GridLogIRL << " -- Iterations  = "<< iter   << "\n";
    std::cout << GridLogIRL << " -- beta(k)     = "<< beta_k << "\n";
    std::cout << GridLogIRL << " -- Nconv       = "<< Nconv  << "\n";
    std::cout << GridLogIRL <<"**************************************************************************"<< std::endl;
  }

 private:
/* Saad PP. 195
1. Choose an initial vector v1 of 2-norm unity. Set β1 ≡ 0, v0 ≡ 0
2. For k = 1,2,...,m Do:
3. wk:=Avk−βkv_{k−1}      
4. αk:=(wk,vk)       // 
5. wk:=wk−αkvk       // wk orthog vk 
6. βk+1 := ∥wk∥2. If βk+1 = 0 then Stop
7. vk+1 := wk/βk+1
8. EndDo
 */
  void step(std::vector<RealD>& lmd,
	    std::vector<RealD>& lme, 
	    BasisFieldVector<Field>& evec,
	    Field& w,int Nm,int k)
  {
    const RealD tiny = 1.0e-20;
    assert( k< Nm );

    GridStopWatch gsw_op,gsw_o;

    Field& evec_k = evec[k];

    _HermOp(evec_k,w);
    std::cout<<GridLogIRL << "_HermOp (poly)" <<std::endl;

    if(k>0) w -= lme[k-1] * evec[k-1];

    ComplexD zalph = innerProduct(evec_k,w); // 4. αk:=(wk,vk)
    RealD     alph = real(zalph);

    w = w - alph * evec_k;// 5. wk:=wk−αkvk

    RealD beta = normalise(w); // 6. βk+1 := ∥wk∥2. If βk+1 = 0 then Stop
    // 7. vk+1 := wk/βk+1

    lmd[k] = alph;
    lme[k] = beta;

    std::cout<<GridLogIRL << "linalg " <<std::endl;

    if (k>0 && k % orth_period == 0) {
      orthogonalize(w,evec,k); // orthonormalise
      std::cout<<GridLogIRL << "orthogonalised " <<std::endl;
    }

    if(k < Nm-1) evec[k+1] = w;

    std::cout<<GridLogIRL << "alpha[" << k << "] = " << zalph << " beta[" << k << "] = "<<beta<<std::endl;
    if ( beta < tiny ) 
      std::cout<<GridLogIRL << " beta is tiny "<<beta<<std::endl;
  }

  void diagonalize_Eigen(std::vector<RealD>& lmd, std::vector<RealD>& lme, 
			 int Nk, int Nm,  
			 Eigen::MatrixXd & Qt, // Nm x Nm
			 GridBase *grid)
  {
    Eigen::MatrixXd TriDiag = Eigen::MatrixXd::Zero(Nk,Nk);

    for(int i=0;i<Nk;i++)   TriDiag(i,i)   = lmd[i];
    for(int i=0;i<Nk-1;i++) TriDiag(i,i+1) = lme[i];
    for(int i=0;i<Nk-1;i++) TriDiag(i+1,i) = lme[i];
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(TriDiag);

    for (int i = 0; i < Nk; i++) {
      lmd[Nk-1-i] = eigensolver.eigenvalues()(i);
    }
    for (int i = 0; i < Nk; i++) {
      for (int j = 0; j < Nk; j++) {
	Qt(Nk-1-i,j) = eigensolver.eigenvectors()(j,i);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // File could end here if settle on Eigen ???
  ///////////////////////////////////////////////////////////////////////////

  void QR_decomp(std::vector<RealD>& lmd,   // Nm 
		 std::vector<RealD>& lme,   // Nm 
		 int Nk, int Nm,            // Nk, Nm
		 Eigen::MatrixXd& Qt,       // Nm x Nm matrix
		 RealD Dsh, int kmin, int kmax)
  {
    int k = kmin-1;
    RealD x;
    
    RealD Fden = 1.0/hypot(lmd[k]-Dsh,lme[k]);
    RealD c = ( lmd[k] -Dsh) *Fden;
    RealD s = -lme[k] *Fden;
      
    RealD tmpa1 = lmd[k];
    RealD tmpa2 = lmd[k+1];
    RealD tmpb  = lme[k];

    lmd[k]   = c*c*tmpa1 +s*s*tmpa2 -2.0*c*s*tmpb;
    lmd[k+1] = s*s*tmpa1 +c*c*tmpa2 +2.0*c*s*tmpb;
    lme[k]   = c*s*(tmpa1-tmpa2) +(c*c-s*s)*tmpb;
    x        =-s*lme[k+1];
    lme[k+1] = c*lme[k+1];
      
    for(int i=0; i<Nk; ++i){
      RealD Qtmp1 = Qt(k,i);
      RealD Qtmp2 = Qt(k+1,i);
      Qt(k,i)  = c*Qtmp1 - s*Qtmp2;
      Qt(k+1,i)= s*Qtmp1 + c*Qtmp2; 
    }

    // Givens transformations
    for(int k = kmin; k < kmax-1; ++k){
      
      RealD Fden = 1.0/hypot(x,lme[k-1]);
      RealD c = lme[k-1]*Fden;
      RealD s = - x*Fden;
	
      RealD tmpa1 = lmd[k];
      RealD tmpa2 = lmd[k+1];
      RealD tmpb  = lme[k];

      lmd[k]   = c*c*tmpa1 +s*s*tmpa2 -2.0*c*s*tmpb;
      lmd[k+1] = s*s*tmpa1 +c*c*tmpa2 +2.0*c*s*tmpb;
      lme[k]   = c*s*(tmpa1-tmpa2) +(c*c-s*s)*tmpb;
      lme[k-1] = c*lme[k-1] -s*x;

      if(k != kmax-2){
	x = -s*lme[k+1];
	lme[k+1] = c*lme[k+1];
      }

      for(int i=0; i<Nk; ++i){
	RealD Qtmp1 = Qt(k,i);
	RealD Qtmp2 = Qt(k+1,i);
	Qt(k,i)     = c*Qtmp1 -s*Qtmp2;
	Qt(k+1,i)   = s*Qtmp1 +c*Qtmp2;
      }
    }
  }

  void diagonalize(std::vector<RealD>& lmd, std::vector<RealD>& lme, 
		   int Nk, int Nm,   
		   Eigen::MatrixXd & Qt,
		   GridBase *grid)
  {
    Qt = Eigen::MatrixXd::Identity(Nm,Nm);
    if ( diagonalisation == IRLdiagonaliseWithDSTEGR ) {
      diagonalize_lapack(lmd,lme,Nk,Nm,Qt,grid);
    } else if ( diagonalisation == IRLdiagonaliseWithQR ) { 
      diagonalize_QR(lmd,lme,Nk,Nm,Qt,grid);
    }  else if ( diagonalisation == IRLdiagonaliseWithEigen ) { 
      diagonalize_Eigen(lmd,lme,Nk,Nm,Qt,grid);
    } else { 
      assert(0);
    }
  }

#ifdef USE_LAPACK
void LAPACK_dstegr(char *jobz, char *range, int *n, double *d, double *e,
                   double *vl, double *vu, int *il, int *iu, double *abstol,
                   int *m, double *w, double *z, int *ldz, int *isuppz,
                   double *work, int *lwork, int *iwork, int *liwork,
                   int *info);
#endif

void diagonalize_lapack(std::vector<RealD>& lmd,
			std::vector<RealD>& lme, 
			int Nk, int Nm,  
			Eigen::MatrixXd& Qt,
			GridBase *grid)
{
#ifdef USE_LAPACK
  const int size = Nm;
  int NN = Nk;
  double evals_tmp[NN];
  double evec_tmp[NN][NN];
  memset(evec_tmp[0],0,sizeof(double)*NN*NN);
  double DD[NN];
  double EE[NN];
  for (int i = 0; i< NN; i++) {
    for (int j = i - 1; j <= i + 1; j++) {
      if ( j < NN && j >= 0 ) {
	if (i==j) DD[i] = lmd[i];
	if (i==j) evals_tmp[i] = lmd[i];
	if (j==(i-1)) EE[j] = lme[j];
      }
    }
  }
  int evals_found;
  int lwork = ( (18*NN) > (1+4*NN+NN*NN)? (18*NN):(1+4*NN+NN*NN)) ;
  int liwork =  3+NN*10 ;
  int iwork[liwork];
  double work[lwork];
  int isuppz[2*NN];
  char jobz = 'V'; // calculate evals & evecs
  char range = 'I'; // calculate all evals
  //    char range = 'A'; // calculate all evals
  char uplo = 'U'; // refer to upper half of original matrix
  char compz = 'I'; // Compute eigenvectors of tridiagonal matrix
  int ifail[NN];
  int info;
  int total = grid->_Nprocessors;
  int node  = grid->_processor;
  int interval = (NN/total)+1;
  double vl = 0.0, vu = 0.0;
  int il = interval*node+1 , iu = interval*(node+1);
  if (iu > NN)  iu=NN;
  double tol = 0.0;
  if (1) {
    memset(evals_tmp,0,sizeof(double)*NN);
    if ( il <= NN){
      LAPACK_dstegr(&jobz, &range, &NN,
		    (double*)DD, (double*)EE,
		    &vl, &vu, &il, &iu, // these four are ignored if second parameteris 'A'
		    &tol, // tolerance
		    &evals_found, evals_tmp, (double*)evec_tmp, &NN,
		    isuppz,
		    work, &lwork, iwork, &liwork,
		    &info);
      for (int i = iu-1; i>= il-1; i--){
	evals_tmp[i] = evals_tmp[i - (il-1)];
	if (il>1) evals_tmp[i-(il-1)]=0.;
	for (int j = 0; j< NN; j++){
	  evec_tmp[i][j] = evec_tmp[i - (il-1)][j];
	  if (il>1) evec_tmp[i-(il-1)][j]=0.;
	}
      }
    }
    {
      grid->GlobalSumVector(evals_tmp,NN);
      grid->GlobalSumVector((double*)evec_tmp,NN*NN);
    }
  } 
  // Safer to sort instead of just reversing it, 
  // but the document of the routine says evals are sorted in increasing order. 
  // qr gives evals in decreasing order.
  for(int i=0;i<NN;i++){
    lmd [NN-1-i]=evals_tmp[i];
    for(int j=0;j<NN;j++){
      Qt((NN-1-i),j)=evec_tmp[i][j];
    }
  }
#else 
  assert(0);
#endif
}

  void diagonalize_QR(std::vector<RealD>& lmd, std::vector<RealD>& lme, 
		      int Nk, int Nm,   
		      Eigen::MatrixXd & Qt,
		      GridBase *grid)
  {
    int QRiter = 100*Nm;
    int kmin = 1;
    int kmax = Nk;

    // (this should be more sophisticated)
    for(int iter=0; iter<QRiter; ++iter){
      
      // determination of 2x2 leading submatrix
      RealD dsub = lmd[kmax-1]-lmd[kmax-2];
      RealD dd = sqrt(dsub*dsub + 4.0*lme[kmax-2]*lme[kmax-2]);
      RealD Dsh = 0.5*(lmd[kmax-2]+lmd[kmax-1] +dd*(dsub/fabs(dsub)));
      // (Dsh: shift)
	
      // transformation
      QR_decomp(lmd,lme,Nk,Nm,Qt,Dsh,kmin,kmax); // Nk, Nm
	
      // Convergence criterion (redef of kmin and kamx)
      for(int j=kmax-1; j>= kmin; --j){
	RealD dds = fabs(lmd[j-1])+fabs(lmd[j]);
	if(fabs(lme[j-1])+dds > dds){
	  kmax = j+1;
	  goto continued;
	}
      }
      QRiter = iter;
      return;

    continued:
      for(int j=0; j<kmax-1; ++j){
	RealD dds = fabs(lmd[j])+fabs(lmd[j+1]);
	if(fabs(lme[j])+dds > dds){
	  kmin = j+1;
	  break;
	}
      }
    }
    std::cout << GridLogError << "[QL method] Error - Too many iteration: "<<QRiter<<"\n";
    abort();
  }

 };
}

#endif
