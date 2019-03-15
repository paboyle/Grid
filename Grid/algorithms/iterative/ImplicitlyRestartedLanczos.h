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

namespace Grid { 

  ////////////////////////////////////////////////////////
  // Move following 100 LOC to lattice/Lattice_basis.h
  ////////////////////////////////////////////////////////
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

    std::vector < vobj , commAllocator<vobj> > B(Nm); // Thread private
       
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

// Extract a single rotated vector
template<class Field>
void basisRotateJ(Field &result,std::vector<Field> &basis,Eigen::MatrixXd& Qt,int j, int k0,int k1,int Nm) 
{
  typedef typename Field::vector_object vobj;
  GridBase* grid = basis[0]._grid;

  result.checkerboard = basis[0].checkerboard;
  parallel_for(int ss=0;ss < grid->oSites();ss++){
    vobj B = zero;
    for(int k=k0; k<k1; ++k){
      B +=Qt(j,k) * basis[k]._odata[ss];
    }
    result._odata[ss] = B;
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

      //////////////////////////////////////
      // idx[i] is a table of desired sources giving a permutation.
      // Swap v[i] with v[idx[i]].
      // Find  j>i for which _vnew[j] = _vold[i],
      // track the move idx[j] => idx[i]
      // track the move idx[i] => i
      //////////////////////////////////////
      size_t j;
      for (j=i;j<idx.size();j++)
	if (idx[j]==i)
	  break;

      assert(idx[i] > i);     assert(j!=idx.size());      assert(idx[j]==i);

      std::swap(_v[i]._odata,_v[idx[i]]._odata); // should use vector move constructor, no data copy
      std::swap(sort_vals[i],sort_vals[idx[i]]);

      idx[j] = idx[i];
      idx[i] = i;
    }
  }
}

inline std::vector<int> basisSortGetIndex(std::vector<RealD>& sort_vals) 
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

/////////////////////////////////////////////////////////////
// Implicitly restarted lanczos
/////////////////////////////////////////////////////////////
template<class Field> class ImplicitlyRestartedLanczosTester 
{
 public:
  virtual int TestConvergence(int j,RealD resid,Field &evec, RealD &eval,RealD evalMaxApprox)=0;
  virtual int ReconstructEval(int j,RealD resid,Field &evec, RealD &eval,RealD evalMaxApprox)=0;
};

enum IRLdiagonalisation { 
  IRLdiagonaliseWithDSTEGR,
  IRLdiagonaliseWithQR,
  IRLdiagonaliseWithEigen
};

template<class Field> class ImplicitlyRestartedLanczosHermOpTester  : public ImplicitlyRestartedLanczosTester<Field>
{
 public:

  LinearFunction<Field>       &_HermOp;
  ImplicitlyRestartedLanczosHermOpTester(LinearFunction<Field> &HermOp) : _HermOp(HermOp)  {  };
  int ReconstructEval(int j,RealD resid,Field &B, RealD &eval,RealD evalMaxApprox)
  {
    return TestConvergence(j,resid,B,eval,evalMaxApprox);
  }
  int TestConvergence(int j,RealD eresid,Field &B, RealD &eval,RealD evalMaxApprox)
  {
    Field v(B);
    RealD eval_poly = eval;
    // Apply operator
    _HermOp(B,v);

    RealD vnum = real(innerProduct(B,v)); // HermOp.
    RealD vden = norm2(B);
    RealD vv0  = norm2(v);
    eval   = vnum/vden;
    v -= eval*B;

    RealD vv = norm2(v) / ::pow(evalMaxApprox,2.0);

    std::cout.precision(13);
    std::cout<<GridLogIRL  << "[" << std::setw(3)<<j<<"] "
	     <<"eval = "<<std::setw(25)<< eval << " (" << eval_poly << ")"
	     <<" |H B[i] - eval[i]B[i]|^2 / evalMaxApprox^2 " << std::setw(25) << vv
	     <<std::endl;

    int conv=0;
    if( (vv<eresid*eresid) ) conv = 1;

    return conv;
  }
};

template<class Field> 
class ImplicitlyRestartedLanczos {
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
  LinearFunction<Field>       &_PolyOp;
  LinearFunction<Field>       &_HermOp;
  ImplicitlyRestartedLanczosTester<Field> &_Tester;
  // Default tester provided (we need a ref to something in default case)
  ImplicitlyRestartedLanczosHermOpTester<Field> SimpleTester;
  /////////////////////////
  // Constructor
  /////////////////////////
  
public:       

  //////////////////////////////////////////////////////////////////
  // PAB:
  //////////////////////////////////////////////////////////////////
  // Too many options  & knobs. 
  // Eliminate:
  //   orth_period
  //   betastp
  //   MinRestart
  //
  // Do we really need orth_period
  // What is the theoretical basis & guarantees of betastp ?
  // Nstop=Nk viable?
  // MinRestart avoidable with new convergence test?
  // Could cut to PolyOp, HermOp, Tester, Nk, Nm, resid, maxiter (+diagonalisation)
  // HermOp could be eliminated if we dropped the Power method for max eval.
  // -- also: The eval, eval2, eval2_copy stuff is still unnecessarily unclear
  //////////////////////////////////////////////////////////////////
 ImplicitlyRestartedLanczos(LinearFunction<Field> & PolyOp,
			    LinearFunction<Field> & HermOp,
			    ImplicitlyRestartedLanczosTester<Field> & Tester,
			    int _Nstop, // sought vecs
			    int _Nk, // sought vecs
			    int _Nm, // spare vecs
			    RealD _eresid, // resid in lmdue deficit 
			    int _MaxIter, // Max iterations
			    RealD _betastp=0.0, // if beta(k) < betastp: converged
			    int _MinRestart=1, int _orth_period = 1,
			    IRLdiagonalisation _diagonalisation= IRLdiagonaliseWithEigen) :
    SimpleTester(HermOp), _PolyOp(PolyOp),      _HermOp(HermOp), _Tester(Tester),
    Nstop(_Nstop)  ,      Nk(_Nk),      Nm(_Nm),
    eresid(_eresid),      betastp(_betastp),
    MaxIter(_MaxIter)  ,      MinRestart(_MinRestart),
    orth_period(_orth_period), diagonalisation(_diagonalisation)  { };

    ImplicitlyRestartedLanczos(LinearFunction<Field> & PolyOp,
			       LinearFunction<Field> & HermOp,
			       int _Nstop, // sought vecs
			       int _Nk, // sought vecs
			       int _Nm, // spare vecs
			       RealD _eresid, // resid in lmdue deficit 
			       int _MaxIter, // Max iterations
			       RealD _betastp=0.0, // if beta(k) < betastp: converged
			       int _MinRestart=1, int _orth_period = 1,
			       IRLdiagonalisation _diagonalisation= IRLdiagonaliseWithEigen) :
    SimpleTester(HermOp),  _PolyOp(PolyOp),      _HermOp(HermOp), _Tester(SimpleTester),
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

  void orthogonalize(Field& w, std::vector<Field>& evec,int k)
  {
    OrthoTime-=usecond()/1e6;
    basisOrthogonalize(evec,w,k);
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
  void calc(std::vector<RealD>& eval, std::vector<Field>& evec,  const Field& src, int& Nconv, bool reverse=false)
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
	normalise(src_n);
	_HermOp(src_n,tmp);
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
      std::partial_sort(eval2.begin(),eval2.begin()+Nm,eval2.end(),std::greater<RealD>());
      std::cout<<GridLogIRL <<" evals sorted "<<std::endl;
      const int chunk=8;
      for(int io=0; io<k2;io+=chunk){
	std::cout<<GridLogIRL << "eval "<< std::setw(3) << io ;
	for(int ii=0;ii<chunk;ii++){
	  if ( (io+ii)<k2 )
	    std::cout<< " "<< std::setw(12)<< eval2[io+ii];
	}
	std::cout << std::endl;
      }

      //////////////////////////////////
      // Implicitly shifted QR transformations
      //////////////////////////////////
      Qt = Eigen::MatrixXd::Identity(Nm,Nm);
      for(int ip=k2; ip<Nm; ++ip){ 
	QR_decomp(eval,lme,Nm,Nm,Qt,eval2[ip],k1,Nm);
      }
      std::cout<<GridLogIRL <<"QR decomposed "<<std::endl;

      assert(k2<Nm);      assert(k2<Nm);      assert(k1>0);

      basisRotate(evec,Qt,k1-1,k2+1,0,Nm,Nm); /// big constraint on the basis
      std::cout<<GridLogIRL <<"basisRotated  by Qt"<<std::endl;
      
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
      }
      Qt = Eigen::MatrixXd::Identity(Nm,Nm);
      diagonalize(eval2,lme2,Nk,Nm,Qt,grid);
      std::cout<<GridLogIRL <<" Diagonalized "<<std::endl;
	  
      Nconv = 0;
      if (iter >= MinRestart) {

	std::cout << GridLogIRL << "Test convergence: rotate subset of vectors to test convergence " << std::endl;

	Field B(grid); B.checkerboard = evec[0].checkerboard;

	//  power of two search pattern;  not every evalue in eval2 is assessed.
	int allconv =1;
	for(int jj = 1; jj<=Nstop; jj*=2){
	  int j = Nstop-jj;
	  RealD e = eval2_copy[j]; // Discard the evalue
	  basisRotateJ(B,evec,Qt,j,0,Nk,Nm);	    
	  if( !_Tester.TestConvergence(j,eresid,B,e,evalMaxApprox) ) {
	    allconv=0;
	  }
	}
	// Do evec[0] for good measure
	{ 
	  int j=0;
	  RealD e = eval2_copy[0]; 
	  basisRotateJ(B,evec,Qt,j,0,Nk,Nm);	    
	  if( !_Tester.TestConvergence(j,eresid,B,e,evalMaxApprox) ) allconv=0;
	}
	if ( allconv ) Nconv = Nstop;

	// test if we converged, if so, terminate
	std::cout<<GridLogIRL<<" #modes converged: >= "<<Nconv<<"/"<<Nstop<<std::endl;
	//	if( Nconv>=Nstop || beta_k < betastp){
	if( Nconv>=Nstop){
	  goto converged;
	}
	  
      } else {
	std::cout << GridLogIRL << "iter < MinRestart: do not yet test for convergence\n";
      } // end of iter loop
    }

    std::cout<<GridLogError<<"\n NOT converged.\n";
    abort();
	
  converged:
    {
      Field B(grid); B.checkerboard = evec[0].checkerboard;
      basisRotate(evec,Qt,0,Nk,0,Nk,Nm);	    
      std::cout << GridLogIRL << " Rotated basis"<<std::endl;
      Nconv=0;
      //////////////////////////////////////////////////////////////////////
      // Full final convergence test; unconditionally applied
      //////////////////////////////////////////////////////////////////////
      for(int j = 0; j<=Nk; j++){
	B=evec[j];
	if( _Tester.ReconstructEval(j,eresid,B,eval2[j],evalMaxApprox) ) {
	  Nconv++;
	}
      }

      if ( Nconv < Nstop )
	std::cout << GridLogIRL << "Nconv ("<<Nconv<<") < Nstop ("<<Nstop<<")"<<std::endl;

      eval=eval2;
      
      //Keep only converged
      eval.resize(Nconv);// Nstop?
      evec.resize(Nconv,grid);// Nstop?
      basisSortInPlace(evec,eval,reverse);
      
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
	    std::vector<Field>& evec,
	    Field& w,int Nm,int k)
  {
    const RealD tiny = 1.0e-20;
    assert( k< Nm );

    GridStopWatch gsw_op,gsw_o;

    Field& evec_k = evec[k];

    _PolyOp(evec_k,w);    std::cout<<GridLogIRL << "PolyOp" <<std::endl;

    if(k>0) w -= lme[k-1] * evec[k-1];

    ComplexD zalph = innerProduct(evec_k,w); // 4. αk:=(wk,vk)
    RealD     alph = real(zalph);

    w = w - alph * evec_k;// 5. wk:=wk−αkvk

    RealD beta = normalise(w); // 6. βk+1 := ∥wk∥2. If βk+1 = 0 then Stop
    // 7. vk+1 := wk/βk+1

    lmd[k] = alph;
    lme[k] = beta;

    if (k>0 && k % orth_period == 0) {
      orthogonalize(w,evec,k); // orthonormalise
      std::cout<<GridLogIRL << "Orthogonalised " <<std::endl;
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
  // File could end here if settle on Eigen ??? !!!
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
