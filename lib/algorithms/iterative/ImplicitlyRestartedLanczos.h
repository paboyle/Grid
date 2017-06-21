    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/ImplicitlyRestartedLanczos.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Chulwoo Jung
Author: Guido Cossu

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
#ifndef GRID_IRL_H
#define GRID_IRL_H

#include <string.h> //memset

namespace Grid {

  enum IRLdiagonalisation { 
    IRLdiagonaliseWithDSTEGR,
    IRLdiagonaliseWithQR,
    IRLdiagonaliseWithEigen
  };
  ////////////////////////////////////////////////////////////////////////////////
  // Helper class for sorting the evalues AND evectors by Field
  // Use pointer swizzle on vectors
  ////////////////////////////////////////////////////////////////////////////////
template<class Field>
class SortEigen {
 private:
  static bool less_lmd(RealD left,RealD right){
    return left > right;
  }  
  static bool less_pair(std::pair<RealD,Field const*>& left,
                        std::pair<RealD,Field const*>& right){
    return left.first > (right.first);
  }  
  
 public:
  void push(std::vector<RealD>& lmd,std::vector<Field>& evec,int N) {
    
    ////////////////////////////////////////////////////////////////////////
    // PAB: FIXME: VERY VERY VERY wasteful: takes a copy of the entire vector set.
    //    : The vector reorder should be done by pointer swizzle somehow
    ////////////////////////////////////////////////////////////////////////
    std::vector<Field> cpy(lmd.size(),evec[0]._grid);
    for(int i=0;i<lmd.size();i++) cpy[i] = evec[i];
    
    std::vector<std::pair<RealD, Field const*> > emod(lmd.size());    

    for(int i=0;i<lmd.size();++i)  emod[i] = std::pair<RealD,Field const*>(lmd[i],&cpy[i]);

    partial_sort(emod.begin(),emod.begin()+N,emod.end(),less_pair);

    typename std::vector<std::pair<RealD, Field const*> >::iterator it = emod.begin();
    for(int i=0;i<N;++i){
      lmd[i]=it->first;
      evec[i]=*(it->second);
      ++it;
    }
  }
  void push(std::vector<RealD>& lmd,int N) {
    std::partial_sort(lmd.begin(),lmd.begin()+N,lmd.end(),less_lmd);
  }
  bool saturated(RealD lmd, RealD thrs) {
    return fabs(lmd) > fabs(thrs);
  }
};

/////////////////////////////////////////////////////////////
// Implicitly restarted lanczos
/////////////////////////////////////////////////////////////
template<class Field> 
class ImplicitlyRestartedLanczos {
private:       
  int MaxIter;   // Max iterations
  int Nstop;     // Number of evecs checked for convergence
  int Nk;        // Number of converged sought
  int Nm;        // Nm -- total number of vectors
  RealD eresid;
  IRLdiagonalisation diagonalisation;
  ////////////////////////////////////
  // Embedded objects
  ////////////////////////////////////
           SortEigen<Field> _sort;
  LinearOperatorBase<Field> &_Linop;
    OperatorFunction<Field> &_poly;

  /////////////////////////
  // Constructor
  /////////////////////////
public:       
 ImplicitlyRestartedLanczos(LinearOperatorBase<Field> &Linop, // op
			    OperatorFunction<Field> & poly,   // polynomial
			    int _Nstop, // really sought vecs
			    int _Nk,    // sought vecs
			    int _Nm,    // total vecs
			    RealD _eresid, // resid in lmd deficit 
			    int _MaxIter,  // Max iterations
			    IRLdiagonalisation _diagonalisation= IRLdiagonaliseWithEigen ) :
    _Linop(Linop),    _poly(poly),
      Nstop(_Nstop), Nk(_Nk), Nm(_Nm),
      eresid(_eresid),  MaxIter(_MaxIter),
      diagonalisation(_diagonalisation)
      { };

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
  void calc(std::vector<RealD>& eval,  std::vector<Field>& evec, const Field& src, int& Nconv)
  {
    
    GridBase *grid = evec[0]._grid;
    assert(grid == src._grid);
    
    std::cout << GridLogMessage <<"**************************************************************************"<< std::endl;
    std::cout << GridLogMessage <<" ImplicitlyRestartedLanczos::calc() starting iteration 0 /  "<< MaxIter<< std::endl;
    std::cout << GridLogMessage <<"**************************************************************************"<< std::endl;
    std::cout << GridLogMessage <<" -- seek   Nk    = " << Nk    <<" vectors"<< std::endl;
    std::cout << GridLogMessage <<" -- accept Nstop = " << Nstop <<" vectors"<< std::endl;
    std::cout << GridLogMessage <<" -- total  Nm    = " << Nm    <<" vectors"<< std::endl;
    std::cout << GridLogMessage <<" -- size of eval = " << eval.size() << std::endl;
    std::cout << GridLogMessage <<" -- size of evec = " << evec.size() << std::endl;
    if ( diagonalisation == IRLdiagonaliseWithDSTEGR ) {
      std::cout << GridLogMessage << "Diagonalisation is DSTEGR "<<std::endl;
    } else if ( diagonalisation == IRLdiagonaliseWithQR ) { 
      std::cout << GridLogMessage << "Diagonalisation is QR "<<std::endl;
    }  else if ( diagonalisation == IRLdiagonaliseWithEigen ) { 
      std::cout << GridLogMessage << "Diagonalisation is Eigen "<<std::endl;
    }
    std::cout << GridLogMessage <<"**************************************************************************"<< std::endl;
    
    assert(Nm == evec.size() && Nm == eval.size());
	
    std::vector<RealD> lme(Nm);  
    std::vector<RealD> lme2(Nm);
    std::vector<RealD> eval2(Nm);
    Eigen::MatrixXd    Qt = Eigen::MatrixXd::Zero(Nm,Nm);
    std::vector<int>   Iconv(Nm);

    std::vector<Field>  B(Nm,grid); // waste of space replicating
    
    Field f(grid);
    Field v(grid);
    
    int k1 = 1;
    int k2 = Nk;
    
    Nconv = 0;
    
    RealD beta_k;
  
    // Set initial vector
    evec[0] = src;
    std::cout << GridLogMessage <<"norm2(src)= " << norm2(src)<<std::endl;
    
    normalise(evec[0]);
    std::cout << GridLogMessage <<"norm2(evec[0])= " << norm2(evec[0]) <<std::endl;
    
    // Initial Nk steps
    for(int k=0; k<Nk; ++k) step(eval,lme,evec,f,Nm,k);
    
    // Restarting loop begins
    int iter;
    for(iter = 0; iter<MaxIter; ++iter){
      
      std::cout<< GridLogMessage <<" **********************"<< std::endl;
      std::cout<< GridLogMessage <<" Restart iteration = "<< iter << std::endl;
      std::cout<< GridLogMessage <<" **********************"<< std::endl;
      
      for(int k=Nk; k<Nm; ++k) step(eval,lme,evec,f,Nm,k);
      
      f *= lme[Nm-1];
      
      // getting eigenvalues
      for(int k=0; k<Nm; ++k){
	eval2[k] = eval[k+k1-1];
	lme2[k] = lme[k+k1-1];
      }
      Qt = Eigen::MatrixXd::Identity(Nm,Nm);
      diagonalize(eval2,lme2,Nm,Nm,Qt,grid);

      // sorting
      _sort.push(eval2,Nm);
      
      // Implicitly shifted QR transformations
      Qt = Eigen::MatrixXd::Identity(Nm,Nm);
      for(int ip=k2; ip<Nm; ++ip){ 
	qr_decomp(eval,lme,Nm,Nm,Qt,eval2[ip],k1,Nm);
      }
    
      for(int i=0; i<(Nk+1); ++i) B[i] = 0.0;
	  
      for(int j=k1-1; j<k2+1; ++j){
	for(int k=0; k<Nm; ++k){
	  B[j].checkerboard = evec[k].checkerboard;
	  B[j] += Qt(j,k) * evec[k];
	}
      }
      for(int j=k1-1; j<k2+1; ++j) evec[j] = B[j];
      
      // Compressed vector f and beta(k2)
      f *= Qt(k2-1,Nm-1);
      f += lme[k2-1] * evec[k2];
      beta_k = norm2(f);
      beta_k = sqrt(beta_k);
      std::cout<< GridLogMessage<<" beta(k) = "<<beta_k<<std::endl;
      
      RealD betar = 1.0/beta_k;
      evec[k2] = betar * f;
      lme[k2-1] = beta_k;
      
      // Convergence test
      for(int k=0; k<Nm; ++k){    
	eval2[k] = eval[k];
	lme2[k] = lme[k];
      }
      Qt = Eigen::MatrixXd::Identity(Nm,Nm);
      diagonalize(eval2,lme2,Nk,Nm,Qt,grid);
      
      for(int k = 0; k<Nk; ++k) B[k]=0.0;
      
      for(int j = 0; j<Nk; ++j){
	for(int k = 0; k<Nk; ++k){
	  B[j].checkerboard = evec[k].checkerboard;
	  B[j] += Qt(j,k) * evec[k];
	}
      }

      Nconv = 0;
      for(int i=0; i<Nk; ++i){
	
	_Linop.HermOp(B[i],v);
	    
	RealD vnum = real(innerProduct(B[i],v)); // HermOp.
	RealD vden = norm2(B[i]);
	eval2[i] = vnum/vden;
	v -= eval2[i]*B[i];
	RealD vv = norm2(v);
	
	std::cout.precision(13);
	std::cout << GridLogMessage << "[" << std::setw(3)<< std::setiosflags(std::ios_base::right) <<i<<"] ";
	std::cout << "eval = "<<std::setw(25)<< std::setiosflags(std::ios_base::left)<< eval2[i];
	std::cout << " |H B[i] - eval[i]B[i]|^2 "<< std::setw(25)<< std::setiosflags(std::ios_base::right)<< vv<< std::endl;
	
	// change the criteria as evals are supposed to be sorted, all evals smaller(larger) than Nstop should have converged
	if((vv<eresid*eresid) && (i == Nconv) ){
	  Iconv[Nconv] = i;
	  ++Nconv;
	}
	
      }  // i-loop end
      
      std::cout<< GridLogMessage <<" #modes converged: "<<Nconv<<std::endl;

      if( Nconv>=Nstop ){
	goto converged;
      }
    } // end of iter loop
    
    std::cout << GridLogMessage <<"**************************************************************************"<< std::endl;
    std::cout<< GridLogError    <<" ImplicitlyRestartedLanczos::calc() NOT converged.";
    std::cout << GridLogMessage <<"**************************************************************************"<< std::endl;
    abort();
	
  converged:
    // Sorting
    eval.resize(Nconv);
    evec.resize(Nconv,grid);
    for(int i=0; i<Nconv; ++i){
      eval[i] = eval2[Iconv[i]];
      evec[i] = B[Iconv[i]];
    }
    _sort.push(eval,evec,Nconv);
    
    std::cout << GridLogMessage <<"**************************************************************************"<< std::endl;
    std::cout << GridLogMessage << "ImplicitlyRestartedLanczos CONVERGED ; Summary :\n";
    std::cout << GridLogMessage <<"**************************************************************************"<< std::endl;
    std::cout << GridLogMessage << " -- Iterations  = "<< iter   << "\n";
    std::cout << GridLogMessage << " -- beta(k)     = "<< beta_k << "\n";
    std::cout << GridLogMessage << " -- Nconv       = "<< Nconv  << "\n";
    std::cout << GridLogMessage <<"**************************************************************************"<< std::endl;
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
    
    _poly(_Linop,evec[k],w);      // 3. wk:=Avk−βkv_{k−1}
    
    if(k>0) w -= lme[k-1] * evec[k-1];
    
    ComplexD zalph = innerProduct(evec[k],w); // 4. αk:=(wk,vk)
    RealD     alph = real(zalph);
    
    w = w - alph * evec[k];// 5. wk:=wk−αkvk
    
    RealD beta = normalise(w); // 6. βk+1 := ∥wk∥2. If βk+1 = 0 then Stop
    // 7. vk+1 := wk/βk+1
    
    lmd[k] = alph;
    lme[k] = beta;
    
    if ( k > 0 ) orthogonalize(w,evec,k); // orthonormalise
    if ( k < Nm-1) evec[k+1] = w;
    
    if ( beta < tiny ) std::cout << GridLogMessage << " beta is tiny "<<beta<<std::endl;
  }
      
  ///////////////////////////////////////////////////////////////////
  // 
  //
  ///////////////////////////////////////////////////////////////////
  void qr_decomp(std::vector<RealD>& lmd,   // Nm 
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
    int Niter = 100*Nm;
    int kmin = 1;
    int kmax = Nk;

    // (this should be more sophisticated)
    for(int iter=0; iter<Niter; ++iter){
      
      // determination of 2x2 leading submatrix
      RealD dsub = lmd[kmax-1]-lmd[kmax-2];
      RealD dd = sqrt(dsub*dsub + 4.0*lme[kmax-2]*lme[kmax-2]);
      RealD Dsh = 0.5*(lmd[kmax-2]+lmd[kmax-1] +dd*(dsub/fabs(dsub)));
      // (Dsh: shift)
	
      // transformation
      qr_decomp(lmd,lme,Nk,Nm,Qt,Dsh,kmin,kmax); // Nk, Nm
	
      // Convergence criterion (redef of kmin and kamx)
      for(int j=kmax-1; j>= kmin; --j){
	RealD dds = fabs(lmd[j-1])+fabs(lmd[j]);
	if(fabs(lme[j-1])+dds > dds){
	  kmax = j+1;
	  goto continued;
	}
      }
      Niter = iter;
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
    std::cout << GridLogError << "[QL method] Error - Too many iteration: "<<Niter<<"\n";
    abort();
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


  static RealD normalise(Field& v) 
  {
    RealD nn = norm2(v);
    nn = sqrt(nn);
    v = v * (1.0/nn);
    return nn;
  }
  
  void orthogonalize(Field& w, std::vector<Field>& evec, int k)
  {
    typedef typename Field::scalar_type MyComplex;
    MyComplex ip;
    
    for(int j=0; j<k; ++j){
      ip = innerProduct(evec[j],w); // are the evecs normalised? ; this assumes so.
      w = w - ip * evec[j];
    }
    normalise(w);
  }

 };
}
#endif
