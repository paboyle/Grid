    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/ImplicitlyRestartedLanczos.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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

#include <algorithms/iterative/DenseMatrix.h>
#include <algorithms/iterative/EigenSort.h>

namespace Grid {

/////////////////////////////////////////////////////////////
// Implicitly restarted lanczos
/////////////////////////////////////////////////////////////


template<class Field> 
    class ImplicitlyRestartedLanczos {

    const RealD small = 1.0e-16;
public:       
    int lock;
    int get;
    int Niter;
    int converged;

    int Nk;      // Number of converged sought
    int Np;      // Np -- Number of spare vecs in kryloc space
    int Nm;      // Nm -- total number of vectors

    RealD eresid;

    SortEigen<Field> _sort;

    LinearOperatorBase<Field> &_Linop;

    OperatorFunction<Field>   &_poly;

    /////////////////////////
    // Constructor
    /////////////////////////
    void init(void){};
    void Abort(int ff, DenseVector<RealD> &evals,  DenseVector<DenseVector<RealD> > &evecs);

    ImplicitlyRestartedLanczos(LinearOperatorBase<Field> &Linop, // op
			       OperatorFunction<Field> & poly,   // polynmial
			       int _Nk, // sought vecs
			       int _Nm, // spare vecs
			       RealD _eresid, // resid in lmdue deficit 
			       int _Niter) : // Max iterations
      _Linop(Linop),
      _poly(poly),
      Nk(_Nk),
      Nm(_Nm),
      eresid(_eresid),
      Niter(_Niter)
    { 
      Np = Nm-Nk; assert(Np>0);
    };

    /////////////////////////
    // Sanity checked this routine (step) against Saad.
    /////////////////////////
    void RitzMatrix(DenseVector<Field>& evec,int k){

      if(1) return;

      GridBase *grid = evec[0]._grid;
      Field w(grid);
      std::cout << "RitzMatrix "<<std::endl;
      for(int i=0;i<k;i++){
	_poly(_Linop,evec[i],w);
	std::cout << "["<<i<<"] ";
	for(int j=0;j<k;j++){
	  ComplexD in = innerProduct(evec[j],w);
	  if ( fabs((double)i-j)>1 ) { 
	    if (abs(in) >1.0e-9 )  { 
	      std::cout<<"oops"<<std::endl;
	      abort();
	    } else 
	      std::cout << " 0 ";
	  } else { 
	    std::cout << " "<<in<<" ";
	  }
	}
	std::cout << std::endl;
      }
    }

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
    void step(DenseVector<RealD>& lmd,
	      DenseVector<RealD>& lme, 
	      DenseVector<Field>& evec,
	      Field& w,int Nm,int k)
    {
      assert( k< Nm );
      
      _poly(_Linop,evec[k],w);      // 3. wk:=Avk−βkv_{k−1}
      if(k>0){
	w -= lme[k-1] * evec[k-1];
      }    

      ComplexD zalph = innerProduct(evec[k],w); // 4. αk:=(wk,vk)
      RealD     alph = real(zalph);

      w = w - alph * evec[k];// 5. wk:=wk−αkvk

      RealD beta = normalise(w); // 6. βk+1 := ∥wk∥2. If βk+1 = 0 then Stop
                                 // 7. vk+1 := wk/βk+1

      const RealD tiny = 1.0e-20;
      if ( beta < tiny ) { 
	std::cout << " beta is tiny "<<beta<<std::endl;
      }
      lmd[k] = alph;
      lme[k]  = beta;

      if (k>0) { 
	orthogonalize(w,evec,k); // orthonormalise
      }
      
      if(k < Nm-1) evec[k+1] = w;
    }

    void qr_decomp(DenseVector<RealD>& lmd,
		   DenseVector<RealD>& lme,
		   int Nk,
		   int Nm,
		   DenseVector<RealD>& Qt,
		   RealD Dsh, 
		   int kmin,
		   int kmax)
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
	RealD Qtmp1 = Qt[i+Nm*k  ];
	RealD Qtmp2 = Qt[i+Nm*(k+1)];
	Qt[i+Nm*k    ] = c*Qtmp1 - s*Qtmp2;
	Qt[i+Nm*(k+1)] = s*Qtmp1 + c*Qtmp2; 
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
	  RealD Qtmp1 = Qt[i+Nm*k    ];
	  RealD Qtmp2 = Qt[i+Nm*(k+1)];
	  Qt[i+Nm*k    ] = c*Qtmp1 -s*Qtmp2;
	  Qt[i+Nm*(k+1)] = s*Qtmp1 +c*Qtmp2;
	}
      }
    }

    void diagonalize(DenseVector<RealD>& lmd,
		     DenseVector<RealD>& lme, 
		     int Nm2,
		     int Nm,
		     DenseVector<RealD>& Qt)
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
	qr_decomp(lmd,lme,Nk,Nm,Qt,Dsh,kmin,kmax);
	
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
      std::cout << "[QL method] Error - Too many iteration: "<<Niter<<"\n";
      abort();
    }

    static RealD normalise(Field& v) 
    {
      RealD nn = norm2(v);
      nn = sqrt(nn);
      v = v * (1.0/nn);
      return nn;
    }

    void orthogonalize(Field& w,
		       DenseVector<Field>& evec,
		       int k)
    {
      typedef typename Field::scalar_type MyComplex;
      MyComplex ip;

      if ( 0 ) {
	for(int j=0; j<k; ++j){
	  normalise(evec[j]);
	  for(int i=0;i<j;i++){
	    ip = innerProduct(evec[i],evec[j]); // are the evecs normalised? ; this assumes so.
	    evec[j] = evec[j] - ip *evec[i];
	  }
	}
      }

      for(int j=0; j<k; ++j){
	ip = innerProduct(evec[j],w); // are the evecs normalised? ; this assumes so.
	w = w - ip * evec[j];
      }
      normalise(w);
    }

    void setUnit_Qt(int Nm, DenseVector<RealD> &Qt) {
      for(int i=0; i<Qt.size(); ++i) Qt[i] = 0.0;
      for(int k=0; k<Nm; ++k) Qt[k + k*Nm] = 1.0;
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
    void calc(DenseVector<RealD>& eval,
	      DenseVector<Field>& evec,
	      const Field& src,
	      int& Nconv)
      {

	GridBase *grid = evec[0]._grid;

	std::cout << " -- Nk = " << Nk << " Np = "<< Np << std::endl;
	std::cout << " -- Nm = " << Nm << std::endl;
	std::cout << " -- size of eval   = " << eval.size() << std::endl;
	std::cout << " -- size of evec  = " << evec.size() << std::endl;
	
	assert(Nm == evec.size() && Nm == eval.size());
	
	DenseVector<RealD> lme(Nm);  
	DenseVector<RealD> lme2(Nm);
	DenseVector<RealD> eval2(Nm);
	DenseVector<RealD> Qt(Nm*Nm);
	DenseVector<int>   Iconv(Nm);

	DenseVector<Field>  B(Nm,grid); // waste of space replicating
	
	Field f(grid);
	Field v(grid);
  
	int k1 = 1;
	int k2 = Nk;

	Nconv = 0;

	RealD beta_k;
  
	// Set initial vector
	// (uniform vector) Why not src??
	//	evec[0] = 1.0;
	evec[0] = src;
	normalise(evec[0]);
	
	// Initial Nk steps
	for(int k=0; k<Nk; ++k) step(eval,lme,evec,f,Nm,k);
	RitzMatrix(evec,Nk);

	// Restarting loop begins
	for(int iter = 0; iter<Niter; ++iter){

	  std::cout<<"\n Restart iteration = "<< iter << std::endl;

	  // 
	  // Rudy does a sort first which looks very different. Getting fed up with sorting out the algo defs.
	  // We loop over 
	  //
	  for(int k=Nk; k<Nm; ++k) step(eval,lme,evec,f,Nm,k);
	  f *= lme[Nm-1];

	  RitzMatrix(evec,k2);
	  
	  // getting eigenvalues
	  for(int k=0; k<Nm; ++k){
	    eval2[k] = eval[k+k1-1];
	    lme2[k] = lme[k+k1-1];
	  }
	  setUnit_Qt(Nm,Qt);
	  diagonalize(eval2,lme2,Nm,Nm,Qt);

	  // sorting
	  _sort.push(eval2,Nm);
	  
	  // Implicitly shifted QR transformations
	  setUnit_Qt(Nm,Qt);
	  for(int ip=k2; ip<Nm; ++ip) 
	    qr_decomp(eval,lme,Nm,Nm,Qt,eval2[ip],k1,Nm);
    
	  for(int i=0; i<(Nk+1); ++i) B[i] = 0.0;
	  
	  for(int j=k1-1; j<k2+1; ++j){
	    for(int k=0; k<Nm; ++k){
	      B[j] += Qt[k+Nm*j] * evec[k];
	    }
	  }
	  for(int j=k1-1; j<k2+1; ++j) evec[j] = B[j];

	  // Compressed vector f and beta(k2)
	  f *= Qt[Nm-1+Nm*(k2-1)];
	  f += lme[k2-1] * evec[k2];
	  beta_k = norm2(f);
	  beta_k = sqrt(beta_k);
	  std::cout<<" beta(k) = "<<beta_k<<std::endl;

	  RealD betar = 1.0/beta_k;
	  evec[k2] = betar * f;
	  lme[k2-1] = beta_k;

	  // Convergence test
	  for(int k=0; k<Nm; ++k){    
	    eval2[k] = eval[k];
	    lme2[k] = lme[k];
	  }
	  setUnit_Qt(Nm,Qt);
	  diagonalize(eval2,lme2,Nk,Nm,Qt);
	  
	  for(int k = 0; k<Nk; ++k) B[k]=0.0;
	  
	  for(int j = 0; j<Nk; ++j){
	    for(int k = 0; k<Nk; ++k){
	      B[j] += Qt[k+j*Nm] * evec[k];
	    }
	  }

	  Nconv = 0;
	  //	  std::cout << std::setiosflags(std::ios_base::scientific);
	  for(int i=0; i<Nk; ++i){

	    _poly(_Linop,B[i],v);
	    
	    RealD vnum = real(innerProduct(B[i],v)); // HermOp.
	    RealD vden = norm2(B[i]);
	    eval2[i] = vnum/vden;
	    v -= eval2[i]*B[i];
	    RealD vv = norm2(v);
	    
	    std::cout << "[" << std::setw(3)<< std::setiosflags(std::ios_base::right) <<i<<"] ";
	    std::cout << "eval = "<<std::setw(25)<< std::setiosflags(std::ios_base::left)<< eval2[i];
	    std::cout <<" |H B[i] - eval[i]B[i]|^2 "<< std::setw(25)<< std::setiosflags(std::ios_base::right)<< vv<< std::endl;
	    
	    if(vv<eresid*eresid){
	      Iconv[Nconv] = i;
	      ++Nconv;
	    }

	  }  // i-loop end
	  //	  std::cout << std::resetiosflags(std::ios_base::scientific);


	  std::cout<<" #modes converged: "<<Nconv<<std::endl;

	  if( Nconv>=Nk ){
	    goto converged;
	  }
	} // end of iter loop
	
	std::cout<<"\n NOT converged.\n";
	abort();
	
      converged:
	// Sorting
	
	eval.clear();
	evec.clear();
	for(int i=0; i<Nconv; ++i){
	  eval.push_back(eval2[Iconv[i]]);
	  evec.push_back(B[Iconv[i]]);
	}
	_sort.push(eval,evec,Nconv);
	
	std::cout << "\n Converged\n Summary :\n";
	std::cout << " -- Iterations  = "<< Nconv  << "\n";
	std::cout << " -- beta(k)     = "<< beta_k << "\n";
	std::cout << " -- Nconv       = "<< Nconv  << "\n";
      }

    /////////////////////////////////////////////////
    // Adapted from Rudy's lanczos factor routine
    /////////////////////////////////////////////////
    int Lanczos_Factor(int start, int end,  int cont,
		       DenseVector<Field> & bq, 
		       Field &bf,
		       DenseMatrix<RealD> &H){
      
      GridBase *grid = bq[0]._grid;

      RealD beta;  
      RealD sqbt;  
      RealD alpha;

      for(int i=start;i<Nm;i++){
	for(int j=start;j<Nm;j++){
	  H[i][j]=0.0;
	}
      }

      std::cout<<"Lanczos_Factor start/end " <<start <<"/"<<end<<std::endl;

      // Starting from scratch, bq[0] contains a random vector and |bq[0]| = 1
      int first;
      if(start == 0){

	std::cout << "start == 0\n"; //TESTING

	_poly(_Linop,bq[0],bf);

	alpha = real(innerProduct(bq[0],bf));//alpha =  bq[0]^dag A bq[0]

	std::cout << "alpha = " << alpha << std::endl;
	
	bf = bf - alpha * bq[0];  //bf =  A bq[0] - alpha bq[0]

	H[0][0]=alpha;

	std::cout << "Set H(0,0) to " << H[0][0] << std::endl;

	first = 1;

      } else {

	first = start;

      }

      // I think start==0 and cont==zero are the same. Test this
      // If so I can drop "cont" parameter?
      if( cont ) assert(start!=0);

      if( start==0 ) assert(cont!=0);

      if( cont){

	beta = 0;sqbt = 0;

	std::cout << "cont is true so setting beta to zero\n";

      }	else {

	beta = norm2(bf);
	sqbt = sqrt(beta);

	std::cout << "beta = " << beta << std::endl;
      }

      for(int j=first;j<end;j++){

	std::cout << "Factor j " << j <<std::endl;

	if(cont){ // switches to factoring; understand start!=0 and initial bf value is right.
	  bq[j] = bf; cont = false;
	}else{
	  bq[j] = (1.0/sqbt)*bf ;

	  H[j][j-1]=H[j-1][j] = sqbt;
	}

	_poly(_Linop,bq[j],bf);

	bf = bf - (1.0/sqbt)*bq[j-1]; 	       //bf = A bq[j] - beta bq[j-1] // PAB this comment was incorrect in beta term??

	alpha = real(innerProduct(bq[j],bf));  //alpha = bq[j]^dag A bq[j]

	bf = bf - alpha*bq[j];                 //bf = A bq[j] - beta bq[j-1] - alpha bq[j]
	RealD fnorm = norm2(bf);

	RealD bck = sqrt( real( conjugate(alpha)*alpha ) + beta );

	beta = fnorm;
	sqbt = sqrt(beta);
	std::cout << "alpha = " << alpha << " fnorm = " << fnorm << '\n';

	///Iterative refinement of orthogonality V = [ bq[0]  bq[1]  ...  bq[M] ]
	int re = 0;
	// FIXME undefined params; how set in Rudy's code
	int ref =0;
	Real rho = 1.0e-8;

	while( re == ref || (sqbt < rho * bck && re < 5) ){

	  Field tmp2(grid);
	  Field tmp1(grid);

	  //bex = V^dag bf
	  DenseVector<ComplexD> bex(j+1);
	  for(int k=0;k<j+1;k++){
	    bex[k] = innerProduct(bq[k],bf);
	  }
	  
	  zero_fermion(tmp2);
	  //tmp2 = V s
	  for(int l=0;l<j+1;l++){
	    RealD nrm = norm2(bq[l]);
	    axpy(tmp1,0.0,bq[l],bq[l]); scale(tmp1,bex[l]); 	//tmp1 = V[j] bex[j]
	    axpy(tmp2,1.0,tmp2,tmp1);					//tmp2 += V[j] bex[j]
	  }

	  //bf = bf - V V^dag bf.   Subtracting off any component in span { V[j] } 
	  RealD btc = axpy_norm(bf,-1.0,tmp2,bf);
	  alpha = alpha + real(bex[j]);	      sqbt = sqrt(real(btc));	      
	  // FIXME is alpha real in RUDY's code?
	  RealD nmbex = 0;for(int k=0;k<j+1;k++){nmbex = nmbex + real( conjugate(bex[k])*bex[k]  );}
	  bck = sqrt( nmbex );
	  re++;
	}
	std::cout << "Iteratively refined orthogonality, changes alpha\n";
	if(re > 1) std::cout << "orthagonality refined " << re << " times" <<std::endl;
	H[j][j]=alpha;
      }

      return end;
    }

    void EigenSort(DenseVector<double> evals,
		   DenseVector<Field>  evecs){
      int N= evals.size();
      _sort.push(evals,evecs, evals.size(),N);
    }

    void ImplicitRestart(int TM, DenseVector<RealD> &evals,  DenseVector<DenseVector<RealD> > &evecs, DenseVector<Field> &bq, Field &bf, int cont)
    {
      std::cout << "ImplicitRestart begin. Eigensort starting\n";

      DenseMatrix<RealD> H; Resize(H,Nm,Nm);

      EigenSort(evals, evecs);

      ///Assign shifts
      int K=Nk;
      int M=Nm;
      int P=Np;
      int converged=0;
      if(K - converged < 4) P = (M - K-1); //one
      //      DenseVector<RealD> shifts(P + shift_extra.size());
      DenseVector<RealD> shifts(P);
      for(int k = 0; k < P; ++k)
	shifts[k] = evals[k]; 

      /// Shift to form a new H and q
      DenseMatrix<RealD> Q; Resize(Q,TM,TM);
      Unity(Q);
      Shift(Q, shifts); // H is implicitly passed in in Rudy's Shift routine

      int ff = K;

      /// Shifted H defines a new K step Arnoldi factorization
      RealD  beta = H[ff][ff-1]; 
      RealD  sig  = Q[TM - 1][ff - 1];
      std::cout << "beta = " << beta << " sig = " << real(sig) <<std::endl;

      std::cout << "TM = " << TM << " ";
      std::cout << norm2(bq[0]) << " -- before" <<std::endl;

      /// q -> q Q
      times_real(bq, Q, TM);

      std::cout << norm2(bq[0]) << " -- after " << ff <<std::endl;
      bf =  beta* bq[ff] + sig* bf;

      /// Do the rest of the factorization
      ff = Lanczos_Factor(ff, M,cont,bq,bf,H);
      
      if(ff < M)
	Abort(ff, evals, evecs);
    }

///Run the Eigensolver
    void Run(int cont, DenseVector<Field> &bq, Field &bf, DenseVector<DenseVector<RealD> > & evecs,DenseVector<RealD> &evals)
    {
      init();

      int M=Nm;

      DenseMatrix<RealD> H; Resize(H,Nm,Nm);
      Resize(evals,Nm);
      Resize(evecs,Nm);

      int ff = Lanczos_Factor(0, M, cont, bq,bf,H); // 0--M to begin with

      if(ff < M) {
	std::cout << "Krylov: aborting ff "<<ff <<" "<<M<<std::endl;
	abort(); // Why would this happen?
      }

      int itcount = 0;
      bool stop = false;

      for(int it = 0; it < Niter && (converged < Nk); ++it) {

	std::cout << "Krylov: Iteration --> " << it << std::endl;
	int lock_num = lock ? converged : 0;
	DenseVector<RealD> tevals(M - lock_num );
	DenseMatrix<RealD> tevecs; Resize(tevecs,M - lock_num,M - lock_num);
	  
	//check residual of polynominal 
	TestConv(H,M, tevals, tevecs);

	if(converged >= Nk)
	    break;

	ImplicitRestart(ff, tevals,tevecs,H);
      }
      Wilkinson<RealD>(H, evals, evecs, small); 
      //      Check();

      std::cout << "Done  "<<std::endl;

    }

   ///H - shift I = QR; H = Q* H Q
    void Shift(DenseMatrix<RealD> & H,DenseMatrix<RealD> &Q, DenseVector<RealD> shifts) {
      
      int P; Size(shifts,P);
      int M; SizeSquare(Q,M);

      Unity(Q);

      int lock_num = lock ? converged : 0;

      RealD t_Househoulder_vector(0.0);
      RealD t_Househoulder_mult(0.0);

      for(int i=0;i<P;i++){

	RealD x, y, z;
	DenseVector<RealD> ck(3), v(3);
	  
	x = H[lock_num+0][lock_num+0]-shifts[i];
	y = H[lock_num+1][lock_num+0];
	ck[0] = x; ck[1] = y; ck[2] = 0; 

	normalise(ck);	///Normalization cancels in PHP anyway
	RealD beta;

	Householder_vector<RealD>(ck, 0, 2, v, beta);
	Householder_mult<RealD>(H,v,beta,0,lock_num+0,lock_num+2,0);
	Householder_mult<RealD>(H,v,beta,0,lock_num+0,lock_num+2,1);
	///Accumulate eigenvector
	Householder_mult<RealD>(Q,v,beta,0,lock_num+0,lock_num+2,1);
	  
	int sw = 0;
	for(int k=lock_num+0;k<M-2;k++){

	  x = H[k+1][k]; 
	  y = H[k+2][k]; 
	  z = (RealD)0.0;
	  if(k+3 <= M-1){
	    z = H[k+3][k];
	  }else{
	    sw = 1; v[2] = 0.0;
	  }

	  ck[0] = x; ck[1] = y; ck[2] = z;

	  normalise(ck);

	  Householder_vector<RealD>(ck, 0, 2-sw, v, beta);
	  Householder_mult<RealD>(H,v, beta,0,k+1,k+3-sw,0);
	  Householder_mult<RealD>(H,v, beta,0,k+1,k+3-sw,1);
	  ///Accumulate eigenvector
	  Householder_mult<RealD>(Q,v, beta,0,k+1,k+3-sw,1);
	}
      }
    }

    void TestConv(DenseMatrix<RealD> & H,int SS, 
		  DenseVector<Field> &bq, Field &bf,
		  DenseVector<RealD> &tevals, DenseVector<DenseVector<RealD> > &tevecs, 
		  int lock, int converged)
    {
      std::cout << "Converged " << converged << " so far." << std::endl;
      int lock_num = lock ? converged : 0;
      int M = Nm;

      ///Active Factorization
      DenseMatrix<RealD> AH; Resize(AH,SS - lock_num,SS - lock_num );

      AH = GetSubMtx(H,lock_num, SS, lock_num, SS);

      int NN=tevals.size();
      int AHsize=SS-lock_num;

      RealD small=1.0e-16;
      Wilkinson<RealD>(AH, tevals, tevecs, small);

      EigenSort(tevals, tevecs);

      RealD resid_nrm=  norm2(bf);

      if(!lock) converged = 0;
#if 0
      for(int i = SS - lock_num - 1; i >= SS - Nk && i >= 0; --i){

	RealD diff = 0;
	diff = abs( tevecs[i][Nm - 1 - lock_num] ) * resid_nrm;

	std::cout << "residual estimate " << SS-1-i << " " << diff << " of (" << tevals[i] << ")" << std::endl;

	if(diff < converged) {

	  if(lock) {
	    
	    DenseMatrix<RealD> Q; Resize(Q,M,M);
	    bool herm = true; 

	    Lock(H, Q, tevals[i], converged, small, SS, herm);

	    times_real(bq, Q, bq.size());
	    bf = Q[M - 1][M - 1]* bf;
	    lock_num++;
	  }
	  converged++;
	  std::cout << " converged on eval " << converged << " of " << Nk << std::endl;
	} else {
	  break;
	}
      }
#endif
      std::cout << "Got " << converged << " so far " <<std::endl;	
    }

    ///Check
    void Check(DenseVector<RealD> &evals,
	       DenseVector<DenseVector<RealD> > &evecs) {

      DenseVector<RealD> goodval(this->get);

      EigenSort(evals,evecs);

      int NM = Nm;

      DenseVector< DenseVector<RealD> > V; Size(V,NM);
      DenseVector<RealD> QZ(NM*NM);

      for(int i = 0; i < NM; i++){
	for(int j = 0; j < NM; j++){
	  // evecs[i][j];
	}
      }
    }


/**
   There is some matrix Q such that for any vector y
   Q.e_1 = y and Q is unitary.
**/
  template<class T>
  static T orthQ(DenseMatrix<T> &Q, DenseVector<T> y){
    int N = y.size();	//Matrix Size
    Fill(Q,0.0);
    T tau;
    for(int i=0;i<N;i++){
      Q[i][0]=y[i];
    }
    T sig = conj(y[0])*y[0];
    T tau0 = abs(sqrt(sig));
    
    for(int j=1;j<N;j++){
      sig += conj(y[j])*y[j]; 
      tau = abs(sqrt(sig) ); 	

      if(abs(tau0) > 0.0){
	
	T gam = conj( (y[j]/tau)/tau0 );
	for(int k=0;k<=j-1;k++){  
	  Q[k][j]=-gam*y[k];
	}
	Q[j][j]=tau0/tau;
      } else {
	Q[j-1][j]=1.0;
      }
      tau0 = tau;
    }
    return tau;
  }

/**
	There is some matrix Q such that for any vector y
	Q.e_k = y and Q is unitary.
**/
  template< class T>
  static T orthU(DenseMatrix<T> &Q, DenseVector<T> y){
    T tau = orthQ(Q,y);
    SL(Q);
    return tau;
  }


/**
	Wind up with a matrix with the first con rows untouched

say con = 2
	Q is such that Qdag H Q has {x, x, val, 0, 0, 0, 0, ...} as 1st colum
	and the matrix is upper hessenberg
	and with f and Q appropriately modidied with Q is the arnoldi factorization

**/

template<class T>
static void Lock(DenseMatrix<T> &H, 	///Hess mtx	
		 DenseMatrix<T> &Q, 	///Lock Transform
		 T val, 		///value to be locked
		 int con, 	///number already locked
		 RealD small,
		 int dfg,
		 bool herm)
{	


  //ForceTridiagonal(H);

  int M = H.dim;
  DenseVector<T> vec; Resize(vec,M-con);

  DenseMatrix<T> AH; Resize(AH,M-con,M-con);
  AH = GetSubMtx(H,con, M, con, M);

  DenseMatrix<T> QQ; Resize(QQ,M-con,M-con);

  Unity(Q);   Unity(QQ);
  
  DenseVector<T> evals; Resize(evals,M-con);
  DenseMatrix<T> evecs; Resize(evecs,M-con,M-con);

  Wilkinson<T>(AH, evals, evecs, small);

  int k=0;
  RealD cold = abs( val - evals[k]); 
  for(int i=1;i<M-con;i++){
    RealD cnew = abs( val - evals[i]);
    if( cnew < cold ){k = i; cold = cnew;}
  }
  vec = evecs[k];

  ComplexD tau;
  orthQ(QQ,vec);
  //orthQM(QQ,AH,vec);

  AH = Hermitian(QQ)*AH;
  AH = AH*QQ;
	

  for(int i=con;i<M;i++){
    for(int j=con;j<M;j++){
      Q[i][j]=QQ[i-con][j-con];
      H[i][j]=AH[i-con][j-con];
    }
  }

  for(int j = M-1; j>con+2; j--){

    DenseMatrix<T> U; Resize(U,j-1-con,j-1-con);
    DenseVector<T> z; Resize(z,j-1-con); 
    T nm = norm(z); 
    for(int k = con+0;k<j-1;k++){
      z[k-con] = conj( H(j,k+1) );
    }
    normalise(z);

    RealD tmp = 0;
    for(int i=0;i<z.size()-1;i++){tmp = tmp + abs(z[i]);}

    if(tmp < small/( (RealD)z.size()-1.0) ){ continue;}	

    tau = orthU(U,z);

    DenseMatrix<T> Hb; Resize(Hb,j-1-con,M);	
	
    for(int a = 0;a<M;a++){
      for(int b = 0;b<j-1-con;b++){
	T sum = 0;
	for(int c = 0;c<j-1-con;c++){
	  sum += H[a][con+1+c]*U[c][b];
	}//sum += H(a,con+1+c)*U(c,b);}
	Hb[b][a] = sum;
      }
    }
	
    for(int k=con+1;k<j;k++){
      for(int l=0;l<M;l++){
	H[l][k] = Hb[k-1-con][l];
      }
    }//H(Hb[k-1-con][l] , l,k);}}

    DenseMatrix<T> Qb; Resize(Qb,M,M);	
	
    for(int a = 0;a<M;a++){
      for(int b = 0;b<j-1-con;b++){
	T sum = 0;
	for(int c = 0;c<j-1-con;c++){
	  sum += Q[a][con+1+c]*U[c][b];
	}//sum += Q(a,con+1+c)*U(c,b);}
	Qb[b][a] = sum;
      }
    }
	
    for(int k=con+1;k<j;k++){
      for(int l=0;l<M;l++){
	Q[l][k] = Qb[k-1-con][l];
      }
    }//Q(Qb[k-1-con][l] , l,k);}}

    DenseMatrix<T> Hc; Resize(Hc,M,M);	
	
    for(int a = 0;a<j-1-con;a++){
      for(int b = 0;b<M;b++){
	T sum = 0;
	for(int c = 0;c<j-1-con;c++){
	  sum += conj( U[c][a] )*H[con+1+c][b];
	}//sum += conj( U(c,a) )*H(con+1+c,b);}
	Hc[b][a] = sum;
      }
    }

    for(int k=0;k<M;k++){
      for(int l=con+1;l<j;l++){
	H[l][k] = Hc[k][l-1-con];
      }
    }//H(Hc[k][l-1-con] , l,k);}}

  }
}


 };

}
#endif

