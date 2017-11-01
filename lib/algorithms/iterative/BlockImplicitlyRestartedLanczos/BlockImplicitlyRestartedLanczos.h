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

#include <zlib.h>
#include <sys/stat.h>

#include <Grid/algorithms/iterative/BlockImplicitlyRestartedLanczos/BlockedGrid.h>
#include <Grid/algorithms/iterative/BlockImplicitlyRestartedLanczos/FieldBasisVector.h>
#include <Grid/algorithms/iterative/BlockImplicitlyRestartedLanczos/BlockProjector.h>

namespace Grid { 

/////////////////////////////////////////////////////////////
// Implicitly restarted lanczos
/////////////////////////////////////////////////////////////

 template<class Field> 
 class BlockImplicitlyRestartedLanczos {

    const RealD small = 1.0e-16;
public:       
    int lock;
    int get;
    int Niter;
    int converged;

    int Nminres; // Minimum number of restarts; only check for convergence after
    int Nstop;   // Number of evecs checked for convergence
    int Nk;      // Number of converged sought
    int Np;      // Np -- Number of spare vecs in kryloc space
    int Nm;      // Nm -- total number of vectors

    int orth_period;

    RealD OrthoTime;

    RealD eresid, betastp;
    SortEigen<Field> _sort;
    LinearFunction<Field> &_HermOp;
    LinearFunction<Field> &_HermOpTest;
    /////////////////////////
    // Constructor
    /////////////////////////

    BlockImplicitlyRestartedLanczos(
			       LinearFunction<Field> & HermOp,
			       LinearFunction<Field> & HermOpTest,
			       int _Nstop, // sought vecs
			       int _Nk, // sought vecs
			       int _Nm, // spare vecs
			       RealD _eresid, // resid in lmdue deficit 
			       RealD _betastp, // if beta(k) < betastp: converged
			       int _Niter, // Max iterations
			       int _Nminres, int _orth_period = 1) :
      _HermOp(HermOp),
      _HermOpTest(HermOpTest),
      Nstop(_Nstop),
      Nk(_Nk),
      Nm(_Nm),
      eresid(_eresid),
      betastp(_betastp),
      Niter(_Niter),
	Nminres(_Nminres),
	orth_period(_orth_period)
    { 
      Np = Nm-Nk; assert(Np>0);
    };

    BlockImplicitlyRestartedLanczos(
			       LinearFunction<Field> & HermOp,
			       LinearFunction<Field> & HermOpTest,
			       int _Nk, // sought vecs
			       int _Nm, // spare vecs
			       RealD _eresid, // resid in lmdue deficit 
			       RealD _betastp, // if beta(k) < betastp: converged
			       int _Niter, // Max iterations
			       int _Nminres,
			       int _orth_period = 1) : 
      _HermOp(HermOp),
      _HermOpTest(HermOpTest),
      Nstop(_Nk),
      Nk(_Nk),
      Nm(_Nm),
      eresid(_eresid),
      betastp(_betastp),
      Niter(_Niter),
	Nminres(_Nminres),
	orth_period(_orth_period)
    { 
      Np = Nm-Nk; assert(Np>0);
    };


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
      assert( k< Nm );

      GridStopWatch gsw_op,gsw_o;

      Field& evec_k = evec[k];

      gsw_op.Start();
      _HermOp(evec_k,w);
      gsw_op.Stop();

      if(k>0){
	w -= lme[k-1] * evec[k-1];
      }    

      ComplexD zalph = innerProduct(evec_k,w); // 4. αk:=(wk,vk)
      RealD     alph = real(zalph);

      w = w - alph * evec_k;// 5. wk:=wk−αkvk

      RealD beta = normalise(w); // 6. βk+1 := ∥wk∥2. If βk+1 = 0 then Stop
                                 // 7. vk+1 := wk/βk+1

      std::cout<<GridLogMessage << "alpha[" << k << "] = " << zalph << " beta[" << k << "] = "<<beta<<std::endl;
      const RealD tiny = 1.0e-20;
      if ( beta < tiny ) { 
	std::cout<<GridLogMessage << " beta is tiny "<<beta<<std::endl;
     }
      lmd[k] = alph;
      lme[k]  = beta;

      gsw_o.Start();
      if (k>0 && k % orth_period == 0) { 
	orthogonalize(w,evec,k); // orthonormalise
      }
      gsw_o.Stop();

      if(k < Nm-1) { 
	evec[k+1] = w;
      }

      std::cout << GridLogMessage << "Timing: operator=" << gsw_op.Elapsed() <<
	" orth=" << gsw_o.Elapsed() << std::endl;

    }

    void qr_decomp(std::vector<RealD>& lmd,
		   std::vector<RealD>& lme,
		   int Nk,
		   int Nm,
		   std::vector<RealD>& Qt,
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

#ifdef USE_LAPACK_IRL
#define LAPACK_INT int
//long long
    void diagonalize_lapack(std::vector<RealD>& lmd,
			    std::vector<RealD>& lme, 
			    int N1,
			    int N2,
			    std::vector<RealD>& Qt,
			    GridBase *grid){

      std::cout << GridLogMessage << "diagonalize_lapack start\n";
      GridStopWatch gsw;

      const int size = Nm;
      //  tevals.resize(size);
      //  tevecs.resize(size);
      LAPACK_INT NN = N1;
      std::vector<double> evals_tmp(NN);
      std::vector<double> evec_tmp(NN*NN);
      memset(&evec_tmp[0],0,sizeof(double)*NN*NN);
      //  double AA[NN][NN];
      std::vector<double> DD(NN);
      std::vector<double> EE(NN);
      for (int i = 0; i< NN; i++)
	for (int j = i - 1; j <= i + 1; j++)
	  if ( j < NN && j >= 0 ) {
	    if (i==j) DD[i] = lmd[i];
	    if (i==j) evals_tmp[i] = lmd[i];
	    if (j==(i-1)) EE[j] = lme[j];
	  }
      LAPACK_INT evals_found;
      LAPACK_INT lwork = ( (18*NN) > (1+4*NN+NN*NN)? (18*NN):(1+4*NN+NN*NN)) ;
      LAPACK_INT liwork =  3+NN*10 ;
      std::vector<LAPACK_INT> iwork(liwork);
      std::vector<double> work(lwork);
      std::vector<LAPACK_INT> isuppz(2*NN);
      char jobz = 'V'; // calculate evals & evecs
      char range = 'I'; // calculate all evals
      //    char range = 'A'; // calculate all evals
      char uplo = 'U'; // refer to upper half of original matrix
      char compz = 'I'; // Compute eigenvectors of tridiagonal matrix
      std::vector<int> ifail(NN);
      LAPACK_INT info;
      //  int total = QMP_get_number_of_nodes();
      //  int node = QMP_get_node_number();
      //  GridBase *grid = evec[0]._grid;
      int total = grid->_Nprocessors;
      int node = grid->_processor;
      int interval = (NN/total)+1;
      double vl = 0.0, vu = 0.0;
      LAPACK_INT il = interval*node+1 , iu = interval*(node+1);
      if (iu > NN)  iu=NN;
      double tol = 0.0;
      if (1) {
	memset(&evals_tmp[0],0,sizeof(double)*NN);
	if ( il <= NN){
	  std::cout << GridLogMessage << "dstegr started" << std::endl; 
	  gsw.Start();
	  dstegr(&jobz, &range, &NN,
		 (double*)&DD[0], (double*)&EE[0],
		 &vl, &vu, &il, &iu, // these four are ignored if second parameteris 'A'
		 &tol, // tolerance
		 &evals_found, &evals_tmp[0], (double*)&evec_tmp[0], &NN,
		 &isuppz[0],
		 &work[0], &lwork, &iwork[0], &liwork,
		 &info);
	  gsw.Stop();
	  std::cout << GridLogMessage << "dstegr completed in " << gsw.Elapsed() << std::endl;
	  for (int i = iu-1; i>= il-1; i--){
	    evals_tmp[i] = evals_tmp[i - (il-1)];
	    if (il>1) evals_tmp[i-(il-1)]=0.;
	    for (int j = 0; j< NN; j++){
	      evec_tmp[i*NN + j] = evec_tmp[(i - (il-1)) * NN + j];
	      if (il>1) evec_tmp[(i-(il-1)) * NN + j]=0.;
	    }
	  }
	}
	{
	  //        QMP_sum_double_array(evals_tmp,NN);
	  //        QMP_sum_double_array((double *)evec_tmp,NN*NN);
	  grid->GlobalSumVector(&evals_tmp[0],NN);
	  grid->GlobalSumVector(&evec_tmp[0],NN*NN);
	}
      } 
      // cheating a bit. It is better to sort instead of just reversing it, but the document of the routine says evals are sorted in increasing order. qr gives evals in decreasing order.
      for(int i=0;i<NN;i++){
	for(int j=0;j<NN;j++)
	  Qt[(NN-1-i)*N2+j]=evec_tmp[i*NN + j];
	lmd [NN-1-i]=evals_tmp[i];
      }

      std::cout << GridLogMessage << "diagonalize_lapack complete\n";
    }
#undef LAPACK_INT 
#endif


    void diagonalize(std::vector<RealD>& lmd,
		     std::vector<RealD>& lme, 
		     int N2,
		     int N1,
		     std::vector<RealD>& Qt,
		     GridBase *grid)
    {

#ifdef USE_LAPACK_IRL
    const int check_lapack=0; // just use lapack if 0, check against lapack if 1

    if(!check_lapack)
	return diagonalize_lapack(lmd,lme,N2,N1,Qt,grid);

	std::vector <RealD> lmd2(N1);
	std::vector <RealD> lme2(N1);
	std::vector<RealD> Qt2(N1*N1);
         for(int k=0; k<N1; ++k){
	    lmd2[k] = lmd[k];
	    lme2[k] = lme[k];
	  }
         for(int k=0; k<N1*N1; ++k)
	Qt2[k] = Qt[k];

//	diagonalize_lapack(lmd2,lme2,Nm2,Nm,Qt,grid);
#endif

      int Niter = 10000*N1;
      int kmin = 1;
      int kmax = N2;
      // (this should be more sophisticated)

      for(int iter=0; ; ++iter){
      if ( (iter+1)%(100*N1)==0) 
      std::cout<<GridLogMessage << "[QL method] Not converged - iteration "<<iter+1<<"\n";

	// determination of 2x2 leading submatrix
	RealD dsub = lmd[kmax-1]-lmd[kmax-2];
	RealD dd = sqrt(dsub*dsub + 4.0*lme[kmax-2]*lme[kmax-2]);
	RealD Dsh = 0.5*(lmd[kmax-2]+lmd[kmax-1] +dd*(dsub/fabs(dsub)));
	// (Dsh: shift)
	
	// transformation
	qr_decomp(lmd,lme,N2,N1,Qt,Dsh,kmin,kmax);
	
	// Convergence criterion (redef of kmin and kamx)
	for(int j=kmax-1; j>= kmin; --j){
	  RealD dds = fabs(lmd[j-1])+fabs(lmd[j]);
	  if(fabs(lme[j-1])+dds > dds){
	    kmax = j+1;
	    goto continued;
	  }
	}
	Niter = iter;
#ifdef USE_LAPACK_IRL
    if(check_lapack){
	const double SMALL=1e-8;
	diagonalize_lapack(lmd2,lme2,N2,N1,Qt2,grid);
	std::vector <RealD> lmd3(N2);
         for(int k=0; k<N2; ++k) lmd3[k]=lmd[k];
        _sort.push(lmd3,N2);
        _sort.push(lmd2,N2);
         for(int k=0; k<N2; ++k){
	    if (fabs(lmd2[k] - lmd3[k]) >SMALL)  std::cout<<GridLogMessage <<"lmd(qr) lmd(lapack) "<< k << ": " << lmd2[k] <<" "<< lmd3[k] <<std::endl;
//	    if (fabs(lme2[k] - lme[k]) >SMALL)  std::cout<<GridLogMessage <<"lme(qr)-lme(lapack) "<< k << ": " << lme2[k] - lme[k] <<std::endl;
	  }
         for(int k=0; k<N1*N1; ++k){
//	    if (fabs(Qt2[k] - Qt[k]) >SMALL)  std::cout<<GridLogMessage <<"Qt(qr)-Qt(lapack) "<< k << ": " << Qt2[k] - Qt[k] <<std::endl;
	}
    }
#endif
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
      std::cout<<GridLogMessage << "[QL method] Error - Too many iteration: "<<Niter<<"\n";
      abort();
    }

#if 1
    template<typename T>
    static RealD normalise(T& v) 
    {
      RealD nn = norm2(v);
      nn = sqrt(nn);
      v = v * (1.0/nn);
      return nn;
    }

    void orthogonalize(Field& w,
		       BasisFieldVector<Field>& evec,
		       int k)
    {
      double t0=-usecond()/1e6;

      evec.orthogonalize(w,k);

      normalise(w);
      t0+=usecond()/1e6;
      OrthoTime +=t0;
    }

    void setUnit_Qt(int Nm, std::vector<RealD> &Qt) {
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

    void calc(std::vector<RealD>& eval,
	      BasisFieldVector<Field>& evec,
	      const Field& src,
	      int& Nconv,
	      bool reverse,
	      int SkipTest)
      {

	GridBase *grid = evec._v[0]._grid;//evec.get(0 + evec_offset)._grid;
	assert(grid == src._grid);

	std::cout<<GridLogMessage << " -- Nk = " << Nk << " Np = "<< Np << std::endl;
	std::cout<<GridLogMessage << " -- Nm = " << Nm << std::endl;
	std::cout<<GridLogMessage << " -- size of eval   = " << eval.size() << std::endl;
	std::cout<<GridLogMessage << " -- size of evec  = " << evec.size() << std::endl;
	
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
	    std::cout << GridLogMessage << " Approximation of largest eigenvalue: " << evalMaxApprox << std::endl;
	    src_n = tmp;
	  }
	}
	
	std::vector<RealD> lme(Nm);  
	std::vector<RealD> lme2(Nm);
	std::vector<RealD> eval2(Nm);
	std::vector<RealD> eval2_copy(Nm);
	std::vector<RealD> Qt(Nm*Nm);


	Field f(grid);
	Field v(grid);
  
	int k1 = 1;
	int k2 = Nk;

	Nconv = 0;

	RealD beta_k;
  
	// Set initial vector
	evec[0] = src;
	normalise(evec[0]);
	std:: cout<<GridLogMessage <<"norm2(evec[0])= " << norm2(evec[0])<<std::endl;
	
	// Initial Nk steps
	OrthoTime=0.;
	double t0=usecond()/1e6;
	for(int k=0; k<Nk; ++k) step(eval,lme,evec,f,Nm,k);
	double t1=usecond()/1e6;
	std::cout<<GridLogMessage <<"IRL::Initial steps: "<<t1-t0<< "seconds"<<std::endl; t0=t1;
	std::cout<<GridLogMessage <<"IRL::Initial steps:OrthoTime "<<OrthoTime<< "seconds"<<std::endl;
	t1=usecond()/1e6;

	// Restarting loop begins
	for(int iter = 0; iter<Niter; ++iter){
	  
	  std::cout<<GridLogMessage<<"\n Restart iteration = "<< iter << std::endl;
	  
	  // 
	  // Rudy does a sort first which looks very different. Getting fed up with sorting out the algo defs.
	  // We loop over 
	  //
	  OrthoTime=0.;
	  for(int k=Nk; k<Nm; ++k) step(eval,lme,evec,f,Nm,k);
	  t1=usecond()/1e6;
	  std::cout<<GridLogMessage <<"IRL:: "<<Np <<" steps: "<<t1-t0<< "seconds"<<std::endl; t0=t1;
	  std::cout<<GridLogMessage <<"IRL::Initial steps:OrthoTime "<<OrthoTime<< "seconds"<<std::endl;
	  f *= lme[Nm-1];
	  
	  t1=usecond()/1e6;

	  
	  // getting eigenvalues
	  for(int k=0; k<Nm; ++k){
	    eval2[k] = eval[k+k1-1];
	    lme2[k] = lme[k+k1-1];
	  }
	  setUnit_Qt(Nm,Qt);
	  diagonalize(eval2,lme2,Nm,Nm,Qt,grid);
	  t1=usecond()/1e6;
	  std::cout<<GridLogMessage <<"IRL:: diagonalize: "<<t1-t0<< "seconds"<<std::endl; t0=t1;
	  
	  // sorting
	  eval2_copy = eval2;

	  _sort.push(eval2,Nm);
	  t1=usecond()/1e6;
	  std::cout<<GridLogMessage <<"IRL:: eval sorting: "<<t1-t0<< "seconds"<<std::endl; t0=t1;
	  
	  // Implicitly shifted QR transformations
	  setUnit_Qt(Nm,Qt);
	  for(int ip=0; ip<k2; ++ip){
	    std::cout<<GridLogMessage << "eval "<< ip << " "<< eval2[ip] << std::endl;
	  }

	  for(int ip=k2; ip<Nm; ++ip){ 
	    std::cout<<GridLogMessage << "qr_decomp "<< ip << " "<< eval2[ip] << std::endl;
	    qr_decomp(eval,lme,Nm,Nm,Qt,eval2[ip],k1,Nm);
	    
	  }
	  t1=usecond()/1e6;
	  std::cout<<GridLogMessage <<"IRL::qr_decomp: "<<t1-t0<< "seconds"<<std::endl; t0=t1;
	  assert(k2<Nm);
	  

	  assert(k2<Nm);
	  assert(k1>0);
	  evec.rotate(Qt,k1-1,k2+1,0,Nm,Nm);
	  
	  t1=usecond()/1e6;
	  std::cout<<GridLogMessage <<"IRL::QR rotation: "<<t1-t0<< "seconds"<<std::endl; t0=t1;
	  fflush(stdout);
	  
	  // Compressed vector f and beta(k2)
	  f *= Qt[Nm-1+Nm*(k2-1)];
	  f += lme[k2-1] * evec[k2];
	  beta_k = norm2(f);
	  beta_k = sqrt(beta_k);
	  std::cout<<GridLogMessage<<" beta(k) = "<<beta_k<<std::endl;
	  
	  RealD betar = 1.0/beta_k;
	  evec[k2] = betar * f;
	  lme[k2-1] = beta_k;
	  
	  // Convergence test
	  for(int k=0; k<Nm; ++k){    
	    eval2[k] = eval[k];
	    lme2[k] = lme[k];

	    std::cout<<GridLogMessage << "eval2[" << k << "] = " << eval2[k] << std::endl;
	  }
	  setUnit_Qt(Nm,Qt);
	  diagonalize(eval2,lme2,Nk,Nm,Qt,grid);
	  t1=usecond()/1e6;
	  std::cout<<GridLogMessage <<"IRL::diagonalize: "<<t1-t0<< "seconds"<<std::endl; t0=t1;
	  
	  
	  Nconv = 0;
	  
	  if (iter >= Nminres) {
	    std::cout << GridLogMessage << "Rotation to test convergence " << std::endl;
	    
	    Field ev0_orig(grid);
	    ev0_orig = evec[0];
	    
	    evec.rotate(Qt,0,Nk,0,Nk,Nm);
	    
	    {
	      std::cout << GridLogMessage << "Test convergence" << std::endl;
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
		std::cout<<GridLogMessage << "[" << std::setw(3)<< std::setiosflags(std::ios_base::right) <<j<<"] "
			 <<"eval = "<<std::setw(25)<< std::setiosflags(std::ios_base::left)<< eval2[j] << " (" << eval2_copy[j] << ")"
			 <<" |H B[i] - eval[i]B[i]|^2 / evalMaxApprox^2 " << std::setw(25)<< std::setiosflags(std::ios_base::right)<< vv
			 <<" "<< vnum/(sqrt(vden)*sqrt(vv0))
			 << " norm(B["<<j<<"])="<< vden <<std::endl;
		
		// change the criteria as evals are supposed to be sorted, all evals smaller(larger) than Nstop should have converged
		if((vv<eresid*eresid) && (j == Nconv) ){
		  Nconv+=SkipTest;
		}
	      }
	      
	      // test if we converged, if so, terminate
	      t1=usecond()/1e6;
	      std::cout<<GridLogMessage <<"IRL::convergence testing: "<<t1-t0<< "seconds"<<std::endl; t0=t1;
	      
	      std::cout<<GridLogMessage<<" #modes converged: "<<Nconv<<std::endl;
	      
	      if( Nconv>=Nstop || beta_k < betastp){
		goto converged;
	      }
	      
	      std::cout << GridLogMessage << "Rotate back" << std::endl;
	      //B[j] +=Qt[k+_Nm*j] * _v[k]._odata[ss];
	      {
		Eigen::MatrixXd qm = Eigen::MatrixXd::Zero(Nk,Nk);
		for (int k=0;k<Nk;k++)
		  for (int j=0;j<Nk;j++)
		    qm(j,k) = Qt[k+Nm*j];
		GridStopWatch timeInv;
		timeInv.Start();
		Eigen::MatrixXd qmI = qm.inverse();
		timeInv.Stop();
		std::vector<RealD> QtI(Nm*Nm);
		for (int k=0;k<Nk;k++)
		  for (int j=0;j<Nk;j++)
		    QtI[k+Nm*j] = qmI(j,k);
		
		RealD res_check_rotate_inverse = (qm*qmI - Eigen::MatrixXd::Identity(Nk,Nk)).norm(); // sqrt( |X|^2 )
		assert(res_check_rotate_inverse < 1e-7);
		evec.rotate(QtI,0,Nk,0,Nk,Nm);
		
		axpy(ev0_orig,-1.0,evec[0],ev0_orig);
		std::cout << GridLogMessage << "Rotation done (in " << timeInv.Elapsed() << " = " << timeInv.useconds() << " us" <<
		  ", error = " << res_check_rotate_inverse << 
		  "); | evec[0] - evec[0]_orig | = " << ::sqrt(norm2(ev0_orig)) << std::endl;
	      }
	    }
	  } else {
	    std::cout << GridLogMessage << "iter < Nminres: do not yet test for convergence\n";
	  } // end of iter loop
	}

	std::cout<<GridLogMessage<<"\n NOT converged.\n";
	abort();
	
      converged:

	if (SkipTest == 1) {
	  eval = eval2;
	} else {

	  // test quickly
	  for (int j=0;j<Nstop;j+=SkipTest) {
	    std::cout<<GridLogMessage << "Eigenvalue[" << j << "] = " << eval2[j] << " (" << eval2_copy[j] << ")" << std::endl;
	  }

	  eval2_copy.resize(eval2.size());
	  eval = eval2_copy;
	}

	evec.sortInPlace(eval,reverse);

	{
	  
	 // test
	 for (int j=0;j<Nstop;j++) {
	   std::cout<<GridLogMessage << " |e[" << j << "]|^2 = " << norm2(evec[j]) << std::endl;
	 }
       }
       
       //_sort.push(eval,evec,Nconv);
       //evec.sort(eval,Nconv);
       
       std::cout<<GridLogMessage << "\n Converged\n Summary :\n";
       std::cout<<GridLogMessage << " -- Iterations  = "<< Nconv  << "\n";
       std::cout<<GridLogMessage << " -- beta(k)     = "<< beta_k << "\n";
       std::cout<<GridLogMessage << " -- Nconv       = "<< Nconv  << "\n";
      }
#endif

    };

}
#endif

