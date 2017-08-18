    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/ImplicitlyRestartedLanczos.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Chulwoo Jung <chulwoo@bnl.gov>

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
#ifndef GRID_IRL_CJ_H
#define GRID_IRL_CJ_H

#include <string.h> //memset

#ifdef USE_LAPACK
#ifdef USE_MKL
#include<mkl_lapack.h>
#else
void LAPACK_dstegr(char *jobz, char *range, int *n, double *d, double *e,
                   double *vl, double *vu, int *il, int *iu, double *abstol,
                   int *m, double *w, double *z, int *ldz, int *isuppz,
                   double *work, int *lwork, int *iwork, int *liwork,
                   int *info);
//#include <lapacke/lapacke.h>
#endif
#endif

//#include <Grid/algorithms/densematrix/DenseMatrix.h>
//#include <Grid/algorithms/iterative/EigenSort.h>

// eliminate temorary vector in calc()
#define MEM_SAVE

namespace Grid {

/////////////////////////////////////////////////////////////
// Implicitly restarted lanczos
/////////////////////////////////////////////////////////////

// creating a seaprate instance to avoid conflicts for the time being

template<class Field> 
    class ImplicitlyRestartedLanczosCJ {

    const RealD small = 1.0e-16;
public:       
    int lock;
    int get;
    int Niter;
    int converged;

    int Nstop;   // Number of evecs checked for convergence
    int Nk;      // Number of converged sought
    int Np;      // Np -- Number of spare vecs in kryloc space
    int Nm;      // Nm -- total number of vectors


    RealD OrthoTime;

    RealD eresid;

    SortEigen<Field> _sort;

    LinearOperatorBase<Field> &_Linop;

    OperatorFunction<Field>   &_poly;

    /////////////////////////
    // Constructor
    /////////////////////////
    void init(void){};
    void Abort(int ff, DenseVector<RealD> &evals,  DenseVector<DenseVector<RealD> > &evecs);

    ImplicitlyRestartedLanczosCJ(
				LinearOperatorBase<Field> &Linop, // op
			       OperatorFunction<Field> & poly,   // polynmial
			       int _Nstop, // sought vecs
			       int _Nk, // sought vecs
			       int _Nm, // spare vecs
			       RealD _eresid, // resid in lmdue deficit 
			       int _Niter) : // Max iterations
      _Linop(Linop),
      _poly(poly),
      Nstop(_Nstop),
      Nk(_Nk),
      Nm(_Nm),
      eresid(_eresid),
      Niter(_Niter)
    { 
      Np = Nm-Nk; assert(Np>0);
    };

    ImplicitlyRestartedLanczosCJ(
				LinearOperatorBase<Field> &Linop, // op
			       OperatorFunction<Field> & poly,   // polynmial
			       int _Nk, // sought vecs
			       int _Nm, // spare vecs
			       RealD _eresid, // resid in lmdue deficit 
			       int _Niter) : // Max iterations
      _Linop(Linop),
      _poly(poly),
      Nstop(_Nk),
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
      std::cout<<GridLogMessage << "RitzMatrix "<<std::endl;
      for(int i=0;i<k;i++){
	_poly(_Linop,evec[i],w);
	std::cout<<GridLogMessage << "["<<i<<"] ";
	for(int j=0;j<k;j++){
	  ComplexD in = innerProduct(evec[j],w);
	  if ( fabs((double)i-j)>1 ) { 
	    if (abs(in) >1.0e-9 )  { 
	      std::cout<<GridLogMessage<<"oops"<<std::endl;
	      abort();
	    } else 
	      std::cout<<GridLogMessage << " 0 ";
	  } else { 
	    std::cout<<GridLogMessage << " "<<in<<" ";
	  }
	}
	std::cout<<GridLogMessage << std::endl;
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
	      Field& w,int Nm,int k, RealD &norm)
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
	 norm=beta;

	std::cout<<GridLogMessage << "alpha = " << zalph << " beta "<<beta<<std::endl;
      const RealD tiny = 1.0e-20;
      if ( beta < tiny ) { 
	std::cout<<GridLogMessage << " beta is tiny "<<beta<<std::endl;
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


#ifdef USE_LAPACK
#ifdef USE_MKL
#define LAPACK_INT MKL_INT
#else
#define LAPACK_INT long long
#endif
    void diagonalize_lapack(DenseVector<RealD>& lmd,
		     DenseVector<RealD>& lme, 
		     int N1,
		     int N2,
		     DenseVector<RealD>& Qt,
		     GridBase *grid){
  const int size = Nm;
//  tevals.resize(size);
//  tevecs.resize(size);
  LAPACK_INT NN = N1;
//  double evals_tmp[NN];
//  double evec_tmp[NN][NN];
  std::vector<double> evals_tmp(NN);
  std::vector<double> evec_tmp(NN*NN);
  memset(evec_tmp.data(),0,sizeof(double)*NN*NN);

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
//  LAPACK_INT lwork = ( (18*NN) > (1+4*NN+NN*NN)? (18*NN):(1+4*NN+NN*NN)) ;
  LAPACK_INT lwork =  1+(18*NN) ;
  LAPACK_INT liwork =  3+NN*10 ;
//  LAPACK_INT iwork[liwork];
//  double work[lwork];
//  LAPACK_INT isuppz[2*NN];
  std::vector<LAPACK_INT> iwork(liwork);
  std::vector<double> work(lwork);
  std::vector<LAPACK_INT> isuppz(2*NN);
  char jobz = 'V'; // calculate evals & evecs
  char range = 'I'; // calculate all evals
  //    char range = 'A'; // calculate all evals
  char uplo = 'U'; // refer to upper half of original matrix
  char compz = 'I'; // Compute eigenvectors of tridiagonal matrix
  int ifail[NN];
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
      memset(evals_tmp.data(),0,sizeof(double)*NN);
      if ( il <= NN){
        printf("total=%d node=%d il=%d iu=%d\n",total,node,il,iu);
#ifdef USE_MKL
        dstegr(&jobz, &range, &NN,
#else
        LAPACK_dstegr(&jobz, &range, &NN,
#endif
            DD.data(), EE.data(),
            &vl, &vu, &il, &iu, // these four are ignored if second parameteris 'A'
            &tol, // tolerance
            &evals_found, evals_tmp.data(), evec_tmp.data(), &NN,
            isuppz.data(),
            work.data(), &lwork, iwork.data(), &liwork,
            &info);
        for (int i = iu-1; i>= il-1; i--){
          printf("node=%d evals_found=%d evals_tmp[%d] = %g\n",node,evals_found, i - (il-1),evals_tmp[i - (il-1)]);
          evals_tmp[i] = evals_tmp[i - (il-1)];
          if (il>1) evals_tmp[i-(il-1)]=0.;
          for (int j = 0; j< NN; j++){
            evec_tmp[i*NN+j] = evec_tmp[(i - (il-1))*NN+j];
            if (il>1) evec_tmp[(i-(il-1))*NN+j]=0.;
          }
        }
      }
      {
         grid->GlobalSumVector(evals_tmp.data(),NN);
         grid->GlobalSumVector(evec_tmp.data(),NN*NN);
      }
    } 
// cheating a bit. It is better to sort instead of just reversing it, but the document of the routine says evals are sorted in increasing order. qr gives evals in decreasing order.
  for(int i=0;i<NN;i++){
    for(int j=0;j<NN;j++)
      Qt[(NN-1-i)*N2+j]=evec_tmp[i*NN+j];
      lmd [NN-1-i]=evals_tmp[i];
  }
}
#undef LAPACK_INT 
#endif
    void diagonalize(DenseVector<RealD>& lmd,
		     DenseVector<RealD>& lme, 
		     int N2,
		     int N1,
		     DenseVector<RealD>& Qt,
		     GridBase *grid)
    {

#ifdef USE_LAPACK
    const int check_lapack=0; // just use lapack if 0, check against lapack if 1

    if(!check_lapack)
	return diagonalize_lapack(lmd,lme,N2,N1,Qt,grid);

	DenseVector <RealD> lmd2(N1);
	DenseVector <RealD> lme2(N1);
	DenseVector<RealD> Qt2(N1*N1);
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
#ifdef USE_LAPACK
    if(check_lapack){
	const double SMALL=1e-8;
	diagonalize_lapack(lmd2,lme2,N2,N1,Qt2,grid);
	DenseVector <RealD> lmd3(N2);
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
      double t0=-usecond()/1e6;
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
      t0+=usecond()/1e6;
      OrthoTime +=t0;
    }

    void setUnit_Qt(int Nm, DenseVector<RealD> &Qt) {
      for(int i=0; i<Qt.size(); ++i) Qt[i] = 0.0;
      for(int k=0; k<Nm; ++k) Qt[k + k*Nm] = 1.0;
    }

// needs more memory
void Rotate0( 
//	int _Nm,
	DenseVector<RealD>& Qt,
	DenseVector<Field>& evec,
	int j0, int j1,
	int _Nk
	)
{
	GridBase *grid = evec[0]._grid;
	DenseVector<Field>  B(Nm,grid); // waste of space replicating
	if (0) {   // old implementation without blocking
	  for(int i=0; i<(Nm); ++i) B[i] = 0.0;
	  
	  for(int j=j0; j<j1; ++j){
	    for(int k=0; k<Nm; ++k){
	    B[j].checkerboard = evec[k].checkerboard;
	      B[j] += Qt[k+Nm*j] * evec[k];
	    }
	  }
	}
	{
	for(int i=0; i<(Nm); ++i) {
		B[i] = 0.0;
	  	B[i].checkerboard = evec[0].checkerboard;
	}
	int j_block = 24; int k_block=24;
PARALLEL_FOR_LOOP
	for(int ss=0;ss < grid->oSites();ss++){
		for(int jj=j0; jj<j1; jj += j_block)
		for(int kk=0; kk<Nm; kk += k_block)
		for(int j=jj; (j<(j1)) && j<(jj+j_block); ++j){
			for(int k=kk; (k<Nm) && k<(kk+k_block) ; ++k){
			    B[j]._odata[ss] +=Qt[k+Nm*j] * evec[k]._odata[ss]; 
			}
		}
	}
	}
	for(int j=j0; j<j1; ++j) evec[j] = B[j];
}

void Rotate( 
//	int _Nm,
	DenseVector<RealD>& Qt,
	DenseVector<Field>& evec,
//int k1, int k2
	int j0, int j1,
	int _Nk
	)
{
	GridBase *grid = evec[0]._grid;
	typedef typename Field::vector_object vobj;
	assert(j0>=0);
	assert(j1<Nm);

#pragma omp parallel
	{
		std::vector < vobj > B(Nm);
#pragma omp for
		for(int ss=0;ss < grid->oSites();ss++){
			for(int j=j0; j<j1; ++j)
				zeroit(B[j]);
			for(int j=j0; j<j1; ++j){
				for(int k=0; k<_Nk ; ++k){
					B[j] +=Qt[k+Nm*j] * evec[k]._odata[ss];
				}
			}
			for(int j=j0; j<j1; ++j){
				evec[j]._odata[ss] = B[j];
			}
		}
	}
}

    void Rotate2( 
//	int Nm,
	DenseVector<RealD>& Qt,
	DenseVector<Field>& evec,
	int k1, int k2
	)
{
	GridBase *grid = evec[0]._grid;
	int j_block = 24; int k_block=24;
	
	assert(k2<Nm);
	assert(k1>0);
	int thr=GridThread::GetThreads();
	int n_field = 1;
	int each = 1;
	if( (Nm*thr)>(grid->oSites()) ) {
		each = (grid->oSites())/Nm ;
		n_field = thr/each + 1;
	}
	std::cout<<GridLogMessage << "thr = " << thr << " n_field= "<< n_field << " each= "<<each << std::endl;
	DenseVector<Field>  B(n_field,grid); 
PARALLEL_FOR_LOOP
	for(int ss=0;ss < grid->oSites();ss++){
		int me = GridThread::ThreadBarrier();
		int i_field = me / each;
		int k_field = me % each;
		assert(i_field < n_field);
		std::cout<<GridLogMessage << "me = " << me << " i_field= "<< i_field << " k_field= "<<k_field << std::endl;
//		printf("thr=%d ss=%d me=%d\n",thr,ss,me);fflush(stdout);
//		assert(Nm*thr<grid->oSites());
		for(int j=0; j<Nm; ++j) B[i_field]._odata[j+Nm*k_field]=0.;
		for(int j=k1-1; j<(k2+1); ++j){
			for(int k=0; k<Nm ; ++k){
			    B[i_field]._odata[j+Nm*k_field] +=Qt[k+Nm*j] * evec[k]._odata[ss]; 
			}
		}
		for(int j=k1-1; j<(k2+1); ++j){
			evec[j]._odata[ss] = B[i_field]._odata[j+Nm*k_field];
		}
	}
}

void ConvCheck0( int _Nk, 
	DenseVector<RealD>& Qt,
	DenseVector<Field>& evec,
	DenseVector<RealD> &eval2,
	DenseVector<int>   &Iconv,
	int &Nconv)
{

	GridBase *grid = evec[0]._grid;
	DenseVector<Field>  B(Nm,grid); // waste of space replicating
	Field v(grid);
if (0) {
	  for(int k = 0; k<_Nk; ++k) B[k]=0.0;
	  
	  for(int j = 0; j<_Nk; ++j){
	    for(int k = 0; k<_Nk; ++k){
	    B[j].checkerboard = evec[k].checkerboard;
	      B[j] += Qt[k+j*Nm] * evec[k];
	    }
	    std::cout<<GridLogMessage << "norm(B["<<j<<"])="<<norm2(B[j])<<std::endl;
	  }
}
if (1) {
	for(int i=0; i<(_Nk+1); ++i) {
		B[i] = 0.0;
	  	B[i].checkerboard = evec[0].checkerboard;
	}

	int j_block = 24; int k_block=24;
PARALLEL_FOR_LOOP
	for(int ss=0;ss < grid->oSites();ss++){
	for(int jj=0; jj<_Nk; jj += j_block)
	for(int kk=0; kk<_Nk; kk += k_block)
	for(int j=jj; (j<_Nk) && j<(jj+j_block); ++j){
	for(int k=kk; (k<_Nk) && k<(kk+k_block) ; ++k){
	    B[j]._odata[ss] +=Qt[k+Nm*j] * evec[k]._odata[ss]; 
	}
	}
	}
}

	  Nconv = 0;
	  //	  std::cout<<GridLogMessage << std::setiosflags(std::ios_base::scientific);
	  for(int i=0; i<_Nk; ++i){

//	    _poly(_Linop,B[i],v);
	    _Linop.HermOp(B[i],v);
	    
	    RealD vnum = real(innerProduct(B[i],v)); // HermOp.
	    RealD vden = norm2(B[i]);
	    RealD vv0 = norm2(v);
	    eval2[i] = vnum/vden;
	    v -= eval2[i]*B[i];
	    RealD vv = norm2(v)/vv0;
	    
	    std::cout.precision(13);
	    std::cout<<GridLogMessage << "[" << std::setw(3)<< std::setiosflags(std::ios_base::right) <<i<<"] ";
	    std::cout<<"eval = "<<std::setw(25)<< std::setiosflags(std::ios_base::left)<< eval2[i];
	    std::cout<<"|H B[i] - eval[i]B[i]|^2/|H B[i]|^2 "<< std::setw(25)<< std::setiosflags(std::ios_base::right)<< vv;
	    std::cout<<" "<< vnum/(sqrt(vden)*sqrt(vv0)) << std::endl;
	    
	// change the criteria as evals are supposed to be sorted, all evals smaller(larger) than Nstop should have converged
	    if((vv<eresid*eresid) && (i == Nconv) ){
	      Iconv[Nconv] = i;
	      ++Nconv;
	    }

	  }  // i-loop end
	  //	  std::cout<<GridLogMessage << std::resetiosflags(std::ios_base::scientific);
}

void FinalCheck( int _Nk,
	DenseVector<RealD>& eval,
	DenseVector<Field>& evec
	)
	{
	GridBase *grid = evec[0]._grid;
	Field v(grid);
	Field B(grid);
		for(int j = 0; j<_Nk; ++j){
		    std::cout<<GridLogMessage << "norm(evec["<<j<<"])="<<norm2(evec[j])<<std::endl;
		    _Linop.HermOp(evec[j],v);
		    
		    RealD vnum = real(innerProduct(evec[j],v)); // HermOp.
		    RealD vden = norm2(evec[j]);
		    RealD vv0 = norm2(v);
		    RealD eval2 = vnum/vden;
		    v -= eval2*evec[j];
		    RealD vv = norm2(v)/vv0;
	    
		    std::cout.precision(13);
		    std::cout<<GridLogMessage << "[" << std::setw(3)<< std::setiosflags(std::ios_base::right) <<j<<"] ";
		    std::cout<<"eval = "<<std::setw(25)<< std::setiosflags(std::ios_base::left)<< eval2;
	    	std::cout<<"|H B[i] - eval[i]B[i]|^2/|H B[i]|^2 "<< std::setw(25)<< std::setiosflags(std::ios_base::right)<< vv;
		    std::cout<<" "<< vnum/(sqrt(vden)*sqrt(vv0)) << std::endl;
			eval[j] = eval2;
	    
		}
	}

void ConvCheck( int _Nk,
	DenseVector<RealD>& Qt,
	DenseVector<Field>& evec,
	DenseVector<RealD> &eval2,
	DenseVector<int>   &Iconv,
	int &Nconv)
	{
	GridBase *grid = evec[0]._grid;
	Field v(grid);
		Field B(grid);
		Nconv = 0;
		for(int j = 0; j<_Nk; ++j){
			B=0.;
			B.checkerboard = evec[0].checkerboard;
		    for(int k = 0; k<_Nk; ++k){
		    	B += Qt[k+j*Nm] * evec[k];
		    }
		    std::cout<<GridLogMessage << "norm(B["<<j<<"])="<<norm2(B)<<std::endl;
//		    _poly(_Linop,B,v);
		    _Linop.HermOp(B,v);
		    
		    RealD vnum = real(innerProduct(B,v)); // HermOp.
		    RealD vden = norm2(B);
		    RealD vv0 = norm2(v);
		    eval2[j] = vnum/vden;
		    v -= eval2[j]*B;
		    RealD vv = norm2(v);
	    
		    std::cout.precision(13);
		    std::cout<<GridLogMessage << "[" << std::setw(3)<< std::setiosflags(std::ios_base::right) <<j<<"] ";
		    std::cout<<"eval = "<<std::setw(25)<< std::setiosflags(std::ios_base::left)<< eval2[j];
		    std::cout<<"|H B[i] - eval[i]B[i]|^2 "<< std::setw(25)<< std::setiosflags(std::ios_base::right)<< vv;
		    std::cout<<" "<< vnum/(sqrt(vden)*sqrt(vv0)) << std::endl;
	    
// change the criteria as evals are supposed to be sorted, all evals smaller(larger) than Nstop should have converged
		    if((vv<eresid*eresid) && (j == Nconv) ){
		      Iconv[Nconv] = j;
		      ++Nconv;
		    }
		}
	}

void ConvRotate2( int _Nk,
	DenseVector<RealD>& Qt,
	DenseVector<Field>& evec,
	DenseVector<RealD> &eval,
	DenseVector<RealD> &eval2,
	DenseVector<int>   &Iconv,
	int &Nconv)
{
	GridBase *grid = evec[0]._grid;
	int thr=GridThread::GetThreads();
       for(int i=0; i<Nconv; ++i)
         eval[i] = eval2[Iconv[i]];
//	int thr=GridThread::GetThreads();
//	printf("thr=%d\n",thr);
	Field B(grid);
PARALLEL_FOR_LOOP
	for(int ss=0;ss < grid->oSites();ss++){
		int me = GridThread::ThreadBarrier();
		printf("thr=%d ss=%d me=%d\n",thr,ss,me);fflush(stdout);
		assert( (Nm*thr)<grid->oSites());
//		auto B2 = evec[0]._odata[0];
//		std::vector < decltype( B2 ) > B(Nm,B2);
		for(int j=0; j<Nconv; ++j) B._odata[Iconv[j]+Nm*me]=0.;
		for(int j=0; j<Nconv; ++j){
			for(int k=0; k<_Nk ; ++k){
			    B._odata[Iconv[j]+Nm*me] +=Qt[k+Nm*Iconv[j]] * evec[k]._odata[ss]; 
			}
		}
		for(int j=0; j<Nconv; ++j){
			evec[j]._odata[ss] = B._odata[Iconv[j]+Nm*me];
		}
	}
}

void ConvRotate( int _Nk,
	DenseVector<RealD>& Qt,
	DenseVector<Field>& evec,
	DenseVector<RealD> &eval,
	DenseVector<RealD> &eval2,
	DenseVector<int>   &Iconv,
	int &Nconv)
{
	GridBase *grid = evec[0]._grid;
	typedef typename Field::vector_object vobj;
#pragma omp parallel
	{
		std::vector < vobj > B(Nm);
#pragma omp for
	for(int ss=0;ss < grid->oSites();ss++){
//		for(int j=0; j<Nconv; ++j) B[j]=0.;
		for(int j=0; j<Nconv; ++j){
			for(int k=0; k<_Nk ; ++k){
			    B[j] +=Qt[k+Nm*Iconv[j]] * evec[k]._odata[ss]; 
			}
		}
		for(int j=0; j<Nconv; ++j){
			evec[j]._odata[ss] = B[j];
		}
	}
	}
	for(int i=0; i<Nconv; ++i) eval[i] = eval2[Iconv[i]];
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

// alternate implementation for minimizing memory usage.  May affect the performance
    void calc(
		DenseVector<RealD>& eval,
	      DenseVector<Field>& evec,
	      const Field& src,
	      int& Nconv)
      {

	GridBase *grid = evec[0]._grid;
	assert(grid == src._grid);

	std::cout<<GridLogMessage << " -- Nk = " << Nk << " Np = "<< Np << std::endl;
	std::cout<<GridLogMessage << " -- Nm = " << Nm << std::endl;
	std::cout<<GridLogMessage << " -- size of eval   = " << eval.size() << std::endl;
	std::cout<<GridLogMessage << " -- size of evec  = " << evec.size() << std::endl;
	
	assert(Nm == evec.size() && Nm == eval.size());
	
	DenseVector<RealD> lme(Nm);  
	DenseVector<RealD> lme2(Nm);
	DenseVector<RealD> eval2(Nm);
	DenseVector<RealD> Qt(Nm*Nm);
	DenseVector<int>   Iconv(Nm);

	
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
	std:: cout<<GridLogMessage <<"norm2(src)= " << norm2(src)<<std::endl;
// << src._grid  << std::endl;
	normalise(evec[0]);
	std:: cout<<GridLogMessage <<"norm2(evec[0])= " << norm2(evec[0]) <<std::endl;
// << evec[0]._grid << std::endl;
	
	// Initial Nk steps
	OrthoTime=0.;
	double t0=usecond()/1e6;
	RealD norm; // sqrt norm of last vector
	for(int k=0; k<Nk; ++k) step(eval,lme,evec,f,Nm,k,norm);
	double t1=usecond()/1e6;
	std::cout<<GridLogMessage <<"IRL::Initial steps: "<<t1-t0<< "seconds"<<std::endl; t0=t1;
	std::cout<<GridLogMessage <<"IRL::Initial steps:OrthoTime "<<OrthoTime<< "seconds"<<std::endl;
//	std:: cout<<GridLogMessage <<"norm2(evec[1])= " << norm2(evec[1]) << std::endl;
//	std:: cout<<GridLogMessage <<"norm2(evec[2])= " << norm2(evec[2]) << std::endl;
	RitzMatrix(evec,Nk);
	t1=usecond()/1e6;
	std::cout<<GridLogMessage <<"IRL::RitzMatrix: "<<t1-t0<< "seconds"<<std::endl; t0=t1;
	for(int k=0; k<Nk; ++k){
//	std:: cout<<GridLogMessage <<"eval " << k << " " <<eval[k] << std::endl;
//	std:: cout<<GridLogMessage <<"lme " << k << " " << lme[k] << std::endl;
	}

	// Restarting loop begins
	for(int iter = 0; iter<Niter; ++iter){

	  std::cout<<GridLogMessage<<"\n Restart iteration = "<< iter << std::endl;

	  // 
	  // Rudy does a sort first which looks very different. Getting fed up with sorting out the algo defs.
	  // We loop over 
	  //
	OrthoTime=0.;
        for(int k=Nk; k<Nm; ++k) step(eval,lme,evec,f,Nm,k,norm);
	t1=usecond()/1e6;
	std::cout<<GridLogMessage <<"IRL:: "<<Np <<" steps: "<<t1-t0<< "seconds"<<std::endl; t0=t1;
	std::cout<<GridLogMessage <<"IRL::Initial steps:OrthoTime "<<OrthoTime<< "seconds"<<std::endl;
	  f *= lme[Nm-1];

	  RitzMatrix(evec,k2);
	t1=usecond()/1e6;
	std::cout<<GridLogMessage <<"IRL:: RitzMatrix: "<<t1-t0<< "seconds"<<std::endl; t0=t1;
	  
	  // getting eigenvalues
	  for(int k=0; k<Nm; ++k){
	    eval2[k] = eval[k+k1-1];
	    lme2[k] = lme[k+k1-1];
	  }
	  setUnit_Qt(Nm,Qt);
	  diagonalize(eval2,lme2,Nm,Nm,Qt,grid);
	t1=usecond()/1e6;
	std::cout<<GridLogMessage <<"IRL:: diagonalize: "<<t1-t0<< "seconds"<<std::endl; t0=t1;
 	int prelNconv=0;
        for(int k=0; k<Nm; ++k){ //check all k's because QR can permutate eigenvectors
	std::cout<<GridLogMessage <<"IRL:: Prel. conv. Test"<<k <<" " << norm<<" "<< fabs(Qt[Nm-1+Nm*k])<<std::endl;
	// Arbitrarily add factor of 10 for now, to be conservative
	if ( norm*fabs(Qt[Nm-1+Nm*k]) < eresid*10 ) prelNconv++;
        }

	  // sorting
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

//	Uses more temorary
//	Rotate0(Qt,evec,k1,k2,Nm);
// 	Uses minimal temporary, possibly with less speed
	Rotate0(Qt,evec,k1-1,k2+1,Nm);
// 	Try if Rotate() doesn't work
//	Rotate2(Qt,evec,k1,k2);
	t1=usecond()/1e6;
	std::cout<<GridLogMessage <<"IRL::QR rotation: "<<t1-t0<< "seconds"<<std::endl; t0=t1;

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
	  }
	  setUnit_Qt(Nm,Qt);
	  diagonalize(eval2,lme2,Nk,Nm,Qt,grid);
	t1=usecond()/1e6;
	std::cout<<GridLogMessage <<"IRL::diagonalize: "<<t1-t0<< "seconds"<<std::endl; t0=t1;
	  
	if(prelNconv < Nstop)
		    std::cout<<GridLogMessage << "Prel. Convergence test ("<<prelNconv<<") failed, skipping a real(and expnesive) one" <<std::endl;
    else
//	ConvCheck0( Nk, Nm, Qt, evec, eval2, Iconv, Nconv);
	ConvCheck( Nk, Qt, evec, eval2, Iconv, Nconv);
	t1=usecond()/1e6;
	std::cout<<GridLogMessage <<"IRL::convergence testing: "<<t1-t0<< "seconds"<<std::endl; t0=t1;


	  std::cout<<GridLogMessage<<" #modes converged: "<<Nconv<<std::endl;

	  if( Nconv>=Nstop ){
	    goto converged;
	  }
	} // end of iter loop
	
	std::cout<<GridLogMessage<<"\n NOT converged.\n";
	abort();
	
      converged:
       // Sorting
       eval.resize(Nconv);
#ifndef MEM_SAVE
       evec.resize(Nconv,grid);
       for(int i=0; i<Nconv; ++i){
         eval[i] = eval2[Iconv[i]];
         evec[i] = B[Iconv[i]];
       }
#else
	Rotate0(Qt,evec,0,Nk,Nm);
	FinalCheck( Nk, eval,evec);
#endif
// Skip sorting, as it doubles the memory usage(!) and can be avoided by diagonalizing "right away"

//      _sort.push(eval,evec,Nconv);

      std::cout<<GridLogMessage << "\n Converged\n Summary :\n";
      std::cout<<GridLogMessage << " -- Iterations  = "<< Nconv  << "\n";
      std::cout<<GridLogMessage << " -- beta(k)     = "<< beta_k << "\n";
      std::cout<<GridLogMessage << " -- Nconv       = "<< Nconv  << "\n";
     }


    void EigenSort(DenseVector<double> evals,
		   DenseVector<Field>  evecs){
      int N= evals.size();
      _sort.push(evals,evecs, evals.size(),N);
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
    T tau0 = fabs(sqrt(sig));
    
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
#endif


 };

}
#endif

