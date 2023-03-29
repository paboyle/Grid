    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/ImplicitlyRestartedBlockLanczos.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Yong-Chull Jang <ypj@quark.phy.bnl.gov> 
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
#ifndef GRID_IRBL_H
#define GRID_IRBL_H

#include <string.h> //memset
#ifdef USE_LAPACK
#include <mkl_lapack.h>
#endif

#undef USE_LAPACK
#define Glog std::cout << GridLogMessage 

#ifdef GRID_CUDA
#include "cublas_v2.h"
#endif

#if 0
#define  CUDA_COMPLEX cuDoubleComplex
#define  CUDA_FLOAT double
#define  MAKE_CUDA_COMPLEX make_cuDoubleComplex
#define  CUDA_GEMM cublasZgemm
#else
#define  CUDA_COMPLEX cuComplex
#define  CUDA_FLOAT float
#define  MAKE_CUDA_COMPLEX make_cuComplex
#define  CUDA_GEMM cublasCgemm
#endif

namespace Grid {

////////////////////////////////////////////////////////////////////////////////
// Helper class for sorting the evalues AND evectors by Field
// Use pointer swizzle on vectors SHOULD GET RID OF IT SOON!
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
    std::vector<Field> cpy(lmd.size(),evec[0].Grid());
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

enum class LanczosType { irbl, rbl };

enum IRBLdiagonalisation { 
  IRBLdiagonaliseWithDSTEGR,
  IRBLdiagonaliseWithQR,
  IRBLdiagonaliseWithEigen
};

/////////////////////////////////////////////////////////////
// Implicitly restarted block lanczos
/////////////////////////////////////////////////////////////
template<class Field> 
class ImplicitlyRestartedBlockLanczos {

private:       
  
  std::string cname = std::string("ImplicitlyRestartedBlockLanczos");
  int MaxIter;   // Max iterations
  int Nstop;     // Number of evecs checked for convergence
  int Nu;        // Number of vecs in the unit block
  int Nk;        // Number of converged sought
  int Nm;        // total number of vectors
  int Nblock_k;    // Nk/Nu
  int Nblock_m;    // Nm/Nu
  int Nconv_test_interval; // Number of skipped vectors when checking a convergence
  RealD eresid;
  IRBLdiagonalisation diagonalisation;
  ////////////////////////////////////
  // Embedded objects
  ////////////////////////////////////
           SortEigen<Field> _sort;
  LinearOperatorBase<Field> &_Linop;
  LinearOperatorBase<Field> &_SLinop;//for split
  OperatorFunction<Field> &_poly;
  GridRedBlackCartesian * f_grid;
  GridRedBlackCartesian * sf_grid;
  int mrhs;
  /////////////////////////
  // BLAS objects
  /////////////////////////
#ifdef GRID_CUDA
  cudaError_t cudaStat;
  CUDA_COMPLEX *w_acc, *evec_acc, *c_acc;
#endif
  int Nevec_acc; // Number of eigenvectors stored in the buffer evec_acc
  
  /////////////////////////
  // Constructor
  /////////////////////////
public:       
 int split_test; //test split in the first iteration
 ImplicitlyRestartedBlockLanczos(LinearOperatorBase<Field> &Linop, // op
 				LinearOperatorBase<Field> &SLinop, // op
				GridRedBlackCartesian * FrbGrid,
				GridRedBlackCartesian * SFrbGrid,
				int _mrhs,
                                 OperatorFunction<Field> & poly,   // polynomial
                                 int _Nstop, // really sought vecs
                                 int _Nconv_test_interval, // conv check interval
                                 int _Nu,    // vecs in the unit block
                                 int _Nk,    // sought vecs
                                 int _Nm,    // total vecs
                                 RealD _eresid, // resid in lmd deficit 
                                 int _MaxIter,  // Max iterations
                                 IRBLdiagonalisation _diagonalisation = IRBLdiagonaliseWithEigen)
   : _Linop(Linop),   _SLinop(SLinop),  _poly(poly),sf_grid(SFrbGrid),f_grid(FrbGrid),
      Nstop(_Nstop), Nconv_test_interval(_Nconv_test_interval), mrhs(_mrhs),
      Nu(_Nu), Nk(_Nk), Nm(_Nm), 
      Nblock_m(_Nm/_Nu), Nblock_k(_Nk/_Nu),
      //eresid(_eresid),  MaxIter(10),
      eresid(_eresid),  MaxIter(_MaxIter),
      diagonalisation(_diagonalisation),split_test(0),
      Nevec_acc(_Nu)
  { assert( (Nk%Nu==0) && (Nm%Nu==0) ); };

  ////////////////////////////////
  // Helpers
  ////////////////////////////////
  static RealD normalize(Field& v, int if_print=0) 
  {
    RealD nn = norm2(v);
    nn = sqrt(nn);
    v = v * (1.0/nn);
    return nn;
  }
  
  void orthogonalize(Field& w, std::vector<Field>& evec, int k, int if_print=0)
  {
    typedef typename Field::scalar_type MyComplex;
//    MyComplex ip;
    ComplexD ip;
    
    for(int j=0; j<k; ++j){
      ip = innerProduct(evec[j],w); 
      if(if_print) 
      if( norm(ip)/norm2(w) > 1e-14)
      Glog<<"orthogonalize before: "<<j<<" of "<<k<<" "<< ip <<std::endl;
      w = w - ip * evec[j];
      if(if_print) {
        ip = innerProduct(evec[j],w); 
        if( norm(ip)/norm2(w) > 1e-14)
          Glog<<"orthogonalize after: "<<j<<" of "<<k<<" "<< ip <<std::endl;
      }
    }
    assert(normalize(w,if_print) != 0);
  }
  void reorthogonalize(Field& w, std::vector<Field>& evec, int k)
  {
     orthogonalize(w, evec, k,1);
  }

  void orthogonalize(std::vector<Field>& w, int _Nu, std::vector<Field>& evec, int k, int if_print=0)
  {
    typedef typename Field::scalar_type MyComplex;
    MyComplex ip;
//    ComplexD ip;
    
    for(int j=0; j<k; ++j){
    for(int i=0; i<_Nu; ++i){
      ip = innerProduct(evec[j],w[i]); 
      w[i] = w[i] - ip * evec[j];
    }}
    for(int i=0; i<_Nu; ++i)
    assert(normalize(w[i],if_print) !=0);
  }


  void orthogonalize_blas(std::vector<Field>& w, std::vector<Field>& evec, int R, int do_print=0) 
  { 
#ifdef GRID_CUDA
    Glog << "cuBLAS orthogonalize" << std::endl;
    
    typedef typename Field::vector_object vobj;
    typedef typename vobj::scalar_type scalar_type;
    typedef typename vobj::vector_type vector_type;
    
    typedef typename Field::scalar_type MyComplex;
    
    GridBase *grid = w[0].Grid();
    const uint64_t sites = grid->lSites();

    int Nbatch = R/Nevec_acc;
    assert( R%Nevec_acc == 0 );
//    Glog << "nBatch, Nevec_acc, R, Nu = " 
//         << Nbatch << "," << Nevec_acc << "," << R << "," << Nu << std::endl;
    
    for (int col=0; col<Nu; ++col) {
//      auto w_v = w[col].View();
      autoView( w_v,w[col], AcceleratorWrite);
      CUDA_COMPLEX *z = reinterpret_cast<CUDA_COMPLEX*>(&w_v[0]);
//      Glog << "col= "<<col <<" z= "<<z <<std::endl;
      for (size_t row=0; row<sites*12; ++row) {
        w_acc[col*sites*12+row] = z[row];
      }
    }
    cublasHandle_t handle;
    cublasStatus_t stat;
    
    stat = cublasCreate(&handle);

    Glog << "cuBLAS Zgemm"<< std::endl;
    
    for (int b=0; b<Nbatch; ++b) {
      for (int col=0; col<Nevec_acc; ++col) {
//        auto evec_v = evec[b*Nevec_acc+col].View();
        autoView( evec_v,evec[b*Nevec_acc+col], AcceleratorWrite);
        CUDA_COMPLEX *z = reinterpret_cast<CUDA_COMPLEX*>(&evec_v[0]);
//        Glog << "col= "<<col <<" z= "<<z <<std::endl;
        for (size_t row=0; row<sites*12; ++row) {
          evec_acc[col*sites*12+row] = z[row];
        }
      }
      CUDA_COMPLEX alpha = MAKE_CUDA_COMPLEX(1.0,0.0);
      CUDA_COMPLEX beta = MAKE_CUDA_COMPLEX(0.0,0.0);
      stat = CUDA_GEMM(handle, CUBLAS_OP_C, CUBLAS_OP_N, Nevec_acc, Nu, 12*sites,
                         &alpha, 
                         evec_acc, 12*sites, w_acc, 12*sites,  
                         &beta, 
                         c_acc, Nevec_acc);
      //Glog << stat << std::endl;
      
      grid->GlobalSumVector((CUDA_FLOAT*)c_acc,2*Nu*Nevec_acc);
      alpha = MAKE_CUDA_COMPLEX(-1.0,0.0);
      beta = MAKE_CUDA_COMPLEX(1.0,0.0);
      stat = CUDA_GEMM(handle, CUBLAS_OP_N, CUBLAS_OP_N, 12*sites, Nu, Nevec_acc,
                         &alpha, 
                         evec_acc, 12*sites, c_acc, Nevec_acc,  
                         &beta, 
                         w_acc, 12*sites);
      //Glog << stat << std::endl;
    }
    for (int col=0; col<Nu; ++col) {
//      auto w_v = w[col].View();
      autoView( w_v,w[col], AcceleratorWrite);
      CUDA_COMPLEX *z = reinterpret_cast<CUDA_COMPLEX*>(&w_v[0]);
      for (size_t row=0; row<sites*12; ++row) {
        z[row] = w_acc[col*sites*12+row];
      }
    }
    for (int i=0; i<Nu; ++i) {
      assert(normalize(w[i],do_print)!=0);
    }
    
    Glog << "cuBLAS Zgemm done"<< std::endl;
    
    cublasDestroy(handle);
    
    Glog << "cuBLAS orthogonalize done" << std::endl;
#else
    Glog << "BLAS wrapper is not implemented" << std::endl;
    exit(1);
#endif
  }


  
  void orthogonalize_blockhead(Field& w, std::vector<Field>& evec, int k, int Nu)
  {
    typedef typename Field::scalar_type MyComplex;
    MyComplex ip;
    
    for(int j=0; j<k; ++j){
      ip = innerProduct(evec[j*Nu],w); 
      w = w - ip * evec[j*Nu];
    }
    normalize(w);
  }
  
  void calc(std::vector<RealD>& eval,  
            std::vector<Field>& evec, 
            const std::vector<Field>& src, int& Nconv, LanczosType Impl)
  {
#ifdef GRID_CUDA
    GridBase *grid = src[0].Grid();
    grid->show_decomposition();

//    printf("GRID_CUDA\n");
    
    // set eigenvector buffers for the cuBLAS calls
    //const uint64_t nsimd = grid->Nsimd();
    const uint64_t sites = grid->lSites();
    
    cudaStat = cudaMallocManaged((void **)&w_acc, Nu*sites*12*sizeof(CUDA_COMPLEX));
//    Glog << "w_acc= "<<w_acc << " "<< cudaStat << std::endl;
cudaStat = cudaMallocManaged((void **)&evec_acc, Nevec_acc*sites*12*sizeof(CUDA_COMPLEX));
//    Glog << "evec_acc= "<<evec_acc << " "<< cudaStat << std::endl;
    cudaStat = cudaMallocManaged((void **)&c_acc, Nu*Nevec_acc*sizeof(CUDA_COMPLEX));
//    Glog << "c_acc= "<<c_acc << " "<< cudaStat << std::endl;
#endif
    switch (Impl) {
      case LanczosType::irbl: 
        calc_irbl(eval,evec,src,Nconv);
        break;
      
      case LanczosType::rbl: 
        calc_rbl(eval,evec,src,Nconv);
        break;
    }
#ifdef GRID_CUDA
    // free eigenvector buffers for the cuBLAS calls
    cudaFree(w_acc);
    cudaFree(evec_acc);
    cudaFree(c_acc);
#endif
  }

  void calc_irbl(std::vector<RealD>& eval,  
                 std::vector<Field>& evec, 
                 const std::vector<Field>& src, int& Nconv)
  {
    std::string fname = std::string(cname+"::calc_irbl()"); 
    GridBase *grid = evec[0].Grid();
    assert(grid == src[0].Grid());
    assert( Nu = src.size() );
    
    Glog << std::string(74,'*') << std::endl;
    Glog << fname + " starting iteration 0 /  "<< MaxIter<< std::endl;
    Glog << std::string(74,'*') << std::endl;
    Glog <<" -- seek   Nk    = "<< Nk    <<" vectors"<< std::endl;
    Glog <<" -- accept Nstop = "<< Nstop <<" vectors"<< std::endl;
    Glog <<" -- total  Nm    = "<< Nm    <<" vectors"<< std::endl;
    Glog <<" -- size of eval = "<< eval.size() << std::endl;
    Glog <<" -- size of evec = "<< evec.size() << std::endl;
    if ( diagonalisation == IRBLdiagonaliseWithEigen ) { 
      Glog << "Diagonalisation is Eigen "<< std::endl;
#ifdef USE_LAPACK
    } else if ( diagonalisation == IRBLdiagonaliseWithLAPACK ) { 
      Glog << "Diagonalisation is LAPACK "<< std::endl;
#endif
    } else {
      abort();
    }
    Glog << std::string(74,'*') << std::endl;
    
    assert(Nm == evec.size() && Nm == eval.size());

    std::vector<std::vector<ComplexD>> lmd(Nu,std::vector<ComplexD>(Nm,0.0));  
    std::vector<std::vector<ComplexD>> lme(Nu,std::vector<ComplexD>(Nm,0.0));  
    std::vector<std::vector<ComplexD>> lmd2(Nu,std::vector<ComplexD>(Nm,0.0));  
    std::vector<std::vector<ComplexD>> lme2(Nu,std::vector<ComplexD>(Nm,0.0));  
    std::vector<RealD> eval2(Nm);
    std::vector<RealD> resid(Nk);

    Eigen::MatrixXcd    Qt = Eigen::MatrixXcd::Zero(Nm,Nm);
    Eigen::MatrixXcd    Q = Eigen::MatrixXcd::Zero(Nm,Nm);

    std::vector<int>   Iconv(Nm);
    std::vector<Field>  B(Nm,grid); // waste of space replicating
    
    std::vector<Field> f(Nu,grid);
    std::vector<Field> f_copy(Nu,grid);
    Field v(grid);
    
    Nconv = 0;
    
    RealD beta_k;
  
    // set initial vector
    for (int i=0; i<Nu; ++i) {
      Glog << "norm2(src[" << i << "])= "<< norm2(src[i]) << std::endl;
      evec[i] = src[i];
      orthogonalize(evec[i],evec,i);
      Glog << "norm2(evec[" << i << "])= "<< norm2(evec[i]) << std::endl;
    }
    
    // initial Nblock_k steps
    for(int b=0; b<Nblock_k; ++b) blockwiseStep(lmd,lme,evec,f,f_copy,b);

    // restarting loop begins
    int iter;
    for(iter = 0; iter<MaxIter; ++iter){
      
      Glog <<"#Restart iteration = "<< iter << std::endl;
      // additional (Nblock_m - Nblock_k) steps
      for(int b=Nblock_k; b<Nblock_m; ++b) blockwiseStep(lmd,lme,evec,f,f_copy,b);
      
      // getting eigenvalues
      for(int u=0; u<Nu; ++u){
        for(int k=0; k<Nm; ++k){
          lmd2[u][k] = lmd[u][k];
          lme2[u][k] = lme[u][k];
        }
      }
      Qt = Eigen::MatrixXcd::Identity(Nm,Nm);
      diagonalize(eval2,lmd2,lme2,Nu,Nm,Nm,Qt,grid);
      _sort.push(eval2,Nm);
      Glog << "#Ritz value before shift: "<< std::endl;
      for(int i=0; i<Nm; ++i){
        std::cout.precision(13);
        std::cout << "[" << std::setw(4)<< std::setiosflags(std::ios_base::right) <<i<<"] ";
        std::cout << "Rval = "<<std::setw(20)<< std::setiosflags(std::ios_base::left)<< eval2[i] << std::endl;
      }
      
      //----------------------------------------------------------------------
      if ( Nm>Nk ) {
        Glog <<" #Apply shifted QR transformations "<<std::endl;
        //int k2 = Nk+Nu;
        int k2 = Nk;
      
        Eigen::MatrixXcd BTDM = Eigen::MatrixXcd::Identity(Nm,Nm);
        Q = Eigen::MatrixXcd::Identity(Nm,Nm);
        
        unpackHermitBlockTriDiagMatToEigen(lmd,lme,Nu,Nblock_m,Nm,Nm,BTDM);

        for(int ip=Nk; ip<Nm; ++ip){ 
          shiftedQRDecompEigen(BTDM,Nu,Nm,eval2[ip],Q);
        }
        
        packHermitBlockTriDiagMatfromEigen(lmd,lme,Nu,Nblock_m,Nm,Nm,BTDM);

        for(int i=0; i<k2; ++i) B[i] = 0.0;
        for(int j=0; j<k2; ++j){
          for(int k=0; k<Nm; ++k){
            B[j].Checkerboard() = evec[k].Checkerboard();
            B[j] += evec[k]*Q(k,j);
          }
        }
        for(int i=0; i<k2; ++i) evec[i] = B[i];

        // reconstruct initial vector for additional pole space
        blockwiseStep(lmd,lme,evec,f,f_copy,Nblock_k-1);

        // getting eigenvalues
        for(int u=0; u<Nu; ++u){
          for(int k=0; k<Nm; ++k){
            lmd2[u][k] = lmd[u][k];
            lme2[u][k] = lme[u][k];
          }
        }
        Qt = Eigen::MatrixXcd::Identity(Nm,Nm);
        diagonalize(eval2,lmd2,lme2,Nu,Nk,Nm,Qt,grid);
        _sort.push(eval2,Nk);
        Glog << "#Ritz value after shift: "<< std::endl;
        for(int i=0; i<Nk; ++i){
          std::cout.precision(13);
          std::cout << "[" << std::setw(4)<< std::setiosflags(std::ios_base::right) <<i<<"] ";
          std::cout << "Rval = "<<std::setw(20)<< std::setiosflags(std::ios_base::left)<< eval2[i] << std::endl;
        }
      }
      //----------------------------------------------------------------------

      // Convergence test
      Glog <<" #Convergence test: "<<std::endl;
      for(int k = 0; k<Nk; ++k) B[k]=0.0;
      for(int j = 0; j<Nk; ++j){
	for(int k = 0; k<Nk; ++k){
	  B[j].Checkerboard() = evec[k].Checkerboard();
	  B[j] += evec[k]*Qt(k,j);
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
        resid[i] = vv;
	
	std::cout.precision(13);
        std::cout << "[" << std::setw(4)<< std::setiosflags(std::ios_base::right) <<i<<"] ";
	std::cout << "eval = "<<std::setw(20)<< std::setiosflags(std::ios_base::left)<< eval2[i];
	std::cout << "   resid^2 = "<< std::setw(20)<< std::setiosflags(std::ios_base::right)<< vv<< std::endl;
	
	// change the criteria as evals are supposed to be sorted, all evals smaller(larger) than Nstop should have converged
	//if( (vv<eresid*eresid) && (i == Nconv) ){
	if (vv<eresid*eresid) {
	  Iconv[Nconv] = i;
	  ++Nconv;
	}
	
      }  // i-loop end
      
      Glog <<" #modes converged: "<<Nconv<<std::endl;
      for(int i=0; i<Nconv; ++i){
	std::cout.precision(13);
        std::cout << "[" << std::setw(4)<< std::setiosflags(std::ios_base::right) <<Iconv[i]<<"] ";
	std::cout << "eval_conv = "<<std::setw(20)<< std::setiosflags(std::ios_base::left)<< eval2[Iconv[i]];
	std::cout << "   resid^2 = "<< std::setw(20)<< std::setiosflags(std::ios_base::right)<< resid[Iconv[i]]<< std::endl;
      } 

      if ( Nconv>=Nstop ) break;

    } // end of iter loop
    
    Glog << std::string(74,'*') << std::endl;
    if ( Nconv<Nstop ) {
      Glog << fname + " NOT converged ; Summary :\n";
    } else {
      Glog << fname + " CONVERGED ; Summary :\n";
      // Sort convered eigenpairs.
      eval.resize(Nconv);
      evec.resize(Nconv,grid);
      for(int i=0; i<Nconv; ++i){
        eval[i] = eval2[Iconv[i]];
        evec[i] = B[Iconv[i]];
      }
      _sort.push(eval,evec,Nconv);
    }
    Glog << std::string(74,'*') << std::endl;
    Glog << " -- Iterations  = "<< iter   << "\n";
    //Glog << " -- beta(k)     = "<< beta_k << "\n";
    Glog << " -- Nconv       = "<< Nconv  << "\n";
    Glog << std::string(74,'*') << std::endl;
  
  }
  
  
  void calc_rbl(std::vector<RealD>& eval,  
                 std::vector<Field>& evec, 
                 const std::vector<Field>& src, int& Nconv)
  {
    std::string fname = std::string(cname+"::calc_rbl()"); 
    GridBase *grid = evec[0].Grid();
    assert(grid == src[0].Grid());
    assert( Nu = src.size() );

    int Np = (Nm-Nk);
    if (Np > 0 && MaxIter > 1) Np /= MaxIter;
    int Nblock_p = Np/Nu;
    for(int i=0;i< evec.size();i++) evec[0].Advise()=AdviseInfrequentUse;
    
    Glog << std::string(74,'*') << std::endl;
    Glog << fname + " starting iteration 0 /  "<< MaxIter<< std::endl;
    Glog << std::string(74,'*') << std::endl;
    Glog <<" -- seek (min) Nk    = "<< Nk    <<" vectors"<< std::endl;
    Glog <<" -- seek (inc) Np    = "<< Np <<" vectors"<< std::endl;
    Glog <<" -- seek (max) Nm    = "<< Nm    <<" vectors"<< std::endl;
    Glog <<" -- accept Nstop     = "<< Nstop <<" vectors"<< std::endl;
    Glog <<" -- size of eval     = "<< eval.size() << std::endl;
    Glog <<" -- size of evec     = "<< evec.size() << std::endl;
    if ( diagonalisation == IRBLdiagonaliseWithEigen ) { 
      Glog << "Diagonalisation is Eigen "<< std::endl;
#ifdef USE_LAPACK
    } else if ( diagonalisation == IRBLdiagonaliseWithLAPACK ) { 
      Glog << "Diagonalisation is LAPACK "<< std::endl;
#endif
    } else {
      abort();
    }
    Glog << std::string(74,'*') << std::endl;
    
    assert(Nm == evec.size() && Nm == eval.size());
	
    std::vector<std::vector<ComplexD>> lmd(Nu,std::vector<ComplexD>(Nm,0.0));  
    std::vector<std::vector<ComplexD>> lme(Nu,std::vector<ComplexD>(Nm,0.0));  
    std::vector<std::vector<ComplexD>> lmd2(Nu,std::vector<ComplexD>(Nm,0.0));  
    std::vector<std::vector<ComplexD>> lme2(Nu,std::vector<ComplexD>(Nm,0.0));  
    std::vector<RealD> eval2(Nk);
    std::vector<RealD> resid(Nm);

    Eigen::MatrixXcd    Qt = Eigen::MatrixXcd::Zero(Nm,Nm);
    Eigen::MatrixXcd    Q = Eigen::MatrixXcd::Zero(Nm,Nm);

    std::vector<int>   Iconv(Nm);
//    int Ntest=Nu;
//    std::vector<Field>  B(Nm,grid); // waste of space replicating
    std::vector<Field>  B(1,grid); // waste of space replicating
    
    std::vector<Field> f(Nu,grid);
    std::vector<Field> f_copy(Nu,grid);
    Field v(grid);
    
    Nconv = 0;
    
//    RealD beta_k;
  
    // set initial vector
    for (int i=0; i<Nu; ++i) {
      Glog << "norm2(src[" << i << "])= "<< norm2(src[i]) << std::endl;
      evec[i] = src[i];
      orthogonalize(evec[i],evec,i);
      Glog << "norm2(evec[" << i << "])= "<< norm2(evec[i]) << std::endl;
    }
//    exit(-43);
    
    // initial Nblock_k steps
    for(int b=0; b<Nblock_k; ++b) blockwiseStep(lmd,lme,evec,f,f_copy,b);

    // restarting loop begins
    int iter;
    int Nblock_l, Nblock_r;
    int Nl, Nr;
    int Nconv_guess = 0;

    for(iter = 0; iter<MaxIter; ++iter){
         
      Glog <<"#Restart iteration = "<< iter << std::endl;
      
      Nblock_l = Nblock_k + iter*Nblock_p;
      Nblock_r = Nblock_l + Nblock_p;
      Nl = Nblock_l*Nu;
      Nr = Nblock_r*Nu;
      eval2.resize(Nr);

      // additional Nblock_p steps
      for(int b=Nblock_l; b<Nblock_r; ++b) blockwiseStep(lmd,lme,evec,f,f_copy,b);
      
      // getting eigenvalues
      for(int u=0; u<Nu; ++u){
        for(int k=0; k<Nr; ++k){
          lmd2[u][k] = lmd[u][k];
          lme2[u][k] = lme[u][k];
        }
      }
      Qt = Eigen::MatrixXcd::Identity(Nr,Nr);
      diagonalize(eval2,lmd2,lme2,Nu,Nr,Nr,Qt,grid);
      _sort.push(eval2,Nr);
      Glog << "#Ritz value: "<< std::endl;
      for(int i=0; i<Nr; ++i){
        std::cout.precision(13);
        std::cout << "[" << std::setw(4)<< std::setiosflags(std::ios_base::right) <<i<<"] ";
        std::cout << "Rval = "<<std::setw(20)<< std::setiosflags(std::ios_base::left)<< eval2[i] << std::endl;
      }
      
      // Convergence test
      Glog <<" #Convergence test: "<<std::endl;
      Nconv = 0;
      for(int j = 0; j<Nr; j+=Nconv_test_interval){
	B[0]=0.0;
        if ( j/Nconv_test_interval == Nconv ) {
          Glog <<" #rotation for next check point evec" 
               << std::setw(4)<< std::setiosflags(std::ios_base::right) 
               << "["<< j <<"]" <<std::endl;
          for(int k = 0; k<Nr; ++k){
            B[0].Checkerboard() = evec[k].Checkerboard();
            B[0] += evec[k]*Qt(k,j);
          }
          
          _Linop.HermOp(B[0],v);
          RealD vnum = real(innerProduct(B[0],v)); // HermOp.
          RealD vden = norm2(B[0]);
          eval2[j] = vnum/vden;
          v -= eval2[j]*B[0];
          RealD vv = norm2(v);
          resid[j] = vv;
          
          std::cout.precision(13);
          std::cout << "[" << std::setw(4)<< std::setiosflags(std::ios_base::right) <<j<<"] ";
          std::cout << "eval = "<<std::setw(20)<< std::setiosflags(std::ios_base::left)<< eval2[j];
          std::cout << "   resid^2 = "<< std::setw(20)<< std::setiosflags(std::ios_base::right)<< vv<< std::endl;
          
          // change the criteria as evals are supposed to be sorted, all evals smaller(larger) than Nstop should have converged
          //if( (vv<eresid*eresid) && (i == Nconv) ){
          if (vv<eresid*eresid) {
            Iconv[Nconv] = j;
            ++Nconv;
          }
        } else {
          break;
        }
      }  // j-loop end
      
      Glog <<" #modes converged: "<<Nconv<<std::endl;
      for(int i=0; i<Nconv; ++i){
	std::cout.precision(13);
        std::cout << "[" << std::setw(4)<< std::setiosflags(std::ios_base::right) <<Iconv[i]<<"] ";
	std::cout << "eval_conv = "<<std::setw(20)<< std::setiosflags(std::ios_base::left)<< eval2[Iconv[i]];
	std::cout << "   resid^2 = "<< std::setw(20)<< std::setiosflags(std::ios_base::right)<< resid[Iconv[i]]<< std::endl;
      } 

      (Nconv > 0 ) ? Nconv_guess = 1 + (Nconv-1)*Nconv_test_interval : Nconv_guess = 0;
      if ( Nconv_guess >= Nstop ) break;

    } // end of iter loop
    
    Glog << std::string(74,'*') << std::endl;
    if ( Nconv_guess < Nstop ) {
      Glog << fname + " NOT converged ; Summary :\n";
    } else {
      Glog << fname + " CONVERGED ; Summary :\n";
      // Sort convered eigenpairs.
      std::vector<Field>  Btmp(Nstop,grid); // waste of space replicating

      for(int i=0; i<Nstop; ++i){
	  Btmp[i]=0.;
          for(int k = 0; k<Nr; ++k){
             Btmp[i].Checkerboard() = evec[k].Checkerboard();
             Btmp[i] += evec[k]*Qt(k,i);
          }
          _Linop.HermOp(Btmp[i],v);
          RealD vnum = real(innerProduct(Btmp[i],v)); // HermOp.
          RealD vden = norm2(Btmp[i]);
//          eval2[j] = vnum/vden;
          v -= vnum/vden*Btmp[i];
          RealD vv = norm2(v);
//          resid[j] = vv;
          
          std::cout.precision(13);
          std::cout << "[" << std::setw(4)<< std::setiosflags(std::ios_base::right) <<i<<"] ";
          std::cout << "eval = "<<std::setw(20)<< std::setiosflags(std::ios_base::left)<< vnum/vden;
          std::cout << "   resid^2 = "<< std::setw(20)<< std::setiosflags(std::ios_base::right)<< vv<< std::endl;
        eval[i] = vnum/vden;
      }
      for(int i=0; i<Nstop; ++i) evec[i] = Btmp[i];
      eval.resize(Nstop);
      evec.resize(Nstop,grid);
      _sort.push(eval,evec,Nstop);
    }
    Glog << std::string(74,'*') << std::endl;
    Glog << " -- Iterations    = "<< iter   << "\n";
    //Glog << " -- beta(k)       = "<< beta_k << "\n";
    Glog << " -- Nconv         = "<< Nconv  << "\n";
    Glog << " -- Nconv (guess) = "<< Nconv_guess  << "\n";
    Glog << std::string(74,'*') << std::endl;
  
  }

private:
  void blockwiseStep(std::vector<std::vector<ComplexD>>& lmd,
	             std::vector<std::vector<ComplexD>>& lme, 
	             std::vector<Field>& evec,
	             std::vector<Field>& w, 
	             std::vector<Field>& w_copy, 
                     int b)
  {
    const RealD tiny = 1.0e-20;
    
    int Nu = w.size();
    int Nm = evec.size();
    assert( b < Nm/Nu );
//    GridCartesian *grid = evec[0]._grid;
    
    // converts block index to full indicies for an interval [L,R)
    int L = Nu*b;
    int R = Nu*(b+1);

    Real beta;

    Glog << "Using split grid"<< std::endl;
//   LatticeGaugeField s_Umu(SGrid);
   assert((Nu%mrhs)==0);
   std::vector<Field>   in(mrhs,f_grid);
     
    Field s_in(sf_grid);
    Field s_out(sf_grid);
   // unnecessary copy. Can or should it be avoided?
    int k_start = 0;
while ( k_start < Nu) {
   Glog << "k_start= "<<k_start<< std::endl;
   for (int u=0; u<mrhs; ++u) in[u] = evec[L+k_start+u];
Glog << "Split "<< std::endl;
   Grid_split(in, s_in);
Glog << "Split done "<< std::endl;
      _poly(_SLinop,s_in,s_out);
Glog << "Unsplit "<< std::endl;
   Grid_unsplit(in,s_out);
Glog << "Unsplit done "<< std::endl;
   for (int u=0; u<mrhs; ++u) w[k_start+u] = in[u];
   k_start +=mrhs;
}
    Glog << "Using split grid done "<< std::endl;
    
// test split in the first iteration
if(split_test){
    Glog << "Split grid testing "<< std::endl;
    // 3. wk:=Avkβkv_{k1}
    for (int k=L, u=0; k<R; ++k, ++u) {
      _poly(_Linop,evec[k],w_copy[u]);      
    }
   for (int u=0; u<Nu; ++u) {
	 w_copy[u] -= w[u];
    Glog << "diff(split - non_split) "<<u<<" " << norm2(w_copy[u]) << std::endl;
    Glog << "Split grid testing done"<< std::endl;
   }
   split_test=0;
}
    Glog << "Poly done"<< std::endl;
    Glog << "LinAlg "<< std::endl;
//    exit(-42);
    
    if (b>0) {
      for (int u=0; u<Nu; ++u) {
        //for (int k=L-Nu; k<L; ++k) {
        for (int k=L-Nu+u; k<L; ++k) {
          w[u] = w[u] - evec[k] * conjugate(lme[u][k]);
        }
      }
    }
    
    // 4. αk:=(vk,wk)
    //for (int u=0; u<Nu; ++u) {
    //  for (int k=L; k<R; ++k) {
    //    lmd[u][k] = innerProduct(evec[k],w[u]);  // lmd = transpose of alpha
    //  }
    //  lmd[u][L+u] = real(lmd[u][L+u]);  // force diagonal to be real
    //}
    for (int u=0; u<Nu; ++u) {
      for (int k=L+u; k<R; ++k) {
        lmd[u][k] = innerProduct(evec[k],w[u]);  // lmd = transpose of alpha
//        Glog <<"lmd "<<u<<" "<<k<<lmd[u][k] -conjugate(innerProduct(evec[u+L],w[k-L]))<<std::endl;
        lmd[k-L][u+L] = conjugate(lmd[u][k]);     // force hermicity
      }
      lmd[u][L+u] = real(lmd[u][L+u]);  // force diagonal to be real
    }
    
    // 5. wk:=wk−αkvk
    for (int u=0; u<Nu; ++u) {
      for (int k=L; k<R; ++k) {
        w[u] = w[u] - evec[k]*lmd[u][k];
      }
      w_copy[u] = w[u];
    }
    Glog << "LinAlg done"<< std::endl;
    
    // In block version, the steps 6 and 7 in Lanczos construction is
    // replaced by the QR decomposition of new basis block.
    // It results block version beta and orthonormal block basis. 
    // Here, QR decomposition is done by using Gram-Schmidt.
    for (int u=0; u<Nu; ++u) {
      for (int k=L; k<R; ++k) {
        lme[u][k] = 0.0;
      }
    }

    // re-orthogonalization for numerical stability
#if 1
    Glog << "Gram Schmidt"<< std::endl;
    orthogonalize(w,Nu,evec,R);
#else
    Glog << "Gram Schmidt using cublas"<< std::endl;
    orthogonalize_blas(w,evec,R);
#endif
    // QR part
    for (int u=1; u<Nu; ++u) {
      orthogonalize(w[u],w,u);
    }
    Glog << "Gram Schmidt done "<< std::endl;
    
    Glog << "LinAlg "<< std::endl;
    for (int u=0; u<Nu; ++u) {
      //for (int v=0; v<Nu; ++v) {
      for (int v=u; v<Nu; ++v) {
        lme[u][L+v] = innerProduct(w[u],w_copy[v]);
      }
      lme[u][L+u] = real(lme[u][L+u]);  // force diagonal to be real
    }
    //lme[0][L] = beta;
    
    for (int u=0; u<Nu; ++u) {
//      Glog << "norm2(w[" << u << "])= "<< norm2(w[u]) << std::endl;
      assert (!isnan(norm2(w[u])));
      for (int k=L+u; k<R; ++k) {
        Glog <<" In block "<< b << "," <<" beta[" << u << "," << k-L << "] = " << lme[u][k] << std::endl;
      }
    }
    Glog << "LinAlg done "<< std::endl;

    if (b < Nm/Nu-1) {
      for (int u=0; u<Nu; ++u) {
        evec[R+u] = w[u];
      }
    }

  }
  
    
  void diagonalize_Eigen(std::vector<RealD>& eval, 
                         std::vector<std::vector<ComplexD>>& lmd,
                         std::vector<std::vector<ComplexD>>& lme, 
			 int Nu, int Nk, int Nm,
			 Eigen::MatrixXcd & Qt, // Nm x Nm
			 GridBase *grid)
  {
    assert( Nk%Nu == 0 && Nm%Nu == 0 );
    assert( Nk <= Nm );
    Eigen::MatrixXcd BlockTriDiag = Eigen::MatrixXcd::Zero(Nk,Nk);
    
    for ( int u=0; u<Nu; ++u ) {
      for (int k=0; k<Nk; ++k ) {
        BlockTriDiag(k,u+(k/Nu)*Nu) = lmd[u][k];
      }
    }
    
    for ( int u=0; u<Nu; ++u ) {
      for (int k=Nu; k<Nk; ++k ) {
        BlockTriDiag(k-Nu,u+(k/Nu)*Nu) = conjugate(lme[u][k-Nu]);
        BlockTriDiag(u+(k/Nu)*Nu,k-Nu) = lme[u][k-Nu];
      }
    }
    //std::cout << BlockTriDiag << std::endl;
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(BlockTriDiag);

    for (int i = 0; i < Nk; i++) {
      eval[Nk-1-i] = eigensolver.eigenvalues()(i);
    }
    for (int i = 0; i < Nk; i++) {
      for (int j = 0; j < Nk; j++) {
	Qt(j,Nk-1-i) = eigensolver.eigenvectors()(j,i);
	//Qt(Nk-1-i,j) = eigensolver.eigenvectors()(i,j);
	//Qt(i,j) = eigensolver.eigenvectors()(i,j);
      }
    }
  }

#ifdef USE_LAPACK
  void diagonalize_lapack(std::vector<RealD>& eval, 
                         std::vector<std::vector<ComplexD>>& lmd,
                         std::vector<std::vector<ComplexD>>& lme, 
			 int Nu, int Nk, int Nm,
			 Eigen::MatrixXcd & Qt, // Nm x Nm
			 GridBase *grid)
  {
    Glog << "diagonalize_lapack: Nu= "<<Nu<<" Nk= "<<Nk<<" Nm= "<<std::endl;
    assert( Nk%Nu == 0 && Nm%Nu == 0 );
    assert( Nk <= Nm );
    Eigen::MatrixXcd BlockTriDiag = Eigen::MatrixXcd::Zero(Nk,Nk);
    
    for ( int u=0; u<Nu; ++u ) {
      for (int k=0; k<Nk; ++k ) {
//        Glog << "lmd "<<u<<" "<<k<<" "<<lmd[u][k] -conjugate(lmd[u][k])<<std::endl;
        BlockTriDiag(k,u+(k/Nu)*Nu) = lmd[u][k];
      }
    }
    
    for ( int u=0; u<Nu; ++u ) {
      for (int k=Nu; k<Nk; ++k ) {
//        Glog << "lme "<<u<<" "<<k<<" "<<lme[u][k] -conjugate(lme[u][k])<<std::endl;
        BlockTriDiag(k-Nu,u+(k/Nu)*Nu) = conjugate(lme[u][k-Nu]);
        BlockTriDiag(u+(k/Nu)*Nu,k-Nu) = lme[u][k-Nu];
      }
    }
    //std::cout << BlockTriDiag << std::endl;
//#ifdef USE_LAPACK
  const int size = Nm;
  MKL_INT NN = Nk;
//  double evals_tmp[NN];
//  double evec_tmp[NN][NN];
  double *evals_tmp = (double *) malloc(NN*sizeof(double));
  MKL_Complex16 *evec_tmp = (MKL_Complex16 *) malloc(NN*NN*sizeof(MKL_Complex16));
  MKL_Complex16 *DD = (MKL_Complex16 *) malloc(NN*NN*sizeof(MKL_Complex16));
  for (int i = 0; i< NN; i++) {
    for (int j = 0; j <NN ; j++) {
        evec_tmp[i*NN+j].real=0.;
        evec_tmp[i*NN+j].imag=0.;
        DD[i*NN+j].real=BlockTriDiag(i,j).real();
        DD[i*NN+j].imag=BlockTriDiag(i,j).imag();
    }
  }
  MKL_INT evals_found;
  MKL_INT lwork = (3*NN);
  MKL_INT lrwork = (24*NN);
  MKL_INT liwork =  NN*10 ;
  MKL_INT iwork[liwork];
  double rwork[lrwork];
//  double work[lwork];
  MKL_Complex16 *work = (MKL_Complex16 *) malloc(lwork*sizeof(MKL_Complex16));
  MKL_INT isuppz[2*NN];
  char jobz = 'V'; // calculate evals & evecs
  char range = 'I'; // calculate all evals
  //    char range = 'A'; // calculate all evals
  char uplo = 'U'; // refer to upper half of original matrix
  char compz = 'I'; // Compute eigenvectors of tridiagonal matrix
  int ifail[NN];
  MKL_INT info;
  int total = grid->_Nprocessors;
  int node  = grid->_processor;
  int interval = (NN/total)+1;
  double vl = 0.0, vu = 0.0;
  MKL_INT il = interval*node+1 , iu = interval*(node+1);
  if (iu > NN)  iu=NN;
  Glog << "node "<<node<<"il "<<il<<"iu "<<iu<<std::endl;
  double tol = 0.0;
  if (1) {
    memset(evals_tmp,0,sizeof(double)*NN);
    if ( il <= NN){
      zheevr(&jobz, &range, &uplo, &NN,
		    DD,  &NN,
		    &vl, &vu, &il, &iu, // these four are ignored if second parameteris 'A'
		    &tol, // tolerance
		    &evals_found, evals_tmp, (MKL_Complex16*)evec_tmp, &NN,
		    isuppz,
		    work, &lwork, 
		    rwork, &lrwork, 
		    iwork, &liwork,
		    &info);
//			(double*)EE,
      for (int i = iu-1; i>= il-1; i--){
	evals_tmp[i] = evals_tmp[i - (il-1)];
	if (il>1) evals_tmp[i-(il-1)]=0.;
	for (int j = 0; j< NN; j++){
	  evec_tmp[i*NN+j] = evec_tmp[(i - (il-1))*NN+j];
	  if (il>1) {
		evec_tmp[(i-(il-1))*NN+j].imag=0.;
		evec_tmp[(i-(il-1))*NN+j].real=0.;
          }
	}
      }
    }
    {
      grid->GlobalSumVector(evals_tmp,NN);
      grid->GlobalSumVector((double*)evec_tmp,2*NN*NN);
    }
  } 
    for (int i = 0; i < Nk; i++) 
      eval[Nk-1-i] = evals_tmp[i];
    for (int i = 0; i < Nk; i++) {
      for (int j = 0; j < Nk; j++) {
//	Qt(j,Nk-1-i) = eigensolver.eigenvectors()(j,i);
	Qt(j,Nk-1-i)=std::complex<double>  
	( evec_tmp[i*Nk+j].real,
	 evec_tmp[i*Nk+j].imag);
//	( evec_tmp[(Nk-1-j)*Nk+Nk-1-i].real,
//	evec_tmp[(Nk-1-j)*Nk+Nk-1-i].imag);
	
      }
    }
    
if (1){
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(BlockTriDiag);

    for (int i = 0; i < Nk; i++) {
      Glog << "eval = "<<i<<" " <<eval[Nk-1-i] <<" "<< eigensolver.eigenvalues()(i) <<std::endl;
//      eval[Nk-1-i] = eigensolver.eigenvalues()(i);
    }
    for (int i = 0; i < Nk; i++) {
      for (int j = 0; j < Nk; j++) {
//	Qt(j,Nk-1-i) = eigensolver.eigenvectors()(j,i);
//        Glog<<"Qt "<<j<<" "<<Nk-1-i<<" = " <<Qt(j,Nk-1-i) <<" "<<eigensolver.eigenvectors()(j,i) <<std::endl;
        MKL_Complex16 tmp = evec_tmp[i*Nk+j];
//        Glog<<"Qt "<<j<<" "<<Nk-1-i<<" = " <<evec_tmp[(Nk-1-j)*Nk+Nk-1-i].real<<" "<<
//evec_tmp[(Nk-1-j)*Nk+Nk-1-i].imag <<" "<<eigensolver.eigenvectors()(j,i) <<std::endl;
	if ( (i<5)&& (j<5))
       Glog<<"Qt "<<j<<" "<<Nk-1-i<<" = " << norm(Qt(j,Nk-1-i))<<" "<<
          norm(eigensolver.eigenvectors()(j,i)) <<std::endl;
      }
    }
}
//  exit(-43);

  free (evals_tmp);
  free (evec_tmp);
  free (DD);
  free (work);
  }
#endif


  void diagonalize(std::vector<RealD>& eval, 
                   std::vector<std::vector<ComplexD>>& lmd, 
                   std::vector<std::vector<ComplexD>>& lme, 
		   int Nu, int Nk, int Nm,   
		   Eigen::MatrixXcd & Qt,
		   GridBase *grid)
  {
    Qt = Eigen::MatrixXcd::Identity(Nm,Nm);
    if ( diagonalisation == IRBLdiagonaliseWithEigen ) { 
      diagonalize_Eigen(eval,lmd,lme,Nu,Nk,Nm,Qt,grid);
#ifdef USE_LAPACK
    } else if ( diagonalisation == IRBLdiagonaliseWithLAPACK ) { 
      diagonalize_lapack(eval,lmd,lme,Nu,Nk,Nm,Qt,grid);
#endif
    } else { 
      assert(0);
    }
  }
  

  void unpackHermitBlockTriDiagMatToEigen(
         std::vector<std::vector<ComplexD>>& lmd,  
         std::vector<std::vector<ComplexD>>& lme,
         int Nu, int Nb, int Nk, int Nm,
         Eigen::MatrixXcd& M)
  {
    //Glog << "unpackHermitBlockTriDiagMatToEigen() begin" << '\n'; 
    assert( Nk%Nu == 0 && Nm%Nu == 0 );
    assert( Nk <= Nm );
    M = Eigen::MatrixXcd::Zero(Nk,Nk);
    
    // rearrange 
    for ( int u=0; u<Nu; ++u ) {
      for (int k=0; k<Nk; ++k ) {
        M(k,u+(k/Nu)*Nu) = lmd[u][k];
      }
    }

    for ( int u=0; u<Nu; ++u ) {
      for (int k=Nu; k<Nk; ++k ) {
        M(k-Nu,u+(k/Nu)*Nu) = conjugate(lme[u][k-Nu]);
        M(u+(k/Nu)*Nu,k-Nu) = lme[u][k-Nu];
      }
    }
    //Glog << "unpackHermitBlockTriDiagMatToEigen() end" << endl; 
  }
 

  void packHermitBlockTriDiagMatfromEigen(
         std::vector<std::vector<ComplexD>>& lmd,
         std::vector<std::vector<ComplexD>>& lme,
         int Nu, int Nb, int Nk, int Nm,
         Eigen::MatrixXcd& M)
  {
    //Glog << "packHermitBlockTriDiagMatfromEigen() begin" << '\n'; 
    assert( Nk%Nu == 0 && Nm%Nu == 0 );
    assert( Nk <= Nm );
    
    // rearrange 
    for ( int u=0; u<Nu; ++u ) {
      for (int k=0; k<Nk; ++k ) {
        lmd[u][k] = M(k,u+(k/Nu)*Nu);
      }
    }

    for ( int u=0; u<Nu; ++u ) {
      for (int k=Nu; k<Nk; ++k ) {
        lme[u][k-Nu] = M(u+(k/Nu)*Nu,k-Nu);
      }
    }
    //Glog << "packHermitBlockTriDiagMatfromEigen() end" << endl; 
  }


  // assume the input matrix M is a band matrix
  void shiftedQRDecompEigen(Eigen::MatrixXcd& M, int Nu, int Nm,
		            RealD Dsh,
		            Eigen::MatrixXcd& Qprod)
  {
    //Glog << "shiftedQRDecompEigen() begin" << '\n'; 
    Eigen::MatrixXcd Q = Eigen::MatrixXcd::Zero(Nm,Nm);
    Eigen::MatrixXcd R = Eigen::MatrixXcd::Zero(Nm,Nm);
    Eigen::MatrixXcd Mtmp = Eigen::MatrixXcd::Zero(Nm,Nm);
    
    Mtmp = M;
    for (int i=0; i<Nm; ++i ) {
      Mtmp(i,i) = M(i,i) - Dsh;
    }
    
    Eigen::HouseholderQR<Eigen::MatrixXcd> QRD(Mtmp);
    Q = QRD.householderQ();
    R = QRD.matrixQR(); // upper triangular part is the R matrix.
                        // lower triangular part used to represent series
                        // of Q sequence.

    // equivalent operation of Qprod *= Q
    //M = Eigen::MatrixXcd::Zero(Nm,Nm);
    
    //for (int i=0; i<Nm; ++i) {
    //  for (int j=0; j<Nm-2*(Nu+1); ++j) {
    //    for (int k=0; k<2*(Nu+1)+j; ++k) {
    //      M(i,j) += Qprod(i,k)*Q(k,j);
    //    }
    //  }
    //}
    //for (int i=0; i<Nm; ++i) {
    //  for (int j=Nm-2*(Nu+1); j<Nm; ++j) {
    //    for (int k=0; k<Nm; ++k) {
    //      M(i,j) += Qprod(i,k)*Q(k,j);
    //    }
    //  }
    //}
    
    Mtmp = Eigen::MatrixXcd::Zero(Nm,Nm);

    for (int i=0; i<Nm; ++i) {
      for (int j=0; j<Nm-(Nu+1); ++j) {
        for (int k=0; k<Nu+1+j; ++k) {
          Mtmp(i,j) += Qprod(i,k)*Q(k,j);
        }
      }
    }
    for (int i=0; i<Nm; ++i) {
      for (int j=Nm-(Nu+1); j<Nm; ++j) {
        for (int k=0; k<Nm; ++k) {
          Mtmp(i,j) += Qprod(i,k)*Q(k,j);
        }
      }
    }
    
    //static int ntimes = 2;
    //for (int j=0; j<Nm-(ntimes*Nu); ++j) {
    //  for (int i=ntimes*Nu+j; i<Nm; ++i) {
    //    Mtmp(i,j) = 0.0;
    //  }
    //}
    //ntimes++;

    Qprod = Mtmp;
     
    // equivalent operation of M = Q.adjoint()*(M*Q)
    Mtmp = Eigen::MatrixXcd::Zero(Nm,Nm);
    
    for (int a=0, i=0, kmax=0; a<Nu+1; ++a) {
      for (int j=0; j<Nm-a; ++j) {
        i = j+a;
        kmax = (Nu+1)+j;
        if (kmax > Nm) kmax = Nm;
        for (int k=i; k<kmax; ++k) { 
          Mtmp(i,j) += R(i,k)*Q(k,j);
        }
        Mtmp(j,i) = conj(Mtmp(i,j));
      }
    }

    for (int i=0; i<Nm; ++i) {
      Mtmp(i,i) = real(Mtmp(i,i)) + Dsh;
    }
    
    M = Mtmp;

    //M = Q.adjoint()*(M*Q);
    //for (int i=0; i<Nm; ++i) {
    //  for (int j=0; j<Nm; ++j) {
    //    if (i==j) M(i,i) = real(M(i,i));
    //    if (j>i)  M(i,j) = conj(M(j,i));
    //    if (i-j > Nu || j-i > Nu) M(i,j) = 0.;
    //  }
    //}
    
    //Glog << "shiftedQRDecompEigen() end" << endl; 
  }

  void exampleQRDecompEigen(void)
  {
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3,3);
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(3,3);
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(3,3);
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(3,3);

    A(0,0) = 12.0;
    A(0,1) = -51.0;
    A(0,2) = 4.0;
    A(1,0) = 6.0;
    A(1,1) = 167.0;
    A(1,2) = -68.0;
    A(2,0) = -4.0;
    A(2,1) = 24.0;
    A(2,2) = -41.0;
    
    Glog << "matrix A before ColPivHouseholder" << std::endl;
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        Glog << "A[" << i << "," << j << "] = " << A(i,j) << '\n';
      }
    }
    Glog << std::endl;

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QRD(A);
    
    Glog << "matrix A after ColPivHouseholder" << std::endl;
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        Glog << "A[" << i << "," << j << "] = " << A(i,j) << '\n';
      }
    }
    Glog << std::endl;
    
    Glog << "HouseholderQ with sequence lenth = nonzeroPiviots" << std::endl;
    Q = QRD.householderQ().setLength(QRD.nonzeroPivots());
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        Glog << "Q[" << i << "," << j << "] = " << Q(i,j) << '\n';
      }
    }
    Glog << std::endl;
    
    Glog << "HouseholderQ with sequence lenth = 1" << std::endl;
    Q = QRD.householderQ().setLength(1);
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        Glog << "Q[" << i << "," << j << "] = " << Q(i,j) << '\n';
      }
    }
    Glog << std::endl;
    
    Glog << "HouseholderQ with sequence lenth = 2" << std::endl;
    Q = QRD.householderQ().setLength(2);
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        Glog << "Q[" << i << "," << j << "] = " << Q(i,j) << '\n';
      }
    }
    Glog << std::endl;
    
    Glog << "matrixR" << std::endl;
    R = QRD.matrixR();
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        Glog << "R[" << i << "," << j << "] = " << R(i,j) << '\n';
      }
    }
    Glog << std::endl;

    Glog << "rank = " << QRD.rank() << std::endl;
    Glog << "threshold = " << QRD.threshold() << std::endl;
    
    Glog << "matrixP" << std::endl;
    P = QRD.colsPermutation();
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        Glog << "P[" << i << "," << j << "] = " << P(i,j) << '\n';
      }
    }
    Glog << std::endl;


    Glog << "QR decomposition without column pivoting" << std::endl;
    
    A(0,0) = 12.0;
    A(0,1) = -51.0;
    A(0,2) = 4.0;
    A(1,0) = 6.0;
    A(1,1) = 167.0;
    A(1,2) = -68.0;
    A(2,0) = -4.0;
    A(2,1) = 24.0;
    A(2,2) = -41.0;
    
    Glog << "matrix A before Householder" << std::endl;
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        Glog << "A[" << i << "," << j << "] = " << A(i,j) << '\n';
      }
    }
    Glog << std::endl;
    
    Eigen::HouseholderQR<Eigen::MatrixXd> QRDplain(A);
    
    Glog << "HouseholderQ" << std::endl;
    Q = QRDplain.householderQ();
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        Glog << "Q[" << i << "," << j << "] = " << Q(i,j) << '\n';
      }
    }
    Glog << std::endl;
    
    Glog << "matrix A after Householder" << std::endl;
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        Glog << "A[" << i << "," << j << "] = " << A(i,j) << '\n';
      }
    }
    Glog << std::endl;
  }

 };
}
#undef Glog
#undef USE_LAPACK
#undef  CUDA_COMPLEX 
#undef  CUDA_FLOAT
#undef  MAKE_CUDA_COMPLEX
#undef  CUDA_GEMM 
#endif
