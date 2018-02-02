    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/ImplicitlyRestartedBlockLanczos.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Chulwoo Jung
Author: Yong-Chull Jang <ypj@quark.phy.bnl.gov> 
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
#ifndef GRID_IRBL_H
#define GRID_IRBL_H

#include <string.h> //memset

#define clog std::cout << GridLogMessage 

namespace Grid {

/////////////////////////////////////////////////////////////
// Implicitly restarted block lanczos
/////////////////////////////////////////////////////////////
template<class Field> 
class ImplicitlyRestartedBlockLanczos {

private:       
  
  std::string cname = std::string("ImplicitlyRestartedBlockLanczos");
  int MaxIter;   // Max iterations
  int Nstop;     // Number of evecs checked for convergence
  int Nu;        // Numbeer of vecs in the unit block
  int Nk;        // Number of converged sought
  int Nm;        // total number of vectors
  int Nblock_k;    // Nk/Nu
  int Nblock_m;    // Nm/Nu
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
 ImplicitlyRestartedBlockLanczos(LinearOperatorBase<Field> &Linop, // op
                                 OperatorFunction<Field> & poly,   // polynomial
                                 int _Nstop, // really sought vecs
                                 int _Nu,    // vecs in the unit block
                                 int _Nk,    // sought vecs
                                 int _Nm,    // total vecs
                                 RealD _eresid, // resid in lmd deficit 
                                 int _MaxIter,  // Max iterations
                                 IRLdiagonalisation _diagonalisation = IRLdiagonaliseWithEigen)
   : _Linop(Linop),    _poly(poly),
      Nstop(_Nstop), Nu(_Nu), Nk(_Nk), Nm(_Nm), 
      Nblock_m(_Nm/_Nu), Nblock_k(_Nk/_Nu),
      //eresid(_eresid),  MaxIter(10),
      eresid(_eresid),  MaxIter(_MaxIter),
      diagonalisation(_diagonalisation)
  { assert( (Nk%Nu==0) && (Nm%Nu==0) ); };

  ////////////////////////////////
  // Helpers
  ////////////////////////////////
  static RealD normalize(Field& v) 
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
      ip = innerProduct(evec[j],w); 
      w = w - ip * evec[j];
    }
    normalize(w);
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
            const std::vector<Field>& src, int& Nconv)
  {
    std::string fname = std::string(cname+"::calc()"); 
    GridBase *grid = evec[0]._grid;
    assert(grid == src[0]._grid);
    assert( Nu = src.size() );
    
    clog << std::string(74,'*') << std::endl;
    clog << fname + " starting iteration 0 /  "<< MaxIter<< std::endl;
    clog << std::string(74,'*') << std::endl;
    clog <<" -- seek   Nk    = "<< Nk    <<" vectors"<< std::endl;
    clog <<" -- accept Nstop = "<< Nstop <<" vectors"<< std::endl;
    clog <<" -- total  Nm    = "<< Nm    <<" vectors"<< std::endl;
    clog <<" -- size of eval = "<< eval.size() << std::endl;
    clog <<" -- size of evec = "<< evec.size() << std::endl;
    if ( diagonalisation == IRLdiagonaliseWithEigen ) { 
      clog << "Diagonalisation is Eigen "<< std::endl;
    } else {
      abort();
    }
    clog << std::string(74,'*') << std::endl;
    
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
      clog << "norm2(src[" << i << "])= "<< norm2(src[i]) << std::endl;
      evec[i] = src[i];
      orthogonalize(evec[i],evec,i);
      clog << "norm2(evec[" << i << "])= "<< norm2(evec[i]) << std::endl;
    }
    
    // initial Nblock_k steps
    for(int b=0; b<Nblock_k; ++b) blockwiseStep(lmd,lme,evec,f,f_copy,b);

    // restarting loop begins
    int iter;
    for(iter = 0; iter<MaxIter; ++iter){
      
      clog <<"#Restart iteration = "<< iter << std::endl;
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
      clog << "#Ritz value before shift: "<< std::endl;
      for(int i=0; i<Nm; ++i){
        std::cout.precision(13);
        std::cout << "[" << std::setw(4)<< std::setiosflags(std::ios_base::right) <<i<<"] ";
        std::cout << "Rval = "<<std::setw(20)<< std::setiosflags(std::ios_base::left)<< eval2[i] << std::endl;
      }
      
      //----------------------------------------------------------------------
      if ( Nm>Nk ) {
        clog <<" #Apply shifted QR transformations "<<std::endl;
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
            B[j].checkerboard = evec[k].checkerboard;
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
        clog << "#Ritz value after shift: "<< std::endl;
        for(int i=0; i<Nk; ++i){
          std::cout.precision(13);
          std::cout << "[" << std::setw(4)<< std::setiosflags(std::ios_base::right) <<i<<"] ";
          std::cout << "Rval = "<<std::setw(20)<< std::setiosflags(std::ios_base::left)<< eval2[i] << std::endl;
        }
      }
      //----------------------------------------------------------------------

      // Convergence test
      clog <<" #Convergence test: "<<std::endl;
      for(int k = 0; k<Nk; ++k) B[k]=0.0;
      for(int j = 0; j<Nk; ++j){
	for(int k = 0; k<Nk; ++k){
	  B[j].checkerboard = evec[k].checkerboard;
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
      
      clog <<" #modes converged: "<<Nconv<<std::endl;
      for(int i=0; i<Nconv; ++i){
	std::cout.precision(13);
        std::cout << "[" << std::setw(4)<< std::setiosflags(std::ios_base::right) <<Iconv[i]<<"] ";
	std::cout << "eval_conv = "<<std::setw(20)<< std::setiosflags(std::ios_base::left)<< eval2[Iconv[i]];
	std::cout << "   resid^2 = "<< std::setw(20)<< std::setiosflags(std::ios_base::right)<< resid[Iconv[i]]<< std::endl;
      } 

      if ( Nconv>=Nstop ) break;

    } // end of iter loop
    
    clog << std::string(74,'*') << std::endl;
    if ( Nconv<Nstop ) {
      clog << fname + " NOT converged ; Summary :\n";
    } else {
      clog << fname + " CONVERGED ; Summary :\n";
      // Sort convered eigenpairs.
      eval.resize(Nconv);
      evec.resize(Nconv,grid);
      for(int i=0; i<Nconv; ++i){
        eval[i] = eval2[Iconv[i]];
        evec[i] = B[Iconv[i]];
      }
      _sort.push(eval,evec,Nconv);
    }
    clog << std::string(74,'*') << std::endl;
    clog << " -- Iterations  = "<< iter   << "\n";
    //clog << " -- beta(k)     = "<< beta_k << "\n";
    clog << " -- Nconv       = "<< Nconv  << "\n";
    clog << std::string(74,'*') << std::endl;
  
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
    
    // converts block index to full indicies for an interval [L,R)
    int L = Nu*b;
    int R = Nu*(b+1);

    Real beta;
    
    // 3. wk:=Avk−βkv_{k−1}
    for (int k=L, u=0; k<R; ++k, ++u) {
      _poly(_Linop,evec[k],w[u]);      
    }
    
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
    
    // In block version, the steps 6 and 7 in Lanczos construction is
    // replaced by the QR decomposition of new basis block.
    // It results block version beta and orthonormal block basis. 
    // Here, QR decomposition is done by using Gram-Schmidt.
    for (int u=0; u<Nu; ++u) {
      for (int k=L; k<R; ++k) {
        lme[u][k] = 0.0;
      }
    }

#if 0
    beta = normalize(w[0]);
    for (int u=1; u<Nu; ++u) {
      //orthogonalize(w[u],w_copy,u);
      orthogonalize(w[u],w,u);
    }
#else
    // re-orthogonalization for numerical stability
    for (int u=0; u<Nu; ++u) {
      orthogonalize(w[u],evec,R);
    }
    // QR part
    for (int u=1; u<Nu; ++u) {
      orthogonalize(w[u],w,u);
    }
#endif
    
    for (int u=0; u<Nu; ++u) {
      //for (int v=0; v<Nu; ++v) {
      for (int v=u; v<Nu; ++v) {
        lme[u][L+v] = innerProduct(w[u],w_copy[v]);
      }
      lme[u][L+u] = real(lme[u][L+u]);  // force diagonal to be real
    }
    //lme[0][L] = beta;
    
    for (int u=0; u<Nu; ++u) {
      clog << "norm2(w[" << u << "])= "<< norm2(w[u]) << std::endl;
      for (int k=L+u; k<R; ++k) {
        clog <<" In block "<< b << ","; 
        std::cout <<" beta[" << u << "," << k-L << "] = ";
        std::cout << lme[u][k] << std::endl;
      }
    }
#if 0    
    // re-orthogonalization for numerical stability
    if (b>0) {
      for (int u=0; u<Nu; ++u) {
        orthogonalize(w[u],evec,R);
      }
      for (int u=1; u<Nu; ++u) {
        orthogonalize(w[u],w,u);
      }
    }
    //if (b>0) {
    //  orthogonalize_blockhead(w[0],evec,b,Nu);
    //  for (int u=1; u<Nu; ++u) {
    //    orthogonalize(w[u],w,u);
    //  }
    //}
#endif

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


  void diagonalize(std::vector<RealD>& eval, 
                   std::vector<std::vector<ComplexD>>& lmd, 
                   std::vector<std::vector<ComplexD>>& lme, 
		   int Nu, int Nk, int Nm,   
		   Eigen::MatrixXcd & Qt,
		   GridBase *grid)
  {
    Qt = Eigen::MatrixXcd::Identity(Nm,Nm);
    if ( diagonalisation == IRLdiagonaliseWithEigen ) { 
      diagonalize_Eigen(eval,lmd,lme,Nu,Nk,Nm,Qt,grid);
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
    //clog << "unpackHermitBlockTriDiagMatToEigen() begin" << '\n'; 
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
    //clog << "unpackHermitBlockTriDiagMatToEigen() end" << endl; 
  }
 

  void packHermitBlockTriDiagMatfromEigen(
         std::vector<std::vector<ComplexD>>& lmd,
         std::vector<std::vector<ComplexD>>& lme,
         int Nu, int Nb, int Nk, int Nm,
         Eigen::MatrixXcd& M)
  {
    //clog << "packHermitBlockTriDiagMatfromEigen() begin" << '\n'; 
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
    //clog << "packHermitBlockTriDiagMatfromEigen() end" << endl; 
  }


  // assume the input matrix M is a band matrix
  void shiftedQRDecompEigen(Eigen::MatrixXcd& M, int Nu, int Nm,
		            RealD Dsh,
		            Eigen::MatrixXcd& Qprod)
  {
    //clog << "shiftedQRDecompEigen() begin" << '\n'; 
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
    
    //clog << "shiftedQRDecompEigen() end" << endl; 
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
    
    clog << "matrix A before ColPivHouseholder" << std::endl;
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        clog << "A[" << i << "," << j << "] = " << A(i,j) << '\n';
      }
    }
    clog << std::endl;

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QRD(A);
    
    clog << "matrix A after ColPivHouseholder" << std::endl;
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        clog << "A[" << i << "," << j << "] = " << A(i,j) << '\n';
      }
    }
    clog << std::endl;
    
    clog << "HouseholderQ with sequence lenth = nonzeroPiviots" << std::endl;
    Q = QRD.householderQ().setLength(QRD.nonzeroPivots());
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        clog << "Q[" << i << "," << j << "] = " << Q(i,j) << '\n';
      }
    }
    clog << std::endl;
    
    clog << "HouseholderQ with sequence lenth = 1" << std::endl;
    Q = QRD.householderQ().setLength(1);
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        clog << "Q[" << i << "," << j << "] = " << Q(i,j) << '\n';
      }
    }
    clog << std::endl;
    
    clog << "HouseholderQ with sequence lenth = 2" << std::endl;
    Q = QRD.householderQ().setLength(2);
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        clog << "Q[" << i << "," << j << "] = " << Q(i,j) << '\n';
      }
    }
    clog << std::endl;
    
    clog << "matrixR" << std::endl;
    R = QRD.matrixR();
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        clog << "R[" << i << "," << j << "] = " << R(i,j) << '\n';
      }
    }
    clog << std::endl;

    clog << "rank = " << QRD.rank() << std::endl;
    clog << "threshold = " << QRD.threshold() << std::endl;
    
    clog << "matrixP" << std::endl;
    P = QRD.colsPermutation();
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        clog << "P[" << i << "," << j << "] = " << P(i,j) << '\n';
      }
    }
    clog << std::endl;


    clog << "QR decomposition without column pivoting" << std::endl;
    
    A(0,0) = 12.0;
    A(0,1) = -51.0;
    A(0,2) = 4.0;
    A(1,0) = 6.0;
    A(1,1) = 167.0;
    A(1,2) = -68.0;
    A(2,0) = -4.0;
    A(2,1) = 24.0;
    A(2,2) = -41.0;
    
    clog << "matrix A before Householder" << std::endl;
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        clog << "A[" << i << "," << j << "] = " << A(i,j) << '\n';
      }
    }
    clog << std::endl;
    
    Eigen::HouseholderQR<Eigen::MatrixXd> QRDplain(A);
    
    clog << "HouseholderQ" << std::endl;
    Q = QRDplain.householderQ();
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        clog << "Q[" << i << "," << j << "] = " << Q(i,j) << '\n';
      }
    }
    clog << std::endl;
    
    clog << "matrix A after Householder" << std::endl;
    for ( int i=0; i<3; i++ ) {
      for ( int j=0; j<3; j++ ) {
        clog << "A[" << i << "," << j << "] = " << A(i,j) << '\n';
      }
    }
    clog << std::endl;
  }

 };
}
#undef clog
#endif
