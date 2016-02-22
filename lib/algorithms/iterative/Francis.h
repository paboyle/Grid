    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/Francis.h

    Copyright (C) 2015

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef FRANCIS_H
#define FRANCIS_H

#include <cstdlib>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <complex>
#include <algorithm>

//#include <timer.h>
//#include <lapacke.h>
//#include <Eigen/Dense>

namespace Grid {

template <class T> int SymmEigensystem(DenseMatrix<T > &Ain, DenseVector<T> &evals, DenseMatrix<T> &evecs, RealD small);
template <class T> int     Eigensystem(DenseMatrix<T > &Ain, DenseVector<T> &evals, DenseMatrix<T> &evecs, RealD small);

/**
  Find the eigenvalues of an upper hessenberg matrix using the Francis QR algorithm.
H =
      x  x  x  x  x  x  x  x  x
      x  x  x  x  x  x  x  x  x
      0  x  x  x  x  x  x  x  x
      0  0  x  x  x  x  x  x  x
      0  0  0  x  x  x  x  x  x
      0  0  0  0  x  x  x  x  x
      0  0  0  0  0  x  x  x  x
      0  0  0  0  0  0  x  x  x
      0  0  0  0  0  0  0  x  x
Factorization is P T P^H where T is upper triangular (mod cc blocks) and P is orthagonal/unitary.
**/
template <class T>
int QReigensystem(DenseMatrix<T> &Hin, DenseVector<T> &evals, DenseMatrix<T> &evecs, RealD small)
{
  DenseMatrix<T> H = Hin; 

  int N ; SizeSquare(H,N);
  int M = N;

  Fill(evals,0);
  Fill(evecs,0);

  T s,t,x=0,y=0,z=0;
  T u,d;
  T apd,amd,bc;
  DenseVector<T> p(N,0);
  T nrm = Norm(H);    ///DenseMatrix Norm
  int n, m;
  int e = 0;
  int it = 0;
  int tot_it = 0;
  int l = 0;
  int r = 0;
  DenseMatrix<T> P; Resize(P,N,N); Unity(P);
  DenseVector<int> trows(N,0);

  /// Check if the matrix is really hessenberg, if not abort
  RealD sth = 0;
  for(int j=0;j<N;j++){
    for(int i=j+2;i<N;i++){
      sth = abs(H[i][j]);
      if(sth > small){
	std::cout << "Non hessenberg H = " << sth << " > " << small << std::endl;
	exit(1);
      }
    }
  }

  do{
    std::cout << "Francis QR Step N = " << N << std::endl;
    /** Check for convergence
      x  x  x  x  x
      0  x  x  x  x
      0  0  x  x  x
      0  0  x  x  x
      0  0  0  0  x
      for this matrix l = 4
     **/
    do{
      l = Chop_subdiag(H,nrm,e,small);
      r = 0;    ///May have converged on more than one eval
      ///Single eval
      if(l == N-1){
        evals[e] = H[l][l];
        N--; e++; r++; it = 0;
      }
      ///RealD eval
      if(l == N-2){
        trows[l+1] = 1;    ///Needed for UTSolve
        apd = H[l][l] + H[l+1][l+1];
        amd = H[l][l] - H[l+1][l+1];
        bc =  (T)4.0*H[l+1][l]*H[l][l+1];
        evals[e]   = (T)0.5*( apd + sqrt(amd*amd + bc) );
        evals[e+1] = (T)0.5*( apd - sqrt(amd*amd + bc) );
        N-=2; e+=2; r++; it = 0;
      }
    } while(r>0);

    if(N ==0) break;

    DenseVector<T > ck; Resize(ck,3);
    DenseVector<T> v;   Resize(v,3);

    for(int m = N-3; m >= l; m--){
      ///Starting vector essentially random shift.
      if(it%10 == 0 && N >= 3 && it > 0){
        s = (T)1.618033989*( abs( H[N-1][N-2] ) + abs( H[N-2][N-3] ) );
        t = (T)0.618033989*( abs( H[N-1][N-2] ) + abs( H[N-2][N-3] ) );
        x = H[m][m]*H[m][m] + H[m][m+1]*H[m+1][m] - s*H[m][m] + t;
        y = H[m+1][m]*(H[m][m] + H[m+1][m+1] - s);
        z = H[m+1][m]*H[m+2][m+1];
      }
      ///Starting vector implicit Q theorem
      else{
        s = (H[N-2][N-2] + H[N-1][N-1]);
        t = (H[N-2][N-2]*H[N-1][N-1] - H[N-2][N-1]*H[N-1][N-2]);
        x = H[m][m]*H[m][m] + H[m][m+1]*H[m+1][m] - s*H[m][m] + t;
        y = H[m+1][m]*(H[m][m] + H[m+1][m+1] - s);
        z = H[m+1][m]*H[m+2][m+1];
      }
      ck[0] = x; ck[1] = y; ck[2] = z;

      if(m == l) break;

      /** Some stupid thing from numerical recipies, seems to work**/
      // PAB.. for heaven's sake quote page, purpose, evidence it works.
      //       what sort of comment is that!?!?!?
      u=abs(H[m][m-1])*(abs(y)+abs(z));
      d=abs(x)*(abs(H[m-1][m-1])+abs(H[m][m])+abs(H[m+1][m+1]));
      if ((T)abs(u+d) == (T)abs(d) ){
	l = m; break;
      }

      //if (u < small){l = m; break;}
    }
    if(it > 100000){
     std::cout << "QReigensystem: bugger it got stuck after 100000 iterations" << std::endl;
     std::cout << "got " << e << " evals " << l << " " << N << std::endl;
      exit(1);
    }
    normalize(ck);    ///Normalization cancels in PHP anyway
    T beta;
    Householder_vector<T >(ck, 0, 2, v, beta);
    Householder_mult<T >(H,v,beta,0,l,l+2,0);
    Householder_mult<T >(H,v,beta,0,l,l+2,1);
    ///Accumulate eigenvector
    Householder_mult<T >(P,v,beta,0,l,l+2,1);
    int sw = 0;      ///Are we on the last row?
    for(int k=l;k<N-2;k++){
      x = H[k+1][k];
      y = H[k+2][k];
      z = (T)0.0;
      if(k+3 <= N-1){
	z = H[k+3][k];
      } else{
	sw = 1; 
	v[2] = (T)0.0;
      }
      ck[0] = x; ck[1] = y; ck[2] = z;
      normalize(ck);
      Householder_vector<T >(ck, 0, 2-sw, v, beta);
      Householder_mult<T >(H,v, beta,0,k+1,k+3-sw,0);
      Householder_mult<T >(H,v, beta,0,k+1,k+3-sw,1);
      ///Accumulate eigenvector
      Householder_mult<T >(P,v, beta,0,k+1,k+3-sw,1);
    }
    it++;
    tot_it++;
  }while(N > 1);
  N = evals.size();
  ///Annoying - UT solves in reverse order;
  DenseVector<T> tmp; Resize(tmp,N);
  for(int i=0;i<N;i++){
    tmp[i] = evals[N-i-1];
  } 
  evals = tmp;
  UTeigenvectors(H, trows, evals, evecs);
  for(int i=0;i<evals.size();i++){evecs[i] = P*evecs[i]; normalize(evecs[i]);}
  return tot_it;
}

template <class T>
int my_Wilkinson(DenseMatrix<T> &Hin, DenseVector<T> &evals, DenseMatrix<T> &evecs, RealD small)
{
  /**
  Find the eigenvalues of an upper Hessenberg matrix using the Wilkinson QR algorithm.
  H =
  x  x  0  0  0  0
  x  x  x  0  0  0
  0  x  x  x  0  0
  0  0  x  x  x  0
  0  0  0  x  x  x
  0  0  0  0  x  x
  Factorization is P T P^H where T is upper triangular (mod cc blocks) and P is orthagonal/unitary.  **/
  return my_Wilkinson(Hin, evals, evecs, small, small);
}

template <class T>
int my_Wilkinson(DenseMatrix<T> &Hin, DenseVector<T> &evals, DenseMatrix<T> &evecs, RealD small, RealD tol)
{
  int N; SizeSquare(Hin,N);
  int M = N;

  ///I don't want to modify the input but matricies must be passed by reference
  //Scale a matrix by its "norm"
  //RealD Hnorm = abs( Hin.LargestDiag() ); H =  H*(1.0/Hnorm);
  DenseMatrix<T> H;  H = Hin;
  
  RealD Hnorm = abs(Norm(Hin));
  H = H * (1.0 / Hnorm);

  // TODO use openmp and memset
  Fill(evals,0);
  Fill(evecs,0);

  T s, t, x = 0, y = 0, z = 0;
  T u, d;
  T apd, amd, bc;
  DenseVector<T> p; Resize(p,N); Fill(p,0);

  T nrm = Norm(H);    ///DenseMatrix Norm
  int n, m;
  int e = 0;
  int it = 0;
  int tot_it = 0;
  int l = 0;
  int r = 0;
  DenseMatrix<T> P; Resize(P,N,N);
  Unity(P);
  DenseVector<int> trows(N, 0);
  /// Check if the matrix is really symm tridiag
  RealD sth = 0;
  for(int j = 0; j < N; ++j)
  {
    for(int i = j + 2; i < N; ++i)
    {
      if(abs(H[i][j]) > tol || abs(H[j][i]) > tol)
      {
	std::cout << "Non Tridiagonal H(" << i << ","<< j << ") = |" << Real( real( H[j][i] ) ) << "| > " << tol << std::endl;
	std::cout << "Warning tridiagonalize and call again" << std::endl;
        // exit(1); // see what is going on
        //return;
      }
    }
  }

  do{
    do{
      //Jasper
      //Check if the subdiagonal term is small enough (<small)
      //if true then it is converged.
      //check start from H.dim - e - 1
      //How to deal with more than 2 are converged?
      //What if Chop_symm_subdiag return something int the middle?
      //--------------
      l = Chop_symm_subdiag(H,nrm, e, small);
      r = 0;    ///May have converged on more than one eval
      //Jasper
      //In this case
      // x  x  0  0  0  0
      // x  x  x  0  0  0
      // 0  x  x  x  0  0
      // 0  0  x  x  x  0
      // 0  0  0  x  x  0
      // 0  0  0  0  0  x  <- l
      //--------------
      ///Single eval
      if(l == N - 1)
      {
        evals[e] = H[l][l];
        N--;
        e++;
        r++;
        it = 0;
      }
      //Jasper
      // x  x  0  0  0  0
      // x  x  x  0  0  0
      // 0  x  x  x  0  0
      // 0  0  x  x  0  0
      // 0  0  0  0  x  x  <- l
      // 0  0  0  0  x  x
      //--------------
      ///RealD eval
      if(l == N - 2)
      {
        trows[l + 1] = 1;    ///Needed for UTSolve
        apd = H[l][l] + H[l + 1][ l + 1];
        amd = H[l][l] - H[l + 1][l + 1];
        bc =  (T) 4.0 * H[l + 1][l] * H[l][l + 1];
        evals[e] = (T) 0.5 * (apd + sqrt(amd * amd + bc));
        evals[e + 1] = (T) 0.5 * (apd - sqrt(amd * amd + bc));
        N -= 2;
        e += 2;
        r++;
        it = 0;
      }
    }while(r > 0);
    //Jasper
    //Already converged
    //--------------
    if(N == 0) break;

    DenseVector<T> ck,v; Resize(ck,2); Resize(v,2);

    for(int m = N - 3; m >= l; m--)
    {
      ///Starting vector essentially random shift.
      if(it%10 == 0 && N >= 3 && it > 0)
      {
        t = abs(H[N - 1][N - 2]) + abs(H[N - 2][N - 3]);
        x = H[m][m] - t;
        z = H[m + 1][m];
      } else {
      ///Starting vector implicit Q theorem
        d = (H[N - 2][N - 2] - H[N - 1][N - 1]) * (T) 0.5;
        t =  H[N - 1][N - 1] - H[N - 1][N - 2] * H[N - 1][N - 2] 
	  / (d + sign(d) * sqrt(d * d + H[N - 1][N - 2] * H[N - 1][N - 2]));
        x = H[m][m] - t;
        z = H[m + 1][m];
      }
      //Jasper
      //why it is here????
      //-----------------------
      if(m == l)
        break;

      u = abs(H[m][m - 1]) * (abs(y) + abs(z));
      d = abs(x) * (abs(H[m - 1][m - 1]) + abs(H[m][m]) + abs(H[m + 1][m + 1]));
      if ((T)abs(u + d) == (T)abs(d))
      {
        l = m;
        break;
      }
    }
    //Jasper
    if(it > 1000000)
    {
      std::cout << "Wilkinson: bugger it got stuck after 100000 iterations" << std::endl;
      std::cout << "got " << e << " evals " << l << " " << N << std::endl;
      exit(1);
    }
    //
    T s, c;
    Givens_calc<T>(x, z, c, s);
    Givens_mult<T>(H, l, l + 1, c, -s, 0);
    Givens_mult<T>(H, l, l + 1, c,  s, 1);
    Givens_mult<T>(P, l, l + 1, c,  s, 1);
    //
    for(int k = l; k < N - 2; ++k)
    {
      x = H.A[k + 1][k];
      z = H.A[k + 2][k];
      Givens_calc<T>(x, z, c, s);
      Givens_mult<T>(H, k + 1, k + 2, c, -s, 0);
      Givens_mult<T>(H, k + 1, k + 2, c,  s, 1);
      Givens_mult<T>(P, k + 1, k + 2, c,  s, 1);
    }
    it++;
    tot_it++;
  }while(N > 1);

  N = evals.size();
  ///Annoying - UT solves in reverse order;
  DenseVector<T> tmp(N);
  for(int i = 0; i < N; ++i)
    tmp[i] = evals[N-i-1];
  evals = tmp;
  //
  UTeigenvectors(H, trows, evals, evecs);
  //UTSymmEigenvectors(H, trows, evals, evecs);
  for(int i = 0; i < evals.size(); ++i)
  {
    evecs[i] = P * evecs[i];
    normalize(evecs[i]);
    evals[i] = evals[i] * Hnorm;
  }
  // // FIXME this is to test
  // Hin.write("evecs3", evecs);
  // Hin.write("evals3", evals);
  // // check rsd
  // for(int i = 0; i < M; i++) {
  //   vector<T> Aevec = Hin * evecs[i];
  //   RealD norm2(0.);
  //   for(int j = 0; j < M; j++) {
  //     norm2 += (Aevec[j] - evals[i] * evecs[i][j]) * (Aevec[j] - evals[i] * evecs[i][j]);
  //   }
  // }
  return tot_it;
}

template <class T>
void Hess(DenseMatrix<T > &A, DenseMatrix<T> &Q, int start){

  /**
  turn a matrix A =
  x  x  x  x  x
  x  x  x  x  x
  x  x  x  x  x
  x  x  x  x  x
  x  x  x  x  x
  into
  x  x  x  x  x
  x  x  x  x  x
  0  x  x  x  x
  0  0  x  x  x
  0  0  0  x  x
  with householder rotations
  Slow.
  */
  int N ; SizeSquare(A,N);
  DenseVector<T > p; Resize(p,N); Fill(p,0);

  for(int k=start;k<N-2;k++){
    //cerr << "hess" << k << std::endl;
    DenseVector<T > ck,v; Resize(ck,N-k-1); Resize(v,N-k-1);
    for(int i=k+1;i<N;i++){ck[i-k-1] = A(i,k);}  ///kth column
    normalize(ck);    ///Normalization cancels in PHP anyway
    T beta;
    Householder_vector<T >(ck, 0, ck.size()-1, v, beta);  ///Householder vector
    Householder_mult<T>(A,v,beta,start,k+1,N-1,0);  ///A -> PA
    Householder_mult<T >(A,v,beta,start,k+1,N-1,1);  ///PA -> PAP^H
    ///Accumulate eigenvector
    Householder_mult<T >(Q,v,beta,start,k+1,N-1,1);  ///Q -> QP^H
  }
  /*for(int l=0;l<N-2;l++){
    for(int k=l+2;k<N;k++){
    A(0,k,l);
    }
    }*/
}

template <class T>
void Tri(DenseMatrix<T > &A, DenseMatrix<T> &Q, int start){
///Tridiagonalize a matrix
  int N; SizeSquare(A,N);
  Hess(A,Q,start);
  /*for(int l=0;l<N-2;l++){
    for(int k=l+2;k<N;k++){
    A(0,l,k);
    }
    }*/
}

template <class T>
void ForceTridiagonal(DenseMatrix<T> &A){
///Tridiagonalize a matrix
  int N ; SizeSquare(A,N);
  for(int l=0;l<N-2;l++){
    for(int k=l+2;k<N;k++){
      A[l][k]=0;
      A[k][l]=0;
    }
  }
}

template <class T>
int my_SymmEigensystem(DenseMatrix<T > &Ain, DenseVector<T> &evals, DenseVector<DenseVector<T> > &evecs, RealD small){
  ///Solve a symmetric eigensystem, not necessarily in tridiagonal form
  int N; SizeSquare(Ain,N);
  DenseMatrix<T > A; A = Ain;
  DenseMatrix<T > Q; Resize(Q,N,N); Unity(Q);
  Tri(A,Q,0);
  int it = my_Wilkinson<T>(A, evals, evecs, small);
  for(int k=0;k<N;k++){evecs[k] = Q*evecs[k];}
  return it;
}


template <class T>
int Wilkinson(DenseMatrix<T> &Ain, DenseVector<T> &evals, DenseVector<DenseVector<T> > &evecs, RealD small){
  return my_Wilkinson(Ain, evals, evecs, small);
}

template <class T>
int SymmEigensystem(DenseMatrix<T> &Ain, DenseVector<T> &evals, DenseVector<DenseVector<T> > &evecs, RealD small){
  return my_SymmEigensystem(Ain, evals, evecs, small);
}

template <class T>
int Eigensystem(DenseMatrix<T > &Ain, DenseVector<T> &evals, DenseVector<DenseVector<T> > &evecs, RealD small){
///Solve a general eigensystem, not necessarily in tridiagonal form
  int N = Ain.dim;
  DenseMatrix<T > A(N); A = Ain;
  DenseMatrix<T > Q(N);Q.Unity();
  Hess(A,Q,0);
  int it = QReigensystem<T>(A, evals, evecs, small);
  for(int k=0;k<N;k++){evecs[k] = Q*evecs[k];}
  return it;
}

}
#endif
