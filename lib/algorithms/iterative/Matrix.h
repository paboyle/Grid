    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/Matrix.h

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
#ifndef MATRIX_H
#define MATRIX_H

#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <complex>
#include <typeinfo>
#include <Grid.h>


/** Sign function **/
template <class T> T sign(T p){return ( p/abs(p) );}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////// Hijack STL containers for our wicked means /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> using Vector = Vector<T>;
template<class T> using Matrix = Vector<Vector<T> >;

template<class T> void Resize(Vector<T > & vec, int N) { vec.resize(N); }

template<class T> void Resize(Matrix<T > & mat, int N, int M) { 
  mat.resize(N);
  for(int i=0;i<N;i++){
    mat[i].resize(M);
  }
}
template<class T> void Size(Vector<T> & vec, int &N) 
{ 
  N= vec.size();
}
template<class T> void Size(Matrix<T> & mat, int &N,int &M) 
{ 
  N= mat.size();
  M= mat[0].size();
}
template<class T> void SizeSquare(Matrix<T> & mat, int &N) 
{ 
  int M; Size(mat,N,M);
  assert(N==M);
}
template<class T> void SizeSame(Matrix<T> & mat1,Matrix<T> &mat2, int &N1,int &M1) 
{ 
  int N2,M2;
  Size(mat1,N1,M1);
  Size(mat2,N2,M2);
  assert(N1==N2);
  assert(M1==M2);
}

//*****************************************
//*	(Complex) Vector operations	*
//*****************************************

/**Conj of a Vector **/
template <class T> Vector<T> conj(Vector<T> p){
	Vector<T> q(p.size());
	for(int i=0;i<p.size();i++){q[i] = conj(p[i]);}
	return q;
}

/** Norm of a Vector**/
template <class T> T norm(Vector<T> p){
	T sum = 0;
	for(int i=0;i<p.size();i++){sum = sum + p[i]*conj(p[i]);}
	return abs(sqrt(sum));
}

/** Norm squared of a Vector **/
template <class T> T norm2(Vector<T> p){
	T sum = 0;
	for(int i=0;i<p.size();i++){sum = sum + p[i]*conj(p[i]);}
	return abs((sum));
}

/** Sum elements of a Vector **/
template <class T> T trace(Vector<T> p){
	T sum = 0;
	for(int i=0;i<p.size();i++){sum = sum + p[i];}
	return sum;
}

/** Fill a Vector with constant c **/
template <class T> void Fill(Vector<T> &p, T c){
	for(int i=0;i<p.size();i++){p[i] = c;}
}
/** Normalize a Vector **/
template <class T> void normalize(Vector<T> &p){
	T m = norm(p);
	if( abs(m) > 0.0) for(int i=0;i<p.size();i++){p[i] /= m;}
}
/** Vector by scalar **/
template <class T, class U> Vector<T> times(Vector<T> p, U s){
	for(int i=0;i<p.size();i++){p[i] *= s;}
	return p;
}
template <class T, class U> Vector<T> times(U s, Vector<T> p){
	for(int i=0;i<p.size();i++){p[i] *= s;}
	return p;
}
/** inner product of a and b = conj(a) . b **/
template <class T> T inner(Vector<T> a, Vector<T> b){
	T m = 0.;
	for(int i=0;i<a.size();i++){m = m + conj(a[i])*b[i];}
	return m;
}
/** sum of a and b = a + b **/
template <class T> Vector<T> add(Vector<T> a, Vector<T> b){
	Vector<T> m(a.size());
	for(int i=0;i<a.size();i++){m[i] = a[i] + b[i];}
	return m;
}
/** sum of a and b = a - b **/
template <class T> Vector<T> sub(Vector<T> a, Vector<T> b){
	Vector<T> m(a.size());
	for(int i=0;i<a.size();i++){m[i] = a[i] - b[i];}
	return m;
}

/** 
 *********************************
 *	Matrices	         *
 *********************************
 **/

template<class T> void Fill(Matrix<T> & mat, T&val) { 
  int N,M;
  Size(mat,N,M);
  for(int i=0;i<N;i++){
  for(int j=0;j<M;j++){
    mat[i][j] = val;
  }}
}

/** Transpose of a matrix **/
Matrix<T> Transpose(Matrix<T> & mat){
  int N,M;
  Size(mat,N,M);
  Matrix C; Resize(C,M,N);
  for(int i=0;i<M;i++){
  for(int j=0;j<N;j++){
    C[i][j] = mat[j][i];
  }} 
  return C;
}
/** Set Matrix to unit matrix **/
template<class T> void Unity(Matrix<T> &mat){
  int N;  SizeSquare(mat,N);
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      if ( i==j ) A[i][j] = 1;
      else        A[i][j] = 0;
    } 
  } 
}
/** Add C * I to matrix **/
template<class T>
void PlusUnit(Matrix<T> & A,T c){
  int dim;  SizeSquare(A,dim);
  for(int i=0;i<dim;i++){A[i][i] = A[i][i] + c;} 
}

/** return the Hermitian conjugate of matrix **/
Matrix<T> HermitianConj(Matrix<T> &mat){

  int dim; SizeSquare(mat,dim);

  Matrix<T> C; Resize(C,dim,dim);

  for(int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      C[i][j] = conj(mat[j][i]);
    } 
  } 
  return C;
}

/** return diagonal entries as a Vector **/
Vector<T> diag(Matrix<T> &A)
{
  int dim; SizeSquare(A,dim);
  Vector<T> d; Resize(d,dim);

  for(int i=0;i<dim;i++){
    d[i] = A[i][i];
  }
  return d;
}

/** Left multiply by a Vector **/
Vector<T> operator *(Vector<T> &B,Matrix<T> &A)
{
  int K,M,N; 
  Size(B,K);
  Size(A,M,N);
  assert(K==M);
  
  Vector<T> C; Resize(C,N);

  for(int j=0;j<N;j++){
    T sum = 0.0;
    for(int i=0;i<M;i++){
      sum += B[i] * A[i][j];
    }
    C[j] =  sum;
  }
  return C; 
}

/** return 1/diagonal entries as a Vector **/
Vector<T> inv_diag(Matrix<T> & A){
  int dim; SizeSquare(A,dim);
  Vector<T> d; Resize(d,dim);
  for(int i=0;i<dim;i++){
    d[i] = 1.0/A[i][i];
  }
  return d;
}
/** Matrix Addition **/
inline Matrix<T> operator + (Matrix<T> &A,Matrix<T> &B)
{
  int N,M  ; SizeSame(A,B,N,M);
  Matrix C; Resize(C,N,M);
  for(int i=0;i<N;i++){
    for(int j=0;j<M;j++){
      C[i][j] = A[i][j] +  B[i][j];
    } 
  } 
  return C;
} 
/** Matrix Subtraction **/
inline Matrix<T> operator- (Matrix<T> & A,Matrix<T> &B){
  int N,M  ; SizeSame(A,B,N,M);
  Matrix C; Resize(C,N,M);
  for(int i=0;i<N;i++){
  for(int j=0;j<M;j++){
    C[i][j] = A[i][j] -  B[i][j];
  }}
  return C;
} 

/** Matrix scalar multiplication **/
inline Matrix<T> operator* (Matrix<T> & A,T c){
  int N,M; Size(A,N,M);
  Matrix C; Resize(C,N,M);
  for(int i=0;i<N;i++){
  for(int j=0;j<M;j++){
    C[i][j] = A[i][j]*c;
  }} 
  return C;
} 
/** Matrix Matrix multiplication **/
inline Matrix<T> operator* (Matrix<T> &A,Matrix<T> &B){
  int K,L,N,M;
  Size(A,K,L);
  Size(B,N,M); assert(L==N);
  Matrix C; Resize(C,K,M);

  for(int i=0;i<K;i++){
    for(int j=0;j<M;j++){
      T sum = 0.0;
      for(int k=0;k<N;k++) sum += A[i][k]*B[k][j];
      C[i][j] =sum;
    }
  }
  return C; 
} 
/** Matrix Vector multiplication **/
inline Vector<T> operator* (Matrix<T> &A,Vector<T> &B){
  int M,N,K;
  Size(A,N,M);
  Size(B,K); assert(K==M);
  Vector<T> C; Resize(C,N);
  for(int i=0;i<N;i++){
    T sum = 0.0;
    for(int j=0;j<M;j++) sum += A[i][j]*B[j];
    C[i] =  sum;
  }
  return C; 
} 

/** Some version of Matrix norm **/
/*
inline T Norm(){ // this is not a usual L2 norm
    T norm = 0;
    for(int i=0;i<dim;i++){
      for(int j=0;j<dim;j++){
	norm += abs(A[i][j]);
    }}
    return norm;
  }
*/

/** Some version of Matrix norm **/
template<class T> T LargestDiag(Matrix<T> &A)
{
  int dim ; SizeSquare(A,dim); 

  T ld = abs(A[0][0]);
  for(int i=1;i<dim;i++){
    T cf = abs(A[i][i]);
    if(abs(cf) > abs(ld) ){ld = cf;}
  }
  return ld;
}

/** Look for entries on the leading subdiagonal that are smaller than 'small' **/
template <class T,class U> int Chop_subdiag(Matrix<T> &A,T norm, int offset, U small)
{
  int dim; SizeSquare(A,dim);
  for(int l = dim - 1 - offset; l >= 1; l--) {             		
    if((U)abs(A[l][l - 1]) < (U)small) {
      A[l][l-1]=(U)0.0;
      return l;
    }
  }
  return 0;
}

/** Look for entries on the leading subdiagonal that are smaller than 'small' **/
template <class T,class U> int Chop_symm_subdiag(Matrix<T> & A,T norm, int offset, U small) 
{
  int dim; SizeSquare(A,dim);
  for(int l = dim - 1 - offset; l >= 1; l--) {
    if((U)abs(A[l][l - 1]) < (U)small) {
      A[l][l - 1] = (U)0.0;
      A[l - 1][l] = (U)0.0;
      return l;
    }
  }
  return 0;
}
/**Assign a submatrix to a larger one**/
template<class T>
void AssignSubMtx(Matrix<T> & A,int row_st, int row_end, int col_st, int col_end, Matrix<T> &S)
{
  for(int i = row_st; i<row_end; i++){
    for(int j = col_st; j<col_end; j++){
      A[i][j] = S[i - row_st][j - col_st];
    }
  }
}

/**Get a square submatrix**/
template <class T>
Matrix<T> GetSubMtx(Matrix<T> &A,int row_st, int row_end, int col_st, int col_end)
{
  Matrix<T> H; Resize(row_end - row_st,col_end-col_st);

  for(int i = row_st; i<row_end; i++){
  for(int j = col_st; j<col_end; j++){
    H[i-row_st][j-col_st]=A[i][j];
  }}
  return H;
}
  
 /**Assign a submatrix to a larger one NB remember Vector Vectors are transposes of the matricies they represent**/
template<class T>
void AssignSubMtx(Matrix<T> & A,int row_st, int row_end, int col_st, int col_end, Matrix<T> &S)
{
  for(int i = row_st; i<row_end; i++){
  for(int j = col_st; j<col_end; j++){
    A[i][j] = S[i - row_st][j - col_st];
  }}
}
  
/** compute b_i A_ij b_j **/ // surprised no Conj
template<class T> T proj(Matrix<T> A, Vector<T> B){
  int dim; SizeSquare(A,dim);
  int dimB; Size(B,dimB);
  assert(dimB==dim);
  T C = 0;
  for(int i=0;i<dim;i++){
    T sum = 0.0;
    for(int j=0;j<dim;j++){
      sum += A[i][j]*B[j];
    }
    C +=  B[i]*sum; // No conj?
  }
  return C; 
}


/*
 *************************************************************
 *
 * Matrix Vector products
 *
 *************************************************************
 */
// Instead make a linop and call my CG;

/// q -> q Q
template <class T,class Fermion> void times(Vector<Fermion> &q, Matrix<T> &Q)
{
  int M; SizeSquare(Q,M);
  int N; Size(q,N); 
  assert(M==N);

  times(q,Q,N);
}

/// q -> q Q
template <class T> void times(multi1d<LatticeFermion> &q, Matrix<T> &Q, int N)
{
  GridBase *grid = q[0]._grid;
  int M; SizeSquare(Q,M);
  int K; Size(q,K); 
  assert(N<M);
  assert(N<K);
  Vector<Fermion> S(N,grid );
  for(int j=0;j<N;j++){
    S[j] = zero;
    for(int k=0;k<N;k++){
      S[j] = S[j] +  q[k]* Q[k][j]; 
    }
  }
  for(int j=0;j<q.size();j++){
    q[j] = S[j];
  }
}
#endif
