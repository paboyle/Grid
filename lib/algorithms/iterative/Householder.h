    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/Householder.h

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
#ifndef HOUSEHOLDER_H
#define HOUSEHOLDER_H

#define TIMER(A) std::cout << GridLogMessage << __FUNC__ << " file "<< __FILE__ <<" line " << __LINE__ << std::endl;
#define ENTER()  std::cout << GridLogMessage << "ENTRY "<<__FUNC__ << " file "<< __FILE__ <<" line " << __LINE__ << std::endl;
#define LEAVE()  std::cout << GridLogMessage << "EXIT  "<<__FUNC__ << " file "<< __FILE__ <<" line " << __LINE__ << std::endl;

#include <cstdlib>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <complex>
#include <algorithm>

namespace Grid {
/** Comparison function for finding the max element in a vector **/
template <class T> bool cf(T i, T j) { 
  return abs(i) < abs(j); 
}

/** 
	Calculate a real Givens angle 
 **/
template <class T> inline void Givens_calc(T y, T z, T &c, T &s){

  RealD mz = (RealD)abs(z);
  
  if(mz==0.0){
    c = 1; s = 0;
  }
  if(mz >= (RealD)abs(y)){
    T t = -y/z;
    s = (T)1.0 / sqrt ((T)1.0 + t * t);
    c = s * t;
  } else {
    T t = -z/y;
    c = (T)1.0 / sqrt ((T)1.0 + t * t);
    s = c * t;
  }
}

template <class T> inline void Givens_mult(DenseMatrix<T> &A,  int i, int k, T c, T s, int dir)
{
  int q ; SizeSquare(A,q);

  if(dir == 0){
    for(int j=0;j<q;j++){
      T nu = A[i][j];
      T w  = A[k][j];
      A[i][j] = (c*nu + s*w);
      A[k][j] = (-s*nu + c*w);
    }
  }

  if(dir == 1){
    for(int j=0;j<q;j++){
      T nu = A[j][i];
      T w  = A[j][k];
      A[j][i] = (c*nu - s*w);
      A[j][k] = (s*nu + c*w);
    }
  }
}

/**
	from input = x;
	Compute the complex Householder vector, v, such that
	P = (I - b v transpose(v) )
	b = 2/v.v

	P | x |    | x | k = 0
	| x |    | 0 | 
	| x | =  | 0 |
	| x |    | 0 | j = 3
	| x |	   | x |

	These are the "Unreduced" Householder vectors.

 **/
template <class T> inline void Householder_vector(DenseVector<T> input, int k, int j, DenseVector<T> &v, T &beta)
{
  int N ; Size(input,N);
  T m = *max_element(input.begin() + k, input.begin() + j + 1, cf<T> );

  if(abs(m) > 0.0){
    T alpha = 0;

    for(int i=k; i<j+1; i++){
      v[i] = input[i]/m;
      alpha = alpha + v[i]*conj(v[i]);
    }
    alpha = sqrt(alpha);
    beta = (T)1.0/(alpha*(alpha + abs(v[k]) ));

    if(abs(v[k]) > 0.0)  v[k] = v[k] + (v[k]/abs(v[k]))*alpha;
    else                 v[k] = -alpha;
  } else{
    for(int i=k; i<j+1; i++){
      v[i] = 0.0;
    } 
  }
}

/**
	from input = x;
	Compute the complex Householder vector, v, such that
	P = (I - b v transpose(v) )
	b = 2/v.v

	Px = alpha*e_dir

	These are the "Unreduced" Householder vectors.

 **/

template <class T> inline void Householder_vector(DenseVector<T> input, int k, int j, int dir, DenseVector<T> &v, T &beta)
{
  int N = input.size();
  T m = *max_element(input.begin() + k, input.begin() + j + 1, cf);
  
  if(abs(m) > 0.0){
    T alpha = 0;

    for(int i=k; i<j+1; i++){
      v[i] = input[i]/m;
      alpha = alpha + v[i]*conj(v[i]);
    }
    
    alpha = sqrt(alpha);
    beta = 1.0/(alpha*(alpha + abs(v[dir]) ));
	
    if(abs(v[dir]) > 0.0) v[dir] = v[dir] + (v[dir]/abs(v[dir]))*alpha;
    else                  v[dir] = -alpha;
  }else{
    for(int i=k; i<j+1; i++){
      v[i] = 0.0;
    } 
  }
}

/**
	Compute the product PA if trans = 0
	AP if trans = 1
	P = (I - b v transpose(v) )
	b = 2/v.v
	start at element l of matrix A
	v is of length j - k + 1 of v are nonzero
 **/

template <class T> inline void Householder_mult(DenseMatrix<T> &A , DenseVector<T> v, T beta, int l, int k, int j, int trans)
{
  int N ; SizeSquare(A,N);

  if(abs(beta) > 0.0){
    for(int p=l; p<N; p++){
      T s = 0;
      if(trans==0){
	for(int i=k;i<j+1;i++) s += conj(v[i-k])*A[i][p];
	s *= beta;
	for(int i=k;i<j+1;i++){ A[i][p] = A[i][p]-s*conj(v[i-k]);}
      } else {
	for(int i=k;i<j+1;i++){ s += conj(v[i-k])*A[p][i];}
	s *= beta;
	for(int i=k;i<j+1;i++){ A[p][i]=A[p][i]-s*conj(v[i-k]);}
      }
    }
  }
}

/**
	Compute the product PA if trans = 0
	AP if trans = 1
	P = (I - b v transpose(v) )
	b = 2/v.v
	start at element l of matrix A
	v is of length j - k + 1 of v are nonzero
	A is tridiagonal
 **/
template <class T> inline void Householder_mult_tri(DenseMatrix<T> &A , DenseVector<T> v, T beta, int l, int M, int k, int j, int trans)
{
  if(abs(beta) > 0.0){

    int N ; SizeSquare(A,N);

    DenseMatrix<T> tmp; Resize(tmp,N,N); Fill(tmp,0); 

    T s;
    for(int p=l; p<M; p++){
      s = 0;
      if(trans==0){
	for(int i=k;i<j+1;i++) s = s + conj(v[i-k])*A[i][p];
      }else{
	for(int i=k;i<j+1;i++) s = s + v[i-k]*A[p][i];
      }
      s = beta*s;
      if(trans==0){
	for(int i=k;i<j+1;i++) tmp[i][p] = tmp(i,p) - s*v[i-k];
      }else{
	for(int i=k;i<j+1;i++) tmp[p][i] = tmp[p][i] - s*conj(v[i-k]);
      }
    }
    for(int p=l; p<M; p++){
      if(trans==0){
	for(int i=k;i<j+1;i++) A[i][p] = A[i][p] + tmp[i][p];
      }else{
	for(int i=k;i<j+1;i++) A[p][i] = A[p][i] + tmp[p][i];
      }
    }
  }
}
}
#endif
