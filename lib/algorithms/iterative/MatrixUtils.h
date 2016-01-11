    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/MatrixUtils.h

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
#ifndef GRID_MATRIX_UTILS_H
#define GRID_MATRIX_UTILS_H

namespace Grid {

  namespace MatrixUtils { 

    template<class T> inline void Size(Matrix<T>& A,int &N,int &M){
      N=A.size(); assert(N>0);
      M=A[0].size();
      for(int i=0;i<N;i++){
	assert(A[i].size()==M);
      }
    }

    template<class T> inline void SizeSquare(Matrix<T>& A,int &N)
    {
      int M;
      Size(A,N,M);
      assert(N==M);
    }

    template<class T> inline void Fill(Matrix<T>& A,T & val)
    { 
      int N,M;
      Size(A,N,M);
      for(int i=0;i<N;i++){
      for(int j=0;j<M;j++){
	A[i][j]=val;
      }}
    }
    template<class T> inline void Diagonal(Matrix<T>& A,T & val)
    { 
      int N;
      SizeSquare(A,N);
      for(int i=0;i<N;i++){
	A[i][i]=val;
      }
    }
    template<class T> inline void Identity(Matrix<T>& A)
    {
      Fill(A,0.0);
      Diagonal(A,1.0);
    }

  };
}
#endif
