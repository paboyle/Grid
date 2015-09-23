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
