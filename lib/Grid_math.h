#ifndef GRID_MATH_H
#define GRID_MATH_H

#include <Grid_math_traits.h>
#include <Grid_math_tensors.h>
#include <Grid_math_arith.h>

//
// Indexing; want to be able to dereference and 
// obtain either an lvalue or an rvalue.
//
namespace Grid {




    ///////////////////////////////////////////////////////////////////////////////////////
    // innerProduct Scalar x Scalar -> Scalar
    // innerProduct Vector x Vector -> Scalar
    // innerProduct Matrix x Matrix -> Scalar
    ///////////////////////////////////////////////////////////////////////////////////////
    template<class l,class r,int N> inline
    auto innerProduct (const iVector<l,N>& lhs,const iVector<r,N>& rhs) -> iScalar<decltype(innerProduct(lhs._internal[0],rhs._internal[0]))>
    {
        typedef decltype(innerProduct(lhs._internal[0],rhs._internal[0])) ret_t;
        iScalar<ret_t> ret=zero;
        for(int c1=0;c1<N;c1++){
            ret._internal += innerProduct(lhs._internal[c1],rhs._internal[c1]);
        }
        return ret;
    }
    template<class l,class r,int N> inline
    auto innerProduct (const iMatrix<l,N>& lhs,const iMatrix<r,N>& rhs) -> iScalar<decltype(innerProduct(lhs._internal[0][0],rhs._internal[0][0]))>
    {
        typedef decltype(innerProduct(lhs._internal[0][0],rhs._internal[0][0])) ret_t;
        iScalar<ret_t> ret=zero;
        iScalar<ret_t> tmp;
        for(int c1=0;c1<N;c1++){
        for(int c2=0;c2<N;c2++){
	  ret._internal+=innerProduct(lhs._internal[c1][c2],rhs._internal[c1][c2]);
        }}
        return ret;
    }
    template<class l,class r> inline
    auto innerProduct (const iScalar<l>& lhs,const iScalar<r>& rhs) -> iScalar<decltype(innerProduct(lhs._internal,rhs._internal))>
    {
        typedef decltype(innerProduct(lhs._internal,rhs._internal)) ret_t;
        iScalar<ret_t> ret;
        ret._internal = innerProduct(lhs._internal,rhs._internal);
        return ret;
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    // outerProduct Scalar x Scalar -> Scalar
    //              Vector x Vector -> Matrix
    ///////////////////////////////////////////////////////////////////////////////////////

template<class l,class r,int N> inline
auto outerProduct (const iVector<l,N>& lhs,const iVector<r,N>& rhs) -> iMatrix<decltype(outerProduct(lhs._internal[0],rhs._internal[0])),N>
{
    typedef decltype(outerProduct(lhs._internal[0],rhs._internal[0])) ret_t;
    iMatrix<ret_t,N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret._internal[c1][c2] = outerProduct(lhs._internal[c1],rhs._internal[c2]);
    }}
    return ret;
}
template<class l,class r> inline
auto outerProduct (const iScalar<l>& lhs,const iScalar<r>& rhs) -> iScalar<decltype(outerProduct(lhs._internal,rhs._internal))>
{
    typedef decltype(outerProduct(lhs._internal,rhs._internal)) ret_t;
    iScalar<ret_t> ret;
    ret._internal = outerProduct(lhs._internal,rhs._internal);
    return ret;
}

inline ComplexF outerProduct(const ComplexF &l, const ComplexF& r)
{
  return l*r;
}
inline ComplexD outerProduct(const ComplexD &l, const ComplexD& r)
{
  return l*r;
}
inline RealF outerProduct(const RealF &l, const RealF& r)
{
  return l*r;
}
inline RealD outerProduct(const RealD &l, const RealD& r)
{
  return l*r;
}
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////// CONJ         ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
 
// Conj function for scalar, vector, matrix
template<class vtype> inline iScalar<vtype> conj(const iScalar<vtype>&r)
{
    iScalar<vtype> ret;
    ret._internal = conj(r._internal);
    return ret;
}
template<class vtype,int N> inline iVector<vtype,N> conj(const iVector<vtype,N>&r)
{
  iVector<vtype,N> ret;
  for(int i=0;i<N;i++){
    ret._internal[i] = conj(r._internal[i]);
  }
  return ret;
}
template<class vtype,int N> inline iMatrix<vtype,N> conj(const iMatrix<vtype,N>&r)
{
  iMatrix<vtype,N> ret;
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    ret._internal[i][j] = conj(r._internal[i][j]);
  }}
  return ret;
}

// Adj function for scalar, vector, matrix
template<class vtype> inline iScalar<vtype> adj(const iScalar<vtype>&r)
{
    iScalar<vtype> ret;
    ret._internal = adj(r._internal);
    return ret;
}
template<class vtype,int N> inline iVector<vtype,N> adj(const iVector<vtype,N>&r)
{
    iVector<vtype,N> ret;
    for(int i=0;i<N;i++){
        ret._internal[i] = adj(r._internal[i]);
    }
    return ret;
}
template<class vtype,int N> inline iMatrix<vtype,N> adj(const iMatrix<vtype,N> &arg)
{
    iMatrix<vtype,N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret._internal[c1][c2]=adj(arg._internal[c2][c1]);
    }}
    return ret;
}



/////////////////////////////////////////////////////////////////
// Transpose all indices
/////////////////////////////////////////////////////////////////

inline ComplexD transpose(ComplexD &rhs){  return rhs;}
inline ComplexF transpose(ComplexF &rhs){  return rhs;}
inline RealD transpose(RealD &rhs){  return rhs;}
inline RealF transpose(RealF &rhs){  return rhs;}

template<class vtype,int N>
  inline typename std::enable_if<isGridTensor<vtype>::value, iMatrix<vtype,N> >::type 
  transpose(iMatrix<vtype,N> arg)
  {
    iMatrix<vtype,N> ret;
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	ret._internal[i][j] = transpose(arg._internal[j][i]); // NB recurses
      }}
    return ret;
  }
template<class vtype,int N>
  inline typename std::enable_if<isGridTensor<vtype>::notvalue, iMatrix<vtype,N> >::type 
  transpose(iMatrix<vtype,N> arg)
  {
    iMatrix<vtype,N> ret;
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	ret._internal[i][j] = arg._internal[j][i]; // Stop recursion if not a tensor type
      }}
    return ret;
  }

template<class vtype>
  inline typename std::enable_if<isGridTensor<vtype>::value, iScalar<vtype> >::type 
  transpose(iScalar<vtype> arg)
  {
    iScalar<vtype> ret;
    ret._internal = transpose(arg._internal); // NB recurses
    return ret;
  }

template<class vtype>
  inline typename std::enable_if<isGridTensor<vtype>::notvalue, iScalar<vtype> >::type 
  transpose(iScalar<vtype> arg)
  {
    iScalar<vtype> ret;
    ret._internal = arg._internal; // NB recursion stops
    return ret;
  }


////////////////////////////////////////////////////////////////////////////////////////////
// Transpose a specific index; instructive to compare this style of recursion termination
// to that of adj; which is easiers?
////////////////////////////////////////////////////////////////////////////////////////////
template<int Level,class vtype,int N> inline 
  typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,Level>::value, iMatrix<vtype,N> >::type 
transposeIndex (const iMatrix<vtype,N> &arg)
{
  iMatrix<vtype,N> ret;
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      ret._internal[i][j] = arg._internal[j][i]; 
  }}
  return ret;
}
// or not
template<int Level,class vtype,int N> inline 
typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,Level>::notvalue, iMatrix<vtype,N> >::type 
transposeIndex (const iMatrix<vtype,N> &arg)
{
  iMatrix<vtype,N> ret;
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      ret._internal[i][j] = transposeIndex<Level>(arg._internal[i][j]); 
  }}
  return ret;
}
template<int Level,class vtype> inline 
typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::notvalue, iScalar<vtype> >::type 
transposeIndex (const iScalar<vtype> &arg)
{
  iScalar<vtype> ret;
  ret._internal=transposeIndex<Level>(arg._internal);
  return ret;
}
template<int Level,class vtype> inline 
typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::value, iScalar<vtype> >::type 
transposeIndex (const iScalar<vtype> &arg)
{
  return arg;
}

//////////////////////////////////////////////////////////////////
// Traces: both all indices and a specific index 
/////////////////////////////////////////////////////////////////

inline ComplexF trace( const ComplexF &arg){    return arg;}
inline ComplexD trace( const ComplexD &arg){    return arg;}
inline RealF trace( const RealF &arg){    return arg;}
inline RealD trace( const RealD &arg){    return arg;}

template<int Level> inline ComplexF traceIndex(const ComplexF arg) { return arg;}
template<int Level> inline ComplexD traceIndex(const ComplexD arg) { return arg;}
template<int Level> inline RealF traceIndex(const RealF arg) { return arg;}
template<int Level> inline RealD traceIndex(const RealD arg) { return arg;}

template<class vtype,int N>
inline auto trace(const iMatrix<vtype,N> &arg) -> iScalar<decltype(trace(arg._internal[0][0]))>
{
    iScalar<decltype( trace(arg._internal[0][0] )) > ret;
    zeroit(ret._internal);
    for(int i=0;i<N;i++){
        ret._internal=ret._internal+trace(arg._internal[i][i]);
    }
    return ret;
}
template<class vtype>
inline auto trace(const iScalar<vtype> &arg) -> iScalar<decltype(trace(arg._internal))>
{
    iScalar<decltype(trace(arg._internal))> ret;
    ret._internal=trace(arg._internal);
    return ret;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Trace Specific indices.
////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
template<int Level,class vtype> inline 
auto traceIndex(const iScalar<vtype> &arg) -> iScalar<decltype(traceIndex<Level>(arg._internal)) >
{
  iScalar<decltype(traceIndex<Level>(arg._internal))> ret;
  ret._internal = traceIndex<Level>(arg._internal);
  return ret;
}
*/
template<int Level,class vtype> inline auto 
traceIndex (const iScalar<vtype> &arg) ->
typename 
std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::notvalue, 
  iScalar<decltype(traceIndex<Level>(arg._internal))> >::type 

{
  iScalar<decltype(traceIndex<Level>(arg._internal))> ret;
  ret._internal=traceIndex<Level>(arg._internal);
  return ret;
}
template<int Level,class vtype> inline auto
traceIndex (const iScalar<vtype> &arg) ->
typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::value, 
                        iScalar<vtype> >::type 
{
  return arg;
}

// If we hit the right index, return scalar and trace it with no further recursion
template<int Level,class vtype,int N> inline 
auto traceIndex(const iMatrix<vtype,N> &arg) ->
  typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::value,  // Index matches
                                                    iScalar<vtype> >::type                              // return scalar
{
  iScalar<vtype> ret;
  zeroit(ret._internal);
  for(int i=0;i<N;i++){
    ret._internal = ret._internal + arg._internal[i][i];
  }
  return ret;
}

// not this level, so recurse
template<int Level,class vtype,int N> inline 
auto traceIndex(const iMatrix<vtype,N> &arg) ->
  typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,Level>::notvalue,// No index match
         iMatrix<decltype(traceIndex<Level>(arg._internal[0][0])),N> >::type     // return matrix
{
  iMatrix<decltype(traceIndex<Level>(arg._internal[0][0])),N> ret;
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    ret._internal[i][j] = traceIndex<Level>(arg._internal[i][j]);
  }}
  return ret;
}

//////////////////////////////////////////////////////////////////////////////
// Peek on a specific index; returns a scalar in that index, tensor inherits rest
//////////////////////////////////////////////////////////////////////////////
// If we hit the right index, return scalar with no further recursion

//template<int Level> inline ComplexF peekIndex(const ComplexF arg) { return arg;}
//template<int Level> inline ComplexD peekIndex(const ComplexD arg) { return arg;}
//template<int Level> inline RealF peekIndex(const RealF arg) { return arg;}
//template<int Level> inline RealD peekIndex(const RealD arg) { return arg;}

// Scalar peek, no indices
template<int Level,class vtype> inline 
  auto peekIndex(const iScalar<vtype> &arg) -> 
  typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::value,  // Index matches
  iScalar<vtype> >::type                              // return scalar
{
  return arg;
}
// Vector peek, one index
template<int Level,class vtype,int N> inline 
  auto peekIndex(const iVector<vtype,N> &arg,int i) -> 
  typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,Level>::value,  // Index matches
  iScalar<vtype> >::type                              // return scalar
{
  iScalar<vtype> ret;                              // return scalar
  ret._internal = arg._internal[i];
  return ret;
}
// Matrix peek, two indices
template<int Level,class vtype,int N> inline 
  auto peekIndex(const iMatrix<vtype,N> &arg,int i,int j) -> 
  typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,Level>::value,  // Index matches
  iScalar<vtype> >::type                              // return scalar
{
  iScalar<vtype> ret;                              // return scalar
  ret._internal = arg._internal[i][j];
  return ret;
}

/////////////
// No match peek for scalar,vector,matrix must forward on either 0,1,2 args. Must have 9 routines with notvalue
/////////////
// scalar
template<int Level,class vtype> inline 
  auto peekIndex(const iScalar<vtype> &arg) ->                                     // Scalar 0 index  
  typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::notvalue,  // Index does NOT match
  iScalar<decltype(peekIndex<Level>(arg._internal))> >::type                       
{
  iScalar<decltype(peekIndex<Level>(arg._internal))> ret;
  ret._internal= peekIndex<Level>(arg._internal);
  return ret;
}
template<int Level,class vtype> inline 
  auto peekIndex(const iScalar<vtype> &arg,int i) ->                             // Scalar 1 index
  typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::notvalue,  // Index does NOT match
  iScalar<decltype(peekIndex<Level>(arg._internal,i))> >::type                       
{
  iScalar<decltype(peekIndex<Level>(arg._internal,i))> ret;
  ret._internal=peekIndex<Level>(arg._internal,i);
  return ret;
}
template<int Level,class vtype> inline 
  auto peekIndex(const iScalar<vtype> &arg,int i,int j) ->                         // Scalar, 2 index
  typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::notvalue,  // Index does NOT match
  iScalar<decltype(peekIndex<Level>(arg._internal,i,j))> >::type                       
{
  iScalar<decltype(peekIndex<Level>(arg._internal,i,j))> ret;
  ret._internal=peekIndex<Level>(arg._internal,i,j);
  return ret;
}
// vector
template<int Level,class vtype,int N> inline 
auto peekIndex(const iVector<vtype,N> &arg) -> 
  typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::notvalue,  // Index does not match
  iVector<decltype(peekIndex<Level>(arg._internal[0])),N> >::type                       
{
  iVector<decltype(peekIndex<Level>(arg._internal[0])),N> ret;
  for(int ii=0;ii<N;ii++){
    ret._internal[ii]=peekIndex<Level>(arg._internal[ii]);
  }
  return ret;
}
template<int Level,class vtype,int N> inline 
  auto peekIndex(const iVector<vtype,N> &arg,int i) -> 
  typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,Level>::notvalue,  // Index does not match
  iVector<decltype(peekIndex<Level>(arg._internal[0],i)),N> >::type                       
{
  iVector<decltype(peekIndex<Level>(arg._internal[0],i)),N> ret;
  for(int ii=0;ii<N;ii++){
    ret._internal[ii]=peekIndex<Level>(arg._internal[ii],i);
  }
  return ret;
}
template<int Level,class vtype,int N> inline 
  auto peekIndex(const iVector<vtype,N> &arg,int i,int j) -> 
  typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,Level>::notvalue,  // Index does not match
  iVector<decltype(peekIndex<Level>(arg._internal[0],i,j)),N> >::type                       
{
  iVector<decltype(peekIndex<Level>(arg._internal[0],i,j)),N> ret;
  for(int ii=0;ii<N;ii++){
    ret._internal[ii]=peekIndex<Level>(arg._internal[ii],i,j);
  }
  return ret;
}
// matrix
template<int Level,class vtype,int N> inline 
auto peekIndex(const iMatrix<vtype,N> &arg) -> 
  typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::notvalue,  // Index does not match
  iMatrix<decltype(peekIndex<Level>(arg._internal[0][0])),N> >::type                       
{
  iMatrix<decltype(peekIndex<Level>(arg._internal[0][0])),N> ret;
  for(int ii=0;ii<N;ii++){
  for(int jj=0;jj<N;jj++){
    ret._internal[ii][jj]=peekIndex<Level>(arg._internal[ii][jj]);// Could avoid this because peeking a scalar is dumb
  }}
  return ret;
}
template<int Level,class vtype,int N> inline 
  auto peekIndex(const iMatrix<vtype,N> &arg,int i) -> 
  typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,Level>::notvalue,  // Index does not match
  iMatrix<decltype(peekIndex<Level>(arg._internal[0],i)),N> >::type                       
{
  iMatrix<decltype(peekIndex<Level>(arg._internal[0],i)),N> ret;
  for(int ii=0;ii<N;ii++){
  for(int jj=0;jj<N;jj++){
    ret._internal[ii][jj]=peekIndex<Level>(arg._internal[ii][jj],i);
  }}
  return ret;
}
template<int Level,class vtype,int N> inline 
  auto peekIndex(const iMatrix<vtype,N> &arg,int i,int j) -> 
  typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,Level>::notvalue,  // Index does not match
  iMatrix<decltype(peekIndex<Level>(arg._internal[0][0],i,j)),N> >::type                       
{
  iMatrix<decltype(peekIndex<Level>(arg._internal[0][0],i,j)),N> ret;
  for(int ii=0;ii<N;ii++){
  for(int jj=0;jj<N;jj++){
    ret._internal[ii][jj]=peekIndex<Level>(arg._internal[ii][jj],i,j);
  }}
  return ret;
}


//////////////////////////////////////////////////////////////////////////////
// Poke a specific index; 
//////////////////////////////////////////////////////////////////////////////

// Scalar poke
template<int Level,class vtype> inline 
  void pokeIndex(iScalar<vtype> &ret, 
		 const typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::value,iScalar<vtype> >::type &arg)
{
  ret._internal = arg._internal;
}
// Vector poke, one index
template<int Level,class vtype,int N> inline 
  void pokeIndex(iVector<vtype,N> &ret, 
		 const typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,Level>::value,iScalar<vtype> >::type &arg,int i)
{
  ret._internal[i] = arg._internal;
}
// Vector poke, two indices
template<int Level,class vtype,int N> inline 
  void pokeIndex(iMatrix<vtype,N> &ret, 
		 const typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,Level>::value,iScalar<vtype> >::type &arg,int i,int j)
{
  ret._internal[i][j] = arg._internal;
}

/////////////
// No match poke for scalar,vector,matrix must forward on either 0,1,2 args. Must have 9 routines with notvalue
/////////////
// scalar
template<int Level,class vtype> inline 
  void pokeIndex(iScalar<vtype> &ret,
		 const typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::notvalue,iScalar<decltype(peekIndex<Level>(ret._internal))> >::type &arg)
		 
{
  pokeIndex<Level>(ret._internal,arg._internal);
}
template<int Level,class vtype> inline 
  void pokeIndex(iScalar<vtype> &ret,
		 const typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::notvalue,iScalar<decltype(peekIndex<Level>(ret._internal,0))> >::type &arg,
		 int i)
		 
{
  pokeIndex<Level>(ret._internal,arg._internal,i);
}
template<int Level,class vtype> inline 
  void pokeIndex(iScalar<vtype> &ret,
		 const typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,Level>::notvalue,iScalar<decltype(peekIndex<Level>(ret._internal,0,0))> >::type &arg,
		 int i,int j)
		 
{
  pokeIndex<Level>(ret._internal,arg._internal,i,j);
}

// Vector
template<int Level,class vtype,int N> inline 
  void pokeIndex(iVector<vtype,N> &ret,
		 const typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,Level>::notvalue,iVector<decltype(peekIndex<Level>(ret._internal)),N> >::type &arg)
		 
{
  for(int ii=0;ii<N;ii++){
    pokeIndex<Level>(ret._internal[ii],arg._internal[ii]);
  }
}
template<int Level,class vtype,int N> inline 
  void pokeIndex(iVector<vtype,N> &ret,
		 const typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,Level>::notvalue,iVector<decltype(peekIndex<Level>(ret._internal,0)),N> >::type &arg,
		 int i)
		 
{
  for(int ii=0;ii<N;ii++){
    pokeIndex<Level>(ret._internal[ii],arg._internal[ii],i);
  }
}
template<int Level,class vtype,int N> inline 
  void pokeIndex(iVector<vtype,N> &ret,
		 const typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,Level>::notvalue,iVector<decltype(peekIndex<Level>(ret._internal,0,0)),N> >::type &arg,
		 int i,int j)
		 
{
  for(int ii=0;ii<N;ii++){
    pokeIndex<Level>(ret._internal[ii],arg._internal[ii],i,j);
  }
}

// Matrix
template<int Level,class vtype,int N> inline 
  void pokeIndex(iMatrix<vtype,N> &ret,
		 const typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,Level>::notvalue,iMatrix<decltype(peekIndex<Level>(ret._internal)),N> >::type &arg)
		 
{
  for(int ii=0;ii<N;ii++){
  for(int jj=0;jj<N;jj++){
    pokeIndex<Level>(ret._internal[ii][jj],arg._internal[ii][jj]);
  }}
}
template<int Level,class vtype,int N> inline 
  void pokeIndex(iMatrix<vtype,N> &ret,
		 const typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,Level>::notvalue,iMatrix<decltype(peekIndex<Level>(ret._internal,0)),N> >::type &arg,
		 int i)
		 
{
  for(int ii=0;ii<N;ii++){
  for(int jj=0;jj<N;jj++){
    pokeIndex<Level>(ret._internal[ii][jj],arg._internal[ii][jj],i);
  }}
}
template<int Level,class vtype,int N> inline 
  void pokeIndex(iMatrix<vtype,N> &ret,
		 const typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,Level>::notvalue,iMatrix<decltype(peekIndex<Level>(ret._internal,0,0)),N> >::type &arg,
		 int i,int j)
		 
{
  for(int ii=0;ii<N;ii++){
  for(int jj=0;jj<N;jj++){
    pokeIndex<Level>(ret._internal[ii][jj],arg._internal[ii][jj],i,j);
  }}
}

/////////////////////////////////////////////////////////////////
// Can only take the real/imag part of scalar objects, since
// lattice objects of different complex nature are non-conformable.
/////////////////////////////////////////////////////////////////
template<class itype> inline auto real(const iScalar<itype> &z) -> iScalar<decltype(real(z._internal))>
{
    iScalar<decltype(real(z._internal))> ret;
    ret._internal = real(z._internal);
    return ret;
}
template<class itype,int N> inline auto real(const iMatrix<itype,N> &z) -> iMatrix<decltype(real(z._internal[0][0])),N>
{
    iMatrix<decltype(real(z._internal[0][0])),N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret._internal[c1][c2] = real(z._internal[c1][c2]);
    }}
    return ret;
}
template<class itype,int N> inline auto real(const iVector<itype,N> &z) -> iVector<decltype(real(z._internal[0])),N>
{
    iVector<decltype(real(z._internal[0])),N> ret;
    for(int c1=0;c1<N;c1++){
        ret._internal[c1] = real(z._internal[c1]);
    }
    return ret;
}
    
template<class itype> inline auto imag(const iScalar<itype> &z) -> iScalar<decltype(imag(z._internal))>
{
    iScalar<decltype(imag(z._internal))> ret;
    ret._internal = imag(z._internal);
    return ret;
}
template<class itype,int N> inline auto imag(const iMatrix<itype,N> &z) -> iMatrix<decltype(imag(z._internal[0][0])),N>
{
    iMatrix<decltype(imag(z._internal[0][0])),N> ret;
    for(int c1=0;c1<N;c1++){
    for(int c2=0;c2<N;c2++){
        ret._internal[c1][c2] = imag(z._internal[c1][c2]);
    }}
    return ret;
}
template<class itype,int N> inline auto imag(const iVector<itype,N> &z) -> iVector<decltype(imag(z._internal[0])),N>
{
    iVector<decltype(imag(z._internal[0])),N> ret;
    for(int c1=0;c1<N;c1++){
        ret._internal[c1] = imag(z._internal[c1]);
    }
    return ret;
}

};
    
#endif
