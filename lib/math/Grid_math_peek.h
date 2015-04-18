#ifndef GRID_MATH_PEEK_H
#define GRID_MATH_PEEK_H
namespace Grid {

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


}
#endif
