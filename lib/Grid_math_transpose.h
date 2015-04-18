#ifndef GRID_MATH_TRANSPOSE_H
#define GRID_MATH_TRANSPOSE_H
namespace Grid {



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

}
#endif
