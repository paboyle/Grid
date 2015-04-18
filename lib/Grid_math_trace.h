#ifndef GRID_MATH_TRACE_H
#define GRID_MATH_TRACE_H
namespace Grid {
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

}
#endif
