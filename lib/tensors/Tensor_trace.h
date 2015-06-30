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

/////////////////////////////////////
// Alternate implementation .. count 
// indices from left. Annoying but
// I'm working around compiler bugs in gcc and earlier icpc
/////////////////////////////////////

// Allow to recurse if vector, but never terminate on a vector
// ltrace of a different index can distribute across the vector index in a replicated way
// but we do not ltrace a vector index.
template<int Level> 
class TensorIndexRecursion {

 public:
  ////////////////////////////////////////////
  // Recursion for tracing a specific index
  ////////////////////////////////////////////
  template<class vtype>
  static auto traceIndex(const iScalar<vtype> arg) ->  iScalar<decltype(TensorIndexRecursion<Level-1>::traceIndex(arg._internal))> 
  {
    iScalar<decltype(TensorIndexRecursion<Level-1>::traceIndex(arg._internal))> ret;
    ret._internal = TensorIndexRecursion<Level-1>::traceIndex(arg._internal);
    return ret;
  }
  template<class vtype,int N>
  static auto traceIndex(const iVector<vtype,N> arg) ->  iVector<decltype(TensorIndexRecursion<Level-1>::traceIndex(arg._internal[0])),N> 
  {
    iVector<decltype(TensorIndexRecursion<Level-1>::traceIndex(arg._internal[0])),N> ret;
    for(int i=0;i<N;i++){
      ret._internal[i] = TensorIndexRecursion<Level-1>::traceIndex(arg._internal[i]);
    }
    return ret;
  }
  template<class vtype,int N>
  static auto traceIndex(const iMatrix<vtype,N> arg) ->  iMatrix<decltype(TensorIndexRecursion<Level-1>::traceIndex(arg._internal[0][0])),N> 
  {
    iMatrix<decltype(TensorIndexRecursion<Level-1>::traceIndex(arg._internal[0][0])),N> ret;
    for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      ret._internal[i][j] = TensorIndexRecursion<Level-1>::traceIndex(arg._internal[i][j]);
    }}
    return ret;
  }
  ////////////////////////////////////////////
  // Recursion for peeking a specific index
  ////////////////////////////////////////////
  template<class vtype>
  static auto peekIndex(const iScalar<vtype> arg,int i) ->  iScalar<decltype(TensorIndexRecursion<Level-1>::peekIndex(arg._internal,0))> 
  {
    iScalar<decltype(TensorIndexRecursion<Level-1>::peekIndex(arg._internal,0))> ret;
    ret._internal = TensorIndexRecursion<Level-1>::peekIndex(arg._internal,i);
    return ret;
  }
  template<class vtype>
  static auto peekIndex(const iScalar<vtype> arg,int i,int j) ->  iScalar<decltype(TensorIndexRecursion<Level-1>::peekIndex(arg._internal,0,0))> 
  {
    iScalar<decltype(TensorIndexRecursion<Level-1>::peekIndex(arg._internal,0,0))> ret;
    ret._internal = TensorIndexRecursion<Level-1>::peekIndex(arg._internal,i,j);
    return ret;
  }

  template<class vtype,int N>
  static auto peekIndex(const iVector<vtype,N> arg,int ii) ->  iVector<decltype(TensorIndexRecursion<Level-1>::peekIndex(arg._internal[0],0)),N> 
  {
    iVector<decltype(TensorIndexRecursion<Level-1>::peekIndex(arg._internal[0],0)),N> ret;
    for(int i=0;i<N;i++){
      ret._internal[i] = TensorIndexRecursion<Level-1>::peekIndex(arg._internal[i],ii);
    }
    return ret;
  }
  template<class vtype,int N>
  static auto peekIndex(const iVector<vtype,N> arg,int ii,int jj) ->  iVector<decltype(TensorIndexRecursion<Level-1>::peekIndex(arg._internal[0],0,0)),N> 
  {
    iVector<decltype(TensorIndexRecursion<Level-1>::peekIndex(arg._internal[0],0,0)),N> ret;
    for(int i=0;i<N;i++){
      ret._internal[i] = TensorIndexRecursion<Level-1>::peekIndex(arg._internal[i],ii,jj);
    }
    return ret;
  }

  template<class vtype,int N>
  static auto peekIndex(const iMatrix<vtype,N> arg,int ii) ->  iMatrix<decltype(TensorIndexRecursion<Level-1>::peekIndex(arg._internal[0][0],0)),N> 
  {
    iMatrix<decltype(TensorIndexRecursion<Level-1>::peekIndex(arg._internal[0][0],0)),N> ret;
    for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      ret._internal[i][j] = TensorIndexRecursion<Level-1>::peekIndex(arg._internal[i][j],ii);
    }}
    return ret;
  }
  template<class vtype,int N>
  static auto peekIndex(const iMatrix<vtype,N> arg,int ii,int jj) ->  iMatrix<decltype(TensorIndexRecursion<Level-1>::peekIndex(arg._internal[0][0],0,0)),N> 
  {
    iMatrix<decltype(TensorIndexRecursion<Level-1>::peekIndex(arg._internal[0][0],0,0)),N> ret;
    for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      ret._internal[i][j] = TensorIndexRecursion<Level-1>::peekIndex(arg._internal[i][j],ii,jj);
    }}
    return ret;
  }
  ////////////////////////////////////////////
  // Recursion for poking a specific index
  ////////////////////////////////////////////

  template<class vtype> inline static 
    void pokeIndex(iScalar<vtype> &ret, const iScalar<decltype(TensorIndexRecursion<Level-1>::peekIndex(ret._internal,0))> &arg, int i)
    {
      TensorIndexRecursion<Level-1>::pokeIndex(ret._internal,arg._internal,i);
    }
  template<class vtype> inline static 
    void pokeIndex(iScalar<vtype> &ret, const iScalar<decltype(TensorIndexRecursion<Level-1>::peekIndex(ret._internal,0,0))> &arg, int i,int j)
    {
      TensorIndexRecursion<Level-1>::pokeIndex(ret._internal,arg._internal,i,j);
    }

  template<class vtype,int N> inline static 
    void pokeIndex(iVector<vtype,N> &ret, const iVector<decltype(TensorIndexRecursion<Level-1>::peekIndex(ret._internal,0)),N> &arg, int i)
    {
      for(int ii=0;ii<N;ii++){
	TensorIndexRecursion<Level-1>::pokeIndex(ret._internal[ii],arg._internal[ii],i);
      }
    }
  template<class vtype,int N> inline static 
    void pokeIndex(iVector<vtype,N> &ret, const iVector<decltype(TensorIndexRecursion<Level-1>::peekIndex(ret._internal,0)),N> &arg, int i,int j)
    {
      for(int ii=0;ii<N;ii++){
	TensorIndexRecursion<Level-1>::pokeIndex(ret._internal[ii],arg._internal[ii],i,j);
      }
    }

  template<class vtype,int N> inline static 
    void pokeIndex(iMatrix<vtype,N> &ret, const iMatrix<decltype(TensorIndexRecursion<Level-1>::peekIndex(ret._internal,0)),N> &arg, int i)
    {
      for(int ii=0;ii<N;ii++){
      for(int jj=0;jj<N;jj++){
	TensorIndexRecursion<Level-1>::pokeIndex(ret._internal[ii][jj],arg._internal[ii][jj],i);
      }}
    }
  template<class vtype,int N> inline static 
    void pokeIndex(iMatrix<vtype,N> &ret, const iMatrix<decltype(TensorIndexRecursion<Level-1>::peekIndex(ret._internal,0)),N> &arg, int i,int j)
    {
      for(int ii=0;ii<N;ii++){
      for(int jj=0;jj<N;jj++){
	TensorIndexRecursion<Level-1>::pokeIndex(ret._internal[ii][jj],arg._internal[ii][jj],i,j);
      }}
    }

  ////////////////////////////////////////////
  // Recursion for transposing a specific index
  ////////////////////////////////////////////
  template<class vtype>
  static auto transposeIndex(const iScalar<vtype> arg) ->  iScalar<vtype> 
  {
    iScalar<vtype> ret;
    ret._internal = TensorIndexRecursion<Level-1>::transposeIndex(arg._internal);
    return ret;
  }
  template<class vtype,int N>
  static auto transposeIndex(const iVector<vtype,N> arg) ->  iVector<vtype,N> 
  {
    iVector<vtype,N> ret;
    for(int i=0;i<N;i++){
      ret._internal[i] = TensorIndexRecursion<Level-1>::transposeIndex(arg._internal[i]);
    }
    return ret;
  }
  template<class vtype,int N>
  static auto transposeIndex(const iMatrix<vtype,N> arg) ->  iMatrix<vtype,N> 
  {
    iMatrix<vtype,N> ret;
    for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      ret._internal[i][j] = TensorIndexRecursion<Level-1>::transposeIndex(arg._internal[i][j]);
    }}
    return ret;
  }
};

////////////////////////////
// strip const & ref quali's
////////////////////////////
#define RemoveCRV(a) typename std::remove_const<typename std::remove_reference<decltype(a)>::type>::type
template<>
class TensorIndexRecursion<0> {
 public:

  /////////////////////////////////////////
  // Ends recursion for trace (scalar/vector/matrix)
  /////////////////////////////////////////
  template<class vtype>
  static auto traceIndex(const iScalar<vtype> arg) ->  iScalar<RemoveCRV(arg._internal)>
  {
    iScalar<RemoveCRV(arg._internal)> ret;
    ret._internal = arg._internal;
    return ret;
  }
  template<class vtype,int N>
  static auto traceIndex(const iVector<vtype,N> arg) ->  iScalar<RemoveCRV(arg._internal[0])>
  {
    iScalar<RemoveCRV(arg._internal[0])> ret;
    ret._internal=zero;
    for(int i=0;i<N;i++){
      ret._internal = ret._internal+ arg._internal[i];
    }
    return ret;
  }
  template<class vtype,int N>
  static auto traceIndex(const iMatrix<vtype,N> arg) ->  iScalar<RemoveCRV(arg._internal[0][0])> 
  {
    iScalar<RemoveCRV(arg._internal[0][0])> ret;
    ret=zero;
    for(int i=0;i<N;i++){
      ret._internal = ret._internal+arg._internal[i][i];
    }
    return ret;
  }
  /////////////////////////////////////////
  // Ends recursion for transpose scalar/vector/matrix
  /////////////////////////////////////////
  template<class vtype>
  static auto transposeIndex(const iScalar<vtype> arg) ->  iScalar<vtype>
  {
    iScalar<vtype> ret;
    ret._internal = arg._internal;
    return ret;
  }
  template<class vtype,int N>
  static auto transposeIndex(const iVector<vtype,N> arg) ->  iVector<vtype,N>
  {
    iScalar<vtype> ret;
    ret._internal=zero;
    for(int i=0;i<N;i++){
      ret._internal = ret._internal+ arg._internal[i];
    }
    return ret;
  }
  template<class vtype,int N>
  static auto transposeIndex(const iMatrix<vtype,N> arg)  ->  iMatrix<vtype,N> 
  {
    iScalar<vtype> ret;
    ret=zero;
    for(int i=0;i<N;i++){
      ret._internal = ret._internal+arg._internal[i][i];
    }
    return ret;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // End recursion for peeking a specific index; single index on vector, double index on matrix
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class vtype,int N>
  static auto peekIndex(const iVector<vtype,N> arg,int ii) ->  iScalar<vtype> 
  {
    iScalar<vtype> ret;
    ret._internal = arg._internal[ii];
    return ret;
  }
  template<class vtype,int N>
  static auto peekIndex(const iMatrix<vtype,N> arg,int ii,int jj) ->  iScalar<vtype>
  {
    iScalar<vtype> ret;
    ret._internal = arg._internal[ii][jj];
    return ret;
  }
  // Vector poke, one index
  template<class vtype,int N> inline static 
    void pokeIndex(iVector<vtype,N> &ret, const iScalar<vtype> &arg,int i)
    {
      ret._internal[i] = arg._internal;
    }
  // Matrix poke two indices
  template<class vtype,int N> inline static 
    void pokeIndex(iMatrix<vtype,N> &ret, const iScalar<vtype> &arg,int i,int j)
    {
      ret._internal[i][j] = arg._internal;
    }

};

////////////////////////////////////////////////////////////////////////////////////////////////////////
// External wrappers
////////////////////////////////////////////////////////////////////////////////////////////////////////

template<int Level,class vtype> inline auto traceIndex (const vtype &arg) -> RemoveCRV(TensorIndexRecursion<Level>::traceIndex(arg))
{
  RemoveCRV(TensorIndexRecursion<Level>::traceIndex(arg)) ret;
  ret=TensorIndexRecursion<Level>::traceIndex(arg);
  return ret;
}
template<int Level,class vtype> inline auto transposeIndex (const vtype &arg) -> RemoveCRV(TensorIndexRecursion<Level>::transposeIndex(arg))
{
  RemoveCRV(TensorIndexRecursion<Level>::transposeIndex(arg)) ret;
  ret=TensorIndexRecursion<Level>::transposeIndex(arg);
  return ret;
}
template<int Level,class vtype> inline auto peekIndex (const vtype &arg,int i) -> RemoveCRV(TensorIndexRecursion<Level>::peekIndex(arg,0))
{
  RemoveCRV(TensorIndexRecursion<Level>::peekIndex(arg,0)) ret;
  ret=TensorIndexRecursion<Level>::peekIndex(arg,i);
  return ret;
}
template<int Level,class vtype> inline auto peekIndex (const vtype &arg,int i,int j) -> RemoveCRV(TensorIndexRecursion<Level>::peekIndex(arg,0,0))
{
  RemoveCRV(TensorIndexRecursion<Level>::peekIndex(arg,0,0)) ret;
  ret=TensorIndexRecursion<Level>::peekIndex(arg,i,j);
  return ret;
}

template<int Level,class vtype> inline 
void pokeIndex (vtype &ret,const decltype(TensorIndexRecursion<Level>::peekIndex(ret,0)) &arg,int i) 
{
  TensorIndexRecursion<Level>::pokeIndex(ret,arg,i);
}

template<int Level,class vtype> inline 
void pokeIndex (vtype &ret,const decltype(TensorIndexRecursion<Level>::peekIndex(ret,0,0)) &arg,int i,int j) 
{
  TensorIndexRecursion<Level>::pokeIndex(ret,arg,i,j);
}


#undef RemoveCRV

}
#endif
