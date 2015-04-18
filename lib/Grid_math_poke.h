#ifndef GRID_MATH_POKE_H
#define GRID_MATH_POKE_H
namespace Grid {

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


}
#endif
