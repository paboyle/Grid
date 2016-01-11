    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_class.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
#ifndef GRID_MATH_TENSORS_H
#define GRID_MATH_TENSORS_H

namespace Grid {

///////////////////////////////////////////////////
// Scalar, Vector, Matrix objects.
// These can be composed to form tensor products of internal indices.
///////////////////////////////////////////////////

// It is useful to NOT have any constructors
// so that these classes assert "is_pod<class> == true"
// because then the standard C++ valarray container eliminates fill overhead on new allocation and 
// non-move copying.
//
// However note that doing this eliminates some syntactical sugar such as 
// calling the constructor explicitly or implicitly
//
class GridTensorBase {};

template<class vtype> class iScalar 
{
public:
  vtype _internal;

  typedef vtype element;
  typedef typename GridTypeMapper<vtype>::scalar_type scalar_type;
  typedef typename GridTypeMapper<vtype>::vector_type vector_type;
  typedef typename GridTypeMapper<vtype>::tensor_reduced tensor_reduced_v;
  typedef iScalar<tensor_reduced_v> tensor_reduced;
  typedef typename GridTypeMapper<vtype>::scalar_object recurse_scalar_object;
  typedef iScalar<recurse_scalar_object> scalar_object;

  // substitutes a real or complex version with same tensor structure
  typedef iScalar<typename GridTypeMapper<vtype>::Complexified > Complexified;
  typedef iScalar<typename GridTypeMapper<vtype>::Realified >    Realified;

  enum { TensorLevel = GridTypeMapper<vtype>::TensorLevel + 1};

  // Scalar no action
  //  template<int Level> using tensor_reduce_level = typename iScalar<GridTypeMapper<vtype>::tensor_reduce_level<Level> >;
  iScalar() = default;
  /*
  iScalar(const iScalar<vtype> &copyme)=default;
  iScalar(iScalar<vtype> &&copyme)=default;
  iScalar<vtype> & operator= (const iScalar<vtype> &copyme) = default;
  iScalar<vtype> & operator= (iScalar<vtype> &&copyme) = default;
  */
  iScalar(scalar_type s) : _internal(s) {};// recurse down and hit the constructor for vector_type
  iScalar(const Zero &z){ *this = zero; };

  iScalar<vtype> & operator= (const Zero &hero){
    zeroit(*this);
    return *this;
  }
  friend strong_inline void vstream(iScalar<vtype> &out,const iScalar<vtype> &in){
    vstream(out._internal,in._internal);
  }
  friend strong_inline void zeroit(iScalar<vtype> &that){
    zeroit(that._internal);
  }
  friend strong_inline void prefetch(iScalar<vtype> &that){
    prefetch(that._internal);
  }
  friend strong_inline void permute(iScalar<vtype> &out,const iScalar<vtype> &in,int permutetype){
    permute(out._internal,in._internal,permutetype);
  }

  // Unary negation
  friend strong_inline iScalar<vtype> operator -(const iScalar<vtype> &r) {
    iScalar<vtype> ret;
    ret._internal= -r._internal;
    return ret;
  }
  // *=,+=,-= operators inherit from corresponding "*,-,+" behaviour
  strong_inline iScalar<vtype> &operator *=(const iScalar<vtype> &r) {
    *this = (*this)*r;
    return *this;
  }
  strong_inline iScalar<vtype> &operator -=(const iScalar<vtype> &r) {
    *this = (*this)-r;
    return *this;
  }
  strong_inline iScalar<vtype> &operator +=(const iScalar<vtype> &r) {
    *this = (*this)+r;
    return *this;
  }
  strong_inline vtype & operator ()(void) {
    return _internal;
  }
  strong_inline const vtype & operator ()(void) const {
    return _internal;
  }

  // Type casts meta programmed, must be pure scalar to match TensorRemove
  template<class U=vtype,class V=scalar_type,IfComplex<V> = 0,IfNotSimd<U> = 0> operator ComplexF () const { return(TensorRemove(_internal)); };
  template<class U=vtype,class V=scalar_type,IfComplex<V> = 0,IfNotSimd<U> = 0> operator ComplexD () const { return(TensorRemove(_internal)); };
  //  template<class U=vtype,class V=scalar_type,IfComplex<V> = 0,IfNotSimd<U> = 0> operator RealD    () const { return(real(TensorRemove(_internal))); }
  template<class U=vtype,class V=scalar_type,IfReal<V>    = 0,IfNotSimd<U> = 0> operator RealD    () const { return TensorRemove(_internal); }
  template<class U=vtype,class V=scalar_type,IfInteger<V> = 0,IfNotSimd<U> = 0> operator Integer  () const { return Integer(TensorRemove(_internal)); }
  
  // convert from a something to a scalar via constructor of something arg
  template<class T,typename std::enable_if<!isGridTensor<T>::value, T>::type* = nullptr > strong_inline iScalar<vtype> operator = (T arg)
    { 
      _internal = arg;
      return *this;
    }

    friend std::ostream& operator<< (std::ostream& stream, const iScalar<vtype> &o){
      stream<< "S {"<<o._internal<<"}";
      return stream;
    };
};
///////////////////////////////////////////////////////////
// Allows to turn scalar<scalar<scalar<double>>>> back to double.
///////////////////////////////////////////////////////////
template<class T>     strong_inline typename std::enable_if<!isGridTensor<T>::value, T>::type TensorRemove(T arg) { return arg;}
template<class vtype> strong_inline auto TensorRemove(iScalar<vtype> arg) -> decltype(TensorRemove(arg._internal))
{
  return TensorRemove(arg._internal);
}
    
template<class vtype,int N> class iVector 
{
public:
  vtype _internal[N];

  typedef vtype element;
  typedef typename GridTypeMapper<vtype>::scalar_type scalar_type;
  typedef typename GridTypeMapper<vtype>::vector_type vector_type;
  typedef typename GridTypeMapper<vtype>::tensor_reduced tensor_reduced_v;
  typedef typename GridTypeMapper<vtype>::scalar_object recurse_scalar_object;
  typedef iScalar<tensor_reduced_v> tensor_reduced;
  typedef iVector<recurse_scalar_object,N> scalar_object;

  // substitutes a real or complex version with same tensor structure
  typedef iVector<typename GridTypeMapper<vtype>::Complexified,N > Complexified;
  typedef iVector<typename GridTypeMapper<vtype>::Realified,N >    Realified;

  template<class T,typename std::enable_if<!isGridTensor<T>::value, T>::type* = nullptr > strong_inline auto operator = (T arg) -> iVector<vtype,N>
    { 
      zeroit(*this);
      for(int i=0;i<N;i++)
	_internal[i] = arg;
      return *this;
    }

  enum { TensorLevel = GridTypeMapper<vtype>::TensorLevel + 1};
  iVector(const Zero &z){ *this = zero; };
  iVector() =default;
  /*
  iVector(const iVector<vtype,N> &copyme)=default;
  iVector(iVector<vtype,N> &&copyme)=default;
  iVector<vtype,N> & operator= (const iVector<vtype,N> &copyme) = default;
  iVector<vtype,N> & operator= (iVector<vtype,N> &&copyme) = default;
  */

  iVector<vtype,N> & operator= (const Zero &hero){
    zeroit(*this);
    return *this;
  }
  friend strong_inline void zeroit(iVector<vtype,N> &that){
    for(int i=0;i<N;i++){
      zeroit(that._internal[i]);
    }
  }
  friend strong_inline void prefetch(iVector<vtype,N> &that){
    for(int i=0;i<N;i++) prefetch(that._internal[i]);
  }
  friend strong_inline void vstream(iVector<vtype,N> &out,const iVector<vtype,N> &in){
    for(int i=0;i<N;i++){
      vstream(out._internal[i],in._internal[i]);
    }
  }
  friend strong_inline void permute(iVector<vtype,N> &out,const iVector<vtype,N> &in,int permutetype){
    for(int i=0;i<N;i++){
      permute(out._internal[i],in._internal[i],permutetype);
    }
  }

  // Unary negation
  friend strong_inline iVector<vtype,N> operator -(const iVector<vtype,N> &r) {
    iVector<vtype,N> ret;
    for(int i=0;i<N;i++) ret._internal[i]= -r._internal[i];
    return ret;
  }
  // *=,+=,-= operators inherit from corresponding "*,-,+" behaviour
  strong_inline iVector<vtype,N> &operator *=(const iScalar<vtype> &r) {
    *this = (*this)*r;
    return *this;
  }
  strong_inline iVector<vtype,N> &operator -=(const iVector<vtype,N> &r) {
    *this = (*this)-r;
    return *this;
  }
  strong_inline iVector<vtype,N> &operator +=(const iVector<vtype,N> &r) {
    *this = (*this)+r;
    return *this;
  }
  strong_inline vtype & operator ()(int i) {
    return _internal[i];
  }
  strong_inline const vtype & operator ()(int i) const {
    return _internal[i];
  }
  friend std::ostream& operator<< (std::ostream& stream, const iVector<vtype,N> &o){
    stream<< "V<"<<N<<">{";
    for(int i=0;i<N;i++) {
      stream<<o._internal[i];
      if (i<N-1)	stream<<",";
    }
    stream<<"}";
    return stream;
  };
  //    strong_inline vtype && operator ()(int i) {
  //      return _internal[i];
  //    }
};
    
template<class vtype,int N> class iMatrix 
{
public:
  vtype _internal[N][N];

  typedef vtype element;
  typedef typename GridTypeMapper<vtype>::scalar_type scalar_type;
  typedef typename GridTypeMapper<vtype>::vector_type vector_type;
  typedef typename GridTypeMapper<vtype>::tensor_reduced tensor_reduced_v;
  typedef typename GridTypeMapper<vtype>::scalar_object recurse_scalar_object;

  // substitutes a real or complex version with same tensor structure
  typedef iMatrix<typename GridTypeMapper<vtype>::Complexified,N > Complexified;
  typedef iMatrix<typename GridTypeMapper<vtype>::Realified,N >    Realified;

  // Tensure removal
  typedef iScalar<tensor_reduced_v> tensor_reduced;
  typedef iMatrix<recurse_scalar_object,N> scalar_object;

  enum { TensorLevel = GridTypeMapper<vtype>::TensorLevel + 1};


  iMatrix(const Zero &z){ *this = zero; };
  iMatrix() =default;
  
  iMatrix& operator=(const iMatrix& rhs){
    for(int i=0;i<N;i++)
      for(int j=0;j<N;j++)
	vstream(_internal[i][j],rhs._internal[i][j]);
    return *this;
  }; 
  
 

  iMatrix(scalar_type s)  { (*this) = s ;};// recurse down and hit the constructor for vector_type

  /*
  iMatrix(const iMatrix<vtype,N> &copyme)=default;
  iMatrix(iMatrix<vtype,N> &&copyme)=default;
  iMatrix<vtype,N> & operator= (const iMatrix<vtype,N> &copyme) = default;
  iMatrix<vtype,N> & operator= (iMatrix<vtype,N> &&copyme) = default;
  */



  iMatrix<vtype,N> & operator= (const Zero &hero){
    zeroit(*this);
    return *this;
  }
  template<class T,typename std::enable_if<!isGridTensor<T>::value, T>::type* = nullptr > strong_inline auto operator = (T arg) -> iMatrix<vtype,N>
    { 
      zeroit(*this);
      for(int i=0;i<N;i++)
	_internal[i][i] = arg;
      return *this;
    }

  friend strong_inline void zeroit(iMatrix<vtype,N> &that){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	zeroit(that._internal[i][j]);
    }}
  }
  friend strong_inline void prefetch(iMatrix<vtype,N> &that){
    for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      prefetch(that._internal[i][j]);
  }
  friend strong_inline void vstream(iMatrix<vtype,N> &out,const iMatrix<vtype,N> &in){
      for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	vstream(out._internal[i][j],in._internal[i][j]);
      }}
    }

  friend strong_inline void permute(iMatrix<vtype,N> &out,const iMatrix<vtype,N> &in,int permutetype){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	permute(out._internal[i][j],in._internal[i][j],permutetype);
    }}
  }


  // Unary negation
  friend strong_inline iMatrix<vtype,N> operator -(const iMatrix<vtype,N> &r) {
    iMatrix<vtype,N> ret;
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	ret._internal[i][j]= -r._internal[i][j];
    }}
    return ret;
  }
  // *=,+=,-= operators inherit from corresponding "*,-,+" behaviour
  template<class T>
  strong_inline iMatrix<vtype,N> &operator *=(const T &r) {
    *this = (*this)*r;
    return *this;
  }
  template<class T>
  strong_inline iMatrix<vtype,N> &operator -=(const T &r) {
    *this = (*this)-r;
    return *this;
  }
  template<class T>
  strong_inline iMatrix<vtype,N> &operator +=(const T &r) {
    *this = (*this)+r;
    return *this;
  }

  // returns an lvalue reference
  strong_inline vtype & operator ()(int i,int j) {
    return _internal[i][j];
  }
  strong_inline const vtype & operator ()(int i,int j) const {
    return _internal[i][j];
  }
  friend std::ostream& operator<< (std::ostream& stream, const iMatrix<vtype,N> &o){
    stream<< "M<"<<N<<">{";
    for(int i=0;i<N;i++) {
      stream<< "{";
      for(int j=0;j<N;j++) {
	stream<<o._internal[i][j];
	if (i<N-1)	stream<<",";
      }
      stream<<"}";
      if(i!=N-1) stream<<"\n\t\t";
    }
    stream<<"}";
    return stream;
  };

  //  strong_inline vtype && operator ()(int i,int j) {
  //    return _internal[i][j];
  //  }

};

template<class v> void vprefetch(const iScalar<v> &vv)
{
  vprefetch(vv._internal);
}
template<class v,int N> void vprefetch(const iVector<v,N> &vv)
{
  for(int i=0;i<N;i++){
    vprefetch(vv._internal[i]);
  }
}
template<class v,int N> void vprefetch(const iMatrix<v,N> &vv)
{
  for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
    vprefetch(vv._internal[i][j]);
  }}
}


}
#endif
