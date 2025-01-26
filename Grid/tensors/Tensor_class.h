/*************************************************************************************
Grid physics library, www.github.com/paboyle/Grid
Source file: ./lib/tensors/Tensor_class.h
Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Michael Marshall <michael.marshall@ed.ac.au>
Author: Christoph Lehner <christoph@lhnr.de>

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
See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
			   /*  END LEGAL */
#pragma once 

NAMESPACE_BEGIN(Grid);

///////////////////////////////////////////////////
// Scalar, Vector, Matrix objects.
// These can be composed to form tensor products of internal indices.
///////////////////////////////////////////////////

// It is useful to NOT have any constructors
// so that these classes assert "is_pod<class> == true"
// because then the standard C++ valarray container eliminates fill overhead on
// new allocation and
// non-move copying.
//
// However note that doing this eliminates some syntactical sugar such as
// calling the constructor explicitly or implicitly
//
class GridTensorBase {};

// Too late to remove these traits from Grid Tensors, so inherit from GridTypeMapper
#define GridVector_CopyTraits \
  using element = vtype; \
  using scalar_type     = typename Traits::scalar_type; \
  using vector_type     = typename Traits::vector_type; \
  using scalar_typeD    = typename Traits::scalar_typeD; \
  using vector_typeD    = typename Traits::vector_typeD; \
  using tensor_reduced  = typename Traits::tensor_reduced; \
  using scalar_object   = typename Traits::scalar_object; \
  using scalar_objectD  = typename Traits::scalar_objectD; \
  using Complexified    = typename Traits::Complexified; \
  using Realified       = typename Traits::Realified; \
  using DoublePrecision = typename Traits::DoublePrecision; \
  using DoublePrecision2= typename Traits::DoublePrecision2; \
  static constexpr int TensorLevel = Traits::TensorLevel

///////////////////////////////////////////////////////////
// Allows to turn scalar<scalar<scalar<double>>>> back to double.
///////////////////////////////////////////////////////////
template <class T>
accelerator_inline typename std::enable_if<!isGridTensor<T>::value, T>::type
TensorRemove(T arg) {
  return arg;
}
template <class vtype>
accelerator_inline auto TensorRemove(iScalar<vtype> arg)
  -> decltype(TensorRemove(arg._internal)) {
  return TensorRemove(arg._internal);
}

template <class vtype>
class iScalar {
public:
  vtype _internal;

  using Traits = GridTypeMapper<iScalar<vtype> >;
  GridVector_CopyTraits;

  static accelerator_inline constexpr int Nsimd(void) { return sizeof(vector_type)/sizeof(scalar_type); } 

  // Scalar no action
  accelerator iScalar() = default;

  friend accelerator_inline void zeroit(iScalar<vtype> &that){
    zeroit(that._internal);
  }

  accelerator_inline iScalar(scalar_type s) : _internal(s){};  // recurse down and hit the constructor for vector_type

  accelerator_inline iScalar(const Zero &z) { zeroit(*this); };

  accelerator_inline iScalar<vtype> &operator=(const Zero &hero) {
    zeroit(*this);  return *this;
  }
  friend accelerator_inline void vstream(iScalar<vtype> &out, const iScalar<vtype> &in) {
    vstream(out._internal, in._internal);
  }
  friend accelerator_inline void vbroadcast(iScalar<vtype> &out,const iScalar<vtype> &in,int lane){
    vbroadcast(out._internal,in._internal,lane);
  }
  friend accelerator_inline void prefetch(iScalar<vtype> &that) {
    prefetch(that._internal);
  }
  friend accelerator_inline void permute(iScalar<vtype> &out, const iScalar<vtype> &in, int permutetype) {
    permute(out._internal, in._internal, permutetype);
  }
  friend accelerator_inline void rotate(iScalar<vtype> &out,const iScalar<vtype> &in,int rot){
    rotate(out._internal,in._internal,rot);
  }
  friend accelerator_inline void exchange(iScalar<vtype> &out1,iScalar<vtype> &out2,
				     const iScalar<vtype> &in1,const iScalar<vtype> &in2,int type)
  {
    exchange(out1._internal,out2._internal,in1._internal, in2._internal,type);
  }

  // Unary negation
  friend accelerator_inline iScalar<vtype> operator-(const iScalar<vtype> &r) {
    iScalar<vtype> ret;
    ret._internal = -r._internal;
    return ret;
  }
  // *=,+=,-= operators inherit from corresponding "*,-,+" behaviour
  accelerator_inline iScalar<vtype> &operator*=(const iScalar<vtype> &r) {
    *this = (*this) * r;
    return *this;
  }
  accelerator_inline iScalar<vtype> &operator-=(const iScalar<vtype> &r) {
    *this = (*this) - r;
    return *this;
  }
  accelerator_inline iScalar<vtype> &operator+=(const iScalar<vtype> &r) {
    *this = (*this) + r;
    return *this;
  }
  accelerator_inline vtype &operator()(void) { return _internal; }
  accelerator_inline const vtype &operator()(void) const { return _internal; }

  // Type casts meta programmed, must be pure scalar to match TensorRemove
  template <class U = vtype, class V = scalar_type, IfComplex<V> = 0, IfNotSimd<U> = 0> accelerator_inline
  operator ComplexF() const {
    return (TensorRemove(_internal));
  }
  template <class U = vtype, class V = scalar_type, IfComplex<V> = 0, IfNotSimd<U> = 0> accelerator_inline
  operator ComplexD() const {
    return (TensorRemove(_internal));
  }
  //             instantiation of "Grid::iScalar<vtype>::operator Grid::RealD() const [with vtype=Grid::Real, U=Grid::Real, V=Grid::RealD, <unnamed>=0, <unnamed>=0U]" 
  template <class U = vtype, class V = scalar_type, IfReal<V> = 0,IfNotSimd<U> = 0> accelerator_inline
  operator RealD() const {
    return (RealD) TensorRemove(_internal);
  }
  template <class U = vtype, class V = scalar_type, IfInteger<V> = 0, IfNotSimd<U> = 0> accelerator_inline
  operator Integer() const {
    return Integer(TensorRemove(_internal));
  }

  // convert from a something to a scalar via constructor of something arg
  template <class T, typename std::enable_if<!isGridTensor<T>::value, T>::type * = nullptr>
  accelerator_inline iScalar<vtype> operator=(T arg) {
    _internal = arg;
    return *this;
  }

  // Convert elements
  template <class ttype>
  accelerator_inline iScalar<vtype> operator=(const iScalar<ttype> &arg) {
    _internal = arg._internal;
    return *this;
  }

  // Host only
  friend std::ostream &operator<<(std::ostream &stream,const iScalar<vtype> &o) {
    stream << "S {" << o._internal << "}";
    return stream;
  };
  // FIXME These will break with change of data layout
  strong_inline const scalar_type * begin() const { return reinterpret_cast<const scalar_type *>(&_internal); }
  strong_inline       scalar_type * begin()       { return reinterpret_cast<      scalar_type *>(&_internal); }
  strong_inline const scalar_type * end()   const { return begin() + Traits::count; }
  strong_inline       scalar_type * end()         { return begin() + Traits::count; }
};

template <class vtype, int N>
class iVector {
public:
  vtype _internal[N];
  
  using Traits = GridTypeMapper<iVector<vtype, N> >;

  GridVector_CopyTraits;

  static accelerator_inline constexpr int Nsimd(void) { return sizeof(vector_type)/sizeof(scalar_type); } 

  template <class T, typename std::enable_if<!isGridTensor<T>::value, T>::type * = nullptr>
  accelerator_inline auto operator=(T arg) -> iVector<vtype, N> {
    zeroit(*this);
    for (int i = 0; i < N; i++) _internal[i] = arg;
    return *this;
  }

  accelerator iVector() = default;
  accelerator_inline iVector(const Zero &z) { zeroit(*this); };

  template<class other>
  accelerator_inline iVector<vtype, N> &operator=(const iVector<other,N> &him)
  {
    for (int i = 0; i < N; i++) {
      _internal[i] = him._internal[i];
    }    
    return *this;
  } 
  accelerator_inline iVector<vtype, N> &operator=(const Zero &hero) {
    zeroit(*this);
    return *this;
  }
  friend accelerator_inline void zeroit(iVector<vtype, N> &that) {
    for (int i = 0; i < N; i++) {
      zeroit(that._internal[i]);
    }
  }
  friend accelerator_inline void prefetch(iVector<vtype, N> &that) {
    for (int i = 0; i < N; i++) prefetch(that._internal[i]);
  }
  friend accelerator_inline void vstream(iVector<vtype, N> &out, const iVector<vtype, N> &in) {
    for (int i = 0; i < N; i++) {
      vstream(out._internal[i], in._internal[i]);
    }
  }
  friend accelerator_inline void vbroadcast(iVector<vtype,N> &out,const iVector<vtype,N> &in,int lane){
    for(int i=0;i<N;i++){
      vbroadcast(out._internal[i],in._internal[i],lane);
    }
  }
  friend accelerator_inline void permute(iVector<vtype,N> &out,const iVector<vtype,N> &in,int permutetype){
    for(int i=0;i<N;i++){
      permute(out._internal[i],in._internal[i],permutetype);
    }
  }
  friend accelerator_inline void rotate(iVector<vtype,N> &out,const iVector<vtype,N> &in,int rot){
    for(int i=0;i<N;i++){
      rotate(out._internal[i],in._internal[i],rot);
    }
  }
  friend accelerator_inline void exchange(iVector<vtype,N> &out1,iVector<vtype,N> &out2,
				     const iVector<vtype,N> &in1,const iVector<vtype,N> &in2,int type){
    for(int i=0;i<N;i++){
      exchange(out1._internal[i],out2._internal[i],in1._internal[i], in2._internal[i],type);
    }
  }

  // Unary negation
  friend accelerator_inline iVector<vtype, N> operator-(const iVector<vtype, N> &r) {
    iVector<vtype, N> ret;
    for (int i = 0; i < N; i++) ret._internal[i] = -r._internal[i];
    return ret;
  }
  // *=,+=,-= operators inherit from corresponding "*,-,+" behaviour
  accelerator_inline iVector<vtype, N> &operator*=(const iScalar<vtype> &r) {
    *this = (*this) * r;
    return *this;
  }
  accelerator_inline iVector<vtype, N> &operator-=(const iVector<vtype, N> &r) {
    *this = (*this) - r;
    return *this;
  }
  accelerator_inline iVector<vtype, N> &operator+=(const iVector<vtype, N> &r) {
    *this = (*this) + r;
    return *this;
  }
  accelerator_inline vtype &operator()(int i) { return _internal[i]; }
  accelerator_inline const vtype &operator()(int i) const { return _internal[i]; }

  // Host
  friend std::ostream &operator<<(std::ostream &stream, const iVector<vtype, N> &o) {
    stream << "V<" << N << ">{";
    for (int i = 0; i < N; i++) {
      stream << o._internal[i];
      if (i < N - 1) stream << ",";
    }
    stream << "}";
    return stream;
  };
  //    strong_inline vtype && operator ()(int i) {
  //      return _internal[i];
  //    }

  // FIXME These will break with change of data layout
  strong_inline const scalar_type * begin() const { return reinterpret_cast<const scalar_type *>(_internal); }
  strong_inline       scalar_type * begin()       { return reinterpret_cast<      scalar_type *>(_internal); }
  strong_inline const scalar_type * end()   const { return begin() + Traits::count; }
  strong_inline       scalar_type * end()         { return begin() + Traits::count; }

};

template <class vtype, int N>
class iMatrix {
public:
  vtype _internal[N][N];

  using Traits = GridTypeMapper<iMatrix<vtype, N> >;

  GridVector_CopyTraits;

  static accelerator_inline constexpr int Nsimd(void) { return sizeof(vector_type)/sizeof(scalar_type); } 

  accelerator_inline iMatrix(const Zero &z) { zeroit(*this); };
  accelerator iMatrix() = default;

  // Allow for type conversion.
  template<class other>
  accelerator_inline iMatrix &operator=(const iMatrix<other,N> &rhs) {
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) 
	_internal[i][j] = rhs._internal[i][j];
    return *this;
  };

  accelerator_inline iMatrix(scalar_type s) {
    (*this) = s;
  };  // recurse down and hit the constructor for vector_type

  accelerator_inline iMatrix<vtype, N> &operator=(const Zero &hero) {
    zeroit(*this);
    return *this;
  }
  template <class T, typename std::enable_if<!isGridTensor<T>::value, T>::type * = nullptr>
  accelerator_inline auto operator=(T arg) -> iMatrix<vtype, N> {
    zeroit(*this);
    for (int i = 0; i < N; i++) _internal[i][i] = arg;
    return *this;
  }

  friend accelerator_inline void zeroit(iMatrix<vtype,N> &that){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	zeroit(that._internal[i][j]);
    }}
  }
  friend accelerator_inline void prefetch(iMatrix<vtype,N> &that){
    for(int i=0;i<N;i++) {
      for(int j=0;j<N;j++) { 
	prefetch(that._internal[i][j]);
    }}
  }
  friend accelerator_inline void vstream(iMatrix<vtype,N> &out,const iMatrix<vtype,N> &in){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	vstream(out._internal[i][j],in._internal[i][j]);
    }}
  }
  friend accelerator_inline void vbroadcast(iMatrix<vtype,N> &out,const iMatrix<vtype,N> &in,int lane){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	vbroadcast(out._internal[i][j],in._internal[i][j],lane);
    }}
  }

  friend accelerator_inline void permute(iMatrix<vtype,N> &out,const iMatrix<vtype,N> &in,int permutetype){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	permute(out._internal[i][j],in._internal[i][j],permutetype);
    }}
  }
  friend accelerator_inline void rotate(iMatrix<vtype,N> &out,const iMatrix<vtype,N> &in,int rot){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
      rotate(out._internal[i][j],in._internal[i][j],rot);
    }}
  }
  friend accelerator_inline void exchange(iMatrix<vtype,N> &out1,iMatrix<vtype,N> &out2,
					  const iMatrix<vtype,N> &in1,const iMatrix<vtype,N> &in2,int type){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	exchange(out1._internal[i][j],out2._internal[i][j],in1._internal[i][j], in2._internal[i][j],type);
    }}
  }
  
  // Unary negation
  friend accelerator_inline iMatrix<vtype, N> operator-(const iMatrix<vtype, N> &r) {
    iMatrix<vtype, N> ret;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
	ret._internal[i][j] = -r._internal[i][j];
    }}
    return ret;
  }
  // *=,+=,-= operators inherit from corresponding "*,-,+" behaviour
  template <class T>
  accelerator_inline iMatrix<vtype, N> &operator*=(const T &r) {
    *this = (*this) * r;
    return *this;
  }
  template <class T>
  accelerator_inline iMatrix<vtype, N> &operator-=(const T &r) {
    *this = (*this) - r;
    return *this;
  }
  template <class T>
  accelerator_inline iMatrix<vtype, N> &operator+=(const T &r) {
    *this = (*this) + r;
    return *this;
  }

  // returns an lvalue reference
  accelerator_inline vtype &operator()(int i, int j) { return _internal[i][j]; }
  accelerator_inline const vtype &operator()(int i, int j) const {
    return _internal[i][j];
  }
  
  // Host function only
  friend std::ostream &operator<<(std::ostream &stream, const iMatrix<vtype, N> &o) {
    stream << "M<" << N << ">{";
    for (int i = 0; i < N; i++) {
      stream << "{";
      for (int j = 0; j < N; j++) {
	stream << o._internal[i][j];
	if (j < N - 1) stream << ",";
      }
      stream << "}";
      if (i != N - 1) stream << "\n\t\t";
    }
    stream << "}";
    return stream;
  };

  //  strong_inline vtype && operator ()(int i,int j) {
  //    return _internal[i][j];
  //  }

  // FIXME These will break with change of data layout
  strong_inline const scalar_type * begin() const { return reinterpret_cast<const scalar_type *>(_internal[0]); }
  strong_inline       scalar_type * begin()       { return reinterpret_cast<      scalar_type *>(_internal[0]); }
  strong_inline const scalar_type * end()   const { return begin() + Traits::count; }
  strong_inline       scalar_type * end()         { return begin() + Traits::count; }
};

template <class v> accelerator_inline
void vprefetch(const iScalar<v> &vv) {
  vprefetch(vv._internal);
}
template <class v, int N> accelerator_inline
void vprefetch(const iVector<v, N> &vv) {
  for (int i = 0; i < N; i++) {
    vprefetch(vv._internal[i]);
  }
}
template <class v, int N> accelerator_inline
void vprefetch(const iMatrix<v, N> &vv) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      vprefetch(vv._internal[i][j]);
    }
  }
}

NAMESPACE_END(Grid);


#ifdef GRID_SYCL
template<class vec> struct sycl::is_device_copyable<Grid::iScalar<vec> > : public std::true_type {};
template<class vec,int N> struct sycl::is_device_copyable<Grid::iVector<vec,N> > : public std::true_type {};
template<class vec,int N> struct sycl::is_device_copyable<Grid::iMatrix<vec,N> > : public std::true_type {};
#endif
