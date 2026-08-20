/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/Simd.h

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: paboyle <paboyle@ph.ed.ac.uk>

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

#include <type_traits>

#include <Grid/simd/Grid_scalar_support.h>

NAMESPACE_BEGIN(Grid);
//////////////////////////////////////
// To take the floating point type of real/complex type
//////////////////////////////////////
template <typename T>
struct RealPart {
  typedef T type;
};
// Grid's own complex alias, so this matches thrust::complex on CUDA/HIP as
// well as std::complex on the host. Qualifying it std:: left Real == the
// complex type itself on device builds, which removes the Grid_simd(Real)
// constructor and so any assignment of a real to a complex lattice.
template <typename T>
struct RealPart<complex<T> > {
  typedef T type;
};


// type alias used to simplify the syntax of std::enable_if
template <typename T> using Invoke = typename T::type;
template <typename Condition, typename ReturnType = void> using EnableIf    = Invoke<std::enable_if<Condition::value, ReturnType> >;
template <typename Condition, typename ReturnType = void> using NotEnableIf = Invoke<std::enable_if<!Condition::value, ReturnType> >;

////////////////////////////////////////////////////////
// Check for complexity with type traits
template <typename T> struct is_complex : public std::false_type {};
template <> struct is_complex<ComplexD> : public std::true_type {};
template <> struct is_complex<ComplexF> : public std::true_type {};

template <typename T> struct is_ComplexD : public std::false_type {};
template <> struct is_ComplexD<ComplexD> : public std::true_type {};

template <typename T> struct is_ComplexF : public std::false_type {};
template <> struct is_ComplexF<ComplexF> : public std::true_type {};

template<typename T, typename V=void> struct is_real : public std::false_type {};
template<typename T> struct is_real<T, typename std::enable_if<std::is_floating_point<T>::value,
  void>::type> : public std::true_type {};

template<typename T, typename V=void> struct is_integer : public std::false_type {};
template<typename T> struct is_integer<T, typename std::enable_if<std::is_integral<T>::value,
  void>::type> : public std::true_type {};

template <typename T>              using IfReal    = Invoke<std::enable_if<is_real<T>::value, int> >;
template <typename T>              using IfComplex = Invoke<std::enable_if<is_complex<T>::value, int> >;
template <typename T>              using IfInteger = Invoke<std::enable_if<is_integer<T>::value, int> >;
template <typename T1,typename T2> using IfSame    = Invoke<std::enable_if<std::is_same<T1,T2>::value, int> >;

template <typename T>              using IfNotReal    = Invoke<std::enable_if<!is_real<T>::value, int> >;
template <typename T>              using IfNotComplex = Invoke<std::enable_if<!is_complex<T>::value, int> >;
template <typename T>              using IfNotInteger = Invoke<std::enable_if<!is_integer<T>::value, int> >;
template <typename T1,typename T2> using IfNotSame    = Invoke<std::enable_if<!std::is_same<T1,T2>::value, int> >;

////////////////////////////////////////////////////////
// Define the operation templates functors
// general forms to allow for vsplat syntax
// need explicit declaration of types when used since
// clang cannot automatically determine the output type sometimes
template <class Out, class Input1, class Input2, class Input3, class Operation>
Out accelerator_inline  trinary(Input1 src_1, Input2 src_2, Input3 src_3, Operation op) {
  return op(src_1, src_2, src_3);
}
template <class Out, class Input1, class Input2, class Operation>
Out accelerator_inline binary(Input1 src_1, Input2 src_2, Operation op) {
  return op(src_1, src_2);
}
template <class Out, class Input, class Operation>
Out accelerator_inline  unary(Input src, Operation op) {
  return op(src);
}
NAMESPACE_END(Grid);
///////////////////////////////////////////////

#include <Grid/simd/Grid_vector_types.h>
#include <Grid/simd/Grid_doubled_vector.h>
#include <Grid/simd/Grid_scalar_types.h>

#ifdef GRID_SYCL
template<> struct sycl::is_device_copyable<Grid::vComplexF> : public std::true_type {};
template<> struct sycl::is_device_copyable<Grid::vComplexD> : public std::true_type {};
template<> struct sycl::is_device_copyable<Grid::vRealF   > : public std::true_type {};
template<> struct sycl::is_device_copyable<Grid::vRealD   > : public std::true_type {};
template<> struct sycl::is_device_copyable<Grid::vInteger > : public std::true_type {};
template<> struct sycl::is_device_copyable<Grid::sComplexF> : public std::true_type {};
template<> struct sycl::is_device_copyable<Grid::sComplexD> : public std::true_type {};
template<> struct sycl::is_device_copyable<Grid::sRealF   > : public std::true_type {};
template<> struct sycl::is_device_copyable<Grid::sRealD   > : public std::true_type {};
template<> struct sycl::is_device_copyable<Grid::sInteger > : public std::true_type {};
#endif

/////////////////////////////////////////
// Detect vector types
/////////////////////////////////////////
NAMESPACE_BEGIN(Grid);
template <typename T>
struct is_simd : public std::false_type {};
template <> struct is_simd<vRealF>     : public std::true_type {};
template <> struct is_simd<vRealD>     : public std::true_type {};
template <> struct is_simd<vRealH>     : public std::true_type {};
template <> struct is_simd<vComplexF>  : public std::true_type {};
template <> struct is_simd<vComplexD>  : public std::true_type {};
template <> struct is_simd<vComplexH>  : public std::true_type {};
template <> struct is_simd<vInteger>   : public std::true_type {};
template <> struct is_simd<sRealF>     : public std::true_type {};
template <> struct is_simd<sRealD>     : public std::true_type {};
template <> struct is_simd<sComplexF>  : public std::true_type {};
template <> struct is_simd<sComplexD>  : public std::true_type {};
template <> struct is_simd<sInteger>   : public std::true_type {};

template <typename T> using IfSimd    = Invoke<std::enable_if<is_simd<T>::value, int> >;
template <typename T> using IfNotSimd = Invoke<std::enable_if<!is_simd<T>::value, unsigned> >;


///////////////////////////////////////////////
// insert / extract with complex support
///////////////////////////////////////////////
template <class S, class V>
accelerator_inline S getlane(const Grid_simd<S, V> &in,int lane) {
  return in.getlane(lane);
}
template <class S, class V>
accelerator_inline void putlane(Grid_simd<S, V> &vec,const S &_S, int lane){
  vec.putlane(_S,lane);
}
template <class S,IfNotSimd<S> = 0 >
accelerator_inline S getlane(const S &in,int lane) {
  return in;
}
template <class S,IfNotSimd<S> = 0 >
accelerator_inline void putlane(S &vec,const S &_S, int lane){
  vec = _S;
}
template <class S, class V>
accelerator_inline S getlane(const Grid_simd2<S, V> &in,int lane) {
  return in.getlane(lane);
}
template <class S, class V>
accelerator_inline void putlane(Grid_simd2<S, V> &vec,const S &_S, int lane){
  vec.putlane(_S,lane);
}
template <class S>
accelerator_inline S getlane(const Grid_simd1<S> &in,int lane) {
  return in.getlane(lane);
}
template <class S>
accelerator_inline void putlane(Grid_simd1<S> &vec,const S &_S, int lane){
  vec.putlane(_S,lane);
}

NAMESPACE_END(Grid);



#include <Grid/simd/Grid_vector_unops.h>
#include <Grid/simd/Grid_scalar_unops.h>

NAMESPACE_BEGIN(Grid);

// Default precision is wired to double
typedef vRealD vReal;
typedef vComplexD vComplex;
 
inline std::ostream& operator<< (std::ostream& stream, const vComplexF &o){
  int nn=vComplexF::Nsimd();
  std::vector<ComplexF,alignedAllocator<ComplexF> > buf(nn);
  vstore(o,&buf[0]);
  stream<<"<";
  for(int i=0;i<nn;i++){
    stream<<buf[i];
    if(i<nn-1) stream<<",";
  }
  stream<<">";
  return stream;
}
 
inline std::ostream& operator<< (std::ostream& stream, const vComplexD &o){
  int nn=vComplexD::Nsimd();
  std::vector<ComplexD,alignedAllocator<ComplexD> > buf(nn);
  vstore(o,&buf[0]);
  stream<<"<";
  for(int i=0;i<nn;i++){
    stream<<buf[i];
    if(i<nn-1) stream<<",";
  }
  stream<<">";
  return stream;
}
inline std::ostream& operator<< (std::ostream& stream, const vComplexD2 &o){
  stream<<"<";
  stream<<o.v[0];
  stream<<o.v[1];
  stream<<">";
  return stream;
}

inline std::ostream& operator<< (std::ostream& stream, const vRealF &o){
  int nn=vRealF::Nsimd();
  std::vector<RealF,alignedAllocator<RealF> > buf(nn);
  vstore(o,&buf[0]);
  stream<<"<";
  for(int i=0;i<nn;i++){
    stream<<buf[i];
    if(i<nn-1) stream<<",";
  }
  stream<<">";
  return stream;
}

inline std::ostream& operator<< (std::ostream& stream, const vRealD &o){
  int nn=vRealD::Nsimd();
  std::vector<RealD,alignedAllocator<RealD> > buf(nn);
  vstore(o,&buf[0]);
  stream<<"<";
  for(int i=0;i<nn;i++){
    stream<<buf[i];
    if(i<nn-1) stream<<",";
  }
  stream<<">";
  return stream;
}
inline std::ostream& operator<< (std::ostream& stream, const vInteger &o){
  int nn=vInteger::Nsimd();
  std::vector<Integer,alignedAllocator<Integer> > buf(nn);
  vstore(o,&buf[0]);
  stream<<"<";
  for(int i=0;i<nn;i++){
    stream<<buf[i];
    if(i<nn-1) stream<<",";
  }
  stream<<">";
  return stream;
}

NAMESPACE_END(Grid)

