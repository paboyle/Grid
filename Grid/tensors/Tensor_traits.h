   /*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid 
    Source file: ./lib/tensors/Tensor_traits.h
    Copyright (C) 2015
Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Christopher Kelly <ckelly@phys.columbia.edu>
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
#ifndef GRID_MATH_TRAITS_H
#define GRID_MATH_TRAITS_H

#include <type_traits>

namespace Grid {

//////////////////////////////////////////////////////////////////////////////////
// Want to recurse: GridTypeMapper<Matrix<vComplexD> >::scalar_type == ComplexD.
// Use of a helper class like this allows us to template specialise and "dress"
// other classes such as RealD == double, ComplexD == std::complex<double> with these
// traits.
//
// It is possible that we could do this more elegantly if I introduced a 
// queryable trait in iScalar, iMatrix and iVector and used the query on vtype in 
// place of the type mapper?
//
// Not sure how to do this, but probably could be done with a research effort
// to study C++11's type_traits.h file. (std::enable_if<isGridTensorType<vtype> >)
//
//////////////////////////////////////////////////////////////////////////////////
  
  template <class T> class GridTypeMapper {
  public:
    typedef typename T::scalar_type scalar_type;
    typedef typename T::vector_type vector_type;
    typedef typename T::vector_typeD vector_typeD;
    typedef typename T::tensor_reduced tensor_reduced;
    typedef typename T::scalar_object scalar_object;
    typedef typename T::Complexified Complexified;
    typedef typename T::Realified Realified;
    typedef typename T::DoublePrecision DoublePrecision;
    enum { TensorLevel = T::TensorLevel };
  };

//////////////////////////////////////////////////////////////////////////////////
// Recursion stops with these template specialisations
//////////////////////////////////////////////////////////////////////////////////
  template<> class GridTypeMapper<RealF> {
  public:
    typedef RealF scalar_type;
    typedef RealF vector_type;
    typedef RealD vector_typeD;
    typedef RealF tensor_reduced ;
    typedef RealF scalar_object;
    typedef ComplexF Complexified;
    typedef RealF Realified;
    typedef RealD DoublePrecision;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<RealD> {
  public:
    typedef RealD scalar_type;
    typedef RealD vector_type;
    typedef RealD vector_typeD;
    typedef RealD tensor_reduced;
    typedef RealD scalar_object;
    typedef ComplexD Complexified;
    typedef RealD Realified;
    typedef RealD DoublePrecision;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<ComplexF> {
  public:
    typedef ComplexF scalar_type;
    typedef ComplexF vector_type;
    typedef ComplexD vector_typeD;
    typedef ComplexF tensor_reduced;
    typedef ComplexF scalar_object;
    typedef ComplexF Complexified;
    typedef RealF Realified;
    typedef ComplexD DoublePrecision;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<ComplexD> {
  public:
    typedef ComplexD scalar_type;
    typedef ComplexD vector_type;
    typedef ComplexD vector_typeD;
    typedef ComplexD tensor_reduced;
    typedef ComplexD scalar_object;
    typedef ComplexD Complexified;
    typedef RealD Realified;
    typedef ComplexD DoublePrecision;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<Integer> {
  public:
    typedef Integer scalar_type;
    typedef Integer vector_type;
    typedef Integer vector_typeD;
    typedef Integer tensor_reduced;
    typedef Integer scalar_object;
    typedef void Complexified;
    typedef void Realified;
    typedef void DoublePrecision;
    enum { TensorLevel = 0 };
  };

  template<> class GridTypeMapper<vRealF> {
  public:
    typedef RealF  scalar_type;
    typedef vRealF vector_type;
    typedef vRealD vector_typeD;
    typedef vRealF tensor_reduced;
    typedef RealF  scalar_object;
    typedef vComplexF Complexified;
    typedef vRealF Realified;
    typedef vRealD DoublePrecision;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<vRealD> {
  public:
    typedef RealD  scalar_type;
    typedef vRealD vector_type;
    typedef vRealD vector_typeD;
    typedef vRealD tensor_reduced;
    typedef RealD  scalar_object;
    typedef vComplexD Complexified;
    typedef vRealD Realified;
    typedef vRealD DoublePrecision;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<vComplexH> {
  public:
    typedef ComplexF  scalar_type;
    typedef vComplexH vector_type;
    typedef vComplexD vector_typeD;
    typedef vComplexH tensor_reduced;
    typedef ComplexF  scalar_object;
    typedef vComplexH Complexified;
    typedef vRealH Realified;
    typedef vComplexD DoublePrecision;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<vComplexF> {
  public:
    typedef ComplexF  scalar_type;
    typedef vComplexF vector_type;
    typedef vComplexD vector_typeD;
    typedef vComplexF tensor_reduced;
    typedef ComplexF  scalar_object;
    typedef vComplexF Complexified;
    typedef vRealF Realified;
    typedef vComplexD DoublePrecision;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<vComplexD> {
  public:
    typedef ComplexD  scalar_type;
    typedef vComplexD vector_type;
    typedef vComplexD vector_typeD;
    typedef vComplexD tensor_reduced;
    typedef ComplexD  scalar_object;
    typedef vComplexD Complexified;
    typedef vRealD Realified;
    typedef vComplexD DoublePrecision;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<vInteger> {
  public:
    typedef  Integer scalar_type;
    typedef vInteger vector_type;
    typedef vInteger vector_typeD;
    typedef vInteger tensor_reduced;
    typedef  Integer scalar_object;
    typedef void Complexified;
    typedef void Realified;
    typedef void DoublePrecision;
    enum { TensorLevel = 0 };
  };

  // First some of my own traits
  template<typename T> struct isGridTensor {
    static const bool value = true;
    static const bool notvalue = false;
  };
  template<> struct isGridTensor<int > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<RealD > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<RealF > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<ComplexD > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<ComplexF > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<Integer > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<vRealD > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<vRealF > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<vComplexD > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<vComplexF > {
    static const bool value = false;
    static const bool notvalue = true;
  };
  template<> struct isGridTensor<vInteger > {
    static const bool value = false;
    static const bool notvalue = true;
  };

  // Match the index
  template<typename T,int Level> struct matchGridTensorIndex {
    static const bool value = (Level==T::TensorLevel);
    static const bool notvalue = (Level!=T::TensorLevel);
  };
  // What is the vtype
  template<typename T> struct isComplex {
    static const bool value = false;
  };
  template<> struct isComplex<ComplexF> {
    static const bool value = true;
  };
  template<> struct isComplex<ComplexD> {
    static const bool value = true;
  };

  //Get the SIMD vector type from a Grid tensor or Lattice<Tensor>
  template<typename T>
  struct getVectorType{
    typedef T type;
  };
  
  //Query if a tensor or Lattice<Tensor> is SIMD vector or scalar
  template<typename T>
  class isSIMDvectorized{
    template<typename U>
    static typename std::enable_if< !std::is_same< typename GridTypeMapper<typename getVectorType<U>::type>::scalar_type,   
      typename GridTypeMapper<typename getVectorType<U>::type>::vector_type>::value, char>::type test(void *);

    template<typename U>
    static double test(...);
  
  public:
    enum {value = sizeof(test<T>(0)) == sizeof(char) };
  };
  
  //Get the precision of a Lattice, tensor or scalar type in units of sizeof(float)
  template<typename T>
  class getPrecision{
  public:
    //get the vector_obj (i.e. a grid Tensor) if its a Lattice<vobj>, do nothing otherwise (i.e. if fundamental or grid Tensor)
    typedef typename getVectorType<T>::type vector_obj; 
    typedef typename GridTypeMapper<vector_obj>::scalar_type scalar_type; //get the associated scalar type. Works on fundamental and tensor types
    typedef typename GridTypeMapper<scalar_type>::Realified real_scalar_type; //remove any std::complex wrapper, should get us to the fundamental type

    enum { value = sizeof(real_scalar_type)/sizeof(float) };
  };
}

#endif

