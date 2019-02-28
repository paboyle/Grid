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

  // Forward declarations
  template<class T>        class iScalar;
  template<class T, int N> class iVector;
  template<class T, int N> class iMatrix;

  // These are the Grid tensors
  template<typename T> struct isGridTensor
  : public std::false_type { static constexpr bool notvalue = true; };
  template<class T> struct isGridTensor<iScalar<T>>
  : public std::true_type { static constexpr bool notvalue = false; };
  template<class T, int N> struct isGridTensor<iVector<T, N>>
  : public std::true_type { static constexpr bool notvalue = false; };
  template<class T, int N> struct isGridTensor<iMatrix<T, N>>
  : public std::true_type { static constexpr bool notvalue = false; };

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

  // This saves repeating common properties for supported Grid Scalar types
  template<typename T, typename V=void> struct GridTypeMapper_Base {};
  // TensorLevel    How many nested grid tensors
  // Rank           Rank of the grid tensor
  // count          Total number of elements, i.e. product of dimensions
  // scalar_size    Size of the underlying fundamental object (tensor_reduced) in bytes
  // size           Total size of all elements in bytes
  // Dimension(dim) Size of dimension dim
  template<typename T> struct GridTypeMapper_Base<T> {
    static constexpr int TensorLevel = 0;
    static constexpr int Rank = 0;
    static constexpr std::size_t count = 1;
    static constexpr std::size_t scalar_size = sizeof(T);
    static constexpr std::size_t size = scalar_size * count;
    static constexpr int Dimension(unsigned int dim) { return 0; }
  };

//////////////////////////////////////////////////////////////////////////////////
// Recursion stops with these template specialisations
//////////////////////////////////////////////////////////////////////////////////

  template<typename T> struct GridTypeMapper {};

  template<> struct GridTypeMapper<RealF> : public GridTypeMapper_Base<RealF> {
    typedef RealF scalar_type;
    typedef RealF vector_type;
    typedef RealD vector_typeD;
    typedef RealF tensor_reduced ;
    typedef RealF scalar_object;
    typedef ComplexF Complexified;
    typedef RealF Realified;
    typedef RealD DoublePrecision;
  };
  template<> struct GridTypeMapper<RealD> : public GridTypeMapper_Base<RealD> {
    typedef RealD scalar_type;
    typedef RealD vector_type;
    typedef RealD vector_typeD;
    typedef RealD tensor_reduced;
    typedef RealD scalar_object;
    typedef ComplexD Complexified;
    typedef RealD Realified;
    typedef RealD DoublePrecision;
  };
  template<> struct GridTypeMapper<ComplexF> : public GridTypeMapper_Base<ComplexF> {
    typedef ComplexF scalar_type;
    typedef ComplexF vector_type;
    typedef ComplexD vector_typeD;
    typedef ComplexF tensor_reduced;
    typedef ComplexF scalar_object;
    typedef ComplexF Complexified;
    typedef RealF Realified;
    typedef ComplexD DoublePrecision;
  };
  template<> struct GridTypeMapper<ComplexD> : public GridTypeMapper_Base<ComplexD> {
    typedef ComplexD scalar_type;
    typedef ComplexD vector_type;
    typedef ComplexD vector_typeD;
    typedef ComplexD tensor_reduced;
    typedef ComplexD scalar_object;
    typedef ComplexD Complexified;
    typedef RealD Realified;
    typedef ComplexD DoublePrecision;
  };
  template<> struct GridTypeMapper<Integer> : public GridTypeMapper_Base<Integer> {
    typedef Integer scalar_type;
    typedef Integer vector_type;
    typedef Integer vector_typeD;
    typedef Integer tensor_reduced;
    typedef Integer scalar_object;
    typedef void Complexified;
    typedef void Realified;
    typedef void DoublePrecision;
  };

  template<> struct GridTypeMapper<vRealF> : public GridTypeMapper_Base<vRealF> {
    typedef RealF  scalar_type;
    typedef vRealF vector_type;
    typedef vRealD vector_typeD;
    typedef vRealF tensor_reduced;
    typedef RealF  scalar_object;
    typedef vComplexF Complexified;
    typedef vRealF Realified;
    typedef vRealD DoublePrecision;
  };
  template<> struct GridTypeMapper<vRealD> : public GridTypeMapper_Base<vRealD> {
    typedef RealD  scalar_type;
    typedef vRealD vector_type;
    typedef vRealD vector_typeD;
    typedef vRealD tensor_reduced;
    typedef RealD  scalar_object;
    typedef vComplexD Complexified;
    typedef vRealD Realified;
    typedef vRealD DoublePrecision;
  };
  template<> struct GridTypeMapper<vComplexH> : public GridTypeMapper_Base<vComplexH> {
    typedef ComplexF  scalar_type;
    typedef vComplexH vector_type;
    typedef vComplexD vector_typeD;
    typedef vComplexH tensor_reduced;
    typedef ComplexF  scalar_object;
    typedef vComplexH Complexified;
    typedef vRealH Realified;
    typedef vComplexD DoublePrecision;
  };
  template<> struct GridTypeMapper<vComplexF> : public GridTypeMapper_Base<vComplexF> {
    typedef ComplexF  scalar_type;
    typedef vComplexF vector_type;
    typedef vComplexD vector_typeD;
    typedef vComplexF tensor_reduced;
    typedef ComplexF  scalar_object;
    typedef vComplexF Complexified;
    typedef vRealF Realified;
    typedef vComplexD DoublePrecision;
  };
  template<> struct GridTypeMapper<vComplexD> : public GridTypeMapper_Base<vComplexD> {
    typedef ComplexD  scalar_type;
    typedef vComplexD vector_type;
    typedef vComplexD vector_typeD;
    typedef vComplexD tensor_reduced;
    typedef ComplexD  scalar_object;
    typedef vComplexD Complexified;
    typedef vRealD Realified;
    typedef vComplexD DoublePrecision;
  };
  template<> struct GridTypeMapper<vInteger> : public GridTypeMapper_Base<vInteger> {
    typedef  Integer scalar_type;
    typedef vInteger vector_type;
    typedef vInteger vector_typeD;
    typedef vInteger tensor_reduced;
    typedef  Integer scalar_object;
    typedef void Complexified;
    typedef void Realified;
    typedef void DoublePrecision;
  };

#define GridTypeMapper_RepeatedTypes \
  typedef typename ObjectTraits::scalar_type scalar_type; \
  typedef typename ObjectTraits::vector_type vector_type; \
  typedef typename ObjectTraits::vector_typeD vector_typeD; \
  typedef typename ObjectTraits::tensor_reduced tensor_reduced; \
  typedef typename ObjectTraits::scalar_object scalar_object; \
  typedef typename ObjectTraits::Complexified Complexified; \
  typedef typename ObjectTraits::Realified Realified; \
  typedef typename ObjectTraits::DoublePrecision DoublePrecision; \
  static constexpr int TensorLevel = BaseTraits::TensorLevel + 1; \
  static constexpr std::size_t scalar_size = BaseTraits::scalar_size; \
  static constexpr std::size_t size = scalar_size * count

  template<typename T> struct GridTypeMapper<iScalar<T>> {
    using ObjectTraits = iScalar<T>;
    using BaseTraits = GridTypeMapper<T>;
    static constexpr int Rank = 1 + BaseTraits::Rank;
    static constexpr std::size_t count = 1 * BaseTraits::count;
    static constexpr int Dimension(unsigned int dim) {
      return ( dim == 0 ) ? 1 : BaseTraits::Dimension(dim - 1); }
    GridTypeMapper_RepeatedTypes;
  };

  template<typename T, int N> struct GridTypeMapper<iVector<T, N>> {
    using ObjectTraits = iVector<T, N>;
    using BaseTraits = GridTypeMapper<T>;
    static constexpr int Rank = 1 + BaseTraits::Rank;
    static constexpr std::size_t count = N * BaseTraits::count;
    static constexpr int Dimension(unsigned int dim) {
      return ( dim == 0 ) ? N : BaseTraits::Dimension(dim - 1); }
    GridTypeMapper_RepeatedTypes;
  };

  template<typename T, int N> struct GridTypeMapper<iMatrix<T, N>> {
    using ObjectTraits = iMatrix<T, N>;
    using BaseTraits = GridTypeMapper<T>;
    static constexpr int Rank = 2 + BaseTraits::Rank;
    static constexpr std::size_t count = N * N * BaseTraits::count;
    static constexpr int Dimension(unsigned int dim) {
      return ( dim == 0 || dim == 1 ) ? N : BaseTraits::Dimension(dim - 2); }
    GridTypeMapper_RepeatedTypes;
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
  
  //Query whether a tensor or Lattice<Tensor> is SIMD vector or scalar
  template<typename T, typename V=void> struct isSIMDvectorized : public std::false_type {};
  template<typename U> struct isSIMDvectorized<U, typename std::enable_if< !std::is_same<
    typename GridTypeMapper<typename getVectorType<U>::type>::scalar_type,
    typename GridTypeMapper<typename getVectorType<U>::type>::vector_type>::value, void>::type>
  : public std::true_type {};

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

