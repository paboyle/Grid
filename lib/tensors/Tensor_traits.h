    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_traits.h

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
    typedef typename T::tensor_reduced tensor_reduced;
    typedef typename T::scalar_object scalar_object;
    typedef typename T::Complexified Complexified;
    typedef typename T::Realified Realified;
    enum { TensorLevel = T::TensorLevel };
  };

//////////////////////////////////////////////////////////////////////////////////
// Recursion stops with these template specialisations
//////////////////////////////////////////////////////////////////////////////////
  template<> class GridTypeMapper<RealF> {
  public:
    typedef RealF scalar_type;
    typedef RealF vector_type;
    typedef RealF tensor_reduced ;
    typedef RealF scalar_object;
    typedef ComplexF Complexified;
    typedef RealF Realified;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<RealD> {
  public:
    typedef RealD scalar_type;
    typedef RealD vector_type;
    typedef RealD tensor_reduced;
    typedef RealD scalar_object;
    typedef ComplexD Complexified;
    typedef RealD Realified;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<ComplexF> {
  public:
    typedef ComplexF scalar_type;
    typedef ComplexF vector_type;
    typedef ComplexF tensor_reduced;
    typedef ComplexF scalar_object;
    typedef ComplexF Complexified;
    typedef RealF Realified;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<ComplexD> {
  public:
    typedef ComplexD scalar_type;
    typedef ComplexD vector_type;
    typedef ComplexD tensor_reduced;
    typedef ComplexD scalar_object;
    typedef ComplexD Complexified;
    typedef RealD Realified;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<Integer> {
  public:
    typedef Integer scalar_type;
    typedef Integer vector_type;
    typedef Integer tensor_reduced;
    typedef Integer scalar_object;
    typedef void Complexified;
    typedef void Realified;
    enum { TensorLevel = 0 };
  };

  template<> class GridTypeMapper<vRealF> {
  public:
    typedef RealF  scalar_type;
    typedef vRealF vector_type;
    typedef vRealF tensor_reduced;
    typedef RealF  scalar_object;
    typedef vComplexF Complexified;
    typedef vRealF Realified;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<vRealD> {
  public:
    typedef RealD  scalar_type;
    typedef vRealD vector_type;
    typedef vRealD tensor_reduced;
    typedef RealD  scalar_object;
    typedef vComplexD Complexified;
    typedef vRealD Realified;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<vComplexF> {
  public:
    typedef ComplexF  scalar_type;
    typedef vComplexF vector_type;
    typedef vComplexF tensor_reduced;
    typedef ComplexF  scalar_object;
    typedef vComplexF Complexified;
    typedef vRealF Realified;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<vComplexD> {
  public:
    typedef ComplexD  scalar_type;
    typedef vComplexD vector_type;
    typedef vComplexD tensor_reduced;
    typedef ComplexD  scalar_object;
    typedef vComplexD Complexified;
    typedef vRealD Realified;
    enum { TensorLevel = 0 };
  };
  template<> class GridTypeMapper<vInteger> {
  public:
    typedef  Integer scalar_type;
    typedef vInteger vector_type;
    typedef vInteger tensor_reduced;
    typedef  Integer scalar_object;
    typedef void Complexified;
    typedef void Realified;
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

}

#endif
