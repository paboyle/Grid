/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Grid/util/Eigen.h
 
 Copyright (C) 2019
 
 Author: Michael Marshall <michael.marshall@ed.ac.uk>
 
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
#ifndef GRID_UTIL_EIGEN_H
#define GRID_UTIL_EIGEN_H
#include <Grid/tensors/Tensor_traits.h>
#include <Grid/Eigen/unsupported/CXX11/Tensor>

namespace Grid {
  // for_all helper function to call the lambda for scalar
  template <typename ETensor, typename Lambda>
  typename std::enable_if<EigenIO::is_tensor_of_scalar<ETensor>::value, void>::type
  for_all_do_lambda( Lambda lambda, typename ETensor::Scalar &scalar, typename ETensor::Index &Seq,
                    std::array<std::size_t, ETensor::NumIndices + GridTypeMapper<typename ETensor::Scalar>::Rank> &MyIndex)
  {
    lambda( scalar, Seq++, MyIndex );
  }
  
  // for_all helper function to call the lambda for container
  template <typename ETensor, typename Lambda>
  typename std::enable_if<EigenIO::is_tensor_of_container<ETensor>::value, void>::type
  for_all_do_lambda( Lambda lambda, typename ETensor::Scalar &container, typename ETensor::Index &Seq,
                    std::array<std::size_t, ETensor::NumIndices + GridTypeMapper<typename ETensor::Scalar>::Rank> &MyIndex)
  {
    using Traits = GridTypeMapper<typename ETensor::Scalar>;
    const auto rank{ETensor::NumIndices};
    const auto InnerRank = Traits::Rank;
    for( typename Traits::scalar_type &Source : container ) {
      lambda(Source, Seq++, MyIndex );
      // Now increment SubIndex
      for( auto i = InnerRank - 1; i != -1 && ++MyIndex[rank + i] == Traits::Dimension(i); i-- )
        MyIndex[rank + i] = 0;
    }
  }
  
  // Calls a lamda (passing index and sequence number) for every member of an Eigen::Tensor
  // For efficiency, iteration proceeds in memory order,
  // ... but parameters guaranteed to be the same regardless of memory order
  template <typename ETensor, typename Lambda>
  typename std::enable_if<EigenIO::is_tensor<ETensor>::value, void>::type
  for_all( ETensor &ET, Lambda lambda )
  {
    using Scalar = typename ETensor::Scalar; // This could be a Container - we'll check later
    const std::size_t NumScalars = ET.size();
    assert( NumScalars > 0 );
    using Index = typename ETensor::Index;
    Index ScalarElementCount{1};
    const auto InnerRank = GridTypeMapper<Scalar>::Rank;
    const auto rank{ETensor::NumIndices};
    std::array<std::size_t, rank + InnerRank> Dims;
    for(auto i = 0; i < rank; i++ ) {
      auto dim = ET.dimension(i);
      assert( dim > 0 );
      Dims[i] = static_cast<std::size_t>(dim);
      assert( Dims[i] == dim ); // check we didn't lose anything in the conversion
      ScalarElementCount *= Dims[i];
    }
    // Check that the number of containers is correct ... and we didn't lose anything in conversions
    assert( NumScalars == ScalarElementCount );
    // If the Scalar is actually a container, add the inner Scalar's dimensions
    size_t InnerScalarCount{1};
    for(auto i = 0; i < InnerRank; i++ ) {
      auto dim = GridTypeMapper<Scalar>::Dimension(i);
      assert( dim > 0 );
      Dims[rank + i] = static_cast<std::size_t>(dim);
      assert( Dims[rank + i] == dim ); // check we didn't lose anything in the conversion
      InnerScalarCount *= dim;
    }
    assert(GridTypeMapper<Scalar>::count == InnerScalarCount);
    assert(GridTypeMapper<Scalar>::size  == sizeof( Scalar ));
    std::array<std::size_t, rank + InnerRank> MyIndex;
    for( auto &idx : MyIndex ) idx = 0;
    Index Seq = 0;
    Scalar * pScalar = ET.data();
    for( std::size_t j = 0; j < NumScalars; j++ ) {
      for_all_do_lambda<ETensor, Lambda>( lambda, * pScalar, Seq, MyIndex );
      // Now increment the index to pass to the lambda (bearing in mind we're walking in memory order)
      if( ETensor::Options & Eigen::RowMajor ) {
        for( auto i = rank - 1; i != -1 && ++MyIndex[i] == Dims[i]; i-- )
          MyIndex[i] = 0;
      } else {
        for( auto i = 0; i < rank && ++MyIndex[i] == Dims[i]; i++ )
          MyIndex[i] = 0;
        size_t NewSeq = 0;
        for( auto i = 0; i < rank + InnerRank ; i++ ) {
          NewSeq *= Dims[i];
          NewSeq += MyIndex[i];
        }
        Seq = static_cast<Index>( NewSeq );
      }
      pScalar++;
    }
  }
  
  // Sequential initialisation of tensors
  // Would have preferred to define template variables for this, but that's c++ 17
  template <typename ETensor>
  typename std::enable_if<EigenIO::is_tensor<ETensor>::value && !is_complex<typename GridTypeMapper<typename ETensor::Scalar>::scalar_type>::value, void>::type
  SequentialInit( ETensor &ET, typename GridTypeMapper<typename ETensor::Scalar>::scalar_type Inc = 1,
                 unsigned short Precision = 0 )
  {
    using Traits = GridTypeMapper<typename ETensor::Scalar>;
    using scalar_type = typename Traits::scalar_type;
    for_all( ET, [&](scalar_type &c, typename ETensor::Index n, const std::array<size_t, ETensor::NumIndices + Traits::Rank> &Dims ) {
      scalar_type x = Inc * static_cast<scalar_type>(n);
      if( Precision ) {
        std::stringstream s;
        s << std::scientific << std::setprecision(Precision) << x;
        s >> x;
      }
      c = x;
    } );
  }
  
  template <typename ETensor>
  typename std::enable_if<EigenIO::is_tensor<ETensor>::value && is_complex<typename GridTypeMapper<typename ETensor::Scalar>::scalar_type>::value, void>::type
  SequentialInit( ETensor &ET, typename GridTypeMapper<typename ETensor::Scalar>::scalar_type Inc={1,-1},
                 unsigned short Precision = 0 )
  {
    using Traits = GridTypeMapper<typename ETensor::Scalar>;
    using scalar_type = typename Traits::scalar_type;
    for_all( ET, [&](scalar_type &c, typename ETensor::Index n, const std::array<size_t, ETensor::NumIndices + Traits::Rank> &Dims ) {
      auto re = Inc.real();
      auto im = Inc.imag();
      re *= n;
      im *= n;
      if( Precision ) {
        std::stringstream s;
        s << std::setprecision(Precision) << re;
        s >> re;
        s.clear();
        s << im;
        s >> im;
      }
      c = scalar_type(re,im);
    } );
  }
  
  // Helper to dump a tensor
#ifdef DEBUG
#define dump_tensor(args...) dump_tensor_func(args)
  template <typename T>
  typename std::enable_if<EigenIO::is_tensor<T>::value, void>::type
  dump_tensor_func(T &t, const char * pName = nullptr)
  {
    using Traits = GridTypeMapper<typename T::Scalar>;
    const auto rank{T::NumIndices};
    const auto &dims = t.dimensions();
    std::cout << "Dumping rank " << rank << ((T::Options & Eigen::RowMajor) ? ", row" : ", column") << "-major tensor ";
    if( pName )
      std::cout << pName;
    for( auto i = 0 ; i < rank; i++ ) std::cout << "[" << dims[i] << "]";
    std::cout << " in memory order:" << std::endl;
    for_all( t, [&](typename Traits::scalar_type &c, typename T::Index index, const std::array<size_t, T::NumIndices + Traits::Rank> &Dims ){
      std::cout << "  ";
      for( auto dim : Dims )
        std::cout << "[" << dim << "]";
      std::cout << " = " << c << std::endl;
    } );
    std::cout << "========================================" << std::endl;
  }
  
  template <typename T>
  typename std::enable_if<!EigenIO::is_tensor<T>::value, void>::type
  dump_tensor_func(T &t, const char * pName = nullptr)
  {
    std::cout << "Dumping non-tensor object ";
    if( pName )
      std::cout << pName;
    std::cout << "=" << t;
  }
  
  // Helper to dump a tensor in memory order
  // Kind of superfluous given the above ... just keeping in case I need to fall back to this
#define DumpMemoryOrder(args...) DumpMemoryOrder_func(args)
  template <typename T>
  typename std::enable_if<EigenIO::is_tensor_of_scalar<T>::value, void>::type
  DumpMemoryOrder_func(T &t, const char * pName = nullptr)
  {
    const auto rank = t.rank();
    const auto &dims = t.dimensions();
    std::cout << "Dumping rank " << rank << ((T::Options & Eigen::RowMajor) ? ", row" : ", column") << "-major tensor ";
    if( pName )
      std::cout << pName;
    for( auto d : dims ) std::cout << "[" << d << "]";
    std::cout << " in memory order:" << std::endl;
    const typename T::Scalar * p = t.data();
    const auto size = t.size();
    const typename T::Scalar * pEnd = p + size;
    if( rank <= 2 ) {
      for( unsigned int i = 0 ; i < t.size() ; i++ )
        std::cout << "[" << i << "]=" << *p++ << " ";
      std::cout << std::endl;
    } else {
      const auto innersize = dims[rank-2] * dims[rank-1];
      using Index = typename T::Index;
      std::vector<Index> idx(rank - 2);
      for( auto &i : idx ) i = 0;
      Index idxCounter = 0;
      while( p < pEnd ) {
        if( T::Options & Eigen::RowMajor ) {
          if( pName )
            std::cout << pName;
          idxCounter = 0;
          for(auto i = 0 ; i < rank - 2 ; i++)
            std::cout << "[" << idx[i] << "]:";
        }
        for( unsigned int i = 0 ; i < innersize ; i++ )
          std::cout << " [" << idxCounter++ << "]=" << *p++;
        if( T::Options & Eigen::RowMajor )
          std::cout << std::endl;
        // Now increment MyIndex
        for( auto i = rank - 3; i != -1 && ++idx[i] == dims[i]; i-- )
          idx[i] = 0;
      }
      if( ! ( T::Options & Eigen::RowMajor ) )
        std::cout << std::endl;
    }
  }
#else
#define dump_tensor(args...)
#define DumpMemoryOrder(args...)
#endif
}
#endif
