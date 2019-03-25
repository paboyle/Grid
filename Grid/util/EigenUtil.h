/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Grid/util/EigenUtil.h
 
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
#ifndef GRID_UTIL_EIGENUTIL_H
#define GRID_UTIL_EIGENUTIL_H
#include <Grid/tensors/Tensor_traits.h>
#include <Grid/Eigen/unsupported/CXX11/Tensor>

namespace Grid {
  // Custom iterator for Eigen tensors
  namespace EigenUtil {
    template <typename ETensor, bool bConst> // Is the tensor constant
    class TensorIterator_raw{
    public:
      using Index  = typename ETensor::Index;
      using Scalar = typename ETensor::Scalar;
      using FullIndex = std::array<Index, ETensor::NumIndices>;
      const ETensor * pET;
      const Index end;      // same as pET->size()
      Index       position; // position (memory order)
      Index       Seq;      // sequence (what our position would be if we were column major)
      FullIndex   indexPos;
      FullIndex   indexSize;

      inline TensorIterator_raw( ETensor & eT, Index pos = 0 ) : pET{&eT}, position{pos}, Seq{pos}, end{pET->size()} {
        for( int i = 0 ; i < ETensor::NumIndices ; i++ ) {
          indexPos[i] = 0;
          indexSize[i] = pET->dimension(i);
        }
      }
      inline TensorIterator_raw<ETensor, bConst> & operator++() {
        auto sz = pET->size();
        if( position < sz ) {
          position++;
          if( ETensor::Options & Eigen::RowMajor ) {
            for( int i = ETensor::NumIndices - 1; i != -1 && ++indexPos[i] == indexSize[i]; i-- )
              indexPos[i] = 0;
            Seq++;
          } else {
            for( int i = 0; i < ETensor::NumIndices && ++indexPos[i] == indexSize[i]; i++ )
              indexPos[i] = 0;
            Seq = 0;
            for( int i = 0; i < ETensor::NumIndices; i++ ) {
              Seq *= indexSize[i];
              Seq += indexPos[i];
            }
          }
        }
        return * this;
      }
      inline typename std::conditional<bConst,const Scalar &,Scalar &>::type operator*() const {
        assert( position >= 0 && position < pET->size() && "Attempt to access Eigen tensor iterator out of range" );
        return ( ( typename std::conditional<bConst,const Scalar *,Scalar*>::type ) pET->data() )[position];
      }
      inline bool operator!=(const TensorIterator_raw<ETensor, bConst> &r)
      { return pET == nullptr || pET != r.pET || position != r.position; }
      // These functions aren't rerquired for iterators, but they make using them easier
      inline bool AtEnd() { return position == end; }
      inline void DumpIndex(void) {
        for( auto dim : indexPos )
          std::cout << "[" << dim << "]";
      }
    };
  }
}

// The only way I could get these iterators to work is to put the begin() and end() functions in the Eigen namespace
// So if Eigen ever defines these, we'll have a conflict and have to change this
namespace Eigen {
  template <typename ETensor> using TensorIterator      = Grid::EigenUtil::TensorIterator_raw<      ETensor, false>;
  template <typename ETensor> using TensorIteratorConst = Grid::EigenUtil::TensorIterator_raw<const ETensor, true>;
  template <typename ETensor>
  inline typename std::enable_if<Grid::EigenIO::is_tensor<ETensor>::value, TensorIterator<ETensor>>::type
    begin( ETensor & ET ) { return TensorIterator<ETensor>(ET); }
  template <typename ETensor>
  inline typename std::enable_if<Grid::EigenIO::is_tensor<ETensor>::value, TensorIterator<ETensor>>::type
    end( ETensor & ET ) { return TensorIterator<ETensor>(ET, ET.size()); }
  template <typename ETensor>
  inline typename std::enable_if<Grid::EigenIO::is_tensor<ETensor>::value, TensorIteratorConst<ETensor>>::type
    begin( const ETensor & ET ) { return TensorIteratorConst<ETensor>(ET); }
  template <typename ETensor>
  inline typename std::enable_if<Grid::EigenIO::is_tensor<ETensor>::value, TensorIteratorConst<ETensor>>::type
    end( const ETensor & ET ) { return TensorIteratorConst<ETensor>(ET, ET.size()); }
}

namespace Grid {
  // for_all helper function to call the lambda for scalar
  template <typename ETensor, typename Lambda>
  typename std::enable_if<EigenIO::is_tensor_of_scalar<ETensor>::value, void>::type
  for_all_do_lambda( Lambda lambda, typename ETensor::Scalar &scalar, typename ETensor::Index &Seq,
      const std::array<typename ETensor::Index, ETensor::NumIndices> &MyIndex,
      const std::array<int, EigenIO::Traits<ETensor>::Rank> &DimGridTensor,
            std::array<int, EigenIO::Traits<ETensor>::Rank> &GridTensorIndex)
  {
    lambda( scalar, Seq++, MyIndex, GridTensorIndex );
  }
  
  // for_all helper function to call the lambda for container
  template <typename ETensor, typename Lambda>
  typename std::enable_if<EigenIO::is_tensor_of_container<ETensor>::value, void>::type
  for_all_do_lambda( Lambda lambda, typename ETensor::Scalar &container, typename ETensor::Index &Seq,
      const std::array<typename ETensor::Index, ETensor::NumIndices> &MyIndex,
      const std::array<int, EigenIO::Traits<ETensor>::Rank> &DimGridTensor,
            std::array<int, EigenIO::Traits<ETensor>::Rank> &GridTensorIndex)
  {
    using Traits = EigenIO::Traits<ETensor>;
    const int InnerRank = Traits::Rank;
    for( typename Traits::scalar_type &Source : container ) {
      lambda(Source, Seq++, MyIndex, GridTensorIndex );
      // Now increment SubIndex
      for( int i = InnerRank - 1; i != -1 && ++GridTensorIndex[i] == DimGridTensor[i]; i-- )
        GridTensorIndex[i] = 0;
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
    using Index = typename ETensor::Index;
    using Traits = EigenIO::Traits<ETensor>;
    // Check that the number of elements in the container is the product of tensor dimensions
    const Index NumScalars = ET.size();
    assert( NumScalars > 0 && "EigenUtil: tensor has no elements" );
    Index ScalarElementCount{1};
    const int rank{ETensor::NumIndices};
    std::array<Index, rank> DimTensor, MyIndex;
    for(int i = 0; i < rank; i++ ) {
      DimTensor[i] = ET.dimension(i);
      ScalarElementCount *= DimTensor[i];
      MyIndex[i] = 0;
    }
    assert( NumScalars == ScalarElementCount && "EigenUtil: tensor size not product of dimensions" );
    // Save the GridTensor dimensions
    const int InnerRank{Traits::Rank};
    std::array<int, InnerRank> DimGridTensor, GridTensorIndex;
    for(int i = 0; i < InnerRank; i++ ) {
      DimGridTensor[i] = Traits::Dimension(i);
      GridTensorIndex[i] = 0;
    }
    // Now walk the tensor in memory order
    Index Seq = 0;
    Scalar * pScalar = ET.data();
    for( Index j = 0; j < NumScalars; j++ ) {
      for_all_do_lambda<ETensor, Lambda>( lambda, * pScalar, Seq, MyIndex, DimGridTensor, GridTensorIndex );
      // Now increment the index to pass to the lambda (bearing in mind we're walking in memory order)
      if( ETensor::Options & Eigen::RowMajor ) {
        for( int i = rank - 1; i != -1 && ++MyIndex[i] == DimTensor[i]; i-- )
          MyIndex[i] = 0;
      } else {
        for( int i = 0; i < rank && ++MyIndex[i] == DimTensor[i]; i++ )
          MyIndex[i] = 0;
        Seq = 0;
        for( int i = 0; i < rank; i++ ) {
          Seq *= DimTensor[i];
          Seq += MyIndex[i];
        }
        Seq *= Traits::count;
      }
      pScalar++;
    }
  }
  
  // Sequential initialisation of tensors up to specified precision
  // Would have preferred to define template variables for this, but that's c++ 17
  template <typename ETensor>
  typename std::enable_if<EigenIO::is_tensor<ETensor>::value && !EigenIO::Traits<ETensor>::is_complex>::type
  SequentialInit( ETensor &ET, typename EigenIO::Traits<ETensor>::scalar_type Inc = 1, unsigned short Precision = 0 )
  {
    using Traits = EigenIO::Traits<ETensor>;
    using scalar_type = typename Traits::scalar_type;
    using Index = typename ETensor::Index;
    for_all( ET, [&](scalar_type &c, Index n, const std::array<Index, ETensor::NumIndices> &TensorIndex,
                     const std::array<int, Traits::Rank> &ScalarIndex ) {
      scalar_type x = Inc * static_cast<scalar_type>(n);
      if( Precision ) {
        std::stringstream s;
        s << std::setprecision(Precision) << x;
        s >> x;
      }
      c = x;
    } );
  }
  
  template <typename ETensor>
  typename std::enable_if<EigenIO::is_tensor<ETensor>::value && EigenIO::Traits<ETensor>::is_complex>::type
  SequentialInit( ETensor &ET, typename EigenIO::Traits<ETensor>::scalar_type Inc={1,-1}, unsigned short Precision = 0 )
  {
    using Traits = EigenIO::Traits<ETensor>;
    using scalar_type = typename Traits::scalar_type;
    using Index = typename ETensor::Index;
    for_all( ET, [&](scalar_type &c, Index n, const std::array<Index, ETensor::NumIndices> &TensorIndex,
                     const std::array<int, Traits::Rank> &ScalarIndex ) {
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
  template <typename T>
  typename std::enable_if<EigenIO::is_tensor<T>::value, void>::type
  dump_tensor(T &t, const char * pName = nullptr)
  {
    using Traits = EigenIO::Traits<T>;
    using scalar_type = typename Traits::scalar_type;
    using Index = typename T::Index;
    const int rank{T::NumIndices};
    const auto &dims = t.dimensions();
    std::cout << "Dumping rank " << rank + Traits::Rank << ((T::Options & Eigen::RowMajor) ? ", row" : ", column") << "-major tensor ";
    if( pName )
      std::cout << pName;
    for( int i = 0 ; i < rank; i++ ) std::cout << "[" << dims[i] << "]";
    for( int i = 0 ; i < Traits::Rank; i++ ) std::cout << "(" << Traits::Dimension(i) << ")";
    std::cout << " in memory order:" << std::endl;
#ifdef OLD_DEFINITION
    for_all( t, [&](scalar_type &c, Index n, const std::array<Index, rank> &TensorIndex,
                    const std::array<int, Traits::Rank> &ScalarIndex ){
      std::cout << "  ";
      for( auto dim : TensorIndex )
        std::cout << "[" << dim << "]";
      for( auto dim : ScalarIndex )
        std::cout << "(" << dim << ")";
      std::cout << " = " << c << std::endl;
    } );
#else
    for( auto it = begin(t); !it.AtEnd(); ++it ) {
      std::cout << "  ";
      it.DumpIndex();
      std::cout << " = " << (const typename T::Scalar)(*it) << std::endl;
    }
#endif
    std::cout << "========================================" << std::endl;
  }
  
  template <typename T>
  typename std::enable_if<!EigenIO::is_tensor<T>::value, void>::type
  dump_tensor(T &t, const char * pName = nullptr)
  {
    std::cout << "Dumping non-tensor object ";
    if( pName ) std::cout << pName;
    std::cout << "=" << t;
  }
}
#endif
