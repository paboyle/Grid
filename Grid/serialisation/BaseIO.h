    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/serialisation/BaseIO.h

    Copyright (C) 2015

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Guido Cossu <guido.cossu@ed.ac.uk>

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
#ifndef GRID_SERIALISATION_ABSTRACT_READER_H
#define GRID_SERIALISATION_ABSTRACT_READER_H

#include <type_traits>
#include <Grid/tensors/Tensors.h>
#include <Grid/serialisation/VectorUtils.h>
#include <Grid/Eigen/unsupported/CXX11/Tensor>

namespace Grid {
  namespace EigenIO {
    template<typename T> struct is_complex : public std::false_type {};
    template<typename T> struct is_complex<std::complex<T>>
        : std::integral_constant<bool, std::is_arithmetic<T>::value> {};

    // Eigen tensors can be composed of arithmetic scalar and complex types
    template<typename T> struct is_scalar : std::integral_constant<bool,
                            std::is_arithmetic<T>::value || is_complex<T>::value> {};

    // Eigen tensors can also be composed of a limited number of containers
    // NB: grid tensors (iScalar, iVector, iMatrix) have stricter limits on complex types
    template <typename T> struct is_container                       : public std::false_type {};
    template <typename T> struct is_container<iScalar<T>>           : public std::true_type {};
    template <typename T, int N> struct is_container<iVector<T, N>> : public std::true_type {};
    template <typename T, int N> struct is_container<iMatrix<T, N>> : public std::true_type {};
    template <typename T, std::size_t N> struct is_container<std::array<T, N>> : public std::true_type {};

    // Is this an Eigen tensor
    template<typename T> struct is_tensor : std::integral_constant<bool,
      std::is_base_of<Eigen::TensorBase<T, Eigen::ReadOnlyAccessors>, T>::value> {};

    // Is this an Eigen tensor of a supported scalar
    template<typename T> struct is_tensor_of_scalar : std::integral_constant<bool,
    is_tensor<T>::value && is_scalar<typename T::Scalar>::value> {};

    // Is this an Eigen tensor of a supported container
    template<typename T> struct is_tensor_of_container : std::integral_constant<bool,
    is_tensor<T>::value && is_container<typename T::Scalar>::value> {};

    // These traits describe the Eigen tensor scalar objects supported for IO
    // Mostly intended for grid tensors, i.e. compositions of iScalar, iVector and iMatrix
    // but also supports fixed size arrays
    template <typename T, typename C = void> struct Traits {}; // C needed for specialisation
    // This defines the bottom level - i.e. it's a description of the underlying scalar
    template <typename T> struct Traits<T, typename std::enable_if<is_scalar<T>::value, void>::type> {
      static constexpr unsigned int depth = 0;  // How many levels of Grid Tensor there are (TensorLevel)
      static constexpr unsigned int rank = 0;   // The rank of the grid tensor (i.e. how many indices used)
      static constexpr unsigned int rank_non_trivial = 0; // As per rank, but excludes those of dimension 1
      static constexpr unsigned int count = 1;  // total number of elements (i.e. product of dimensions)
      using scalar_type = T;                              // Type of the underlying scalar
      static constexpr std::size_t scalar_size = sizeof(T);    // Size of the underlying scalar in bytes
      static constexpr std::size_t size = scalar_size * count; // total size of elements in bytes
      static constexpr std::size_t Dimension(unsigned int dim) { return 0; } // Dimension size
      static constexpr std::size_t DimensionNT(unsigned int dim) { return 0; } // non-trivial dim size
      // e.g. iScalar<iVector<Complex,1>>
      //      depth = 2
      //      rank  = 1
      //      rank_non_trivial = 0
      //      count  = 1
      // e.g. iVector<iMatrix<Complex,3>,4>
      //      depth = 2
      //      rank  = 3
      //      rank_non_trivial = 3
      //      count  = 36
      // e.g. iScalar<iVector<iMatrix<Complex,4>,3>>
      //      depth = 3
      //      rank  = 3
      //      rank_non_trivial = 3
      //      count  = 48
    };
    template <typename T> struct Traits<iScalar<T>> {
      static constexpr unsigned int depth = 1 + Traits<T>::depth;
      static constexpr unsigned int rank = 0 + Traits<T>::rank;
      static constexpr unsigned int rank_non_trivial = 0 + Traits<T>::rank_non_trivial;
      static constexpr unsigned int count = 1 * Traits<T>::count;
      using scalar_type = typename Traits<T>::scalar_type;
      static constexpr std::size_t scalar_size = Traits<T>::scalar_size;
      static constexpr std::size_t size = scalar_size * count;
      static constexpr std::size_t Dimension(unsigned int dim) {
        return ( dim == 0 ) ? 1 : Traits<T>::Dimension(dim - 1); }
      static constexpr std::size_t DimensionNT(unsigned int dim) {
        return Traits<T>::DimensionNT(dim); }
    };
    template <typename T, int N> struct Traits<iVector<T, N>> {
      static constexpr unsigned int depth = 1 + Traits<T>::depth;
      static constexpr unsigned int rank = 1 + Traits<T>::rank;
      static constexpr unsigned int rank_non_trivial = (N>1 ? 1 : 0) + Traits<T>::rank_non_trivial;
      static constexpr unsigned int count = N * Traits<T>::count;
      using scalar_type = typename Traits<T>::scalar_type;
      static constexpr std::size_t scalar_size = Traits<T>::scalar_size;
      static constexpr std::size_t size = scalar_size * count;
      static constexpr std::size_t Dimension(unsigned int dim) {
        return ( dim == 0 ) ? N : Traits<T>::Dimension(dim - 1); }
      static constexpr std::size_t DimensionNT(unsigned int dim) {
        return ( N == 1 ) ? Traits<T>::DimensionNT(dim) : ( dim == 0 ) ? N : Traits<T>::DimensionNT(dim - 1);
      }
    };
    template <typename T, int N> struct Traits<iMatrix<T, N>> {
      static constexpr unsigned int depth = 1 + Traits<T>::depth;
      static constexpr unsigned int rank = 2 + Traits<T>::rank;
      static constexpr unsigned int rank_non_trivial = (N>1 ? 2 : 0) + Traits<T>::rank_non_trivial;
      static constexpr unsigned int count = N * N * Traits<T>::count;
      using scalar_type = typename Traits<T>::scalar_type;
      static constexpr std::size_t scalar_size = Traits<T>::scalar_size;
      static constexpr std::size_t size = scalar_size * count;
      static constexpr std::size_t Dimension(unsigned int dim) {
        return ( dim == 0 || dim == 1 ) ? N : Traits<T>::Dimension(dim - 2); }
      static constexpr std::size_t DimensionNT(unsigned int dim) {
        return ( N == 1 ) ? Traits<T>::DimensionNT(dim) : ( dim == 0 || dim == 1 ) ? N : Traits<T>::DimensionNT(dim - 2);
      }
    };
    template <typename T, int N> struct Traits<std::array<T, N>> : Traits<iVector<T, N>> {};
    // Tensors have the same traits as their top-level scalar
    // Shouldn't be necessary ... but I make the mistake of getting traits of the tensor so often
    // that I am tempted to define this.
    // HOWEVER, Eigen tensors have a dynamic size flavour, but the scalars are (currently) all fixed size
    //template <typename T> struct Traits<T, typename std::enable_if<is_tensor<T>::value, void>::type> : Traits<T> {};
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
    const auto InnerRank = EigenIO::Traits<Scalar>::rank_non_trivial;
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
    // If the Scalar is actually a container, add the inner Scalar's non-trivial dimensions
    size_t InnerScalarCount{1};
    for(auto i = 0; i < InnerRank; i++ ) {
      auto dim = EigenIO::Traits<Scalar>::DimensionNT(i);
      assert( dim > 1 );
      Dims[rank + i] = static_cast<std::size_t>(dim);
      assert( Dims[rank + i] == dim ); // check we didn't lose anything in the conversion
      InnerScalarCount *= dim;
    }
    assert(EigenIO::Traits<Scalar>::count == InnerScalarCount);
    assert(EigenIO::Traits<Scalar>::size  == sizeof( Scalar ));
    std::array<std::size_t, rank + InnerRank> MyIndex;
    for( auto &idx : MyIndex ) idx = 0;
    Index Seq = 0;
    Scalar * pScalar = ET.data();
    for( std::size_t j = 0; j < NumScalars; j++ ) {
      // if constexpr is C++ 17 ... but otherwise need two specialisations (Container vs Scalar)
      if constexpr ( InnerRank == 0 ) {
        lambda( * pScalar, Seq++, MyIndex );
      } else {
        for( typename Scalar::scalar_type &Source : * pScalar ) {
          lambda(Source, Seq++, MyIndex );
          // Now increment SubIndex
          for( auto i = rank + InnerRank - 1; i != rank - 1 && ++MyIndex[i] == Dims[i]; i-- )
            MyIndex[i] = 0;
        }
      }
      // Now increment the index to pass to the lambda (bearing in mind we're walking in memory order)
      if( ETensor::Options & Eigen::RowMajor ) {
        for( auto i = rank - 1; i != -1 && ++MyIndex[i] == Dims[i]; i-- )
          MyIndex[i] = 0;
      } else {
        for( auto i = 0; i < rank && ++MyIndex[i] == Dims[i]; i++ )
          MyIndex[i] = 0;
        Seq = 0;
        for( auto i = 0; i < rank + InnerRank ; i++ ) {
          Seq *= Dims[i];
          Seq += MyIndex[i];
        }
      }
      pScalar++;
    }
  }

  // Helper to dump a tensor in memory order
  template <typename T>
  typename std::enable_if<EigenIO::is_tensor_of_scalar<T>::value, void>::type
  DumpMemoryOrder(T t, const char * pName = nullptr)
  {
    const auto dims = t.dimensions();
    const auto rank = t.rank();
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

  // Abstract writer/reader classes ////////////////////////////////////////////
  // static polymorphism implemented using CRTP idiom
  class Serializable;

  // Static abstract writer
  template <typename T>
  class Writer
  {
  public:
    Writer(void);
    virtual ~Writer(void) = default;
    void push(const std::string &s);
    void pop(void);
    template <typename U>
    typename std::enable_if<std::is_base_of<Serializable, U>::value, void>::type
    write(const std::string& s, const U &output);
    template <typename U>
    typename std::enable_if<!std::is_base_of<Serializable, U>::value
                         && !std::is_base_of<Eigen::TensorBase<U, Eigen::ReadOnlyAccessors>, U>::value, void>::type
    write(const std::string& s, const U &output);
    template <typename U>
    void write(const std::string &s, const iScalar<U> &output);
    template <typename U, int N>
    void write(const std::string &s, const iVector<U, N> &output);
    template <typename U, int N>
    void write(const std::string &s, const iMatrix<U, N> &output);
    template <typename ETensor>
    typename std::enable_if<std::is_base_of<Eigen::TensorBase<ETensor, Eigen::ReadOnlyAccessors>, ETensor>::value && EigenIO::is_scalar<typename ETensor::Scalar>::value, void>::type
    write(const std::string &s, const ETensor &output);
    template <typename ETensor>
    typename std::enable_if<std::is_base_of<Eigen::TensorBase<ETensor, Eigen::ReadOnlyAccessors>, ETensor>::value && EigenIO::is_container<typename ETensor::Scalar>::value, void>::type
    write(const std::string &s, const ETensor &output);

    void         scientificFormat(const bool set);
    bool         isScientific(void);
    void         setPrecision(const unsigned int prec);
    unsigned int getPrecision(void);
  private:
    T            *upcast;
    bool         scientific_{false};
    unsigned int prec_{0};
  };
  
  // Static abstract reader
  template <typename T>
  class Reader
  {
  public:
    Reader(void);
    virtual ~Reader(void) = default;
    bool push(const std::string &s);
    void pop(void);
    template <typename U>
    typename std::enable_if<std::is_base_of<Serializable, U>::value, void>::type
    read(const std::string& s, U &output);
    template <typename U>
    typename std::enable_if<!std::is_base_of<Serializable, U>::value, void>::type
    read(const std::string& s, U &output);
    template <typename U>
    void read(const std::string &s, iScalar<U> &output);
    template <typename U, int N>
    void read(const std::string &s, iVector<U, N> &output);
    template <typename U, int N>
    void read(const std::string &s, iMatrix<U, N> &output);
  protected:
    template <typename U>
    void fromString(U &output, const std::string &s);
  private:
    T *upcast;
  };

   // What is the vtype
  template<typename T> struct isReader {
    static const bool value = false;
  };
  template<typename T> struct isWriter {
    static const bool value = false;
  };

  // Writer template implementation
  template <typename T>
  Writer<T>::Writer(void)
  {
    upcast = static_cast<T *>(this);
  }
  
  template <typename T>
  void Writer<T>::push(const std::string &s)
  {
    upcast->push(s);
  }
  
  template <typename T>
  void Writer<T>::pop(void)
  {
    upcast->pop();
  }
  
  template <typename T>
  template <typename U>
  typename std::enable_if<std::is_base_of<Serializable, U>::value, void>::type
  Writer<T>::write(const std::string &s, const U &output)
  {
    U::write(*this, s, output);
  }
  
  template <typename T>
  template <typename U>
  typename std::enable_if<!std::is_base_of<Serializable, U>::value
                       && !std::is_base_of<Eigen::TensorBase<U, Eigen::ReadOnlyAccessors>, U>::value, void>::type
  Writer<T>::write(const std::string &s, const U &output)
  {
    upcast->writeDefault(s, output);
  }


  template <typename T>
  template <typename U>
  void Writer<T>::write(const std::string &s, const iScalar<U> &output)
  {
    upcast->writeDefault(s, tensorToVec(output));
  }

  template <typename T>
  template <typename U, int N>
  void Writer<T>::write(const std::string &s, const iVector<U, N> &output)
  {
    upcast->writeDefault(s, tensorToVec(output));
  }

  template <typename T>
  template <typename U, int N>
  void Writer<T>::write(const std::string &s, const iMatrix<U, N> &output)
  {
    upcast->writeDefault(s, tensorToVec(output));
  }

  // Eigen::Tensors of arithmetic/complex base type
  template <typename T>
  template <typename ETensor>
  typename std::enable_if<std::is_base_of<Eigen::TensorBase<ETensor, Eigen::ReadOnlyAccessors>, ETensor>::value && EigenIO::is_scalar<typename ETensor::Scalar>::value, void>::type
  Writer<T>::write(const std::string &s, const ETensor &output)
  {
    const typename ETensor::Index NumElements{output.size()};
    assert( NumElements > 0 );
    if( NumElements == 1 )
      upcast->writeDefault(s, * output.data());
    else {
      // We're only interested in non-trivial dimensions (i.e. dimensions > 1)
      unsigned int TrivialDimCount{0};
      std::vector<size_t> NonTrivialDims;
      NonTrivialDims.reserve(output.NumDimensions); // Make sure we only do one malloc
      for(auto i = 0; i < output.NumDimensions; i++ ) {
        auto dim = output.dimension(i);
        if( dim <= 1 ) {
          TrivialDimCount++;
          assert( dim == 1 ); // Not expecting dimension to be <= 0
        } else {
          size_t s = static_cast<size_t>(dim);
          assert( s == dim ); // check we didn't lose anything in the conversion
          NonTrivialDims.push_back(s);
        }
      }
      // NB: NumElements > 1 implies this is not a scalar, so some dims should be left
      assert( output.NumDimensions > TrivialDimCount );
      // If the Tensor isn't in Row-Major order, then we'll need to copy it's data
      const bool CopyData{ETensor::Layout != Eigen::StorageOptions::RowMajor};
      using Scalar = typename ETensor::Scalar;
      const Scalar * pWriteBuffer;
      Scalar * pCopyBuffer = nullptr;
      if( !CopyData )
        pWriteBuffer = output.data();
      else {
        // Regardless of the Eigen::Tensor storage order, the copy will be Row Major
        pCopyBuffer = static_cast<Scalar *>(malloc(sizeof(Scalar) * NumElements));
        pWriteBuffer = pCopyBuffer;
        std::array<typename ETensor::Index, ETensor::NumIndices> MyIndex;
        for( auto &idx : MyIndex ) idx = 0;
        for( typename ETensor::Index n = 0; n < NumElements; n++ ) {
          pCopyBuffer[n] = output( MyIndex );
          // Now increment the index
          for( int i = output.NumDimensions - 1; i >= 0 && ++MyIndex[i] == output.dimension(i); i-- )
            MyIndex[i] = 0;
        }
      }
      upcast->template writeMultiDim<typename ETensor::Scalar>(s, NonTrivialDims, pWriteBuffer, NumElements);
      if( pCopyBuffer ) free( pCopyBuffer );
    }
  }
  
  // Eigen::Tensors of Grid tensors (iScalar, iVector, iMatrix)
  template <typename T>
  template <typename ETensor/*, typename U, int N*/>
  typename std::enable_if<std::is_base_of<Eigen::TensorBase<ETensor, Eigen::ReadOnlyAccessors>, ETensor>::value && EigenIO::is_container<typename ETensor::Scalar>::value, void>::type
  Writer<T>::write(const std::string &s, const ETensor &output)
  {
    const typename ETensor::Index NumElements{output.size()};
    assert( NumElements > 0 );
    if( NumElements == 1 )
      upcast->writeDefault(s, tensorToVec(* output.data()));
    else {
      // We're only interested in non-trivial dimensions (i.e. dimensions > 1)
      unsigned int TrivialDimCount{0};
      std::vector<size_t> NonTrivialDims;
      NonTrivialDims.reserve(output.NumDimensions + EigenIO::Traits<typename ETensor::Scalar>::rank_non_trivial); // Make sure we only do one malloc
      for(auto i = 0; i < output.NumDimensions; i++ ) {
        auto dim = output.dimension(i);
        if( dim <= 1 ) {
          TrivialDimCount++;
          assert( dim == 1 ); // Not expecting dimension to be <= 0
        } else {
          size_t s = static_cast<size_t>(dim);
          assert( s == dim ); // check we didn't lose anything in the conversion
          NonTrivialDims.push_back(s);
        }
      }
      // NB: NumElements > 1 implies this is not a scalar, so some dims should be left
      assert( output.NumDimensions > TrivialDimCount );
      // Now add the extra dimensions, based on object zero
      typename TensorToVec<typename ETensor::Scalar>::type ttv = tensorToVec(* output.data());
      Flatten<typename TensorToVec<typename ETensor::Scalar>::type> f(ttv);
      const std::vector<size_t> & ExtraDims{f.getDim()};
      assert(ExtraDims.size() == EigenIO::Traits<typename ETensor::Scalar>::rank_non_trivial);
      size_t ExtraCount{1};
      for( auto i : ExtraDims ) {
        assert( i > 0 );
        ExtraCount *= i;
        NonTrivialDims.push_back(i);
      }
      assert(EigenIO::Traits<typename ETensor::Scalar>::count == ExtraCount);
      assert(EigenIO::Traits<typename ETensor::Scalar>::size  == sizeof( typename ETensor::Scalar ));
      // If the Tensor isn't in Row-Major order, then we'll need to copy it's data
      const bool CopyData{ETensor::Layout != Eigen::StorageOptions::RowMajor};
      using Scalar = typename ETensor::Scalar::scalar_type;
      const Scalar * pWriteBuffer;
      Scalar * pCopyBuffer = nullptr;
      const typename ETensor::Index TotalNumElements = NumElements * ExtraCount;
      if( !CopyData )
        pWriteBuffer = output.data()->begin();
      else {
        // Regardless of the Eigen::Tensor storage order, the copy will be Row Major
        pCopyBuffer = static_cast<Scalar *>(malloc(TotalNumElements * sizeof(Scalar)));
        pWriteBuffer = pCopyBuffer;
        Scalar * pCopy = pCopyBuffer;
        std::array<typename ETensor::Index, ETensor::NumIndices> MyIndex;
        for( auto &idx : MyIndex ) idx = 0;
        for( typename ETensor::Index n = 0; n < NumElements; n++ ) {
          // Copy the grid tensor
          for( const Scalar &Source : output( MyIndex ) )
            * pCopy ++ =  Source;
          // Now increment the index
          for( int i = output.NumDimensions - 1; i >= 0 && ++MyIndex[i] == output.dimension(i); i-- )
            MyIndex[i] = 0;
        }
      }
      upcast->template writeMultiDim<Scalar>(s, NonTrivialDims, pWriteBuffer, TotalNumElements);
      if( pCopyBuffer ) free( pCopyBuffer );
    }
  }

  template <typename T>
  void Writer<T>::scientificFormat(const bool set)
  {
    scientific_ = set;
  }

  template <typename T>
  bool Writer<T>::isScientific(void)
  {
    return scientific_;
  }

  template <typename T>
  void Writer<T>::setPrecision(const unsigned int prec)
  {
    prec_ = prec;
  }

  template <typename T>
  unsigned int Writer<T>::getPrecision(void)
  {
    return prec_;
  }
  
  // Reader template implementation
  template <typename T>
  Reader<T>::Reader(void)
  {
    upcast = static_cast<T *>(this);
  }
  
  template <typename T>
  bool Reader<T>::push(const std::string &s)
  {
    return upcast->push(s);
  }
  
  template <typename T>
  void Reader<T>::pop(void)
  {
    upcast->pop();
  }
  
  template <typename T>
  template <typename U>
  typename std::enable_if<std::is_base_of<Serializable, U>::value, void>::type
  Reader<T>::read(const std::string &s, U &output)
  {
    U::read(*this, s, output);
  }
  
  template <typename T>
  template <typename U>
  typename std::enable_if<!std::is_base_of<Serializable, U>::value, void>::type
  Reader<T>::read(const std::string &s, U &output)
  {
    upcast->readDefault(s, output);
  }

  template <typename T>
  template <typename U>
  void Reader<T>::read(const std::string &s, iScalar<U> &output)
  {
    typename TensorToVec<iScalar<U>>::type v;

    upcast->readDefault(s, v);
    vecToTensor(output, v);
  }

  template <typename T>
  template <typename U, int N>
  void Reader<T>::read(const std::string &s, iVector<U, N> &output)
  {
    typename TensorToVec<iVector<U, N>>::type v;
    
    upcast->readDefault(s, v);
    vecToTensor(output, v);
  }
  
  template <typename T>
  template <typename U, int N>
  void Reader<T>::read(const std::string &s, iMatrix<U, N> &output)
  {
    typename TensorToVec<iMatrix<U, N>>::type v;
    
    upcast->readDefault(s, v);
    vecToTensor(output, v);
  }

  template <typename T>
  template <typename U>
  void Reader<T>::fromString(U &output, const std::string &s)
  {
    std::istringstream is(s);
    
    is.exceptions(std::ios::failbit);
    try
    {
      is >> std::boolalpha >> output;
    }
    catch(std::ios_base::failure &e)
    {
      std::cerr << "numerical conversion failure on '" << s << "' ";
      std::cerr << "(typeid: " << typeid(U).name() << ")" << std::endl;
      abort();
    }
  }

  // serializable base class ///////////////////////////////////////////////////
  class Serializable
  {
  public:
    template <typename T>
    static inline void write(Writer<T> &WR,const std::string &s,
                             const Serializable &obj)
    {}
    
    template <typename T>
    static inline void read(Reader<T> &RD,const std::string &s,
                            Serializable &obj)
    {}
    
    friend inline std::ostream & operator<<(std::ostream &os,
                                            const Serializable &obj)
    {
      return os;
    }

    template <typename T>
    static inline typename std::enable_if<!std::is_base_of<Eigen::TensorBase<T, Eigen::ReadOnlyAccessors>, T>::value, bool>::type
    CompareMember(const T &lhs, const T &rhs) {
      return lhs == rhs;
    }

    template <typename T>
    static inline typename std::enable_if<std::is_base_of<Eigen::TensorBase<T, Eigen::ReadOnlyAccessors>, T>::value, bool>::type
    CompareMember(const T &lhs, const T &rhs) {
      Eigen::Tensor<bool, 0, T::Options> bResult = (lhs == rhs).all();
      return bResult(0);
    }

    template <typename T>
    static inline typename std::enable_if<std::is_base_of<Eigen::TensorBase<T, Eigen::ReadOnlyAccessors>, T>::value, bool>::type
    CompareMember(const std::vector<T> &lhs, const std::vector<T> &rhs) {
      const auto NumElements{lhs.size()};
      bool bResult = ( NumElements == rhs.size() );
      for( auto i = 0 ; i < NumElements && bResult ; i++ ) {
        Eigen::Tensor<bool, 0, T::Options> b = (lhs[i] == rhs[i]).all();
        bResult = b(0);
      }
      return bResult;
    }

    template <typename T>
    static inline typename std::enable_if<!std::is_base_of<Eigen::TensorBase<T, Eigen::ReadOnlyAccessors>, T>::value, void>::type
    WriteMember(std::ostream &os, const T &object) {
      os << object;
    }
    
    template <typename T>
    static inline typename std::enable_if<std::is_base_of<Eigen::TensorBase<T, Eigen::ReadOnlyAccessors>, T>::value, void>::type
    WriteMember(std::ostream &os, const T &object) {
      os << "Eigen::Tensor";
    }
  };

  // Generic writer interface //////////////////////////////////////////////////
  template <typename T>
  inline void push(Writer<T> &w, const std::string &s) {
    w.push(s);
  }
  
  template <typename T>
  inline void push(Writer<T> &w, const char *s)
  {
    w.push(std::string(s));
  }
  
  template <typename T>
  inline void pop(Writer<T> &w)
  {
    w.pop();
  }
  
  template <typename T, typename U>
  inline void write(Writer<T> &w, const std::string& s, const U &output)
  {
    w.write(s, output);
  }
  
  // Generic reader interface //////////////////////////////////////////////////
  template <typename T>
  inline bool push(Reader<T> &r, const std::string &s)
  {
    return r.push(s);
  }
  
  template <typename T>
  inline bool push(Reader<T> &r, const char *s)
  {
    return r.push(std::string(s));
  }
  
  template <typename T>
  inline void pop(Reader<T> &r)
  {
    r.pop();
  }
  
  template <typename T, typename U>
  inline void read(Reader<T> &r, const std::string &s, U &output)
  {
    r.read(s, output);
  }
}

#endif
