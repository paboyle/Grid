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
  // TODO Support Grid::complex from GPU port
  template<typename T> using Grid_complex = std::complex<T>;

  // Returns original type, except for Grid_complex, where it returns the underlying type
  template<typename T> struct RealType { using type = T; };
  template<typename T> struct RealType<Grid_complex<T>> { using type = T; };

  namespace EigenIO {
    template<typename T> struct is_complex : public std::false_type {};
    template<typename T> struct is_complex<Grid_complex<T>>
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
    template<typename T, typename C = void> struct is_tensor_of_scalar : public std::false_type {};
    template<typename T> struct is_tensor_of_scalar<T, typename std::enable_if<is_tensor<T>::value && is_scalar<typename T::Scalar>::value, void>::type> : public std::true_type {};

    // Is this an Eigen tensor of a supported container
    template<typename T, typename C = void> struct is_tensor_of_container : public std::false_type {};
    template<typename T> struct is_tensor_of_container<T, typename std::enable_if<is_tensor<T>::value && is_container<typename T::Scalar>::value, void>::type> : public std::true_type {};

    // Is this a fixed-size Eigen tensor
    template<typename T, typename C = void> struct is_tensor_fixed : public std::false_type {};
    template<typename Scalar_, typename Dimensions_, int Options_, typename IndexType>
    struct is_tensor_fixed<Eigen::TensorFixedSize<Scalar_, Dimensions_, Options_, IndexType>>
        : public std::true_type {};
    template<typename Scalar_, typename Dimensions_, int Options_, typename IndexType,
              int MapOptions_, template <class> class MapPointer_>
    struct is_tensor_fixed<Eigen::TensorMap<Eigen::TensorFixedSize<Scalar_, Dimensions_,
                                            Options_, IndexType>, MapOptions_, MapPointer_>>
        : public std::true_type {};

    // Is this a variable-size Eigen tensor
    template<typename T, typename C = void> struct is_tensor_variable : public std::false_type {};
    template<typename T> struct is_tensor_variable<T, typename std::enable_if<is_tensor<T>::value
        && !is_tensor_fixed<T>::value, void>::type> : public std::true_type {};

    // These traits describe the Eigen tensor scalar and container objects supported for IO
    // Containers arbitrarily deeply nested compositions of fixed size objects:
    // ... grid tensors (iScalar, iVector, and iMatrix) and fixed size arrays
    template <typename T, typename C = void> struct Traits {}; // C needed for specialisation
    // This defines the bottom level - i.e. it's a description of the underlying scalar
    template <typename T> struct Traits<T, typename std::enable_if<is_scalar<T>::value, void>::type> {
      using scalar_type = T; // Type of the underlying scalar
      using scalar_real = typename RealType<scalar_type>::type; // real type underlying scalar_type
      static constexpr unsigned int rank = 0;   // The rank of the grid tensor (i.e. how many indices used)
      static constexpr unsigned int rank_non_trivial = 0; // As per rank, but excludes those of dimension 1
      static constexpr unsigned int count = 1;  // total number of elements (i.e. product of dimensions)
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
      using scalar_type = typename Traits<T>::scalar_type;
      using scalar_real = typename RealType<scalar_type>::type;
      static constexpr unsigned int rank = 1 + Traits<T>::rank;
      static constexpr unsigned int rank_non_trivial = 0 + Traits<T>::rank_non_trivial;
      static constexpr unsigned int count = 1 * Traits<T>::count;
      static constexpr std::size_t scalar_size = Traits<T>::scalar_size;
      static constexpr std::size_t size = scalar_size * count;
      static constexpr std::size_t Dimension(unsigned int dim) {
        return ( dim == 0 ) ? 1 : Traits<T>::Dimension(dim - 1); }
      static constexpr std::size_t DimensionNT(unsigned int dim) {
        return Traits<T>::DimensionNT(dim); }
    };
    template <typename T, int N> struct Traits<iVector<T, N>> {
      using scalar_type = typename Traits<T>::scalar_type;
      using scalar_real = typename RealType<scalar_type>::type;
      static constexpr unsigned int rank = 1 + Traits<T>::rank;
      static constexpr unsigned int rank_non_trivial = (N>1 ? 1 : 0) + Traits<T>::rank_non_trivial;
      static constexpr unsigned int count = N * Traits<T>::count;
      static constexpr std::size_t scalar_size = Traits<T>::scalar_size;
      static constexpr std::size_t size = scalar_size * count;
      static constexpr std::size_t Dimension(unsigned int dim) {
        return ( dim == 0 ) ? N : Traits<T>::Dimension(dim - 1); }
      static constexpr std::size_t DimensionNT(unsigned int dim) {
        return ( N == 1 ) ? Traits<T>::DimensionNT(dim) : ( dim == 0 ) ? N : Traits<T>::DimensionNT(dim - 1);
      }
    };
    template <typename T, int N> struct Traits<iMatrix<T, N>> {
      using scalar_type = typename Traits<T>::scalar_type;
      using scalar_real = typename RealType<scalar_type>::type;
      static constexpr unsigned int rank = 2 + Traits<T>::rank;
      static constexpr unsigned int rank_non_trivial = (N>1 ? 2 : 0) + Traits<T>::rank_non_trivial;
      static constexpr unsigned int count = N * N * Traits<T>::count;
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
      if constexpr ( EigenIO::is_scalar<Scalar>::value ) {
        lambda( * pScalar, Seq++, MyIndex );
      } else {
        for( typename EigenIO::Traits<Scalar>::scalar_type &Source : * pScalar ) {
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

  // Sequential initialisation of tensors
  // Would have preferred to define template variables for this, but that's c++ 17
  template <typename ETensor>
  typename std::enable_if<EigenIO::is_tensor<ETensor>::value && !EigenIO::is_complex<typename EigenIO::Traits<typename ETensor::Scalar>::scalar_type>::value, void>::type
  SequentialInit( ETensor &ET, typename EigenIO::Traits<typename ETensor::Scalar>::scalar_type Inc = 1 )
  {
    using Traits = EigenIO::Traits<typename ETensor::Scalar>;
    using scalar_type = typename Traits::scalar_type;
    for_all( ET, [&](scalar_type &c, typename ETensor::Index n, const std::array<size_t, ETensor::NumIndices + Traits::rank_non_trivial> &Dims ) {
      c = Inc * static_cast<scalar_type>(n);
    } );
  }

  template <typename ETensor>
  typename std::enable_if<EigenIO::is_tensor<ETensor>::value && EigenIO::is_complex<typename EigenIO::Traits<typename ETensor::Scalar>::scalar_type>::value, void>::type
  SequentialInit( ETensor &ET, typename EigenIO::Traits<typename ETensor::Scalar>::scalar_type Inc={1,-1})
  {
    using Traits = EigenIO::Traits<typename ETensor::Scalar>;
    using scalar_type = typename Traits::scalar_type;
    for_all( ET, [&](scalar_type &c, typename ETensor::Index n, const std::array<size_t, ETensor::NumIndices + Traits::rank_non_trivial> &Dims ) {
      c = Inc * static_cast<typename RealType<scalar_type>::type>(n);
    } );
  }
  
  // Helper to dump a tensor
#ifdef DEBUG
#define dump_tensor(args...) dump_tensor_func(args)
  template <typename T>
  typename std::enable_if<EigenIO::is_tensor<T>::value, void>::type
  dump_tensor_func(T &t, const char * pName = nullptr)
  {
    using Traits = typename EigenIO::Traits<typename T::Scalar>;
    const auto rank{T::NumIndices};
    const auto &dims = t.dimensions();
    std::cout << "Dumping rank " << rank << ((T::Options & Eigen::RowMajor) ? ", row" : ", column") << "-major tensor ";
    if( pName )
      std::cout << pName;
    for( auto i = 0 ; i < rank; i++ ) std::cout << "[" << dims[i] << "]";
    std::cout << " in memory order:" << std::endl;
    for_all( t, [&](typename Traits::scalar_type &c, typename T::Index index, const std::array<size_t, T::NumIndices + Traits::rank_non_trivial> &Dims ){
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
                         && !EigenIO::is_tensor<U>::value, void>::type
    write(const std::string& s, const U &output);
    template <typename U>
    void write(const std::string &s, const iScalar<U> &output);
    template <typename U, int N>
    void write(const std::string &s, const iVector<U, N> &output);
    template <typename U, int N>
    void write(const std::string &s, const iMatrix<U, N> &output);
    template <typename ETensor>
    typename std::enable_if<EigenIO::is_tensor<ETensor>::value, void>::type
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
    typename std::enable_if<!std::is_base_of<Serializable, U>::value
                         && !EigenIO::is_tensor<U>::value, void>::type
    read(const std::string& s, U &output);
    template <typename U>
    void read(const std::string &s, iScalar<U> &output);
    template <typename U, int N>
    void read(const std::string &s, iVector<U, N> &output);
    template <typename U, int N>
    void read(const std::string &s, iMatrix<U, N> &output);
    template <typename ETensor>
    typename std::enable_if<EigenIO::is_tensor<ETensor>::value, void>::type
    read(const std::string &s, ETensor &output);
    template <typename ETensor>
    typename std::enable_if<EigenIO::is_tensor_fixed<ETensor>::value, void>::type
    Reshape(ETensor &t, const std::array<typename ETensor::Index, ETensor::NumDimensions> &dims );
    template <typename ETensor>
    typename std::enable_if<EigenIO::is_tensor_variable<ETensor>::value, void>::type
    Reshape(ETensor &t, const std::array<typename ETensor::Index, ETensor::NumDimensions> &dims );
    /*template <typename ETensor>
    typename std::enable_if<EigenIO::is_tensor_fixed<ETensor>::value, std::size_t>::type
    DimSize(ETensor &t, std::size_t dim );
    template <typename ETensor>
    typename std::enable_if<EigenIO::is_tensor_variable<ETensor>::value, std::size_t>::type
    DimSize(ETensor &t, std::size_t dim );*/
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
                       && !EigenIO::is_tensor<U>::value, void>::type
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
  
  // Eigen::Tensors of Grid tensors (iScalar, iVector, iMatrix)
  template <typename T>
  template <typename ETensor/*, typename U, int N*/>
  typename std::enable_if<EigenIO::is_tensor<ETensor>::value, void>::type
  Writer<T>::write(const std::string &s, const ETensor &output)
  {
    using Index = typename ETensor::Index;
    using Container = typename ETensor::Scalar; // NB: could be same as Scalar
    using Traits = EigenIO::Traits<Container>;
    using Scalar = typename Traits::scalar_type; // type of the underlying scalar
    constexpr unsigned int TensorRank{ETensor::NumIndices};
    constexpr unsigned int ContainerRank{Traits::rank}; // Only non-zero for containers
    constexpr unsigned int TotalRank{TensorRank + ContainerRank};
    const Index NumElements{output.size()};
    assert( NumElements > 0 );

    // Get the dimensionality of the tensor
    std::vector<std::size_t>  TotalDims(TotalRank);
    for(auto i = 0; i < TensorRank; i++ ) {
      auto dim = output.dimension(i);
      TotalDims[i] = static_cast<size_t>(dim);
      assert( TotalDims[i] == dim ); // check we didn't lose anything in the conversion
    }
    for(auto i = 0; i < ContainerRank; i++ )
      TotalDims[TensorRank + i] = Traits::Dimension(i);

    // If the Tensor isn't in Row-Major order, then we'll need to copy it's data
    const bool CopyData{NumElements > 1 && ETensor::Layout != Eigen::StorageOptions::RowMajor};
    const Scalar * pWriteBuffer;
    Scalar * pCopyBuffer = nullptr;
    const Index TotalNumElements = NumElements * Traits::count;
    if( !CopyData ) {
      if constexpr ( ContainerRank == 0 )
        pWriteBuffer = output.data();
      else
        pWriteBuffer = output.data()->begin();
    } else {
      // Regardless of the Eigen::Tensor storage order, the copy will be Row Major
      pCopyBuffer = static_cast<Scalar *>(malloc(TotalNumElements * sizeof(Scalar)));
      pWriteBuffer = pCopyBuffer;
      Scalar * pCopy = pCopyBuffer;
      std::array<Index, TensorRank> MyIndex;
      for( auto &idx : MyIndex ) idx = 0;
      for( auto n = 0; n < NumElements; n++ ) {
        if constexpr ( ContainerRank == 0 )
          * pCopy ++ = output( MyIndex );
        else {
          for( const Scalar &Source : output( MyIndex ) )
            * pCopy ++ = Source;
        }
        // Now increment the index
        for( int i = output.NumDimensions - 1; i >= 0 && ++MyIndex[i] == output.dimension(i); i-- )
          MyIndex[i] = 0;
      }
    }
    upcast->template writeMultiDim<Scalar>(s, TotalDims, pWriteBuffer, TotalNumElements);
    if( pCopyBuffer ) free( pCopyBuffer );
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
  typename std::enable_if<!std::is_base_of<Serializable, U>::value
                       && !EigenIO::is_tensor<U>::value, void>::type
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
  template <typename ETensor>
  typename std::enable_if<EigenIO::is_tensor<ETensor>::value, void>::type
  Reader<T>::read(const std::string &s, ETensor &output)
  {
    // alias to element type
    using Container = typename ETensor::Scalar;
    using Traits = EigenIO::Traits<Container>;
    using Scalar = typename Traits::scalar_type;

    // read the (flat) data and dimensionality
    std::vector<std::size_t>   dimData;
    std::vector<Scalar>  buf;
    upcast->readMultiDim( s, buf, dimData );
    // Make sure that the number of elements read matches dimensions read
    std::size_t NumElements = 1;
    std::size_t RankRequired = 0;
    std::vector<typename ETensor::Index>   dimNonTrivial;
    dimNonTrivial.reserve(dimData.size());
    for( auto d : dimData ) {
      NumElements *= d;
      if( d > 1 ) {
        RankRequired++;
        dimNonTrivial.push_back(d);
      }
    }
    assert( NumElements == buf.size() && "Number of elements read back <> product of dimensions" );
    // If our scalar object is a Container, make sure it's dimensions match what we read back
    const auto InnerRank{Traits::rank_non_trivial};
    if ( InnerRank > 0 ) {
      assert( RankRequired >= InnerRank && "Tensor Container too complex for data" );
      for( auto i = InnerRank - 1 ; i != -1 ; i-- ) {
        auto d = dimNonTrivial[--RankRequired];
        assert( d == Traits::DimensionNT(i) && "Tensor Container dimensions don't match data" );
        NumElements /= d;
        dimNonTrivial.pop_back();
      }
    }
    // Make sure our object has the right rank
    assert( ETensor::NumDimensions >= RankRequired );
    bool bShapeOK = true;
    std::size_t RankNonTrivial = 0;
    const auto & dims{output.dimensions()};
    using ETDims = std::array<typename ETensor::Index, ETensor::NumDimensions>;
    ETDims dimsNew;
    // Make sure fixed dimension objects have allocated memory
    /*if constexpr( EigenIO::is_tensor_fixed<ETensor>::value ) {
      for( auto &d : dimsNew ) d = 0;
      output( dimsNew ) = 0;
    }*/
    for( auto i = 0, j = 0 ; bShapeOK && i < ETensor::NumDimensions ; i++ ) {
      auto d = dims[i];
      if( d < 1 )
        bShapeOK = false;
      else if( d > 1 ) {
        RankNonTrivial++;
        if( d != dimNonTrivial[j] )
          bShapeOK = false;
        j++;
      }
      dimsNew[i] = d;
    }
    //if( RankNonTrivial == 0 ) RankNonTrivial++;
    // Make the tensor the same size as the data read
    if ( !bShapeOK || RankNonTrivial != RankRequired ) {
      for( auto i = 0 ; i < ETensor::NumDimensions ; i++ )
        dimsNew[i] = ( i < RankRequired ) ? dimNonTrivial[i] : 1;
      Reshape(output, dimsNew);
    }
    // Copy the data into the tensor
    ETDims MyIndex;
    for( auto &d : MyIndex ) d = 0;
    std::size_t idx = 0;
    for( auto n = 0 ; n < NumElements ; n++ ) {
      Container & c = output( MyIndex );
      if constexpr ( EigenIO::is_scalar<Container>::value ) {
        c = buf[idx++];
      } else {
        for( Scalar & s : c )
          s = buf[idx++];
      }
      // Now increment the index
      for( int i = ETensor::NumDimensions - 1; i >= 0 && ++MyIndex[i] == dims[i]; i-- )
        MyIndex[i] = 0;
    }
  }

  template <typename T>
  template <typename ETensor>
  typename std::enable_if<EigenIO::is_tensor_fixed<ETensor>::value, void>::type
  Reader<T>::Reshape(ETensor &t, const std::array<typename ETensor::Index, ETensor::NumDimensions> &dims )
  {
    assert( 0 && "EigenIO: Fixed tensor dimensions can't be changed" );
  }

  template <typename T>
  template <typename ETensor>
  typename std::enable_if<EigenIO::is_tensor_variable<ETensor>::value, void>::type
  Reader<T>::Reshape(ETensor &t, const std::array<typename ETensor::Index, ETensor::NumDimensions> &dims )
  {
    //t.reshape( dims );
    t.resize( dims );
  }

  /*template <typename T>
  template <typename ETensor>
  typename std::enable_if<EigenIO::is_tensor_fixed<ETensor>::value, std::size_t>::type
  Reader<T>::DimSize(ETensor &t, std::size_t dim )
  {
    return 0;//ETensor::Dimension[dim];
  }

  template <typename T>
  template <typename ETensor>
  typename std::enable_if<EigenIO::is_tensor_variable<ETensor>::value, std::size_t>::type
  Reader<T>::DimSize(ETensor &t, std::size_t dim )
  {
    return t.dimension(dim);
  }*/

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

    template <typename T1, typename T2>
    static inline typename std::enable_if<!EigenIO::is_tensor<T1>::value || !EigenIO::is_tensor<T2>::value, bool>::type
    CompareMember(const T1 &lhs, const T2 &rhs) {
      return lhs == rhs;
    }

    template <typename T1, typename T2>
    static inline typename std::enable_if<EigenIO::is_tensor<T1>::value && EigenIO::is_tensor<T2>::value, bool>::type
    CompareMember(const T1 &lhs, const T2 &rhs) {
      // First check whether dimensions match (Eigen tensor library will assert if they don't match)
      bool bReturnValue = (T1::NumIndices == T2::NumIndices);
      for( auto i = 0 ; bReturnValue && i < T1::NumIndices ; i++ )
          bReturnValue = ( lhs.dimension(i)) == rhs.dimension(i);
      if( bReturnValue ) {
        Eigen::Tensor<bool, 0, T1::Options> bResult = (lhs == rhs).all();
        bReturnValue = bResult(0);
      }
      return bReturnValue;
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
