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
  // Abstract writer/reader classes ////////////////////////////////////////////
  // static polymorphism implemented using CRTP idiom
  class Serializable;

  // which types are supported scalar types for Eigen::Tensor
  template<typename T> struct is_eigen_tensor_scalar : std::integral_constant<bool,
  std::is_arithmetic<T>::value || Grid::is_complex<T>::value> {};

  // Helper to allow iteration through an Eigen::Tensor (using a lambda)
  template <typename ETensor, typename Lambda>
  typename std::enable_if<std::is_base_of<Eigen::TensorBase<ETensor, Eigen::ReadOnlyAccessors>, ETensor>::value && is_grid_tensor<typename ETensor::Scalar>::value, void>::type
  for_all( ETensor &ET, Lambda lambda )
  {
    using Scalar = typename ETensor::Scalar::scalar_type;
    const std::size_t NumElements = ET.size();
    assert( NumElements > 0 );
    if( NumElements == 1 ) {
      const auto MyRank{grid_tensor_att<typename ETensor::Scalar>::rank_non_trivial};
      std::vector<std::size_t> SubIndex(MyRank);
      for( auto &idx : SubIndex ) idx = 0;
      typename ETensor::Index n = 0;
      for( Scalar &Source : * ET.data() ) {
        lambda(Source, n++, &SubIndex[0] );
        // Now increment SubIndex
        for( auto i = MyRank - 1; i >= 0 && ++SubIndex[i] == 11/*ReducedDims[i]*/; i-- )
          SubIndex[i] = 0;
      }
    }
    else {
      // We're only interested in non-trivial dimensions (i.e. dimensions > 1)
      unsigned int TrivialDimCount{0};
      std::vector<size_t> ReducedDims;
      ReducedDims.reserve(ET.NumDimensions + grid_tensor_att<typename ETensor::Scalar>::rank_non_trivial); // Make sure we only do one malloc
      for(auto i = 0; i < ET.NumDimensions; i++ ) {
        auto dim = ET.dimension(i);
        if( dim <= 1 ) {
          TrivialDimCount++;
          assert( dim == 1 ); // Not expecting dimension to be <= 0
        } else {
          size_t s = static_cast<size_t>(dim);
          assert( s == dim ); // check we didn't lose anything in the conversion
          ReducedDims.push_back(s);
        }
      }
      // NB: NumElements > 1 implies this is not a scalar, so some dims should be left
      assert( ET.NumDimensions > TrivialDimCount );
      // Now add the extra dimensions, based on object zero
      typename TensorToVec<typename ETensor::Scalar>::type ttv = tensorToVec(* ET.data());
      Flatten<typename TensorToVec<typename ETensor::Scalar>::type> f(ttv);
      const std::vector<size_t> & ExtraDims{f.getDim()};
      assert(ExtraDims.size() == grid_tensor_att<typename ETensor::Scalar>::rank_non_trivial);
      size_t ExtraCount{1};
      for( auto i : ExtraDims ) {
        assert( i > 0 );
        ExtraCount *= i;
        ReducedDims.push_back(i);
      }
      assert(grid_tensor_att<typename ETensor::Scalar>::count == ExtraCount);
      assert(grid_tensor_att<typename ETensor::Scalar>::size  == sizeof( typename ETensor::Scalar ));
      const unsigned int ReducedDimsSize = static_cast<unsigned int>(ReducedDims.size());
      assert( ReducedDimsSize == ReducedDims.size() );
      const typename ETensor::Index TotalNumElements = NumElements * ExtraCount;
      std::array<typename ETensor::Index, ETensor::NumIndices> MyIndex;
      for( auto &idx : MyIndex ) idx = 0;
      std::vector<std::size_t> SubIndex(ReducedDimsSize);
      for( auto &idx : SubIndex ) idx = 0;
      for( typename ETensor::Index n = 0; n < TotalNumElements; ) {
        for( Scalar &Source : ET( MyIndex ) ) {
          lambda(Source, n++, &SubIndex[0] );
          // Now increment MyIndex
          for( auto i = ET.NumDimensions - 1; i >= 0 && ++MyIndex[i] == ET.dimension(i); i-- )
            MyIndex[i] = 0;
          // Now increment SubIndex
          for( auto i = ReducedDimsSize - 1; i >= 0 && ++SubIndex[i] == ReducedDims[i]; i-- )
            SubIndex[i] = 0;
        }
      }
    }
  }

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
    typename std::enable_if<std::is_base_of<Eigen::TensorBase<ETensor, Eigen::ReadOnlyAccessors>, ETensor>::value && is_eigen_tensor_scalar<typename ETensor::Scalar>::value, void>::type
    write(const std::string &s, const ETensor &output);
    template <typename ETensor>
    typename std::enable_if<std::is_base_of<Eigen::TensorBase<ETensor, Eigen::ReadOnlyAccessors>, ETensor>::value && is_grid_tensor<typename ETensor::Scalar>::value, void>::type
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
  typename std::enable_if<std::is_base_of<Eigen::TensorBase<ETensor, Eigen::ReadOnlyAccessors>, ETensor>::value && is_eigen_tensor_scalar<typename ETensor::Scalar>::value, void>::type
  Writer<T>::write(const std::string &s, const ETensor &output)
  {
    const typename ETensor::Index NumElements{output.size()};
    assert( NumElements > 0 );
    if( NumElements == 1 )
      upcast->writeDefault(s, * output.data());
    else {
      // We're only interested in non-trivial dimensions (i.e. dimensions > 1)
      unsigned int TrivialDimCount{0};
      std::vector<size_t> ReducedDims;
      ReducedDims.reserve(output.NumDimensions); // Make sure we only do one malloc
      for(auto i = 0; i < output.NumDimensions; i++ ) {
        auto dim = output.dimension(i);
        if( dim <= 1 ) {
          TrivialDimCount++;
          assert( dim == 1 ); // Not expecting dimension to be <= 0
        } else {
          size_t s = static_cast<size_t>(dim);
          assert( s == dim ); // check we didn't lose anything in the conversion
          ReducedDims.push_back(s);
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
      upcast->template writeMultiDim<typename ETensor::Scalar>(s, ReducedDims, pWriteBuffer, NumElements);
      if( pCopyBuffer ) free( pCopyBuffer );
    }
  }
  
  // Eigen::Tensors of Grid tensors (iScalar, iVector, iMatrix)
  template <typename T>
  template <typename ETensor/*, typename U, int N*/>
  typename std::enable_if<std::is_base_of<Eigen::TensorBase<ETensor, Eigen::ReadOnlyAccessors>, ETensor>::value && is_grid_tensor<typename ETensor::Scalar>::value, void>::type
  Writer<T>::write(const std::string &s, const ETensor &output)
  {
    const typename ETensor::Index NumElements{output.size()};
    assert( NumElements > 0 );
    if( NumElements == 1 )
      upcast->writeDefault(s, tensorToVec(* output.data()));
    else {
      // We're only interested in non-trivial dimensions (i.e. dimensions > 1)
      unsigned int TrivialDimCount{0};
      std::vector<size_t> ReducedDims;
      ReducedDims.reserve(output.NumDimensions + grid_tensor_att<typename ETensor::Scalar>::rank_non_trivial); // Make sure we only do one malloc
      for(auto i = 0; i < output.NumDimensions; i++ ) {
        auto dim = output.dimension(i);
        if( dim <= 1 ) {
          TrivialDimCount++;
          assert( dim == 1 ); // Not expecting dimension to be <= 0
        } else {
          size_t s = static_cast<size_t>(dim);
          assert( s == dim ); // check we didn't lose anything in the conversion
          ReducedDims.push_back(s);
        }
      }
      // NB: NumElements > 1 implies this is not a scalar, so some dims should be left
      assert( output.NumDimensions > TrivialDimCount );
      // Now add the extra dimensions, based on object zero
      typename TensorToVec<typename ETensor::Scalar>::type ttv = tensorToVec(* output.data());
      Flatten<typename TensorToVec<typename ETensor::Scalar>::type> f(ttv);
      const std::vector<size_t> & ExtraDims{f.getDim()};
      assert(ExtraDims.size() == grid_tensor_att<typename ETensor::Scalar>::rank_non_trivial);
      size_t ExtraCount{1};
      for( auto i : ExtraDims ) {
        assert( i > 0 );
        ExtraCount *= i;
        ReducedDims.push_back(i);
      }
      assert(grid_tensor_att<typename ETensor::Scalar>::count == ExtraCount);
      assert(grid_tensor_att<typename ETensor::Scalar>::size  == sizeof( typename ETensor::Scalar ));
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
      upcast->template writeMultiDim<Scalar>(s, ReducedDims, pWriteBuffer, TotalNumElements);
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
