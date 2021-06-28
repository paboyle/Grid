/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: ./Grid/serialisation/VectorUtils.h
 
 Copyright (C) 2015
 
 Author: Antonin Portelli <antonin.portelli@me.com>
 Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 Author: paboyle <paboyle@ph.ed.ac.uk>
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
#ifndef GRID_SERIALISATION_VECTORUTILS_H
#define GRID_SERIALISATION_VECTORUTILS_H

#include <type_traits>
#include <Grid/tensors/Tensors.h>

namespace Grid {
  // Pair IO utilities /////////////////////////////////////////////////////////
  // helper function to parse input in the format "<obj1 obj2>"
  template <typename T1, typename T2>
  inline std::istream & operator>>(std::istream &is, std::pair<T1, T2> &buf)
  {
    T1 buf1;
    T2 buf2;
    char c;

    // Search for "pair" delimiters.
    do
    {
      is.get(c);
    } while (c != '(' && !is.eof());
    if (c == '(')
    {
      int start = is.tellg();
      do
      {
        is.get(c);
      } while (c != ')' && !is.eof());
      if (c == ')')
      {
        int end = is.tellg();
        int psize = end - start - 1;

        // Only read data between pair limiters.
        is.seekg(start);
        std::string tmpstr(psize, ' ');
        is.read(&tmpstr[0], psize);
        std::istringstream temp(tmpstr);
        temp >> buf1 >> buf2;
        buf = std::make_pair(buf1, buf2);
        is.seekg(end);
      }
    }
    is.peek();
    return is;
  }
  
  // output to streams for pairs
  template <class T1, class T2>
  inline std::ostream & operator<<(std::ostream &os, const std::pair<T1, T2> &p)
  {
    os << "(" << p.first << " " << p.second << ")";
    return os;
  }
  
  // std::vector<std:vector<...>> nested to specified Rank //////////////////////////////////
  template<typename T, unsigned int Rank>
  struct NestedStdVector {
    typedef typename std::vector<typename NestedStdVector<T, Rank - 1>::type> type;
  };
  
  template<typename T>
  struct NestedStdVector<T,0> {
    typedef T type;
  };
  
  // Grid scalar tensors to nested std::vectors //////////////////////////////////
  template <typename T>
  struct TensorToVec
  {
    typedef T type;
  };

  template <typename T>
  struct TensorToVec<iScalar<T>>
  {
    typedef typename TensorToVec<T>::type type;
  };

  template <typename T, int N>
  struct TensorToVec<iVector<T, N>>
  {
    typedef typename std::vector<typename TensorToVec<T>::type> type;
  };

  template <typename T, int N>
  struct TensorToVec<iMatrix<T, N>>
  {
    typedef typename std::vector<std::vector<typename TensorToVec<T>::type>> type;
  };

  template <typename T>
  void tensorDim(std::vector<size_t> &dim, const T &t, const bool wipe = true)
  {
    if (wipe)
    {
      dim.clear();
    }
  }

  template <typename T>
  void tensorDim(std::vector<size_t> &dim, const iScalar<T> &t, const bool wipe = true)
  {
    if (wipe)
    {
      dim.clear();
    }
    tensorDim(dim, t._internal, false);
  }

  template <typename T, int N>
  void tensorDim(std::vector<size_t> &dim, const iVector<T, N> &t, const bool wipe = true)
  {
    if (wipe)
    {
      dim.clear();
    }
    dim.push_back(N);
    tensorDim(dim, t._internal[0], false);
  }

  template <typename T, int N>
  void tensorDim(std::vector<size_t> &dim, const iMatrix<T, N> &t, const bool wipe = true)
  {
    if (wipe)
    {
      dim.clear();
    }
    dim.push_back(N);
    dim.push_back(N);
    tensorDim(dim, t._internal[0][0], false);
  }

  template <typename T>
  typename TensorToVec<T>::type tensorToVec(const T &t)
  {
    return t;
  }

  template <typename T>
  typename TensorToVec<iScalar<T>>::type tensorToVec(const iScalar<T>& t)
  {
    return tensorToVec(t._internal);
  }

  template <typename T, int N>
  typename TensorToVec<iVector<T, N>>::type tensorToVec(const iVector<T, N>& t)
  {
    typename TensorToVec<iVector<T, N>>::type v;

    v.resize(N);
    for (unsigned int i = 0; i < N; i++) 
    {
      v[i] = tensorToVec(t._internal[i]);
    }

    return v;
  }

  template <typename T, int N>
  typename TensorToVec<iMatrix<T, N>>::type tensorToVec(const iMatrix<T, N>& t)
  {
    typename TensorToVec<iMatrix<T, N>>::type v;

    v.resize(N);
    for (unsigned int i = 0; i < N; i++)
    {
      v[i].resize(N);
      for (unsigned int j = 0; j < N; j++) 
      {
        v[i][j] = tensorToVec(t._internal[i][j]);
      }
    }

    return v;
  }

  template <typename T>
  void vecToTensor(T &t, const typename TensorToVec<T>::type &v)
  {
    t = v;
  }


  template <typename T>
  void vecToTensor(iScalar<T> &t, const typename TensorToVec<iScalar<T>>::type &v)
  {
    vecToTensor(t._internal, v);
  }

  template <typename T, int N>
  void vecToTensor(iVector<T, N> &t, const typename TensorToVec<iVector<T, N>>::type &v)
  {
    for (unsigned int i = 0; i < N; i++) 
    {
      vecToTensor(t._internal[i], v[i]);
    }
  }

  template <typename T, int N>
  void vecToTensor(iMatrix<T, N> &t, const typename TensorToVec<iMatrix<T, N>>::type &v)
  {
    for (unsigned int i = 0; i < N; i++)
    for (unsigned int j = 0; j < N; j++)
    {
      vecToTensor(t._internal[i][j], v[i][j]);
    }
  }

  // is_flattenable<T>::value is true if T is a std::vector<> which can be flattened //////////////////////
  template <typename T, typename V = void>
  struct is_flattenable : std::false_type
  {
    using type      = T;
    using grid_type = T;
    static constexpr int vecRank = 0;
    static constexpr bool isGridTensor = false;
    static constexpr bool children_flattenable = std::is_arithmetic<T>::value or is_complex<T>::value;
  };

  template <typename T>
  struct is_flattenable<T, typename std::enable_if<isGridTensor<T>::value>::type> : std::false_type
  {
    using type      = typename GridTypeMapper<T>::scalar_type;
    using grid_type = T;
    static constexpr int vecRank = 0;
    static constexpr bool isGridTensor = true;
    static constexpr bool children_flattenable = true;
  };

  template <typename T>
  struct is_flattenable<std::vector<T>, typename std::enable_if<is_flattenable<T>::children_flattenable>::type>
  : std::true_type
  {
    using type      = typename is_flattenable<T>::type;
    using grid_type = typename is_flattenable<T>::grid_type;
    static constexpr bool isGridTensor = is_flattenable<T>::isGridTensor;
    static constexpr int vecRank = is_flattenable<T>::vecRank + 1;
    static constexpr bool children_flattenable = true;
  };
  
  // Vector flattening utility class ////////////////////////////////////////////
  // Class to flatten a multidimensional std::vector
  template <typename V>
  class Flatten
  {
  public:
    using Scalar  = typename is_flattenable<V>::type;
    static constexpr bool isGridTensor = is_flattenable<V>::isGridTensor;
  public:
    explicit                    Flatten(const V &vector);
    const V &                   getVector(void)     const { return vector_; }
    const std::vector<Scalar> & getFlatVector(void) const { return flatVector_; }
    const std::vector<size_t> & getDim(void)        const { return dim_; }
  private:
    template <typename W> typename std::enable_if<!is_flattenable<W>::value && !is_flattenable<W>::isGridTensor>::type
    accumulate(const W &e);
    template <typename W> typename std::enable_if<!is_flattenable<W>::value &&  is_flattenable<W>::isGridTensor>::type
    accumulate(const W &e);
    template <typename W> typename std::enable_if< is_flattenable<W>::value>::type
    accumulate(const W &v);
    template <typename W> typename std::enable_if<!is_flattenable<W>::value && !is_flattenable<W>::isGridTensor>::type
    accumulateDim(const W &e) {} // Innermost is a scalar - do nothing
    template <typename W> typename std::enable_if<!is_flattenable<W>::value &&  is_flattenable<W>::isGridTensor>::type
    accumulateDim(const W &e);
    template <typename W> typename std::enable_if< is_flattenable<W>::value>::type
    accumulateDim(const W &v);
  private:
    const V             &vector_;
    std::vector<Scalar> flatVector_;
    std::vector<size_t> dim_;
  };
  
  // Class to reconstruct a multidimensional std::vector
  template <typename V>
  class Reconstruct
  {
  public:
    using Scalar  = typename is_flattenable<V>::type;
    static constexpr bool isGridTensor = is_flattenable<V>::isGridTensor;
  public:
    Reconstruct(const std::vector<Scalar> &flatVector,
                const std::vector<size_t> &dim);
    const V &                   getVector(void)     const { return vector_; }
    const std::vector<Scalar> & getFlatVector(void) const { return flatVector_; }
    const std::vector<size_t> & getDim(void)        const { return dim_; }
  private:
    template <typename W> typename std::enable_if<!is_flattenable<W>::value && !is_flattenable<W>::isGridTensor>::type
    fill(W &v);
    template <typename W> typename std::enable_if<!is_flattenable<W>::value &&  is_flattenable<W>::isGridTensor>::type
    fill(W &v);
    template <typename W> typename std::enable_if< is_flattenable<W>::value>::type
    fill(W &v);
    template <typename W> typename std::enable_if< is_flattenable<W>::value &&  is_flattenable<W>::vecRank==1>::type
    resize(W &v, const unsigned int dim);
    template <typename W> typename std::enable_if< is_flattenable<W>::value && (is_flattenable<W>::vecRank>1)>::type
    resize(W &v, const unsigned int dim);
    template <typename W> typename std::enable_if<!is_flattenable<W>::isGridTensor>::type
    checkInnermost(const W &e) {} // Innermost is a scalar - do nothing
    template <typename W> typename std::enable_if< is_flattenable<W>::isGridTensor>::type
    checkInnermost(const W &e);
  private:
    V                         vector_;
    const std::vector<Scalar> &flatVector_;
    std::vector<size_t>       dim_;
    size_t                    ind_{0};
    unsigned int              dimInd_{0};
  };

  // Flatten class template implementation
  template <typename V>
  template <typename W> typename std::enable_if<!is_flattenable<W>::value && !is_flattenable<W>::isGridTensor>::type
  Flatten<V>::accumulate(const W &e)
  {
    flatVector_.push_back(e);
  }
  
  template <typename V>
  template <typename W> typename std::enable_if<!is_flattenable<W>::value && is_flattenable<W>::isGridTensor>::type
  Flatten<V>::accumulate(const W &e)
  {
    for (const Scalar &x: e) {
      flatVector_.push_back(x);
    }
  }

  template <typename V>
  template <typename W> typename std::enable_if<is_flattenable<W>::value>::type
  Flatten<V>::accumulate(const W &v)
  {
    for (auto &e: v)
    {
      accumulate(e);
    }
  }
  
  template <typename V>
  template <typename W> typename std::enable_if<!is_flattenable<W>::value && is_flattenable<W>::isGridTensor>::type
  Flatten<V>::accumulateDim(const W &e)
  {
    using Traits = GridTypeMapper<typename is_flattenable<W>::grid_type>;
    for (int rank=0; rank < Traits::Rank; ++rank)
      dim_.push_back(Traits::Dimension(rank));
  }
  
  template <typename V>
  template <typename W> typename std::enable_if<is_flattenable<W>::value>::type
  Flatten<V>::accumulateDim(const W &v)
  {
    dim_.push_back(v.size());
    accumulateDim(v[0]);
  }
  
  template <typename V>
  Flatten<V>::Flatten(const V &vector)
  : vector_(vector)
  {
    accumulateDim(vector_);
    std::size_t TotalSize{ dim_[0] };
    for (int i = 1; i < dim_.size(); ++i) {
      TotalSize *= dim_[i];
    }
    flatVector_.reserve(TotalSize);
    accumulate(vector_);
  }
  
  // Reconstruct class template implementation
  template <typename V>
  template <typename W> typename std::enable_if<!is_flattenable<W>::value && !is_flattenable<W>::isGridTensor>::type
  Reconstruct<V>::fill(W &v)
  {
    v = flatVector_[ind_++];
  }
  
  template <typename V>
  template <typename W> typename std::enable_if<!is_flattenable<W>::value &&  is_flattenable<W>::isGridTensor>::type
  Reconstruct<V>::fill(W &v)
  {
    for (auto &e: v)
    {
      e = flatVector_[ind_++];
    }
  }

  template <typename V>
  template <typename W> typename std::enable_if<is_flattenable<W>::value>::type
  Reconstruct<V>::fill(W &v)
  {
    for (auto &e: v)
    {
      fill(e);
    }
  }
  
  template <typename V>
  template <typename W> typename std::enable_if<is_flattenable<W>::value && is_flattenable<W>::vecRank==1>::type
  Reconstruct<V>::resize(W &v, const unsigned int dim)
  {
    v.resize(dim_[dim]);
  }
  
  template <typename V>
  template <typename W> typename std::enable_if<is_flattenable<W>::value && (is_flattenable<W>::vecRank>1)>::type
  Reconstruct<V>::resize(W &v, const unsigned int dim)
  {
    v.resize(dim_[dim]);
    for (auto &e: v)
    {
      resize(e, dim + 1);
    }
  }
  
  template <typename V>
  template <typename W> typename std::enable_if<is_flattenable<W>::isGridTensor>::type
  Reconstruct<V>::checkInnermost(const W &)
  {
    using Traits = GridTypeMapper<typename is_flattenable<W>::grid_type>;
    const int gridRank{Traits::Rank};
    const int dimRank{static_cast<int>(dim_.size())};
    assert(dimRank >= gridRank && "Tensor rank too low for Grid tensor");
    for (int i=0; i<gridRank; ++i) {
      assert(dim_[dimRank - gridRank + i] == Traits::Dimension(i) && "Tensor dimension doesn't match Grid tensor");
    }
    dim_.resize(dimRank - gridRank);
  }

  template <typename V>
  Reconstruct<V>::Reconstruct(const std::vector<Scalar> &flatVector,
                              const std::vector<size_t> &dim)
  : flatVector_(flatVector)
  , dim_(dim)
  {
    checkInnermost(vector_);
    assert(dim_.size() == is_flattenable<V>::vecRank && "Tensor rank doesn't match nested std::vector rank");
    resize(vector_, 0);
    fill(vector_);
  }
  
  // Vector IO utilities ///////////////////////////////////////////////////////
  // helper function to read space-separated values
  template <typename T>
  std::vector<T> strToVec(const std::string s)
  {
    std::istringstream sstr(s);
    std::vector<T>     v;
    
    for(T buf; sstr >> buf;)
    {
      v.push_back(buf);
    }
    
    return v;
  }
  
  // output to streams for vectors
  template < class T >
  inline std::ostream & operator<<(std::ostream &os, const std::vector<T> &v)
  {
    os << "[";
    for (unsigned int i = 0; i < v.size(); ++i)
    {
      os << v[i];
      if (i < v.size() - 1)
      {
        os << " ";
      }
    }
    os << "]";
    
    return os;
  }

  // In general, scalar types are considered "flattenable" (regularly shaped)
  template <typename T>
  bool isRegularShapeHelper(const std::vector<T> &, std::vector<std::size_t> &, int, bool)
  {
    return true;
  }

  template <typename T>
  bool isRegularShapeHelper(const std::vector<std::vector<T>> &v, std::vector<std::size_t> &Dims, int Depth, bool bFirst)
  {
    if( bFirst)
    {
      assert( Dims.size() == Depth     && "Bug: Delete this message after testing" );
      Dims.push_back(v[0].size());
      if (!Dims[Depth])
        return false;
    }
    else
    {
      assert( Dims.size() >= Depth + 1 && "Bug: Delete this message after testing" );
    }
    for (std::size_t i = 0; i < v.size(); ++i)
    {
      if (v[i].size() != Dims[Depth] || !isRegularShapeHelper(v[i], Dims, Depth + 1, bFirst && i==0))
      {
        return false;
      }
    }
    return true;
  }

  template <typename T>
  bool isRegularShape(const T &t) { return true; }

  template <typename T>
  bool isRegularShape(const std::vector<T> &v) { return !v.empty(); }

  // Return non-zero if all dimensions of this std::vector<std::vector<T>> are regularly shaped
  template <typename T>
  bool isRegularShape(const std::vector<std::vector<T>> &v)
  {
    if (v.empty() || v[0].empty())
      return false;
    // Make sure all of my rows are the same size
    std::vector<std::size_t> Dims;
    Dims.reserve(is_flattenable<T>::vecRank);
    Dims.push_back(v.size());
    Dims.push_back(v[0].size());
    for (std::size_t i = 0; i < Dims[0]; ++i)
    {
      if (v[i].size() != Dims[1] || !isRegularShapeHelper(v[i], Dims, 2, i==0))
      {
        return false;
      }
    }
    return true;
  }
}

// helper function to read space-separated values
template <typename T>
std::string vecToStr(const std::vector<T> &v)
{
  using Grid::operator<<;
  
  std::ostringstream sstr;

  sstr << v;

  return sstr.str();
}

#endif
