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

  // Vector element trait //////////////////////////////////////////////////////  
  template <typename T>
  struct element
  {
    typedef T type;
    static constexpr bool is_number = false;
  };
  
  template <typename T>
  struct element<std::vector<T>>
  {
    typedef typename element<T>::type type;
    static constexpr bool is_number = std::is_arithmetic<T>::value
                                      or is_complex<T>::value
                                      or element<T>::is_number;
  };
  
  // Vector flattening utility class ////////////////////////////////////////////
  // Class to flatten a multidimensional std::vector
  template <typename V>
  class Flatten
  {
  public:
    typedef typename element<V>::type Element;
  public:
    explicit                     Flatten(const V &vector);
    const V &                    getVector(void);
    const std::vector<Element> & getFlatVector(void);
    const std::vector<size_t>  & getDim(void);
  private:
    void accumulate(const Element &e);
    template <typename W>
    void accumulate(const W &v);
    void accumulateDim(const Element &e);
    template <typename W>
    void accumulateDim(const W &v);
  private:
    const V              &vector_;
    std::vector<Element> flatVector_;
    std::vector<size_t>  dim_;
  };
  
  // Class to reconstruct a multidimensional std::vector
  template <typename V>
  class Reconstruct
  {
  public:
    typedef typename element<V>::type Element;
  public:
    Reconstruct(const std::vector<Element> &flatVector,
                const std::vector<size_t> &dim);
    const V &                    getVector(void);
    const std::vector<Element> & getFlatVector(void);
    const std::vector<size_t>  & getDim(void);
  private:
    void fill(std::vector<Element> &v);
    template <typename W>
    void fill(W &v);
    void resize(std::vector<Element> &v, const unsigned int dim);
    template <typename W>
    void resize(W &v, const unsigned int dim);
  private:
    V                          vector_;
    const std::vector<Element> &flatVector_;
    std::vector<size_t>        dim_;
    size_t                     ind_{0};
    unsigned int               dimInd_{0};
  };

  // Flatten class template implementation
  template <typename V>
  void Flatten<V>::accumulate(const Element &e)
  {
    flatVector_.push_back(e);
  }
  
  template <typename V>
  template <typename W>
  void Flatten<V>::accumulate(const W &v)
  {
    for (auto &e: v)
    {
      accumulate(e);
    }
  }
  
  template <typename V>
  void Flatten<V>::accumulateDim(const Element &e) {};
  
  template <typename V>
  template <typename W>
  void Flatten<V>::accumulateDim(const W &v)
  {
    dim_.push_back(v.size());
    accumulateDim(v[0]);
  }
  
  template <typename V>
  Flatten<V>::Flatten(const V &vector)
  : vector_(vector)
  {
    accumulate(vector_);
    accumulateDim(vector_);
  }
  
  template <typename V>
  const V & Flatten<V>::getVector(void)
  {
    return vector_;
  }
  
  template <typename V>
  const std::vector<typename Flatten<V>::Element> &
  Flatten<V>::getFlatVector(void)
  {
    return flatVector_;
  }
  
  template <typename V>
  const std::vector<size_t> & Flatten<V>::getDim(void)
  {
    return dim_;
  }
  
  // Reconstruct class template implementation
  template <typename V>
  void Reconstruct<V>::fill(std::vector<Element> &v)
  {
    for (auto &e: v)
    {
      e = flatVector_[ind_++];
    }
  }
  
  template <typename V>
  template <typename W>
  void Reconstruct<V>::fill(W &v)
  {
    for (auto &e: v)
    {
      fill(e);
    }
  }
  
  template <typename V>
  void Reconstruct<V>::resize(std::vector<Element> &v, const unsigned int dim)
  {
    v.resize(dim_[dim]);
  }
  
  template <typename V>
  template <typename W>
  void Reconstruct<V>::resize(W &v, const unsigned int dim)
  {
    v.resize(dim_[dim]);
    for (auto &e: v)
    {
      resize(e, dim + 1);
    }
  }
  
  template <typename V>
  Reconstruct<V>::Reconstruct(const std::vector<Element> &flatVector,
                              const std::vector<size_t> &dim)
  : flatVector_(flatVector)
  , dim_(dim)
  {
    resize(vector_, 0);
    fill(vector_);
  }
  
  template <typename V>
  const V & Reconstruct<V>::getVector(void)
  {
    return vector_;
  }
  
  template <typename V>
  const std::vector<typename Reconstruct<V>::Element> &
  Reconstruct<V>::getFlatVector(void)
  {
    return flatVector_;
  }
  
  template <typename V>
  const std::vector<size_t> & Reconstruct<V>::getDim(void)
  {
    return dim_;
  }

  // Vector IO utilities ///////////////////////////////////////////////////////
  // helper function to read space-separated values
  template <typename T>
  std::vector<T> strToVec(const std::string s)
  {
    std::istringstream sstr(s);
    T                  buf;
    std::vector<T>     v;
    
    while(!sstr.eof())
    {
      sstr >> buf;
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