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

namespace Grid {
  // Grid scalar tensors to nested std::vectors //////////////////////////////////
  template <typename T, typename V>
  void tensorToVec(V &v, const T& t)
  {
    v = t;
  }

  template <typename T, typename V>
  void tensorToVec(V &v, const iScalar<T>& t)
  {
    tensorToVec(v, t._internal);
  }

  template <typename T, typename V, int N>
  void tensorToVec(std::vector<V> &v, const iVector<T, N>& t)
  {
    v.resize(N);
    for (unsigned int i = 0; i < N; i++) 
    {
      tensorToVec(v[i], t._internal[i]);
    }
  }

  template <typename T, typename V, int N>
  void tensorToVec(std::vector<std::vector<V>> &v, const iMatrix<T, N>& t)
  {
    v.resize(N, std::vector<V>(N));
    for (unsigned int i = 0; i < N; i++)
    for (unsigned int j = 0; j < N; j++) 
    {
      tensorToVec(v[i][j], t._internal[i][j]);
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
    } while (c != '{' && !is.eof());
    if (c == '{')
    {
      int start = is.tellg();
      do
      {
        is.get(c);
      } while (c != '}' && !is.eof());
      if (c == '}')
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
    os << "{" << p.first << " " << p.second << "}";
    return os;
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
    for (auto &x: v)
    {
      os << x << " ";
    }
    if (v.size() > 0)
    {
      os << "\b";
    }
    os << "]";
    
    return os;
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
    typename std::enable_if<!std::is_base_of<Serializable, U>::value, void>::type
    write(const std::string& s, const U &output);
  private:
    T *upcast;
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



  // Generic writer interface
  // serializable base class
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
  };
  
  // Flatten class template implementation /////////////////////////////////////
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
  
  // Reconstruct class template implementation /////////////////////////////////
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
  
  // Generic reader interface
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
  
  // Writer template implementation ////////////////////////////////////////////
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
  typename std::enable_if<!std::is_base_of<Serializable, U>::value, void>::type
  Writer<T>::write(const std::string &s, const U &output)
  {
    upcast->writeDefault(s, output);
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
  void Reader<T>::fromString(U &output, const std::string &s)
  {
    std::istringstream is(s);
    
    is.exceptions(std::ios::failbit);
    try
    {
      is >> std::boolalpha >> output;
    }
    catch(std::istringstream::failure &e)
    {
      std::cerr << "numerical conversion failure on '" << s << "' ";
      std::cerr << "(typeid: " << typeid(U).name() << ")" << std::endl;
      abort();
    }
  }
}

#endif
