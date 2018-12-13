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

namespace Grid {
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
    template <typename U>
    void write(const std::string &s, const iScalar<U> &output);
    template <typename U, int N>
    void write(const std::string &s, const iVector<U, N> &output);
    template <typename U, int N>
    void write(const std::string &s, const iMatrix<U, N> &output);
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
  typename std::enable_if<!std::is_base_of<Serializable, U>::value, void>::type
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
