    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/serialisation/BaseIO.h

    Copyright (C) 2015

Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef GRID_SERIALISATION_ABSTRACT_READER_H
#define GRID_SERIALISATION_ABSTRACT_READER_H

#include <type_traits>

namespace Grid {
  
  class Serializable {};
  
  // static polymorphism implemented using CRTP idiom
  
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
    typename std::enable_if<std::is_enum<U>::value, void>::type
    write(const std::string& s, const U &output);
    template <typename U>
    typename std::enable_if<
      !(std::is_base_of<Serializable, U>::value or std::is_enum<U>::value),
      void>::type
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
    void push(const std::string &s);
    void pop(void);
    template <typename U>
    typename std::enable_if<std::is_base_of<Serializable, U>::value, void>::type
    read(const std::string& s, U &output);
    template <typename U>
    typename std::enable_if<std::is_enum<U>::value, void>::type
    read(const std::string& s, U &output);
    template <typename U>
    typename std::enable_if<
      !(std::is_base_of<Serializable, U>::value or std::is_enum<U>::value),
      void>::type
    read(const std::string& s, U &output);
  protected:
    template <typename U>
    void fromString(U &output, const std::string &s);
  private:
    T *upcast;
  };
  
  // Generic writer interface
  template <typename T>
  inline void push(Writer<T> &w, const std::string &s)
  {
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
  inline void push(Reader<T> &r, const std::string &s)
  {
    r.push(s);
  }
  
  template <typename T>
  inline void push(Reader<T> &r, const char *s)
  {
    r.push(std::string(s));
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
  
  template < class T >
  inline std::ostream& operator << (std::ostream& os, const std::vector<T>& v)
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
  typename std::enable_if<std::is_enum<U>::value, void>::type
  Writer<T>::write(const std::string &s, const U &output)
  {
    EnumIO<U>::write(*this, s, output);
  }
  
  template <typename T>
  template <typename U>
  typename std::enable_if<
    !(std::is_base_of<Serializable, U>::value or std::is_enum<U>::value),
    void>::type
  Writer<T>::write(const std::string &s, const U &output)
  {
    upcast->writeDefault(s, output);
  }
  
  // Reader template implementation ////////////////////////////////////////////
  template <typename T>
  Reader<T>::Reader(void)
  {
    upcast = static_cast<T *>(this);
  }
  
  template <typename T>
  void Reader<T>::push(const std::string &s)
  {
    upcast->push(s);
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
  typename std::enable_if<std::is_enum<U>::value, void>::type
  Reader<T>::read(const std::string &s, U &output)
  {
    EnumIO<U>::read(*this, s, output);
  }
  
  template <typename T>
  template <typename U>
  typename std::enable_if<
    !(std::is_base_of<Serializable, U>::value or std::is_enum<U>::value),
    void>::type
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
