    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/serialisation/BinaryIO.h

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
#ifndef GRID_SERIALISATION_BINARY_READER_H
#define GRID_SERIALISATION_BINARY_READER_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <cassert>

namespace Grid {
  
  class BinaryWriter: public Writer<BinaryWriter>
  {
  public:
    BinaryWriter(const std::string &fileName);
    virtual ~BinaryWriter(void) = default;
    void push(const std::string &s) {};
    void pop(void) {};
    template <typename U>
    void writeDefault(const std::string &s, const U &x);
    template <typename U>
    void writeDefault(const std::string &s, const std::vector<U> &x);
    void writeDefault(const std::string &s, const char *x);
    template <typename U>
    void writeMultiDim(const std::string &s, const std::vector<size_t> & Dimensions, const U * pDataRowMajor, size_t NumElements);
  private:
    std::ofstream file_;
  };
  
  class BinaryReader: public Reader<BinaryReader>
  {
  public:
    BinaryReader(const std::string &fileName);
    virtual ~BinaryReader(void) = default;
    bool push(const std::string &s) {return true;}
    void pop(void) {};
    template <typename U>
    void readDefault(const std::string &s, U &output);
    template <typename U>
    void readDefault(const std::string &s, std::vector<U> &output);
    template <typename U>
    void readMultiDim(const std::string &s, std::vector<U> &buf, std::vector<size_t> &dim);
  private:
    std::ifstream file_;
  };
  
  // Writer template implementation ////////////////////////////////////////////
  template <typename U>
  void BinaryWriter::writeDefault(const std::string &s, const U &x)
  {
    file_.write((char *)&x, sizeof(U));
  }
  
  template <>
  void BinaryWriter::writeDefault(const std::string &s, const std::string &x);
  
  template <typename U>
  void BinaryWriter::writeDefault(const std::string &s, const std::vector<U> &x)
  {
    uint64_t sz = x.size();
    
    write("", sz);
    for (uint64_t i = 0; i < sz; ++i)
    {
      write("", x[i]);
    }
  }
  
  template <typename U>
  void BinaryWriter::writeMultiDim(const std::string &s, const std::vector<size_t> & Dimensions, const U * pDataRowMajor, size_t NumElements)
  {
    uint64_t rank = static_cast<uint64_t>( Dimensions.size() );
    uint64_t tmp = 1;
    for( auto i = 0 ; i < rank ; i++ )
      tmp *= Dimensions[i];
    assert( tmp == NumElements && "Dimensions don't match size of data being written" );
    // Total number of elements
    write("", tmp);
    // Number of dimensions
    write("", rank);
    // Followed by each dimension
    for( auto i = 0 ; i < rank ; i++ ) {
      tmp = Dimensions[i];
      write("", tmp);
    }
    for( auto i = 0; i < NumElements; ++i)
      write("", pDataRowMajor[i]);
  }

  // Reader template implementation ////////////////////////////////////////////
  template <typename U>
  void BinaryReader::readDefault(const std::string &s, U &output)
  {
    file_.read((char *)&output, sizeof(U));
  }
  
  template <>
  void BinaryReader::readDefault(const std::string &s, std::string &output);
  
  template <typename U>
  void BinaryReader::readDefault(const std::string &s, std::vector<U> &output)
  {
    uint64_t sz;
    
    read("", sz);
    output.resize(sz);
    for (uint64_t i = 0; i < sz; ++i)
    {
      read("", output[i]);
    }
  }

  template <typename U>
  void BinaryReader::readMultiDim(const std::string &s, std::vector<U> &buf, std::vector<size_t> &dim)
  {
    // Number of elements
    uint64_t NumElements;
    read("", NumElements);
    // Number of dimensions
    uint64_t rank;
    read("", rank);
    // Followed by each dimension
    uint64_t count = 1;
    dim.resize(rank);
    uint64_t tmp;
    for( auto i = 0 ; i < rank ; i++ ) {
      read("", tmp);
      dim[i] = tmp;
      count *= tmp;
    }
    assert( count == NumElements && "Dimensions don't match size of data being read" );
    buf.resize(count);
    for( auto i = 0; i < count; ++i)
      read("", buf[i]);
  }
}

#endif
