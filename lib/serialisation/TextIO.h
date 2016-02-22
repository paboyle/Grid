    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/serialisation/TextIO.h

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
#ifndef GRID_SERIALISATION_TEXT_READER_H
#define GRID_SERIALISATION_TEXT_READER_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <cassert>

namespace Grid
{
  
  class TextWriter: public Writer<TextWriter>
  {
  public:
    TextWriter(const std::string &fileName);
    virtual ~TextWriter(void) = default;
    void push(const std::string &s);
    void pop(void);
    template <typename U>
    void writeDefault(const std::string &s, const U &x);
    template <typename U>
    void writeDefault(const std::string &s, const std::vector<U> &x);
  private:
    void indent(void);
  private:
    std::ofstream file_;
    int           level_{0};
  };
  
  class TextReader: public Reader<TextReader>
  {
  public:
    TextReader(const std::string &fileName);
    virtual ~TextReader(void) = default;
    void push(const std::string &s);
    void pop(void);
    template <typename U>
    void readDefault(const std::string &s, U &output);
    template <typename U>
    void readDefault(const std::string &s, std::vector<U> &output);
  private:
    void checkIndent(void);
  private:
    std::ifstream file_;
    int           level_{0};
  };
  
  // Writer template implementation ////////////////////////////////////////////
  template <typename U>
  void TextWriter::writeDefault(const std::string &s, const U &x)
  {
    indent();
    file_ << std::boolalpha << x << std::endl;
  }
  
  template <typename U>
  void TextWriter::writeDefault(const std::string &s, const std::vector<U> &x)
  {
    uint64_t sz = x.size();
    
    write(s, sz);
    for (uint64_t i = 0; i < sz; ++i)
    {
      write(s, x[i]);
    }
  }
  
  // Reader template implementation ////////////////////////////////////////////
  template <typename U>
  void TextReader::readDefault(const std::string &s, U &output)
  {
    std::string buf;
    
    readDefault(s, buf);
    fromString(output, buf);
  }
  
  template <>
  void TextReader::readDefault(const std::string &s, std::string &output);
  
  template <typename U>
  void TextReader::readDefault(const std::string &s, std::vector<U> &output)
  {
    uint64_t sz;
    
    read("", sz);
    output.resize(sz);
    for (uint64_t i = 0; i < sz; ++i)
    {
      read("", output[i]);
    }
  }
}

#endif

