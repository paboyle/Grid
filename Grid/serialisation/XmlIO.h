    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/serialisation/XmlIO.h

    Copyright (C) 2015

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#ifndef GRID_SERIALISATION_XML_READER_H
#define GRID_SERIALISATION_XML_READER_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <cassert>

#include <Grid/pugixml/pugixml.h>
#include <Grid/GridCore.h>

namespace Grid
{
  void xmlCheckParse(const pugi::xml_parse_result &result, const std::string name);
  
  class XmlWriter: public Writer<XmlWriter>
  {    
  public:
    XmlWriter(const std::string &fileName, std::string toplev = std::string("grid") );
    virtual ~XmlWriter(void);
    void push(const std::string &s);
    void pushXmlString(const std::string &s);
    void pop(void);
    template <typename U>
    void writeDefault(const std::string &s, const U &x);
    template <typename U>
    void writeDefault(const std::string &s, const std::vector<U> &x);
    template <typename U>
    void writeMultiDim(const std::string &s, const std::vector<size_t> & Dimensions, const U * pDataRowMajor, size_t NumElements);
    std::string docString(void);
    std::string string(void);
  private:
    const std::string  indent_{"  "};
    pugi::xml_document doc_;
    pugi::xml_node     node_;
    std::string        fileName_;
  };
  
  class XmlReader: public Reader<XmlReader>
  {
  public:
    XmlReader(const std::string &fileName, const bool isBuffer = false, 
              std::string toplev = std::string("grid") );
    virtual ~XmlReader(void) = default;
    bool push(const std::string &s = "");
    void pop(void);
    bool nextElement(const std::string &s = "");
    template <typename U>
    void readDefault(const std::string &s, U &output);
    template <typename U>
    void readDefault(const std::string &s, std::vector<U> &output);
    template <typename U>
    void readMultiDim(const std::string &s, std::vector<U> &buf, std::vector<size_t> &dim);
    void readCurrentSubtree(std::string &s);
  private:
    void checkParse(const pugi::xml_parse_result &result, const std::string name);
  private:
    const std::string  indent_{"  "};
    pugi::xml_document doc_;
    pugi::xml_node     node_;
    std::string        fileName_;
  };

  template <>
  struct isReader< XmlReader > {
    static const bool value = true;
  };

  template <>
  struct isWriter<XmlWriter > {
    static const bool value = true;
  };
  
  // Writer template implementation ////////////////////////////////////////////
  template <typename U>
  void XmlWriter::writeDefault(const std::string &s, const U &x)
  {
    std::ostringstream os;
    
    if (getPrecision())
    {
      os.precision(getPrecision());
    }
    if (isScientific())
    {
      os << std::scientific;
    }
    os << std::boolalpha << x;
    pugi::xml_node leaf = node_.append_child(s.c_str());
    leaf.append_child(pugi::node_pcdata).set_value(os.str().c_str());
  }
  
  template <typename U>
  void XmlWriter::writeDefault(const std::string &s, const std::vector<U> &x)
  {
    std::vector<size_t> dims(1);
    dims[0] = x.size();
    writeMultiDim(s, dims, &x[0], dims[0]);
  }

  template <typename U>
  void XmlWriter::writeMultiDim(const std::string &s, const std::vector<size_t> & Dimensions, const U * pDataRowMajor, size_t NumElements)
  {
    push(s);
    if( Dimensions.size() > 1 )
    {
      for( auto d : Dimensions )
        write("dim", d);
    }
    while (NumElements--)
    {
      write("elem", *pDataRowMajor++);
    }
    pop();
  }

  // Reader template implementation ////////////////////////////////////////////
  template <typename U>
  void XmlReader::readDefault(const std::string &s, U &output)
  {
    std::string buf;
    
    readDefault(s, buf);
    fromString(output, buf);
  }
  
  template <>
  void XmlReader::readDefault(const std::string &s, std::string &output);
  
  template <typename U>
  void XmlReader::readDefault(const std::string &s, std::vector<U> &output)
  {
    std::vector<size_t> dims;
    readMultiDim(s, output, dims);
    assert( dims.size() == 1 && dims[0] == output.size() && "XmlIO: Expected 1D vector" );
  }

  template <typename U>
  void XmlReader::readMultiDim(const std::string &s, std::vector<U> &buf, std::vector<size_t> &dim)
  {
    unsigned int i = 0;
    unsigned int Rank = 0;
    if (!push(s))
    {
      std::cout << GridLogWarning << "XML: cannot open node '" << s << "'";
      std::cout << std::endl;
    } else {
      while (node_.child("dim"))
      {
        dim.resize(Rank + 1);
        read("dim", dim[Rank]);
        node_.child("dim").set_name("dim-done");
        Rank++;
      }
      while (node_.child("elem"))
      {
        buf.resize(i + 1);
        read("elem", buf[i]);
        node_.child("elem").set_name("elem-done");
        i++;
      }
      pop();
      if( Rank == 0 )
        dim.push_back(i);
    }
  }
}
#endif
