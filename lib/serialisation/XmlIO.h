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

#include "pugixml/pugixml.h"

namespace Grid
{
  
  class XmlWriter: public Writer<XmlWriter>
  {
    
  public:
    XmlWriter(const std::string &fileName);
    virtual ~XmlWriter(void);
    void push(const std::string &s);
    void pop(void);
    template <typename U>
    void writeDefault(const std::string &s, const U &x);
    template <typename U>
    void writeDefault(const std::string &s, const std::vector<U> &x);
  private:
    pugi::xml_document doc_;
    pugi::xml_node     node_;
    std::string        fileName_;
  };
  
  class XmlReader: public Reader<XmlReader>
  {
  public:
    XmlReader(const std::string &fileName);
    virtual ~XmlReader(void) = default;
    void push(const std::string &s);
    void pop(void);
    bool nextElement(const std::string &s);
    template <typename U>
    void readDefault(const std::string &s, U &output);
    template <typename U>
    void readDefault(const std::string &s, std::vector<U> &output);
  private:
    pugi::xml_document doc_;
    pugi::xml_node     node_;
    std::string        fileName_;
  };
  
  // Writer template implementation ////////////////////////////////////////////
  template <typename U>
  void XmlWriter::writeDefault(const std::string &s, const U &x)
  {
    std::ostringstream os;
    
    os << std::boolalpha << x;
    pugi::xml_node leaf = node_.append_child(s.c_str());
    leaf.append_child(pugi::node_pcdata).set_value(os.str().c_str());
  }
  
  template <typename U>
  void XmlWriter::writeDefault(const std::string &s, const std::vector<U> &x)
  {
    push(s);
    for (auto &x_i: x)
    {
      write("elem", x_i);
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
    std::string    buf;
    unsigned int   i = 0;
    
    push(s);
    while (node_.child("elem"))
    {
      output.resize(i + 1);
      read("elem", output[i]);
      node_.child("elem").set_name("elem-done");
      i++;
    }
    pop();
  }
  
}
#endif
