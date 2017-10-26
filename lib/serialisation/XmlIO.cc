    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/serialisation/XmlIO.cc

    Copyright (C) 2015

Author: Antonin Portelli <antonin.portelli@me.com>
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
#include <Grid/GridCore.h>

using namespace Grid;
using namespace std;

// Writer implementation ///////////////////////////////////////////////////////
XmlWriter::XmlWriter(const string &fileName, string toplev) : fileName_(fileName)
{
  if ( toplev == std::string("") ) {
    node_=doc_;
  } else { 
    node_=doc_.append_child();
    node_.set_name(toplev.c_str());
  }
}

XmlWriter::~XmlWriter(void)
{
  if ( fileName_ != std::string("") ) { 
    doc_.save_file(fileName_.c_str(), "  ");
  }
}

void XmlWriter::push(const string &s)
{
  node_ = node_.append_child(s.c_str());
}

void XmlWriter::pop(void)
{
  node_ = node_.parent();
}
std::string XmlWriter::XmlString(void)
{
  std::ostringstream oss; 
  doc_.save(oss);
  return oss.str();
}

XmlReader::XmlReader(const char *xmlstring,string toplev) : fileName_("")
{
  pugi::xml_parse_result result;
  result = doc_.load_string(xmlstring);
  if ( !result ) {
    cerr << "XML error description (from char *): " << result.description() << "\nXML\n"<< xmlstring << "\n";
    cerr << "XML error offset      (from char *) " << result.offset         << "\nXML\n"<< xmlstring <<"\n";
    abort();
  }
  if ( toplev == std::string("") ) {
    node_ = doc_;
  } else { 
    node_ = doc_.child(toplev.c_str());
  }
}

// Reader implementation ///////////////////////////////////////////////////////
XmlReader::XmlReader(const string &fileName,string toplev) : fileName_(fileName)
{
  pugi::xml_parse_result result;
  result = doc_.load_file(fileName_.c_str());
  if ( !result ) {
    cerr << "XML error description: " << result.description() <<" "<< fileName_ <<"\n";
    cerr << "XML error offset     : " << result.offset        <<" "<< fileName_ <<"\n";
    abort();
  }
  if ( toplev == std::string("") ) {
    node_ = doc_;
  } else { 
    node_ = doc_.child(toplev.c_str());
  }
}

bool XmlReader::push(const string &s)
{

  if (node_.child(s.c_str()) == NULL )
    return false;

  node_ = node_.child(s.c_str());
  return true;
    
}

void XmlReader::pop(void)
{
  node_ = node_.parent();
}

bool XmlReader::nextElement(const std::string &s)
{
  if (node_.next_sibling(s.c_str()))
    {
      node_ = node_.next_sibling(s.c_str());
      
      return true;
    }
  else
    {
      return false;
    }

}

template <>
void XmlReader::readDefault(const string &s, string &output)
{
  output = node_.child(s.c_str()).first_child().value();
}
