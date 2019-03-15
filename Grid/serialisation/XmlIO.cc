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

void Grid::xmlCheckParse(const pugi::xml_parse_result &result, const std::string name)
{
  if (!result) 
  {
    std::cerr << "XML parsing error for " << name << std::endl;
    std::cerr << "XML error description: " << result.description() << std::endl;
    std::cerr << "XML error offset     : " << result.offset << std::endl;
    abort();
  }
}

// Writer implementation ///////////////////////////////////////////////////////
XmlWriter::XmlWriter(const std::string &fileName, std::string toplev) : fileName_(fileName)
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
    doc_.save_file(fileName_.c_str(), indent_.c_str());
  }
}

void XmlWriter::push(const std::string &s)
{
  node_ = node_.append_child(s.c_str());
}

void XmlWriter::pushXmlString(const std::string &s)
{
  pugi::xml_document doc;
  auto               result = doc.load_buffer(s.c_str(), s.size());

  xmlCheckParse(result, "fragment\n'" + s +"'");
  for (pugi::xml_node child = doc.first_child(); child; child = child.next_sibling())
  {
      node_ = node_.append_copy(child);
  }
  pop();
}

void XmlWriter::pop(void)
{
  node_ = node_.parent();
}

std::string XmlWriter::docString(void)
{
  std::ostringstream oss; 
  doc_.save(oss, indent_.c_str());
  return oss.str();
}

std::string XmlWriter::string(void)
{
  std::ostringstream oss; 
  doc_.save(oss, indent_.c_str(), pugi::format_default | pugi::format_no_declaration);
  return oss.str();
}

// Reader implementation ///////////////////////////////////////////////////////
XmlReader::XmlReader(const std::string &s,  const bool isBuffer, 
                     std::string toplev) 
{
  pugi::xml_parse_result result;
  
  if (isBuffer)
  {
    fileName_ = "<string>";
    result    = doc_.load_string(s.c_str());
    xmlCheckParse(result, "string\n'" + s + "'");
  }
  else
  {
    fileName_ = s;
    result    = doc_.load_file(s.c_str());
    xmlCheckParse(result, "file '" + fileName_ + "'");
  }
  if ( toplev == std::string("") ) {
  node_ = doc_;
  } else { 
    node_ = doc_.child(toplev.c_str());
  }
}

#define XML_SAFE_NODE(expr)\
if (expr)\
{\
  node_ = expr;\
  return true;\
}\
else\
{\
  return false;\
}

bool XmlReader::push(const std::string &s)
{
  if (s.empty())
  {
    XML_SAFE_NODE(node_.first_child());
  }
  else
  {
    XML_SAFE_NODE(node_.child(s.c_str()));
  }
}

void XmlReader::pop(void)
{
  node_ = node_.parent();
}

bool XmlReader::nextElement(const std::string &s)
{
  if (s.empty())
  {
    XML_SAFE_NODE(node_.next_sibling());
  }
  else
  {
    XML_SAFE_NODE(node_.next_sibling(s.c_str()));
  }
}

void XmlReader::readCurrentSubtree(std::string &s)
{
  std::ostringstream oss; 
  pugi::xml_document doc;

  doc.append_copy(node_);
  doc.save(oss, indent_.c_str(), pugi::format_default | pugi::format_no_declaration);
  s = oss.str();
}

template <>
void XmlReader::readDefault(const std::string &s, std::string &output)
{
  if (node_.child(s.c_str()))
  {
    output = node_.child(s.c_str()).first_child().value();
  }
  else
  {
    std::cout << GridLogWarning << "XML: cannot open node '" << s << "'";
    std::cout << std::endl;

    output = ""; 
  }
}
