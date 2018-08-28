    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/serialisation/JSON_IO.cc

    Copyright (C) 2016

    Author: Guido Cossu<guido.cossu@ed.ac.uk>

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
#include <Grid/Grid.h>

using namespace Grid;

// Writer implementation ///////////////////////////////////////////////////////
JSONWriter::JSONWriter(const std::string &fileName)
: fileName_(fileName), ss_("{ ", std::ostringstream::ate){}

JSONWriter::~JSONWriter(void)
{
  // close
  delete_comma();
  ss_ << "}";  

  // write prettified JSON to file
  std::ofstream os(fileName_);
  //std::cout << "JSONWriter::~JSONWriter" << std::endl;
  os << std::setw(2) << json::parse(ss_.str()) << std::endl;
}

void JSONWriter::push(const std::string &s)
{
  // adding a nested object
  if (s.size())
    ss_ << " \""<<s<<"\" : {";
  else
    ss_ << " {";
}

void JSONWriter::pop(void)
{
  //std::cout << "JSONWriter::pop" << std::endl;
  delete_comma();
  ss_ << "},";
}

void JSONWriter::delete_comma()
{
  std::string dlast = ss_.str();
  dlast.pop_back(); // deletes the last comma
  ss_.str(dlast);
}


// here we are hitting a g++ bug (Bug 56480)
// compiles fine with clang
// have to wrap in the Grid namespace
// annoying, but necessary for TravisCI
namespace Grid
{
  void JSONWriter::writeDefault(const std::string &s,	const std::string &x)
  {
    //std::cout << "JSONWriter::writeDefault(string) : " << s <<  std::endl;
    std::ostringstream os;
    os << std::boolalpha << x;
    if (s.size())
      ss_ << "\""<< s << "\" : \"" << os.str() << "\" ," ;
    else
     ss_ << os.str() << " ," ;
  }
}// namespace Grid 


// Reader implementation ///////////////////////////////////////////////////////
JSONReader::JSONReader(const std::string &fileName)
: fileName_(fileName)
{
  std::ifstream file(fileName_);
  file >> jobject_;

  // test
  // serialize to standard output
  //std::cout << "JSONReader::JSONReader : " << jobject_ << endl; 
  jcur_ = jobject_;
}

bool JSONReader::push(const std::string &s)
{
  if (s.size()){
    jold_.push_back(jcur_);
    do_pop.push_back(true);
    try
    {
      jcur_ = jcur_[s]; 
    }
    catch (std::out_of_range& e)
    {
      std::cout << "out of range: " << e.what() << '\n';
      return false;
    }
    //cout << "JSONReader::push : " << s << " : "<< jcur_ << endl;
  }
  else
  {
    do_pop.push_back(false);
  }


  return true;
}

void JSONReader::pop(void)
{
  if (do_pop.back()){
    jcur_ = jold_.back();
    jold_.pop_back();
    do_pop.pop_back();
  }
  else
    do_pop.pop_back();

  //cout << "JSONReader::pop : " << jcur_ << endl;
}

bool JSONReader::nextElement(const std::string &s)
{
  // Work in progress
  // JSON dictionaries do not support multiple names 
  // Same name objects must be packed in vectors
  ++it_;
  
  //if (it_ == it_end_){
  //  return false;
  //}

  jcur_ = *it_; 
  //cout << "JSONReader::nextElement(string) : " << s << " : "<< jcur_ << endl;
  //return true;

    return false;
}

template <>
void JSONReader::readDefault(const std::string &s, std::string &output)
{
  //cout << "JSONReader::readDefault(string) : " << s<< " " << jcur_ << endl;
  if (s.size()){
    //std::cout << "String: "<< jcur_[s] << std::endl;
    output = jcur_[s];
  }
  else
  {
    //std::cout << "String: "<< jcur_ << std::endl;
    output = jcur_;    
  }
}
