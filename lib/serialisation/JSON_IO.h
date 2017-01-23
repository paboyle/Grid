    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/serialisation/JSON_IO.h

    Copyright (C) 2015

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
#ifndef GRID_SERIALISATION_JSON_IO_H
#define GRID_SERIALISATION_JSON_IO_H

#include <iostream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>

#include <Grid/json/json.hpp>

// for convenience
using json = nlohmann::json;

namespace Grid
{
  
  class JSONWriter: public Writer<JSONWriter>
  {
    
  public:
    JSONWriter(const std::string &fileName);
    virtual ~JSONWriter(void);
    void push(const std::string &s);
    void pop(void);
    template <typename U>
    void writeDefault(const std::string &s, const U &x);
    template <typename U>
    void writeDefault(const std::string &s, const std::vector<U> &x);

    template<std::size_t N>
    void writeDefault(const std::string &s, const char(&x)[N]);

  private:
    void delete_comma();
    std::string         fileName_;
    std::ostringstream  ss_;
  };
  
  class JSONReader: public Reader<JSONReader>
  {
  public:
    JSONReader(const std::string &fileName);
    virtual ~JSONReader(void) = default;
    bool push(const std::string &s);
    void pop(void);
    bool nextElement(const std::string &s);
    template <typename U>
    void readDefault(const std::string &s, U &output);
    template <typename U>
    void readDefault(const std::string &s, std::vector<U> &output);
  private:
    json                jobject_; // main object
    json                jcur_;  // current json object
    std::vector<json>   jold_;  // previous json object
    std::string         fileName_;
    std::vector<bool>   do_pop;
  };
  
  // Writer template implementation ////////////////////////////////////////////
  template <typename U>
  void JSONWriter::writeDefault(const std::string &s, const U &x)
  {
    std::ostringstream os;
    os << std::boolalpha << x;
    if (s.size())
      ss_ << "\""<< s << "\" : " << os.str() << " ," ;
    else
     ss_ << os.str() << " ," ;
  }

  template <typename U>
  void JSONWriter::writeDefault(const std::string &s, const std::vector<U> &x)
  {
    if (s.size())
      ss_ << " \""<<s<<"\" : [";
    else
      ss_ << " [";
    for (auto &x_i: x)
    {
      write("", x_i);
    }
    delete_comma();
    ss_<< "],";
  }
  
  template<std::size_t N>
  void JSONWriter::writeDefault(const std::string &s, const char(&x)[N]){
    if (s.size())
    ss_ << "\""<< s << "\" : \"" << x << "\" ," ;
    else
    ss_ << "\"" << x << "\" ," ; 
  }

  // Reader template implementation ////////////////////////////////////////////
  template <typename U>
  void JSONReader::readDefault(const std::string &s, U &output)
  {
    //std::string buf;
    std::cout << "JSONReader::readDefault(U) : " << s << "  :  "<< jcur_ << std::endl;
    //readDefault(s, output);

    //std::cout << s << "   " << buf << std::endl;
    //fromString(output, buf);
  
    if (s.size()){
      std::cout << "String: "<< jcur_[s] << std::endl;
      output = jcur_[s];
    }
    else
    {
      std::cout << "String: "<< jcur_ << std::endl;
      output = jcur_;    
    }


  }
  
  template <>
  void JSONReader::readDefault(const std::string &s, std::string &output);
  
  template <typename U>
  void JSONReader::readDefault(const std::string &s, std::vector<U> &output)
  {
    std::string    buf;
    unsigned int   i = 0;
    std::cout << "JSONReader::readDefault(vec) : " << jcur_ << std::endl;
    if (s.size())
      push(s);
    
    json j = jcur_;
    for (json::iterator it = j.begin(); it != j.end(); ++it) {
      jcur_ = *it;
      std::cout << "Value: " << it.value() << "\n";
      output.resize(i + 1);
      read("", output[i++]);
    }

    jcur_ = j;
    if (s.size())
      pop();
  }
  
}
#endif
