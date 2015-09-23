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


namespace Grid {

class XMLWriter  : public Writer {
private:

  pugi::xml_document doc;
  pugi::xml_node     node;
  std::string file;

public:

  XMLWriter(const std::string &_file) : file(_file)
  {
    node=doc.append_child();
    node.set_name("document");
  }

  ~XMLWriter()
  {
    //    simple_walker walker;
    //    doc.traverse(walker);

    doc.save_file(file.c_str(),"  ");
  }

  void push(const std::string &s)
  {
    node = node.append_child(s.c_str());
  }
  void pop(void) {
    node = node.parent();
  }

  void write( const std::string& s,const std::string &output      ) { 
    pugi::xml_node leaf=node.append_child(s.c_str());
    leaf.append_child(pugi::node_pcdata).set_value(output.c_str());
  };

  void write( const std::string& s,const  int16_t    output      ) { writeInternal(s,output); };
  void write( const std::string& s,const uint16_t    output      ) { writeInternal(s,output); };
  void write( const std::string& s,const  int32_t    output      ) { writeInternal(s,output); };
  void write( const std::string& s,const uint32_t    output      ) { writeInternal(s,output); };
  void write( const std::string& s,const  int64_t    output      ) { writeInternal(s,output); };
  void write( const std::string& s,const uint64_t    output      ) { writeInternal(s,output); };
  void write( const std::string& s,const  float      output      ) { writeInternal(s,output); };
  void write( const std::string& s,const double      output      ) { writeInternal(s,output); };
  void write( const std::string& s,const bool        output      ) { writeInternal(s,output); };

private:

  template<class T> void writeInternal( const std::string& s,const T output ){
    std::ostringstream os;
    os << std::boolalpha << output;
    write(s,os.str());
  }
  
};


class XMLReader : public Reader {
private:

  pugi::xml_document doc;
  pugi::xml_node     node;

public:


  XMLReader(const std::string &_file) 
  {
    pugi::xml_parse_result result = doc.load_file(_file.c_str());

    if ( !result ) { 
      std::cout << "XML error description: " << result.description() << "\n";
      std::cout << "XML error offset:      " << result.offset << "\n";
    }

    assert(result);

    //    simple_walker walker;
    //    doc.traverse(walker);

    node= doc.child("document");
  }

  ~XMLReader()  {  }

  void read( const std::string& s,std::string &output      ) { 
    output=node.child(s.c_str()).first_child().value();
  };
  void push(const std::string &s)
  {
    node = node.child(s.c_str());
  }
  void pop(void) {
    node = node.parent();
  }

  void read( const std::string& s,  int16_t  &output      ) { readInternal(s,output); };
  void read( const std::string& s, uint16_t  &output      ) { readInternal(s,output); };
  void read( const std::string& s,  int32_t  &output      ) { readInternal(s,output); };
  void read( const std::string& s, uint32_t  &output      ) { readInternal(s,output); };
  void read( const std::string& s,  int64_t  &output      ) { readInternal(s,output); };
  void read( const std::string& s, uint64_t  &output      ) { readInternal(s,output); };
  void read( const std::string& s,  float    &output      ) { readInternal(s,output); };
  void read( const std::string& s, double    &output      ) { readInternal(s,output); };
  void read( const std::string& s, bool      &output      ) { readInternal(s,output); };


private:

  template<class T> void readInternal( const std::string& path, T &output ){
    std::string asString;
    read(path,asString);
    convert(asString,output);
  }

  template<class T> void convert(const std::string &asString,T &output)
  {
    std::istringstream is(asString);  is.exceptions(std::ios::failbit);
    try {
      is >> std::boolalpha >> output;
    } catch(std::istringstream::failure e) {
      std::cerr << "XML read failure on "<<" "<<asString<<" "<<typeid(T).name()<<std::endl;
    }
    assert( is.tellg()==-1);
  }
  
};

}
#endif
