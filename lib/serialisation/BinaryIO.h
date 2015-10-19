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

class BinaryWriter : public Writer {
private:
  
  std::ofstream file;

public:

  BinaryWriter(const std::string &_file) : file(_file,std::ios::binary|std::ios::out) {}

  ~BinaryWriter()  {}

  // Binary is scopeless
  void push(const std::string &s) {} 
  void pop(void) {}

  void write(const std::string& s,const std::string &output) { 
    uint32_t sz = output.size();
    write(s,sz);
    const char * cstr = output.c_str();
    for(int c=0;c<output.size();c++){
      write(s,cstr[c]);
    }
  };
  void write( const std::string& s,const  char       output      ) { writeInternal(s,output); };
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
    // FIXME --- htons, htonl, htno64 etc..
    file.write((char *)&output,sizeof(T));
  }
  
};


class BinaryReader : public Reader{
private:

  std::ifstream file;

public:


  BinaryReader(const std::string &_file) : file(_file,std::ios::binary|std::ios::in) {}

  ~BinaryReader() {}

  // Binary is scopeless
  void push(const std::string &s) { }
  void pop(void) { }

  void read( const std::string& s,std::string &output ) { 

    output.clear();

    uint32_t sz;
    file.read((char *)&sz,sizeof(sz));

    for(int c=0;c<sz;c++){
      char ch;
      file.read(&ch,sizeof(ch));
      output.push_back(ch);
    }

  };

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
    file.read((char *)&output,sizeof(output)); // byte order??
  }

};

}
#endif
