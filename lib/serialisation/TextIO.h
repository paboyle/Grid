#ifndef GRID_SERIALISATION_TEXT_READER_H
#define GRID_SERIALISATION_TEXT_READER_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <cassert>

namespace Grid {

class TextWriter  : public Writer {
private:

  std::ofstream file;
  int level;
  void indent(void) {
    for(int i=0;i<level;i++){
      file <<"\t";
    }
  }
public:

  TextWriter(const std::string &_file) : file(_file,std::ios::out) {
    level=0;
  }

  ~TextWriter()  {  }

  void push(const std::string &s)
  {
    //    std::string tmp = s;
    //    write(s,tmp);
    level++;
  }
  void pop(void) {
    level--;
  }

  void write( const std::string& s,const std::string &output      ) { 
    indent();
    file<<output<<std::endl;
  };
  void write( const std::string& s,  int16_t    output      ) { writeInternal(s,output); };
  void write( const std::string& s, uint16_t    output      ) { writeInternal(s,output); };
  void write( const std::string& s,  int32_t    output      ) { writeInternal(s,output); };
  void write( const std::string& s, uint32_t    output      ) { writeInternal(s,output); };
  void write( const std::string& s,  int64_t    output      ) { writeInternal(s,output); };
  void write( const std::string& s, uint64_t    output      ) { writeInternal(s,output); };
  void write( const std::string& s,  float      output      ) { writeInternal(s,output); };
  void write( const std::string& s, double      output      ) { writeInternal(s,output); };
  void write( const std::string& s, bool        output      ) { writeInternal(s,output); };

private:

  template<class T> void writeInternal( const std::string& s, T output ){
    indent();
    file << std::boolalpha << output<<std::endl;
  }
  
};

class TextReader : public Reader {
private:

  std::ifstream file;
  int level;

public:


 TextReader(const std::string &_file) : file(_file,std::ios::in) { level = 0;};

  ~TextReader()  {  }

  void read( const std::string& s,std::string &output      ) { 
    char c='a';
    for(int i=0;i<level;i++){
      file.get(c);
      if ( c != '\t' ) 
	std::cout << "mismatch on tab "<<c<<" level "<< level<< " i "<< i<<std::endl;
    }
    output.clear();
    std::getline(file,output);
  };
  void push(const std::string &s)  { 
    //  std::string tmp; read(s,tmp); 
    level++;
  }
  void pop(void) { level--; }

  template<class T>
  void read( const std::string& s, std::vector<T>  &output      ) { 
    
    push(s);

    uint64_t n; read("N",n);

    // skip the vector length
    T tmp;
    output.resize(0);
    for(int i=0;i<n;i++){
      std::ostringstream oss;      oss << "elem" << i;
      read(*this,oss.str(),tmp);
      output.push_back(tmp);
    }

    pop();

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

