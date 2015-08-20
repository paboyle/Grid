#ifndef GRID_SERIALISATION_READER_H
#define GRID_SERIALISATION_READER_H

#include <serialisation/MacroMagic.h>

namespace Grid {

// Generic reader writer interface

template< class Writer> void push(Writer & WR,const std::string &s) { WR.push(s);}
template< class Writer> void push(Writer & WR,const char *s)        { WR.push(std::string(s));}
template< class Writer> void pop (Writer & WR)                      { WR.pop();}

template< class Writer> void write(Writer& wr, const std::string& s,const char * output      ) { wr.iwrite(s,std::string(output)); };
template< class Writer> void write(Writer& wr, const std::string& s,const std::string &output) { wr.iwrite(s,output); };
template< class Writer> void write(Writer& wr, const std::string& s,  int16_t    output      ) { wr.iwrite(s,output); };
template< class Writer> void write(Writer& wr, const std::string& s, uint16_t    output      ) { wr.iwrite(s,output); };
template< class Writer> void write(Writer& wr, const std::string& s,  int32_t    output      ) { wr.iwrite(s,output); };
template< class Writer> void write(Writer& wr, const std::string& s, uint32_t    output      ) { wr.iwrite(s,output); };
template< class Writer> void write(Writer& wr, const std::string& s,  int64_t    output      ) { wr.iwrite(s,output); };
template< class Writer> void write(Writer& wr, const std::string& s, uint64_t    output      ) { wr.iwrite(s,output); };
template< class Writer> void write(Writer& wr, const std::string& s,  float      output      ) { wr.iwrite(s,output); };
template< class Writer> void write(Writer& wr, const std::string& s, double      output      ) { wr.iwrite(s,output); };
template< class Writer> void write(Writer& wr, const std::string& s, bool        output      ) { wr.iwrite(s,output); };

template< class Reader> void read(Reader& rd, const std::string& s,std::string &output)       { rd.iread(s,output); };
template< class Reader> void read(Reader& rd, const std::string& s,  int16_t   &output      ) { rd.iread(s,output); };
template< class Reader> void read(Reader& rd, const std::string& s, uint16_t   &output      ) { rd.iread(s,output); };
template< class Reader> void read(Reader& rd, const std::string& s,  int32_t   &output      ) { rd.iread(s,output); };
template< class Reader> void read(Reader& rd, const std::string& s, uint32_t   &output      ) { rd.iread(s,output); };
template< class Reader> void read(Reader& rd, const std::string& s,  int64_t   &output      ) { rd.iread(s,output); };
template< class Reader> void read(Reader& rd, const std::string& s, uint64_t   &output      ) { rd.iread(s,output); };
template< class Reader> void read(Reader& rd, const std::string& s,  float     &output      ) { rd.iread(s,output); };
template< class Reader> void read(Reader& rd, const std::string& s, double     &output      ) { rd.iread(s,output); };
template< class Reader> void read(Reader& rd, const std::string& s, bool       &output      ) { rd.iread(s,output); };


template<class Writer, class T>
void write(Writer& wr, const std::string& s,const std::vector<T> output ) { 
  push(wr,s);
  for(int i=0;i<output.size();i++){
    std::ostringstream oss;      oss << "elem" << i;
    write(wr,oss.str(),output[i]);
  }
  pop(wr);
};

template<class Reader, class T>
void read(Reader& rd, const std::string& s,std::vector<T> &output ) { 
  rd.iread(s,output);
};

template < class T >
inline std::ostream& operator << (std::ostream& os, const std::vector<T>& v) 
{
  os << "[";
  for (typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    {
      os << " " << *ii;
    }
  os << " ]";
  return os;
}

}

// Todo:
//#include <serialisation/CoutReader.h>
//#include <serialisation/TextReader.h>
//#include <serialisation/JSONReader.h>
//#include <serialisation/YAMLReader.h>
#include <serialisation/XMLReader.h>

namespace Grid {

using XMLPolicy::Reader;
using XMLPolicy::Writer;

}

#endif
