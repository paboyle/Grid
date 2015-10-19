#ifndef GRID_SERIALISATION_READER_H
#define GRID_SERIALISATION_READER_H

#include <serialisation/MacroMagic.h>
#include <serialisation/BaseIO.h>
#include <stdint.h>

namespace Grid {


  // Generic reader writer interface
  inline void push(Writer & WR,const std::string &s) { WR.push(s);}
  inline void push(Writer & WR,const char *s)        { WR.push(std::string(s));}
  inline void pop (Writer & WR)                      { WR.pop();}

  //  inline void write(Writer& wr, const std::string& s,const char * output      ) { wr.write(s,std::string(output)); };
  inline void write(Writer& wr, const std::string& s,const std::string &output) { wr.write(s,output); };
  inline void write(Writer& wr, const std::string& s,const  int16_t    output      ) { wr.write(s,output); };
  inline void write(Writer& wr, const std::string& s,const uint16_t    output      ) { wr.write(s,output); };
  inline void write(Writer& wr, const std::string& s,const  int32_t    output      ) { wr.write(s,output); };
  inline void write(Writer& wr, const std::string& s,const uint32_t    output      ) { wr.write(s,output); };
  inline void write(Writer& wr, const std::string& s,const  int64_t    output      ) { wr.write(s,output); };
  inline void write(Writer& wr, const std::string& s,const uint64_t    output      ) { wr.write(s,output); };
  inline void write(Writer& wr, const std::string& s,const  float      output      ) { wr.write(s,output); };
  inline void write(Writer& wr, const std::string& s,const double      output      ) { wr.write(s,output); };
  inline void write(Writer& wr, const std::string& s,const bool        output      ) { wr.write(s,output); };
  
  inline void push(Reader & WR,const std::string &s) { WR.push(s);}
  inline void push(Reader & WR,const char *s)        { WR.push(std::string(s));}
  inline void pop (Reader & WR)                      { WR.pop();}
  
  inline void read(Reader& rd, const std::string& s,std::string &output)       { rd.read(s,output); };
  inline void read(Reader& rd, const std::string& s,  int16_t   &output      ) { rd.read(s,output); };
  inline void read(Reader& rd, const std::string& s, uint16_t   &output      ) { rd.read(s,output); };
  inline void read(Reader& rd, const std::string& s,  int32_t   &output      ) { rd.read(s,output); };
  inline void read(Reader& rd, const std::string& s, uint32_t   &output      ) { rd.read(s,output); };
  inline void read(Reader& rd, const std::string& s,  int64_t   &output      ) { rd.read(s,output); };
  inline void read(Reader& rd, const std::string& s, uint64_t   &output      ) { rd.read(s,output); };
  inline void read(Reader& rd, const std::string& s,  float     &output      ) { rd.read(s,output); };
  inline void read(Reader& rd, const std::string& s, double     &output      ) { rd.read(s,output); };
  inline void read(Reader& rd, const std::string& s, bool       &output      ) { rd.read(s,output); };


  template<class T> void write(Writer& wr, const std::string& s,const std::vector<T> output ) { 
    push(wr,s);
    uint64_t sz =output.size();
    write(wr,"N",sz);
    for(int i=0;i<output.size();i++){
      std::ostringstream oss;      oss << "elem" << i;
      write(wr,oss.str(),output[i]);
    }
    pop(wr);
  };


  
  template<class T>
  void read(Reader& rd, const std::string& s,std::vector<T> &output ) { 

    push(rd,s);

    uint64_t sz; read(rd,"N",sz);
    // skip the vector length
    T tmp;
    output.resize(0);
    for(int i=0;i<sz;i++){
      std::ostringstream oss;      oss << "elem" << i;
      read(rd,oss.str(),tmp);
      output.push_back(tmp);
    }
    pop(rd);
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

//////////////////////////////////////////
// Todo:
//////////////////////////////////////////
#include <serialisation/BinaryIO.h>
#include <serialisation/TextIO.h>
//#include <serialisation/JsonIO.h>
//#include <serialisation/YamlIO.h>
#include <serialisation/XmlIO.h>

//////////////////////////////////////////
// Select the default serialiser use ifdef's
//////////////////////////////////////////
namespace Grid {
  typedef XMLReader DefaultReader;
  typedef XMLWriter DefaultWriter;
}
#endif
