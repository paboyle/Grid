#ifndef GRID_SERIALISATION_ABSTRACT_READER_H
#define GRID_SERIALISATION_ABSTRACT_READER_H

namespace Grid {
 class Writer  {
public:

  virtual ~Writer() {};

  virtual void push(const std::string &s) = 0;
  virtual void pop(void) =0;
  virtual void write( const std::string& s,const std::string &output      ) =0;
  virtual void write( const std::string& s,  int16_t    output      ) =0;
  virtual void write( const std::string& s, uint16_t    output      ) =0;
  virtual void write( const std::string& s,  int32_t    output      ) =0;
  virtual void write( const std::string& s, uint32_t    output      ) =0;
  virtual void write( const std::string& s,  int64_t    output      ) =0;
  virtual void write( const std::string& s, uint64_t    output      ) =0;
  virtual void write( const std::string& s,  float      output      ) =0;
  virtual void write( const std::string& s, double      output      ) =0;
  virtual void write( const std::string& s, bool        output      ) =0;
  
};

class Reader {
public:



  virtual ~Reader() {};

  virtual void read( const std::string& s,std::string &output      ) =0;
  virtual void push(const std::string &s) =0;
  virtual void pop(void) = 0;

  virtual void read( const std::string& s,  int16_t  &output      ) =0;
  virtual void read( const std::string& s, uint16_t  &output      ) =0;
  virtual void read( const std::string& s,  int32_t  &output      ) =0;
  virtual void read( const std::string& s, uint32_t  &output      ) =0;
  virtual void read( const std::string& s,  int64_t  &output      ) =0;
  virtual void read( const std::string& s, uint64_t  &output      ) =0;
  virtual void read( const std::string& s,  float    &output      ) =0;
  virtual void read( const std::string& s, double    &output      ) =0;
  virtual void read( const std::string& s, bool      &output      ) =0;
  
};

}
#endif
