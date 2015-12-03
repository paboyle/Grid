#ifndef GRID_SERIALISATION_TEXT_READER_H
#define GRID_SERIALISATION_TEXT_READER_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <cassert>

namespace Grid
{
  
  class TextWriter: public Writer<TextWriter>
  {
  public:
    TextWriter(const std::string &fileName);
    virtual ~TextWriter(void) = default;
    void push(const std::string &s);
    void pop(void);
    template <typename U>
    void writeDefault(const std::string &s, const U &output);
    template <typename U>
    void writeDefault(const std::string &s, const std::vector<U> &output);
  private:
    void indent(void);
  private:
    std::ofstream file_;
    int           level_{0};
  };
  
  class TextReader: public Reader<TextReader>
  {
  public:
    TextReader(const std::string &fileName);
    virtual ~TextReader(void) = default;
    void push(const std::string &s);
    void pop(void);
    template <typename U>
    void readDefault(const std::string &s, U &output);
    template <typename U>
    void readDefault(const std::string &s, std::vector<U> &output);
  private:
    void checkIndent(void);
  private:
    std::ifstream file_;
    int           level_{0};
  };
  
  // Writer template implementation ////////////////////////////////////////////
  template <typename U>
  void TextWriter::writeDefault(const std::string &s, const U &output)
  {
    indent();
    file_ << std::boolalpha << output << std::endl;
  }
  
  template <typename U>
  void TextWriter::writeDefault(const std::string &s, const std::vector<U> &output)
  {
    uint64_t sz = output.size();
    
    write(s, sz);
    for (uint64_t i = 0; i < sz; ++i)
    {
      write(s, output[i]);
    }
  }
  
  // Reader template implementation ////////////////////////////////////////////
  template <typename U>
  void TextReader::readDefault(const std::string &s, U &output)
  {
    std::string buf;
    
    readDefault(s, buf);
    fromString(output, buf);
  }
  
  template <typename U>
  void TextReader::readDefault(const std::string &s, std::vector<U> &output)
  {
    uint64_t sz;
    
    read("", sz);
    output.resize(sz);
    for (uint64_t i = 0; i < sz; ++i)
    {
      read("", output[i]);
    }
  }
}

#endif

