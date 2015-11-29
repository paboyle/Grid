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
  
  class BinaryWriter: public Writer<BinaryWriter>
  {
  public:
    BinaryWriter(const std::string &fileName);
    virtual ~BinaryWriter(void) = default;
    void push(const std::string &s) {};
    void pop(void) {};
    template <typename U>
    void writeDefault(const std::string &s, const U &x);
    template <typename U>
    void writeDefault(const std::string &s, const std::vector<U> &x);
  private:
    std::ofstream file_;
  };
  
  class BinaryReader: public Reader<BinaryReader>
  {
  public:
    BinaryReader(const std::string &fileName);
    virtual ~BinaryReader(void) = default;
    void push(const std::string &s) {};
    void pop(void) {};
    template <typename U>
    void readDefault(const std::string &s, U &output);
    template <typename U>
    void readDefault(const std::string &s, std::vector<U> &output);
  private:
    std::ifstream file_;
  };
  
  // Writer template implementation ////////////////////////////////////////////
  template <typename U>
  void BinaryWriter::writeDefault(const std::string &s, const U &x)
  {
    file_.write((char *)&x, sizeof(U));
  }
  
  template <typename U>
  void BinaryWriter::writeDefault(const std::string &s, const std::vector<U> &x)
  {
    uint64_t sz = x.size();
    
    write("", sz);
    for (uint64_t i = 0; i < sz; ++i)
    {
      write("", x[i]);
    }
  }
  
  // Reader template implementation ////////////////////////////////////////////
  template <typename U>
  void BinaryReader::readDefault(const std::string &s, U &output)
  {
    file_.read((char *)&output, sizeof(U));
  }
  
  template <typename U>
  void BinaryReader::readDefault(const std::string &s, std::vector<U> &output)
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
