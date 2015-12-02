#include <Grid.h>


namespace Grid {
// Writer implementation ///////////////////////////////////////////////////////
BinaryWriter::BinaryWriter(const std::string &fileName)
: file_(fileName, std::ios::binary|std::ios::out)
{}

template <>
void BinaryWriter::writeDefault(const std::string &s, const std::string &output)
{
  uint64_t sz = output.size();
  
  write("", sz);
  for (uint64_t i = 0; i < sz; ++i)
  {
    write("", output[i]);
  }
}

// Reader implementation ///////////////////////////////////////////////////////
BinaryReader::BinaryReader(const std::string &fileName)
: file_(fileName, std::ios::binary|std::ios::in)
{}

template <>
void BinaryReader::readDefault(const std::string &s, std::string &output)
{
  uint64_t sz;
  
  read("", sz);
  output.reserve(sz);
  file_.read((char *)output.data(), sz);
}
}
