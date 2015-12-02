#include <Grid.h>

using namespace Grid;
using namespace std;

// Writer implementation ///////////////////////////////////////////////////////
BinaryWriter::BinaryWriter(const string &fileName)
: file_(fileName, ios::binary|ios::out)
{}

template <>
void BinaryWriter::writeDefault(const string &s, const string &output)
{
  uint64_t sz = output.size();
  
  write("", sz);
  for (uint64_t i = 0; i < sz; ++i)
  {
    write("", output[i]);
  }
}

// Reader implementation ///////////////////////////////////////////////////////
BinaryReader::BinaryReader(const string &fileName)
: file_(fileName, ios::binary|ios::in)
{}

template <>
void BinaryReader::readDefault(const string &s, string &output)
{
  uint64_t sz;
  
  read("", sz);
  output.reserve(sz);
  file_.read((char *)output.data(), sz);
}
