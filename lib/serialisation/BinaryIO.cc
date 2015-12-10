#include <Grid.h>

using namespace Grid;
using namespace std;

// Writer implementation ///////////////////////////////////////////////////////
BinaryWriter::BinaryWriter(const string &fileName)
: file_(fileName, ios::binary|ios::out)
{}

template <>
void BinaryWriter::writeDefault(const string &s, const string &x)
{
    uint64_t sz = x.size();
    
    write("", sz);
    for (uint64_t i = 0; i < sz; ++i)
    {
        write("", x[i]);
    }
}

void BinaryWriter::writeDefault(const string &s, const char *x)
{
  string sx(x);
  
  writeDefault(s, sx);
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
    output.resize(sz);
    file_.read((char *)output.data(), sz);
}
