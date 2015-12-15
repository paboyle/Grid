#include <Grid.h>

using namespace Grid;
using namespace std;

// Writer implementation ///////////////////////////////////////////////////////
TextWriter::TextWriter(const string &fileName)
: file_(fileName, ios::out)
{}

void TextWriter::push(const string &s)
{
  level_++;
};

void TextWriter::pop(void)
{
  level_--;
};

void TextWriter::indent(void)
{
  for (int i = 0; i < level_; ++i)
  {
    file_ << '\t';
  }
};

// Reader implementation ///////////////////////////////////////////////////////
TextReader::TextReader(const string &fileName)
: file_(fileName, ios::in)
{}

void TextReader::push(const string &s)
{
  level_++;
};

void TextReader::pop(void)
{
  level_--;
};

void TextReader::checkIndent(void)
{
  char c;
  
  for (int i = 0; i < level_; ++i)
  {
    file_.get(c);
    if (c != '\t')
    {
      cerr << "mismatch on tab " << c << " level " << level_;
      cerr << " i "<< i << endl;
      abort();
    }
  }
}

template <>
void TextReader::readDefault(const std::string &s, std::string &output)
{
    checkIndent();
    output.clear();
    getline(file_, output);
}
