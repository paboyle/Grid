#include <Grid.h>

namespace Grid {
// Writer implementation ///////////////////////////////////////////////////////
TextWriter::TextWriter(const std::string &fileName)
: file_(fileName, std::ios::out)
{}

void TextWriter::push(const std::string &s)
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
TextReader::TextReader(const std::string &fileName)
: file_(fileName, std::ios::in)
{}

void TextReader::push(const std::string &s)
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
      std::cerr << "mismatch on tab " << c << " level " << level_;
      std::cerr << " i "<< i <<std::endl;
      std::abort();
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
}
