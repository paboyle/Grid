#include <Grid.h>

namespace Grid {
// Writer implementation ///////////////////////////////////////////////////////
XmlWriter::XmlWriter(const std::string &fileName)
: fileName_(fileName)
{
  node_ = doc_.append_child();
  node_.set_name("grid");
}

XmlWriter::~XmlWriter(void)
{
  doc_.save_file(fileName_.c_str(), "  ");
}

void XmlWriter::push(const std::string &s)
{
  node_ = node_.append_child(s.c_str());
}

void XmlWriter::pop(void)
{
  node_ = node_.parent();
}

// Reader implementation ///////////////////////////////////////////////////////
XmlReader::XmlReader(const std::string &fileName)
: fileName_(fileName)
{
  pugi::xml_parse_result result = doc_.load_file(fileName_.c_str());
  
  if ( !result )
  {
    std::cerr << "XML error description: " << result.description() << "\n";
    std::cerr << "XML error offset     : " << result.offset        << "\n";
    std::abort();
  }
  
  node_ = doc_.child("grid");
}

void XmlReader::push(const std::string &s)
{
  node_ = node_.child(s.c_str());
}

void XmlReader::pop(void)
{
  node_ = node_.parent();
}

template <>
void XmlReader::readDefault(const std::string &s, std::string &output)
{
  output = node_.child(s.c_str()).first_child().value();
}
}
