#include <Grid.h>

using namespace Grid;
using namespace std;

// Writer implementation ///////////////////////////////////////////////////////
XmlWriter::XmlWriter(const string &fileName)
: fileName_(fileName)
{
  node_ = doc_.append_child();
  node_.set_name("grid");
}

XmlWriter::~XmlWriter(void)
{
  doc_.save_file(fileName_.c_str(), "  ");
}

void XmlWriter::push(const string &s)
{
  node_ = node_.append_child(s.c_str());
}

void XmlWriter::pop(void)
{
  node_ = node_.parent();
}

// Reader implementation ///////////////////////////////////////////////////////
XmlReader::XmlReader(const string &fileName)
: fileName_(fileName)
{
  pugi::xml_parse_result result = doc_.load_file(fileName_.c_str());
  
  if ( !result )
  {
    cerr << "XML error description: " << result.description() << "\n";
    cerr << "XML error offset     : " << result.offset        << "\n";
    abort();
  }
  
  node_ = doc_.child("grid");
}

void XmlReader::push(const string &s)
{
  node_ = node_.child(s.c_str());
}

void XmlReader::pop(void)
{
  node_ = node_.parent();
}

bool XmlReader::nextElement(const std::string &s)
{
  if (node_.next_sibling(s.c_str()))
  {
    node_ = node_.next_sibling(s.c_str());
    
    return true;
  }
  else
  {
    return false;
  }
}

template <>
void XmlReader::readDefault(const string &s, string &output)
{
  output = node_.child(s.c_str()).first_child().value();
}
