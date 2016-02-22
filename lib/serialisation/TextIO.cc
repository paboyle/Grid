    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/serialisation/TextIO.cc

    Copyright (C) 2015

Author: Antonin Portelli <antonin.portelli@me.com>
Author: paboyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
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
