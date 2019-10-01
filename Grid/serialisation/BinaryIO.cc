    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/serialisation/BinaryIO.cc

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
#include <Grid/GridCore.h>

NAMESPACE_BEGIN(Grid);

// Writer implementation ///////////////////////////////////////////////////////
BinaryWriter::BinaryWriter(const std::string &fileName)
: file_(fileName, std::ios::binary|std::ios::out)
{}

template <>
void BinaryWriter::writeDefault(const std::string &s, const std::string &x)
{
    uint64_t sz = x.size();
    
    write("", sz);
    for (uint64_t i = 0; i < sz; ++i)
    {
        write("", x[i]);
    }
}

void BinaryWriter::writeDefault(const std::string &s, const char *x)
{
  std::string sx(x);
  
  writeDefault(s, sx);
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
    output.resize(sz);
    file_.read((char *)output.data(), sz);
}

NAMESPACE_END(Grid);
