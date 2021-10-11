/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: ./Grid/serialisation/VectorUtils.h
 
 Copyright (C) 2015
 
 Author: Antonin Portelli <antonin.portelli@me.com>
 Author: Peter Boyle <paboyle@ed.ac.uk>
 Author: Guido Cossu <guido.cossu@ed.ac.uk>
 Author: Michael Marshall <michael.marshall@ed.ac.uk>

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

#include <Grid/Grid.h>

using namespace Grid;
#ifndef H5_NO_NAMESPACE
using namespace H5NS; // Compile error here? Try adding --enable-cxx to hdf5 configure
#endif

// Writer implementation ///////////////////////////////////////////////////////
Hdf5Writer::Hdf5Writer(const std::string &fileName)
: fileName_(fileName)
, file_(fileName.c_str(), H5F_ACC_TRUNC)
{
  group_ = file_.openGroup("/");
  writeSingleAttribute(dataSetThres_, HDF5_GRID_GUARD "dataset_threshold",
                       Hdf5Type<unsigned int>::type());
}

void Hdf5Writer::push(const std::string &s)
{
  group_ = group_.createGroup(s);
  path_.push_back(s);
}

void Hdf5Writer::pop(void)
{
  path_.pop_back();
  if (path_.empty())
  {
    group_ = file_.openGroup("/");
  }
  else
  {
    auto binOp = [](const std::string &a, const std::string &b)->std::string
    {
      return a + "/" + b;
    };
    
    group_ = group_.openGroup(std::accumulate(path_.begin(), path_.end(),
                                              std::string(""), binOp));
  }
}

template <>
void Hdf5Writer::writeDefault(const std::string &s, const std::string &x)
{
  StrType     strType(PredType::C_S1, x.size());
  
  writeSingleAttribute(*(x.data()), s, strType);
}

void Hdf5Writer::writeDefault(const std::string &s, const char *x)
{
  std::string sx(x);
  
  writeDefault(s, sx);
}

Group & Hdf5Writer::getGroup(void)
{
  return group_;
}

// Reader implementation ///////////////////////////////////////////////////////
Hdf5Reader::Hdf5Reader(const std::string &fileName, const bool readOnly)
: fileName_(fileName)
, file_(fileName.c_str(), readOnly ? H5F_ACC_RDONLY : H5F_ACC_RDWR)
{
  group_ = file_.openGroup("/");
  readSingleAttribute(dataSetThres_, HDF5_GRID_GUARD "dataset_threshold",
                      Hdf5Type<unsigned int>::type());
}

bool Hdf5Reader::push(const std::string &s)
{
  group_ = group_.openGroup(s);
  path_.push_back(s);
  
  return true;
}

void Hdf5Reader::pop(void)
{
  path_.pop_back();
  if (path_.empty())
  {
    group_ = file_.openGroup("/");
  }
  else
  {
    auto binOp = [](const std::string &a, const std::string &b)->std::string
    {
      return a + "/" + b;
    };
    
    group_ = group_.openGroup(std::accumulate(path_.begin(), path_.end(),
                                              std::string(""), binOp));
  }
}

template <>
void Hdf5Reader::readDefault(const std::string &s, std::string &x)
{
  Attribute attribute;
  
  attribute       = group_.openAttribute(s);
  StrType strType = attribute.getStrType();
  
  x.resize(strType.getSize());
  attribute.read(strType, &(x[0]));
}

Group & Hdf5Reader::getGroup(void)
{
  return group_;
}
