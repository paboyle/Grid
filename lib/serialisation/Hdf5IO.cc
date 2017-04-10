#include <Grid/Grid.h>

using namespace Grid;
#ifndef H5_NO_NAMESPACE
using namespace H5NS;
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

// Reader implementation ///////////////////////////////////////////////////////
Hdf5Reader::Hdf5Reader(const std::string &fileName)
: fileName_(fileName)
, file_(fileName.c_str(), H5F_ACC_RDONLY)
{
  group_ = file_.openGroup("/");
  readSingleAttribute(dataSetThres_, HDF5_GRID_GUARD "dataset_threshold",
                      Hdf5Type<unsigned int>::type());
}

void Hdf5Reader::push(const std::string &s)
{
  group_ = group_.openGroup(s);
  path_.push_back(s);
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
