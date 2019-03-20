#ifndef GRID_SERIALISATION_HDF5_H
#define GRID_SERIALISATION_HDF5_H

#include <stack>
#include <string>
#include <vector>
#include <H5Cpp.h>
#include <Grid/tensors/Tensors.h>
#include "Hdf5Type.h"

#ifndef H5_NO_NAMESPACE
#define H5NS H5
#endif

// default thresold above which datasets are used instead of attributes
#ifndef HDF5_DEF_DATASET_THRES
#define HDF5_DEF_DATASET_THRES 6u
#endif

// name guard for Grid metadata
#define HDF5_GRID_GUARD "_Grid_"

namespace Grid
{
  class Hdf5Writer: public Writer<Hdf5Writer>
  {
  public:
    Hdf5Writer(const std::string &fileName);
    virtual ~Hdf5Writer(void) = default;
    void push(const std::string &s);
    void pop(void);
    void writeDefault(const std::string &s, const char *x);
    template <typename U>
    void writeDefault(const std::string &s, const U &x);
    template <typename U>
    typename std::enable_if<element<std::vector<U>>::is_number, void>::type
    writeDefault(const std::string &s, const std::vector<U> &x);
    template <typename U>
    typename std::enable_if<!element<std::vector<U>>::is_number, void>::type
    writeDefault(const std::string &s, const std::vector<U> &x);
    H5NS::Group & getGroup(void);
  private:
    template <typename U>
    void writeSingleAttribute(const U &x, const std::string &name,
                              const H5NS::DataType &type);
  private:
    std::string              fileName_;
    std::vector<std::string> path_;
    H5NS::H5File             file_;
    H5NS::Group              group_;
    unsigned int             dataSetThres_{HDF5_DEF_DATASET_THRES};
  };
  
  class Hdf5Reader: public Reader<Hdf5Reader>
  {
  public:
    Hdf5Reader(const std::string &fileName, const bool readOnly = true);
    virtual ~Hdf5Reader(void) = default;
    bool push(const std::string &s);
    void pop(void);
    template <typename U>
    void readDefault(const std::string &s, U &output);
    template <typename U>
    typename std::enable_if<element<std::vector<U>>::is_number, void>::type
    readDefault(const std::string &s, std::vector<U> &x);
    template <typename U>
    typename std::enable_if<!element<std::vector<U>>::is_number, void>::type
    readDefault(const std::string &s, std::vector<U> &x);
    H5NS::Group & getGroup(void);
  private:
    template <typename U>
    void readSingleAttribute(U &x, const std::string &name,
                             const H5NS::DataType &type);
  private:
    std::string              fileName_;
    std::vector<std::string> path_;
    H5NS::H5File             file_;
    H5NS::Group              group_;
    unsigned int             dataSetThres_;
  };
  
  // Writer template implementation ////////////////////////////////////////////
  template <typename U>
  void Hdf5Writer::writeSingleAttribute(const U &x, const std::string &name,
                                        const H5NS::DataType &type)
  {
    H5NS::Attribute attribute;
    hsize_t         attrDim = 1;
    H5NS::DataSpace attrSpace(1, &attrDim);
    
    attribute = group_.createAttribute(name, type, attrSpace);
    attribute.write(type, &x);
  }
  
  template <typename U>
  void Hdf5Writer::writeDefault(const std::string &s, const U &x)
  {
    writeSingleAttribute(x, s, Hdf5Type<U>::type());
  }
  
  template <>
  void Hdf5Writer::writeDefault(const std::string &s, const std::string &x);
  
  template <typename U>
  typename std::enable_if<element<std::vector<U>>::is_number, void>::type
  Hdf5Writer::writeDefault(const std::string &s, const std::vector<U> &x)
  {
    // alias to element type
    typedef typename element<std::vector<U>>::type Element;
    
    // flatten the vector and getting dimensions
    Flatten<std::vector<U>> flat(x);
    std::vector<hsize_t> dim;
    const auto           &flatx = flat.getFlatVector();
    
    for (auto &d: flat.getDim())
    {
      dim.push_back(d);
    }
    
    // write to file
    H5NS::DataSpace dataSpace(dim.size(), dim.data());
    
    if (flatx.size() > dataSetThres_)
    {
      H5NS::DataSet           dataSet;
      H5NS::DSetCreatPropList plist;
      
      plist.setChunk(dim.size(), dim.data());
      plist.setFletcher32();
      dataSet = group_.createDataSet(s, Hdf5Type<Element>::type(), dataSpace, plist);
      dataSet.write(flatx.data(), Hdf5Type<Element>::type());
    }
    else
    {
      H5NS::Attribute attribute;
      
      attribute = group_.createAttribute(s, Hdf5Type<Element>::type(), dataSpace);
      attribute.write(Hdf5Type<Element>::type(), flatx.data());
    }
  }
  
  template <typename U>
  typename std::enable_if<!element<std::vector<U>>::is_number, void>::type
  Hdf5Writer::writeDefault(const std::string &s, const std::vector<U> &x)
  {
    push(s);
    writeSingleAttribute(x.size(), HDF5_GRID_GUARD "vector_size",
                         Hdf5Type<uint64_t>::type());
    for (hsize_t i = 0; i < x.size(); ++i)
    {
      write(s + "_" + std::to_string(i), x[i]);
    }
    pop();
  }
  
  // Reader template implementation ////////////////////////////////////////////
  template <typename U>
  void Hdf5Reader::readSingleAttribute(U &x, const std::string &name,
                                       const H5NS::DataType &type)
  {
    H5NS::Attribute attribute;
    
    attribute = group_.openAttribute(name);
    attribute.read(type, &x);
  }
  
  template <typename U>
  void Hdf5Reader::readDefault(const std::string &s, U &output)
  {
    readSingleAttribute(output, s, Hdf5Type<U>::type());
  }
  
  template <>
  void Hdf5Reader::readDefault(const std::string &s, std::string &x);
  
  template <typename U>
  typename std::enable_if<element<std::vector<U>>::is_number, void>::type
  Hdf5Reader::readDefault(const std::string &s, std::vector<U> &x)
  {
    // alias to element type
    typedef typename element<std::vector<U>>::type Element;
    
    // read the dimensions
    H5NS::DataSpace       dataSpace;
    std::vector<hsize_t>  hdim;
    std::vector<size_t>   dim;
    hsize_t               size = 1;
    
    if (group_.attrExists(s))
    {
      dataSpace = group_.openAttribute(s).getSpace();
    }
    else
    {
      dataSpace = group_.openDataSet(s).getSpace();
    }
    hdim.resize(dataSpace.getSimpleExtentNdims());
    dataSpace.getSimpleExtentDims(hdim.data());
    for (auto &d: hdim)
    {
      dim.push_back(d);
      size *= d;
    }
    
    // read the flat vector
    std::vector<Element> buf(size);

    if (size > dataSetThres_)
    {
      H5NS::DataSet dataSet;
      
      dataSet = group_.openDataSet(s);
      dataSet.read(buf.data(), Hdf5Type<Element>::type());
    }
    else
    {
      H5NS::Attribute attribute;
      
      attribute = group_.openAttribute(s);
      attribute.read(Hdf5Type<Element>::type(), buf.data());
    }
    
    // reconstruct the multidimensional vector
    Reconstruct<std::vector<U>> r(buf, dim);
    
    x = r.getVector();
  }
  
  template <typename U>
  typename std::enable_if<!element<std::vector<U>>::is_number, void>::type
  Hdf5Reader::readDefault(const std::string &s, std::vector<U> &x)
  {
    uint64_t size;
    
    push(s);
    readSingleAttribute(size, HDF5_GRID_GUARD "vector_size",
                        Hdf5Type<uint64_t>::type());
    x.resize(size);
    for (hsize_t i = 0; i < x.size(); ++i)
    {
      read(s + "_" + std::to_string(i), x[i]);
    }
    pop();
  }
}

#endif
