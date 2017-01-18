#ifndef GRID_SERIALISATION_HDF5_H
#define GRID_SERIALISATION_HDF5_H

#include <stack>
#include <string>
#include <vector>
#include <H5Cpp.h>
#include "Hdf5Type.h"

#ifndef H5_NO_NAMESPACE
#define H5NS H5
#endif

// default thresold above which datasets are used instead of attributes
#ifndef H5_DEF_DATASET_THRES
#define H5_DEF_DATASET_THRES 6u
#endif

namespace Grid
{
  template <typename T>
  struct is_arithmetic_vector
  {
    static constexpr bool value = false;
  };
  
  template <typename T>
  struct is_arithmetic_vector<std::vector<T>>
  {
    static constexpr bool value = std::is_arithmetic<T>::value
                                  or is_arithmetic_vector<T>::value;
  };
  
  class Hdf5Writer: public Writer<Hdf5Writer>
  {
  public:
    Hdf5Writer(const std::string &fileName);
    virtual ~Hdf5Writer(void);
    void push(const std::string &s);
    void pop(void);
    void writeDefault(const std::string &s, const char *x);
    template <typename U>
    void writeDefault(const std::string &s, const U &x);
    template <typename U>
    typename std::enable_if<is_arithmetic_vector<std::vector<U>>::value
                            and std::is_arithmetic<U>::value, void>::type
    writeDefault(const std::string &s, const std::vector<U> &x);
    template <typename U>
    typename std::enable_if<is_arithmetic_vector<std::vector<U>>::value
                            and !std::is_arithmetic<U>::value, void>::type
    writeDefault(const std::string &s, const std::vector<U> &x);
    template <typename U>
    typename std::enable_if<!is_arithmetic_vector<std::vector<U>>::value, void>::type
    writeDefault(const std::string &s, const std::vector<U> &x);
  private:
    std::string              fileName_;
    std::vector<std::string> path_;
    std::vector<hsize_t>     dim_;
    bool                     multiDim_{true};
    H5NS::H5File             file_;
    H5NS::Group              group_;
    unsigned int             datasetThres_{H5_DEF_DATASET_THRES};
  };
  
  class Hdf5Reader: public Reader<Hdf5Reader>
  {
  public:
    Hdf5Reader(const std::string &fileName);
    virtual ~Hdf5Reader(void);
    void push(const std::string &s);
    void pop(void);
    template <typename U>
    void readDefault(const std::string &s, U &output);
    template <typename U>
    void readDefault(const std::string &s, std::vector<U> &output);
  private:
  };
  
  // Writer template implementation ////////////////////////////////////////////
  template <typename U>
  void Hdf5Writer::writeDefault(const std::string &s, const U &x)
  {
    H5NS::Attribute attribute;
    hsize_t         attrDim = 1;
    H5NS::DataSpace attrSpace(1, &attrDim);
    
    attribute = group_.createAttribute(s, *Hdf5Type<U>::type, attrSpace);
    attribute.write(*Hdf5Type<U>::type, &x);
  }
  
  template <>
  void Hdf5Writer::writeDefault(const std::string &s, const std::string &x);
  
  template <typename U>
  typename std::enable_if<is_arithmetic_vector<std::vector<U>>::value
                          and std::is_arithmetic<U>::value, void>::type
  Hdf5Writer::writeDefault(const std::string &s, const std::vector<U> &x)
  {
    hsize_t size = 1;
    
    dim_.push_back(x.size());
    for (auto d: dim_)
    {
      size *= d;
    }
    
    H5NS::DataSpace dataspace(dim_.size(), dim_.data());
    
    if (size > datasetThres_)
    {
      H5NS::DataSet dataset;
      
      dataset = group_.createDataSet(s, *Hdf5Type<U>::type, dataspace);
      dataset.write(x.data(), *Hdf5Type<U>::type);
    }
    else
    {
      H5NS::Attribute attribute;
      
      attribute = group_.createAttribute(s, *Hdf5Type<U>::type, dataspace);
      attribute.write(*Hdf5Type<U>::type, x.data());
    }
    dim_.clear();
    multiDim_ = true;
  }
  
  template <typename U>
  typename std::enable_if<is_arithmetic_vector<std::vector<U>>::value
                          and !std::is_arithmetic<U>::value, void>::type
  Hdf5Writer::writeDefault(const std::string &s, const std::vector<U> &x)
  {
    hsize_t firstSize = x[0].size();
    
    for (auto &v: x)
    {
      multiDim_ = (multiDim_ and (v.size() == firstSize));
    }
    assert(multiDim_);
    dim_.push_back(x.size());
    writeDefault(s, x[0]);
  }
  
  template <typename U>
  typename std::enable_if<!is_arithmetic_vector<std::vector<U>>::value, void>::type
  Hdf5Writer::writeDefault(const std::string &s, const std::vector<U> &x)
  {
    push(s);
    for (hsize_t i = 0; i < x.size(); ++i)
    {
      write(s + "_" + std::to_string(i), x[i]);
    }
    pop();
  }
  
  // Reader template implementation ////////////////////////////////////////////
  template <typename U>
  void Hdf5Reader::readDefault(const std::string &s, U &output)
  {
    
  }
  
  template <typename U>
  void Hdf5Reader::readDefault(const std::string &s, std::vector<U> &output)
  {
    
  }
}

#endif
