/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: ./Grid/serialisation/VectorUtils.h
 
 Copyright (C) 2015
 
 Author: Peter Boyle <paboyle@ed.ac.uk>
 Author: Antonin Portelli <antonin.portelli@me.com>
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

#ifndef GRID_SERIALISATION_HDF5_H
#define GRID_SERIALISATION_HDF5_H

#include <stack>
#include <string>
#include <list>
#include <vector>
#include <H5Cpp.h>
#include <Grid/tensors/Tensors.h>
#include "Hdf5Type.h"

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
    void writeRagged(const std::string &s, const std::vector<U> &x);
    template <typename U>
    typename std::enable_if<is_flattenable<std::vector<U>>::value>::type
    writeDefault(const std::string &s, const std::vector<U> &x);
    template <typename U>
    typename std::enable_if<!is_flattenable<std::vector<U>>::value>::type
    writeDefault(const std::string &s, const std::vector<U> &x) { writeRagged(s, x); }
    template <typename U>
    void writeMultiDim(const std::string &s, const std::vector<size_t> & Dimensions, const U * pDataRowMajor, size_t NumElements);
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
    const unsigned int       dataSetThres_{HDF5_DEF_DATASET_THRES};
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
    void readRagged(const std::string &s, std::vector<U> &x);
    template <typename U>
    typename std::enable_if<is_flattenable<std::vector<U>>::value>::type
    readDefault(const std::string &s, std::vector<U> &x);
    template <typename U>
    typename std::enable_if<!is_flattenable<std::vector<U>>::value>::type
    readDefault(const std::string &s, std::vector<U> &x) { readRagged(s, x); }
    template <typename U>
    void readMultiDim(const std::string &s, std::vector<U> &buf, std::vector<size_t> &dim);
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
  void Hdf5Writer::writeMultiDim(const std::string &s, const std::vector<size_t> & Dimensions, const U * pDataRowMajor, size_t NumElements)
  {
    // Hdf5 needs the dimensions as hsize_t
    const int rank = static_cast<int>(Dimensions.size());
    std::vector<hsize_t> dim(rank);
    for(int i = 0; i < rank; i++)
      dim[i] = Dimensions[i];
    // write the entire dataset to file
    H5NS::DataSpace dataSpace(rank, dim.data());

    if (NumElements > dataSetThres_)
    {
      // Make sure 1) each dimension; and 2) chunk size is < 4GB
      const hsize_t MaxElements = ( sizeof( U ) == 1 ) ? 0xffffffff : 0x100000000 / sizeof( U );
      hsize_t ElementsPerChunk = 1;
      bool bTooBig = false;
      for( int i = rank - 1 ; i != -1 ; i-- ) {
        auto &d = dim[i];
        if( bTooBig )
          d = 1; // Chunk size is already as big as can be - remaining dimensions = 1
        else {
          // If individual dimension too big, reduce by prime factors if possible
          while( d > MaxElements && ( d & 1 ) == 0 )
            d >>= 1;
          const char ErrorMsg[] = " dimension > 4GB and not divisible by 2^n. "
                                  "Hdf5IO chunk size will be inefficient. NB Serialisation is not intended for large datasets - please consider alternatives.";
          if( d > MaxElements ) {
            std::cout << GridLogWarning << "Individual" << ErrorMsg << std::endl;
            hsize_t quotient = d / MaxElements;
            if( d % MaxElements )
              quotient++;
            d /= quotient;
          }
          // Now make sure overall size is not too big
          hsize_t OverflowCheck = ElementsPerChunk;
          ElementsPerChunk *= d;
          assert( OverflowCheck == ElementsPerChunk / d && "Product of dimensions overflowed hsize_t" );
          // If product of dimensions too big, reduce by prime factors
          while( ElementsPerChunk > MaxElements && ( ElementsPerChunk & 1 ) == 0 ) {
            bTooBig = true;
            d >>= 1;
            ElementsPerChunk >>= 1;
          }
          if( ElementsPerChunk > MaxElements ) {
            std::cout << GridLogWarning << "Product of" << ErrorMsg << std::endl;
            hsize_t quotient = ElementsPerChunk / MaxElements;
            if( ElementsPerChunk % MaxElements )
              quotient++;
            d /= quotient;
            ElementsPerChunk /= quotient;
          }
        }
      }
      H5NS::DataSet           dataSet;
      H5NS::DSetCreatPropList plist;
      plist.setChunk(rank, dim.data());
      plist.setFletcher32();
      dataSet = group_.createDataSet(s, Hdf5Type<U>::type(), dataSpace, plist);
      dataSet.write(pDataRowMajor, Hdf5Type<U>::type());
    }
    else
    {
      H5NS::Attribute attribute;
      attribute = group_.createAttribute(s, Hdf5Type<U>::type(), dataSpace);
      attribute.write(Hdf5Type<U>::type(), pDataRowMajor);
    }
  }

  template <typename U>
  typename std::enable_if<is_flattenable<std::vector<U>>::value>::type
  Hdf5Writer::writeDefault(const std::string &s, const std::vector<U> &x)
  {
    if (isRegularShape(x))
    {
      // alias to element type
      using Scalar = typename is_flattenable<std::vector<U>>::type;
      
      // flatten the vector and getting dimensions
      Flatten<std::vector<U>> flat(x);
      std::vector<size_t> dim;
      const auto           &flatx = flat.getFlatVector();
      for (auto &d: flat.getDim())
        dim.push_back(d);
      writeMultiDim<Scalar>(s, dim, &flatx[0], flatx.size());
    }
    else
    {
      writeRagged(s, x);
    }
  }
  
  template <typename U>
  void Hdf5Writer::writeRagged(const std::string &s, const std::vector<U> &x)
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
  void Hdf5Reader::readMultiDim(const std::string &s, std::vector<U> &buf, std::vector<size_t> &dim)
  {
    // alias to element type
    using Scalar = typename is_flattenable<std::vector<U>>::type;
    
    // read the dimensions
    H5NS::DataSpace       dataSpace;
    std::vector<hsize_t>  hdim;
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
    buf.resize(size);
    
    if (size > dataSetThres_)
    {
      H5NS::DataSet dataSet;
      
      dataSet = group_.openDataSet(s);
      dataSet.read(buf.data(), Hdf5Type<Scalar>::type());
    }
    else
    {
      H5NS::Attribute attribute;
      
      attribute = group_.openAttribute(s);
      attribute.read(Hdf5Type<Scalar>::type(), buf.data());
    }
  }

  template <typename U>
  typename std::enable_if<is_flattenable<std::vector<U>>::value>::type
  Hdf5Reader::readDefault(const std::string &s, std::vector<U> &x)
  {
    if (H5Lexists        (group_.getId(), s.c_str(), H5P_DEFAULT) > 0
     && H5Aexists_by_name(group_.getId(), s.c_str(), HDF5_GRID_GUARD "vector_size", H5P_DEFAULT ) > 0)
    {
      readRagged(s, x);
    }
    else
    {
      // alias to element type
      using Scalar = typename is_flattenable<std::vector<U>>::type;

      std::vector<size_t>   dim;
      std::vector<Scalar>   buf;
      readMultiDim( s, buf, dim );

      // reconstruct the multidimensional vector
      Reconstruct<std::vector<U>> r(buf, dim);

      x = r.getVector();
    }
  }
  
  template <typename U>
  void Hdf5Reader::readRagged(const std::string &s, std::vector<U> &x)
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
