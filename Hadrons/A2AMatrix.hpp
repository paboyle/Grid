/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/A2AMatrix.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>

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
#ifndef A2A_Matrix_hpp_
#define A2A_Matrix_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

template <typename T, typename MetadataType>
class A2AMatrixIo
{
public:
    A2AMatrixIo(void) = default;
    A2AMatrixIo(std::string filename, std::string dataname, 
                const unsigned int nt, const unsigned int ni,
                const unsigned int nj, const unsigned int blockSize);
    ~A2AMatrixIo(void) = default;
    void initFile(const MetadataType &d);
    void saveBlock(const T *data, const unsigned int i, const unsigned int j);
private:
    std::string  filename_, dataname_;
    unsigned int nt_, ni_, nj_, blockSize_;
};

template <typename T, typename MetadataType>
A2AMatrixIo<T, MetadataType>::A2AMatrixIo(std::string filename, 
                                          std::string dataname, 
                                          const unsigned int nt, 
                                          const unsigned int ni,
                                          const unsigned int nj,
                                          const unsigned int blockSize)
: filename_(filename), dataname_(dataname)
, nt_(nt), ni_(ni), nj_(nj), blockSize_(blockSize)
{}

template <typename T, typename MetadataType>
void A2AMatrixIo<T, MetadataType>::initFile(const MetadataType &d)
{
#ifdef HAVE_HDF5
    std::vector<hsize_t>    dim = {static_cast<hsize_t>(nt_), 
                                   static_cast<hsize_t>(ni_), 
                                   static_cast<hsize_t>(nj_)},
                            chunk = {static_cast<hsize_t>(nt_), 
                                     static_cast<hsize_t>(blockSize_), 
                                     static_cast<hsize_t>(blockSize_)};
    H5NS::DataSpace         dataspace(dim.size(), dim.data());
    H5NS::DataSet           dataset;
    H5NS::DSetCreatPropList plist;
    
    // create empty file just with metadata
    {
        Hdf5Writer writer(filename_);
        write(writer, dataname_, d);
    }

    // create the dataset
    Hdf5Reader reader(filename_);

    push(reader, dataname_);
    auto &group = reader.getGroup();
    plist.setChunk(chunk.size(), chunk.data());
    dataset = group.createDataSet("data", Hdf5Type<T>::type(), dataspace, plist);
#else
    HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}

template <typename T, typename MetadataType>
void A2AMatrixIo<T, MetadataType>::saveBlock(const T *data, 
                                             const unsigned int i, 
                                             const unsigned int j)
{
#ifdef HAVE_HDF5
    Hdf5Reader           reader(filename_);
    std::vector<hsize_t> count = {nt_, blockSize_, blockSize_},
                         offset = {0, static_cast<hsize_t>(i),
                                   static_cast<hsize_t>(j)},
                         stride = {1, 1, 1},
                         block  = {1, 1, 1}; 
    H5NS::DataSpace      memspace(count.size(), count.data()), dataspace;
    H5NS::DataSet        dataset;
    size_t               shift;

    push(reader, dataname_);
    auto &group = reader.getGroup();
    dataset     = group.openDataSet("data");
    dataspace   = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                              stride.data(), block.data());
    dataset.write(data, Hdf5Type<T>::type(), memspace, dataspace);
#else
    HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}

END_HADRONS_NAMESPACE

#endif // A2A_Matrix_hpp_
