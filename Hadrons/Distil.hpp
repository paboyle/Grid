/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Distil.hpp
 
 Copyright (C) 2015-2019
 
 Author: Felix Erben <ferben@ed.ac.uk>
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 
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

#ifndef Hadrons_Distil_hpp_
#define Hadrons_Distil_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 NamedTensor
 This is an Eigen::Tensor of type Scalar_ and rank NumIndices_ (row-major order)
 They can be persisted to disk in tag Name_, and IndexNames are validated on load.
 TODO: WHAT TO SAVE / VALIDATE ON LOAD (Override to warn instead of assert on load)
 Ensemble string
 Configuration number
 Noise unique string
 Distillation parameters
 ******************************************************************************/

extern const std::string NamedTensorFileExtension;

template<typename Scalar_, int NumIndices_, const std::string &Name_, const std::array<std::string,NumIndices_> &IndexNames_>
class NamedTensor : Serializable
{
public:
    using Scalar = Scalar_;
    static constexpr int NumIndices = NumIndices_;
    using ET = Eigen::Tensor<Scalar_, NumIndices_, Eigen::RowMajor>;
    using Index = typename ET::Index;
    GRID_SERIALIZABLE_CLASS_MEMBERS(NamedTensor,
                                    ET,                       tensor,
                                    std::vector<std::string>, IndexNames );
    
    // Get the default index names as std::vector
    std::vector<std::string> DefaultIndexNames()
    {
        std::vector<std::string> names{NumIndices_};
        for (std::size_t i = 0; i < NumIndices_; i++)
            names[i] = IndexNames_[i];
        return names;
    }
    
    // Default constructor (assumes tensor will be loaded from file)
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE NamedTensor() : IndexNames{DefaultIndexNames()} {}
    
    // Construct a named tensor explicitly specifying size of each dimension
    template<typename... IndexTypes>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE NamedTensor(Eigen::Index firstDimension, IndexTypes... otherDimensions)
    : tensor(firstDimension, otherDimensions...), IndexNames{DefaultIndexNames()}
    {
        assert(sizeof...(otherDimensions) + 1 == NumIndices_ && "NamedTensor: dimensions != tensor rank");
    }
    
    // Do my index names match the default for my type?
    bool ValidateIndexNames() const
    {
        bool bSame{ IndexNames.size() == NumIndices_ };
        for( std::size_t i = 0; bSame && i < NumIndices_; i++ )
        {
            bSame = IndexNames[i].size() == IndexNames_[i].size()
            && std::equal( IndexNames[i].begin(), IndexNames[i].end(), IndexNames_[i].begin(),
                          [](const char & c1, const char & c2){ return c1 == c2 || std::toupper(c1) == std::toupper(c2); });
        }
        return bSame;
    }
    
#ifdef HAVE_HDF5
    using Default_Reader = Grid::Hdf5Reader;
    using Default_Writer = Grid::Hdf5Writer;
#else
    using Default_Reader = Grid::BinaryReader;
    using Default_Writer = Grid::BinaryWriter;
#endif
    
    template<typename Writer> void write(Writer &w, const std::string &Tag = Name_) const
    { write(w, Tag, *this); }
    
    inline void write(const std::string &filename, const std::string &Tag = Name_) const
    {
        std::string sFileName{filename};
        sFileName.append( NamedTensorFileExtension );
        LOG(Message) << "Writing " << Name_ << " to file " << sFileName << " tag " << Tag << std::endl;
        Default_Writer w( sFileName );
        write( w, Tag );
    }
    
    // Read and validate index names
    template<typename Reader> void read(Reader &r, bool bValidate = true, const std::string &Tag = Name_)
    {
        // Grab index names and dimensions
        std::vector<std::string> OldIndexNames{std::move(IndexNames)};
        typename ET::Dimensions OldDimensions{tensor.dimensions()};
        read(r, Tag, *this);
        const typename ET::Dimensions & NewDimensions{tensor.dimensions()};
        for (int i = 0; i < NumIndices_; i++)
            assert(OldDimensions[i] == 0 || OldDimensions[i] == NewDimensions[i] && "NamedTensor::read dimension size");
        if (bValidate)
            assert(ValidateIndexNames() && "NamedTensor::read dimension name");
    }
    
    inline void read (const std::string &filename, bool bValidate = true, const std::string &Tag = Name_)
    {
        std::string sFileName{filename};
        sFileName.append( NamedTensorFileExtension );
        LOG(Message) << "Reading " << Name_ << " from file " << sFileName << " tag " << Tag << std::endl;
        Default_Reader r(sFileName);
        read(r, bValidate, Tag);
    }
};

/******************************************************************************
 Common elements for distillation
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(MDistil)

//Eigenvectors of the Laplacian
using LapEvecs = Grid::Hadrons::EigenPack<LatticeColourVector>;

// Noise vector (index order: nnoise, nt, nvec, ns)
using NoiseTensor = Eigen::Tensor<Complex, 4, Eigen::RowMajor>;

extern const std::string PerambTensorName;
extern const std::array<std::string, 6> PerambIndexNames;
using PerambTensor = NamedTensor<SpinVector, 6, PerambTensorName, PerambIndexNames>;

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif // Hadrons_Distil_hpp_
