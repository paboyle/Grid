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
 Eigen::Tensor of type Scalar_ and rank NumIndices_ (row-major order), together with a name for each index.
 Index names are mutable, but tensor dimensionality is not (size of each dimension is mutable).
 They can be persisted to / restored from disk, by default using tag Name.
 During restore from disk, these validations are performed:
   1) Tensor dimensionality must match
   2) IndexNames are validated against current values
   3) If the tensor has non-zero size, the tensor being loaded must have same extent in each dimension
 ******************************************************************************/

extern const std::string NamedTensorFileExtension;

template<typename Scalar_, int NumIndices_>
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

    // Name of the object and Index names as set in the constructor
    const std::string              &Name;
    const std::vector<std::string> &DefaultIndexNames;

    virtual ~NamedTensor(){};
    // Default constructor (assumes tensor will be loaded from file)
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE NamedTensor(const std::string &Name_, const std::vector<std::string> &IndexNames_)
    : IndexNames{IndexNames_}, Name{Name_}, DefaultIndexNames{IndexNames_} {}
    
    // Construct a named tensor explicitly specifying size of each dimension
    template<typename... IndexTypes>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE NamedTensor(const std::string &Name_, const std::vector<std::string> &IndexNames_,
                                                      Eigen::Index firstDimension, IndexTypes... otherDimensions)
    : tensor(firstDimension, otherDimensions...), IndexNames{IndexNames_}, Name{Name_}, DefaultIndexNames{IndexNames_}
    {
        assert(sizeof...(otherDimensions) + 1 == NumIndices_ && "NamedTensor: dimensions != tensor rank");
    }

    // Do my index names match the default for my type?
    bool ValidateIndexNames( const std::vector<std::string> &CheckNames ) const
    {
        assert( CheckNames.size() == NumIndices_ && "Bug: CheckNames don't match NumIndices_" );
        bool bSame{ IndexNames.size() == NumIndices_ };
        for( std::size_t i = 0; bSame && i < NumIndices_; i++ )
        {
            bSame = IndexNames[i].size() == CheckNames[i].size()
            && std::equal( IndexNames[i].begin(), IndexNames[i].end(), CheckNames[i].begin(),
                          [](const char & c1, const char & c2){ return c1 == c2 || std::toupper(c1) == std::toupper(c2); });
        }
        return bSame;
    }
    bool ValidateIndexNames() const { return ValidateIndexNames(DefaultIndexNames); }

#ifdef HAVE_HDF5
    using Default_Reader = Grid::Hdf5Reader;
    using Default_Writer = Grid::Hdf5Writer;
#else
    using Default_Reader = Grid::BinaryReader;
    using Default_Writer = Grid::BinaryWriter;
#endif
    
    void write(const std::string &FileName, const std::string &Tag) const
    {
        std::string FileName_{FileName};
        FileName_.append( NamedTensorFileExtension );
        LOG(Message) << "Writing " << Name << " to file " << FileName_ << " tag " << Tag << std::endl;
        Default_Writer w( FileName_ );
        write( w, Tag, *this );
    }
    void write(const std::string &FileName) const { return write(FileName, Name); }

    // Read tensor.
    // Validate:
    //  1) index names (if requested)
    //  2) index dimensions (if they are non-zero when called)
    template<typename Reader> void read(Reader &r, bool bValidate, const std::string &Tag)
    {
        // Grab index names and dimensions
        std::vector<std::string> OldIndexNames{std::move(IndexNames)};
        const typename ET::Dimensions OldDimensions{tensor.dimensions()};
        read(r, Tag, *this);
        const typename ET::Dimensions & NewDimensions{tensor.dimensions()};
        for (int i = 0; i < NumIndices_; i++)
            assert(OldDimensions[i] == 0 || OldDimensions[i] == NewDimensions[i] && "NamedTensor::read dimension size");
        if (bValidate)
            assert(ValidateIndexNames(OldIndexNames) && "NamedTensor::read dimension name");
    }
    template<typename Reader> void read(Reader &r, bool bValidate = true) { read(r, bValidate, Name); }

    inline void read (const std::string &FileName, bool bValidate, const std::string &Tag)
    {
        Default_Reader r(FileName + NamedTensorFileExtension);
        read(r, bValidate, Tag);
    }
    inline void read (const std::string &FileName, bool bValidate= true) { return read(FileName, bValidate, Name); }
};

/******************************************************************************
 Common elements for distillation
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(MDistil)

//Eigenvectors of the Laplacian
using LapEvecs = Grid::Hadrons::EigenPack<LatticeColourVector>;

// Noise vector (index order: nnoise, nt, nvec, ns)

class NoiseTensor : public NamedTensor<Complex, 4>
{
    static const std::string               Name_;
    static const std::vector<std::string>  DefaultIndexNames_;
    public:
    // Default constructor (assumes tensor will be loaded from file)
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE NoiseTensor() : NamedTensor{Name_, DefaultIndexNames_} {}

    // Construct a named tensor explicitly specifying size of each dimension
    template<typename... IndexTypes>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE NoiseTensor(Eigen::Index nNoise, Eigen::Index nT, Eigen::Index nVec, Eigen::Index nS)
    : NamedTensor{Name_, DefaultIndexNames_, nNoise, nT, nVec, nS} {}
};

class PerambTensor : public NamedTensor<SpinVector, 6>
{
    static const std::string               Name_;
    static const std::vector<std::string>  DefaultIndexNames_;
    public:
    // Default constructor (assumes tensor will be loaded from file)
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE PerambTensor() : NamedTensor{Name_, DefaultIndexNames_} {}

    // Construct a named tensor explicitly specifying size of each dimension
    template<typename... IndexTypes>
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE PerambTensor(Eigen::Index nT, Eigen::Index nVec, Eigen::Index LI, Eigen::Index nNoise, Eigen::Index nT_inv, Eigen::Index SI)
    : NamedTensor{Name_, DefaultIndexNames_, nT, nVec, LI, nNoise, nT_inv, SI} {}
};

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif // Hadrons_Distil_hpp_
