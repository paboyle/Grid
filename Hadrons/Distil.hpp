/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/Distil.hpp

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

#ifndef Hadrons_MDistil_Distil_hpp_
#define Hadrons_MDistil_Distil_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 Common elements for distillation
 ******************************************************************************/

using LapEvecs = Grid::Hadrons::EigenPack<LatticeColourVector>;

// Noise vector index order: nnoise, nt, nvec, ns
using NoiseTensor = Eigen::Tensor<Complex, 4, Eigen::RowMajor>;

struct DistilParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(DistilParameters,
                                  int, nnoise,
                                  int, tsrc,
                                  std::string, TI,
                                  std::string, LI,
                                  std::string, SI )
  DistilParameters() = default;
  template <class ReaderClass> DistilParameters(Reader<ReaderClass>& Reader){read(Reader,"Distil",*this);}
  
  // Numeric parameter is allowed to be empty (in which case it = Default),
  // but assert during setup() if specified but not numeric
  
  static int ParameterDefault( const std::string & s, int Default, bool bCalledFromSetup )
  {
    int i = Default;
    if( s.length() > 0 ) {
      std::istringstream ss( s );
      ss >> i;
      if( bCalledFromSetup )
        assert( !ss.fail() && "Parameter should either be empty or integer" );
    }
    return i;
  }
};

#define DISTIL_PARAMETERS_DEFINE( inSetup ) \
const int Nt{env().getDim(Tdir)}; \
const int nvec{par().nvec}; \
const int nnoise{par().Distil.nnoise}; \
const int tsrc{par().Distil.tsrc}; \
const int TI{Hadrons::MDistil::DistilParameters::ParameterDefault(par().Distil.TI, Nt,   inSetup)}; \
const int LI{Hadrons::MDistil::DistilParameters::ParameterDefault(par().Distil.LI, nvec, inSetup)}; \
const int SI{Hadrons::MDistil::DistilParameters::ParameterDefault(par().Distil.SI, Ns,   inSetup)}; \
const bool full_tdil{ TI == Nt }; \
const bool exact_distillation{ full_tdil && LI == nvec }; \
const int Nt_inv{ full_tdil ? 1 : TI }

/******************************************************************************
 Default for distillation file operations. For now only used by NamedTensor
 ******************************************************************************/

#ifdef HAVE_HDF5
using Default_Reader = Grid::Hdf5Reader;
using Default_Writer = Grid::Hdf5Writer;
static const char * FileExtension = ".h5";
#else
using Default_Reader = Grid::BinaryReader;
using Default_Writer = Grid::BinaryWriter;
static const char * FileExtension = ".dat";
#endif

/******************************************************************************
 NamedTensor object
 This is an Eigen::Tensor of type Scalar_ and rank NumIndices_ (row-major order)
 They can be persisted to disk
 IndexNames contains one name for each index, and IndexNames are validated on load.
 WHAT TO SAVE / VALIDATE ON LOAD (Override to warn instead of assert on load)
 Ensemble string
 Configuration number
 Noise unique string
 Distillation parameters
 
 ******************************************************************************/

template<typename Scalar_, int NumIndices_>
class NamedTensor : Serializable
{
public:
  using Scalar = Scalar_;
  static constexpr int NumIndices = NumIndices_;
  using ET = Eigen::Tensor<Scalar_, NumIndices_, Eigen::RowMajor>;
  using Index = typename ET::Index;
  GRID_SERIALIZABLE_CLASS_MEMBERS(NamedTensor
                                  , ET, tensor
                                  , std::vector<std::string>, IndexNames
                                  );
  // Named tensors are intended to be a superset of Eigen tensor
  inline operator ET&() { return tensor; }
  template<typename... IndexTypes>
  inline const Scalar_& operator()(const std::array<Eigen::Index, NumIndices_> &Indices) const
  { return tensor.operator()(Indices); }
  inline Scalar_& operator()(const std::array<Eigen::Index, NumIndices_> &Indices)
  { return tensor.operator()(Indices); }
  template<typename... IndexTypes>
  inline const Scalar_& operator()(Eigen::Index firstDimension, IndexTypes... otherDimensions) const
  {
    // The number of indices used to access a tensor coefficient must be equal to the rank of the tensor.
    assert(sizeof...(otherDimensions) + 1 == NumIndices_ && "NamedTensor: dimensions != tensor rank");
    return tensor.operator()(std::array<Eigen::Index, NumIndices_>{{firstDimension, otherDimensions...}});
  }
  template<typename... IndexTypes>
  inline Scalar_& operator()(Eigen::Index firstDimension, IndexTypes... otherDimensions)
  {
    // The number of indices used to access a tensor coefficient must be equal to the rank of the tensor.
    assert(sizeof...(otherDimensions) + 1 == NumIndices_ && "NamedTensor: dimensions != tensor rank");
    return tensor.operator()(std::array<Eigen::Index, NumIndices_>{{firstDimension, otherDimensions...}});
  }
  
  // Construct a named tensor explicitly specifying size of each dimension
  template<typename... IndexTypes>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE NamedTensor(const std::array<std::string,NumIndices_> &IndexNames_, Eigen::Index firstDimension, IndexTypes... otherDimensions)
  : tensor(firstDimension, otherDimensions...), IndexNames{NumIndices}
  {
    // The number of dimensions used to construct a tensor must be equal to the rank of the tensor.
    assert(sizeof...(otherDimensions) + 1 == NumIndices_ && "NamedTensor: dimensions != tensor rank");
    for( int i = 0; i < NumIndices_; i++ )
      IndexNames[i] = IndexNames_[i];
      }
  
  // Default constructor (assumes tensor will be loaded from file)
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE NamedTensor() : IndexNames{NumIndices_} {}
  
  // Construct a named tensor without specifying size of each dimension (because it will be loaded from file)
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE NamedTensor(const std::array<std::string,NumIndices_> &IndexNames_)
  : IndexNames{NumIndices_}
  {
    for( int i = 0; i < NumIndices_; i++ )
      IndexNames[i] = IndexNames_[i];
      }
  
  bool ValidateIndexNames( std::size_t NumNames, const std::string * MatchNames ) const;
  
  // Read/Write in any format
  template<typename Reader> inline void read (Reader &r, const char * pszTag = nullptr);
  template<typename Writer> inline void write(Writer &w, const char * pszTag = nullptr) const;
  // Read/Write in default format, i.e. HDF5 if present, else binary
  inline void read (const char * filename, const char * pszTag = nullptr);
  inline void write(const char * filename, const char * pszTag = nullptr) const;
};

// Is this a named tensor
template<typename T, typename V = void> struct is_named_tensor : public std::false_type {};
template<typename Scalar_, int NumIndices_> struct is_named_tensor<NamedTensor<Scalar_, NumIndices_>> : public std::true_type {};
template<typename T> struct is_named_tensor<T, typename std::enable_if<std::is_base_of<NamedTensor<typename T::Scalar, T::NumIndices>, T>::value>::type> : public std::true_type {};

/******************************************************************************
 PerambTensor object
 ******************************************************************************/

using PerambTensor = NamedTensor<SpinVector, 6>;
static const std::array<std::string, 6> PerambIndexNames{"nT", "nVec", "LI", "nNoise", "nT_inv", "SI"};

/******************************************************************************
 Write NamedTensor
 ******************************************************************************/

template<typename Scalar_, int NumIndices_>
template<typename Writer>
void NamedTensor<Scalar_, NumIndices_>::write(Writer &w, const char * pszTag)const{
  if( pszTag == nullptr )
    pszTag = "NamedTensor";
  LOG(Message) << "Writing NamedTensor to tag " << pszTag << std::endl;
  write(w, pszTag, *this);
}

template<typename Scalar_, int NumIndices_>
void NamedTensor<Scalar_, NumIndices_>::write(const char * filename, const char * pszTag)const{
  std::string sFileName{filename};
  sFileName.append( MDistil::FileExtension );
  LOG(Message) << "Writing NamedTensor to file " << sFileName << std::endl;
  MDistil::Default_Writer w( sFileName );
  write( w, pszTag );
}

/******************************************************************************
 Validate named tensor index names
 ******************************************************************************/

template<typename Scalar_, int NumIndices_>
bool NamedTensor<Scalar_, NumIndices_>::ValidateIndexNames( std::size_t NumNames, const std::string * MatchNames ) const {
  bool bSame{ NumNames == NumIndices_ && IndexNames.size() == NumIndices_ };
  for( std::size_t i = 0; bSame && i < NumIndices_; i++ )
  {
    bSame = MatchNames[i].size() == IndexNames[i].size()
    && std::equal( MatchNames[i].begin(), MatchNames[i].end(), IndexNames[i].begin(),
                  [](const char & c1, const char & c2){ return c1 == c2 || std::toupper(c1) == std::toupper(c2); });
  }
  return bSame;
}

/******************************************************************************
 Read NamedTensor
 ******************************************************************************/

template<typename Scalar_, int NumIndices_>
template<typename Reader>
void NamedTensor<Scalar_, NumIndices_>::read(Reader &r, const char * pszTag) {
  if( pszTag == nullptr )
    pszTag = "NamedTensor";
  // Grab index names and dimensions
  std::vector<std::string> OldIndexNames{std::move(IndexNames)};
  typename ET::Dimensions OldDimensions{tensor.dimensions()};
  LOG(Message) << "Reading NamedTensor from tag  " << pszTag << std::endl;
  read(r, pszTag, *this);
  const typename ET::Dimensions & NewDimensions{tensor.dimensions()};
  for( int i = 0; i < NumIndices_; i++ )
    assert(OldDimensions[i] == 0 || OldDimensions[i] == NewDimensions[i] && "NamedTensor::load dimension size");
  assert( ValidateIndexNames( OldIndexNames.size(), &OldIndexNames[0] ) && "NamedTensor::load dimension name" );
}

template<typename Scalar_, int NumIndices_>
void NamedTensor<Scalar_, NumIndices_>::read(const char * filename, const char * pszTag) {
  std::string sFileName{filename};
  sFileName.append( MDistil::FileExtension );
  LOG(Message) << "Reading NamedTensor from file " << sFileName << std::endl;
  MDistil::Default_Reader r( sFileName );
  read( r, pszTag );
}

/******************************************************************************
 Make a lower dimensional grid in preparation for local slice operations
 ******************************************************************************/

inline GridCartesian * MakeLowerDimGrid( GridCartesian * gridHD )
{
  int nd{static_cast<int>(gridHD->_ndimension)};
  Coordinate latt_size   = gridHD->_gdimensions;
  latt_size[nd-1] = 1;
  Coordinate simd_layout = GridDefaultSimd(nd-1, vComplex::Nsimd());
  simd_layout.push_back( 1 );
  Coordinate mpi_layout  = gridHD->_processors;
  mpi_layout[nd-1] = 1;
  GridCartesian * gridLD = new GridCartesian(latt_size,simd_layout,mpi_layout,*gridHD);
  return gridLD;
}

/*************************************************************************************
 Rotate eigenvectors into our phase convention
 First component of first eigenvector is real and positive
 *************************************************************************************/

inline void RotateEigen(std::vector<LatticeColourVector> & evec)
{
  ColourVector cv0;
  auto grid = evec[0].Grid();
  Coordinate siteFirst(grid->Nd(),0);
  peekSite(cv0, evec[0], siteFirst);
  Grid::Complex cplx0 = cv0()()(0);
  if( cplx0.imag() == 0 )
    std::cout << GridLogMessage << "RotateEigen() : Site 0 : " << cplx0 << " => already meets phase convention" << std::endl;
  else {
    const Real cplx0_mag = Grid::abs(cplx0);
#ifdef GRID_NVCC
    const Grid::Complex phase = thrust::conj(cplx0 / cplx0_mag);
    const Real argphase = thrust::arg(phase);
#else
    const Grid::Complex phase = std::conj(cplx0 / cplx0_mag);
    const Real argphase = std::arg(phase);
#endif
    std::cout << GridLogMessage << "RotateEigen() : Site 0 : |" << cplx0 << "|=" << cplx0_mag << " => phase=" << (argphase / 3.14159265) << " pi" << std::endl;
    {
      // TODO: Only really needed on the master slice
      for( int k = 0 ; k < evec.size() ; k++ )
        evec[k] *= phase;
      if(grid->IsBoss()){
        for( int c = 0 ; c < Nc ; c++ )
          cv0()()(c) *= phase;
        cplx0.imag(0); // This assumes phase convention is real, positive (so I get rid of rounding error)
        //pokeSite(cv0, evec[0], siteFirst);
        pokeLocalSite(cv0, evec[0], siteFirst);
      }
    }
  }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif // Hadrons_MDistil_Distil_hpp_
