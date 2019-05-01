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

/******************************************************************************
 A consistent set of cross-platform methods for big endian <-> host byte ordering
 I imagine this exists already?
 This can be removed once the (deprecated) NamedTensor::ReadBinary & WriteBinary methods are deleted
 ******************************************************************************/

#if defined(__linux__)
#  include <endian.h>
#elif defined(__FreeBSD__) || defined(__NetBSD__)
#  include <sys/endian.h>
#elif defined(__OpenBSD__)
#  include <sys/types.h>
#  define be16toh(x) betoh16(x)
#  define be32toh(x) betoh32(x)
#  define be64toh(x) betoh64(x)
#elif defined(__APPLE__)
#include <libkern/OSByteOrder.h>
#define htobe16(x) OSSwapHostToBigInt16(x)
#define htole16(x) OSSwapHostToLittleInt16(x)
#define be16toh(x) OSSwapBigToHostInt16(x)
#define le16toh(x) OSSwapLittleToHostInt16(x)

#define htobe32(x) OSSwapHostToBigInt32(x)
#define htole32(x) OSSwapHostToLittleInt32(x)
#define be32toh(x) OSSwapBigToHostInt32(x)
#define le32toh(x) OSSwapLittleToHostInt32(x)

#define htobe64(x) OSSwapHostToBigInt64(x)
#define htole64(x) OSSwapHostToLittleInt64(x)
#define be64toh(x) OSSwapBigToHostInt64(x)
#define le64toh(x) OSSwapLittleToHostInt64(x)
#endif

/******************************************************************************
 This potentially belongs in CartesianCommunicator
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(Grid)
inline void SliceShare( GridBase * gridLowDim, GridBase * gridHighDim, void * Buffer, int BufferSize )
{
  // Work out which dimension is the spread-out dimension
  assert(gridLowDim);
  assert(gridHighDim);
  const int iNumDims{(const int)gridHighDim->_gdimensions.size()};
  assert(iNumDims == gridLowDim->_gdimensions.size());
  int dimSpreadOut = -1;
  std::vector<int> coor(iNumDims);
  for( int i = 0 ; i < iNumDims ; i++ ) {
    coor[i] = gridHighDim->_processor_coor[i];
    if( gridLowDim->_gdimensions[i] != gridHighDim->_gdimensions[i] ) {
      assert( dimSpreadOut == -1 );
      assert( gridLowDim->_processors[i] == 1 ); // easiest assumption to make for now
      dimSpreadOut = i;
    }
  }
  if( dimSpreadOut != -1 && gridHighDim->_processors[dimSpreadOut] != gridLowDim->_processors[dimSpreadOut] ) {
    // Make sure the same number of data elements exist on each slice
    const int NumSlices{gridHighDim->_processors[dimSpreadOut] / gridLowDim->_processors[dimSpreadOut]};
    assert(gridHighDim->_processors[dimSpreadOut] == gridLowDim->_processors[dimSpreadOut] * NumSlices);
    const int SliceSize{BufferSize/NumSlices};
    //CCC_DEBUG_DUMP(Buffer, NumSlices, SliceSize);
    assert(BufferSize == SliceSize * NumSlices);
//#ifndef USE_LOCAL_SLICES
//    assert(0); // Can't do this without MPI (should really test whether MPI is defined)
//#else
    const auto MyRank{gridHighDim->ThisRank()};
    std::vector<CommsRequest_t> reqs(0);
    int MySlice{coor[dimSpreadOut]};
    char * const _buffer{(char *)Buffer};
    char * const MyData{_buffer + MySlice * SliceSize};
    for(int i = 1; i < NumSlices ; i++ ){
      int SendSlice = ( MySlice + i ) % NumSlices;
      int RecvSlice = ( MySlice - i + NumSlices ) % NumSlices;
      char * const RecvData{_buffer + RecvSlice * SliceSize};
      coor[dimSpreadOut] = SendSlice;
      const auto SendRank{gridHighDim->RankFromProcessorCoor(coor)};
      coor[dimSpreadOut] = RecvSlice;
      const auto RecvRank{gridHighDim->RankFromProcessorCoor(coor)};
      std::cout << GridLogMessage << "Send slice " << MySlice << " (" << MyRank << ") to " << SendSlice << " (" << SendRank
      << "), receive slice from " << RecvSlice << " (" << RecvRank << ")" << std::endl;
      gridHighDim->SendToRecvFromBegin(reqs,MyData,SendRank,RecvData,RecvRank,SliceSize);
      //memcpy(RecvData,MyData,SliceSize); // Debug
    }
    gridHighDim->SendToRecvFromComplete(reqs);
    std::cout << GridLogMessage << "Slice data shared." << std::endl;
    //CCC_DEBUG_DUMP(Buffer, NumSlices, SliceSize);
//#endif
  }
}

/*************************************************************************************

 Not sure where the right home for this is? But presumably in Grid
 
 -Grad^2 (Peardon, 2009, pg 2, equation 3, https://arxiv.org/abs/0905.2160)
 Field      Type of field the operator will be applied to
 GaugeField Gauge field the operator will smear using
 
 *************************************************************************************/

template<typename Field, typename GaugeField=LatticeGaugeField>
class LinOpPeardonNabla : public LinearOperatorBase<Field>, public LinearFunction<Field> {
  typedef typename GaugeField::vector_type vCoeff_t;
protected: // I don't really mind if _gf is messed with ... so make this public?
  //GaugeField & _gf;
  int          nd; // number of spatial dimensions
  std::vector<Lattice<iColourMatrix<vCoeff_t> > > U;
public:
  // Construct this operator given a gauge field and the number of dimensions it should act on
  LinOpPeardonNabla( GaugeField& gf, int dimSpatial = Grid::QCD::Tdir ) : /*_gf(gf),*/ nd{dimSpatial} {
    assert(dimSpatial>=1);
    for( int mu = 0 ; mu < nd ; mu++ )
      U.push_back(PeekIndex<LorentzIndex>(gf,mu));
      }
  
  // Apply this operator to "in", return result in "out"
  void operator()(const Field& in, Field& out) {
    assert( nd <= in._grid->Nd() );
    conformable( in, out );
    out = ( ( Real ) ( 2 * nd ) ) * in;
    Field _tmp(in._grid);
    typedef typename GaugeField::vector_type vCoeff_t;
    //Lattice<iColourMatrix<vCoeff_t> > U(in._grid);
    for( int mu = 0 ; mu < nd ; mu++ ) {
      //U = PeekIndex<LorentzIndex>(_gf,mu);
      out -= U[mu] * Cshift( in, mu, 1);
      _tmp = adj( U[mu] ) * in;
      out -= Cshift(_tmp,mu,-1);
    }
  }
  
  void OpDiag (const Field &in, Field &out) { assert(0); };
  void OpDir  (const Field &in, Field &out,int dir,int disp) { assert(0); };
  void Op     (const Field &in, Field &out) { assert(0); };
  void AdjOp  (const Field &in, Field &out) { assert(0); };
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2) { assert(0); };
  void HermOp(const Field &in, Field &out) { operator()(in,out); };
};

template<typename Field>
class LinOpPeardonNablaHerm : public LinearFunction<Field> {
public:
  OperatorFunction<Field>   & _poly;
  LinearOperatorBase<Field> &_Linop;
  
  LinOpPeardonNablaHerm(OperatorFunction<Field> & poly,LinearOperatorBase<Field>& linop)
  : _poly{poly}, _Linop{linop} {}
  
  void operator()(const Field& in, Field& out) {
    _poly(_Linop,in,out);
  }
};

END_MODULE_NAMESPACE // Grid

/******************************************************************************
 Common elements for distillation
 ******************************************************************************/

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MDistil)

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

  int getTI( const Environment & env, bool bCalledFromSetup ) const {
    return ParameterDefault( TI, env.getDim(Tdir), bCalledFromSetup ); }
};

#define DISTIL_PARAMETERS_DEFINE( inSetup ) \
const int Nt{env().getDim(Tdir)}; \
const int nvec{par().nvec}; \
const int Ns{Grid::QCD::Ns}; \
const int nnoise{par().Distil.nnoise}; \
const int tsrc{par().Distil.tsrc}; \
const int TI{par().Distil.getTI(env(), inSetup)}; \
const int LI{Hadrons::MDistil::DistilParameters::ParameterDefault(par().Distil.LI, nvec, inSetup)}; \
const int SI{Hadrons::MDistil::DistilParameters::ParameterDefault(par().Distil.SI, Ns, inSetup)}; \
const bool full_tdil{ TI == Nt }; \
const bool exact_distillation{ full_tdil && LI == nvec }; \
const int Nt_inv{ full_tdil ? 1 : TI }

class BFieldIO: Serializable{
public:
  using BaryonTensorSet = Eigen::Tensor<Complex, 7>;
  GRID_SERIALIZABLE_CLASS_MEMBERS(BFieldIO, BaryonTensorSet, BField );
};

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
 Scalar_ objects are assumed to be composite objects of size Endian_Scalar_Size.
 (Disable big-endian by setting Endian_Scalar_Size=1).
 NB: Endian_Scalar_Size will disappear when ReadBinary & WriteBinary retired
 IndexNames contains one name for each index, and IndexNames are validated on load.
 WHAT TO SAVE / VALIDATE ON LOAD (Override to warn instead of assert on load)
 Ensemble string
 Configuration number
 Noise unique string
 Distillation parameters

 ******************************************************************************/

template<typename Scalar_, int NumIndices_, uint16_t Endian_Scalar_Size_ = sizeof(Scalar_)>
class NamedTensor : Serializable
{
public:
  using Scalar = Scalar_;
  static constexpr int NumIndices = NumIndices_;
  static constexpr uint16_t Endian_Scalar_Size = Endian_Scalar_Size_;
  using ET = Eigen::Tensor<Scalar_, NumIndices_, Eigen::RowMajor>;
  using Index = typename ET::Index;
  GRID_SERIALIZABLE_CLASS_MEMBERS(NamedTensor
                                  , ET, tensor
                                  , std::vector<std::string>, IndexNames
                                  );
public:
  // Named tensors are intended to be a superset of Eigen tensor
  inline operator ET&() const { return tensor; }
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
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE NamedTensor(std::array<std::string,NumIndices_> &IndexNames_)
  : IndexNames{NumIndices_}
  {
    for( int i = 0; i < NumIndices_; i++ )
      IndexNames[i] = IndexNames_[i];
  }
  
  // Share data for timeslices we calculated with other nodes
  inline void SliceShare( GridCartesian * gridLowDim, GridCartesian * gridHighDim ) {
    Grid::SliceShare( gridLowDim, gridHighDim, tensor.data(), (int) (tensor.size() * sizeof(Scalar_)));
  }

  bool ValidateIndexNames( int iNumNames, const std::string * MatchNames ) const;

  // Read/Write in any format
  template<typename Reader> inline void read (Reader &r, const char * pszTag = nullptr);
  template<typename Writer> inline void write(Writer &w, const char * pszTag = nullptr) const;
  // Read/Write in default format, i.e. HDF5 if present, else binary
  inline void read (const char * filename, const char * pszTag = nullptr);
  inline void write(const char * filename, const char * pszTag = nullptr) const;
  // Original I/O implementation. This will be removed when we're sure it's no longer needed
  EIGEN_DEPRECATED inline void ReadBinary (const std::string filename); // To be removed
  EIGEN_DEPRECATED inline void WriteBinary(const std::string filename); // To be removed
};

// Is this a named tensor
template<typename T, typename V = void> struct is_named_tensor : public std::false_type {};
template<typename Scalar_, int NumIndices_, uint16_t Endian_Scalar_Size_> struct is_named_tensor<NamedTensor<Scalar_, NumIndices_, Endian_Scalar_Size_>> : public std::true_type {};
template<typename T> struct is_named_tensor<T, typename std::enable_if<std::is_base_of<NamedTensor<typename T::Scalar, T::NumIndices, T::Endian_Scalar_Size_>, T>::value>::type> : public std::true_type {};

/******************************************************************************
 PerambTensor object
 Endian_Scalar_Size can be removed once (deprecated) NamedTensor::ReadBinary & WriteBinary methods deleted
 ******************************************************************************/

//template<typename Scalar_, int NumIndices_, uint16_t Endian_Scalar_Size = sizeof(Scalar_)>
using PerambTensor = NamedTensor<SpinVector, 6, sizeof(Real)>;
static const std::array<std::string, 6> PerambIndexNames{"nT", "nVec", "LI", "nNoise", "nT_inv", "SI"};

/******************************************************************************
 Save NamedTensor binary format (NB: On-disk format is Big Endian)
 Assumes the Scalar_ objects are contiguous (no padding)
 ******************************************************************************/

template<typename Scalar_, int NumIndices_, uint16_t Endian_Scalar_Size>
void NamedTensor<Scalar_, NumIndices_, Endian_Scalar_Size>::WriteBinary(const std::string filename) {
  LOG(Message) << "Writing NamedTensor to \"" << filename << "\"" << std::endl;
  std::ofstream w(filename, std::ios::binary);
  // Enforce assumption that the scalar is composed of fundamental elements of size Endian_Scalar_Size
  assert((Endian_Scalar_Size == 1 || Endian_Scalar_Size == 2 || Endian_Scalar_Size == 4 || Endian_Scalar_Size == 8 )
         && "NamedTensor error: Endian_Scalar_Size should be 1, 2, 4 or 8");
  assert((sizeof(Scalar_) % Endian_Scalar_Size) == 0 && "NamedTensor error: Scalar_ is not composed of Endian_Scalar_Size" );
  // Size of the data (in bytes)
  const uint32_t Scalar_Size{sizeof(Scalar_)};
  const auto NumElements{tensor.size()};
  const std::streamsize TotalDataSize{static_cast<std::streamsize>(NumElements * Scalar_Size)};
  uint64_t u64 = htobe64(static_cast<uint64_t>(TotalDataSize));
  w.write(reinterpret_cast<const char *>(&u64), sizeof(u64));
  // Size of a Scalar_
  uint32_t u32{htobe32(Scalar_Size)};
  w.write(reinterpret_cast<const char *>(&u32), sizeof(u32));
  // Endian_Scalar_Size
  uint16_t u16{htobe16(Endian_Scalar_Size)};
  w.write(reinterpret_cast<const char *>(&u16), sizeof(u16));
  // number of dimensions which aren't 1
  u16 = static_cast<uint16_t>(this->NumIndices);
  for( auto dim : tensor.dimensions() )
    if( dim == 1 )
      u16--;
  u16 = htobe16( u16 );
  w.write(reinterpret_cast<const char *>(&u16), sizeof(u16));
  // dimensions together with names
  int d = 0;
  for( auto dim : tensor.dimensions() ) {
    if( dim != 1 ) {
      // size of this dimension
      u16 = htobe16( static_cast<uint16_t>( dim ) );
      w.write(reinterpret_cast<const char *>(&u16), sizeof(u16));
      // length of this dimension name
      u16 = htobe16( static_cast<uint16_t>( IndexNames[d].size() ) );
      w.write(reinterpret_cast<const char *>(&u16), sizeof(u16));
      // dimension name
      w.write(IndexNames[d].c_str(), IndexNames[d].size());
    }
    d++;
  }
  // Actual data
  char * const pStart{reinterpret_cast<char *>(tensor.data())};
  // Swap to network byte order in place (alternative is to copy memory - still slow)
  void * const pEnd{pStart + TotalDataSize};
  if(Endian_Scalar_Size == 8)
    for(uint64_t * p = reinterpret_cast<uint64_t *>(pStart) ; p < pEnd ; p++ )
      * p = htobe64( * p );
  else if(Endian_Scalar_Size == 4)
    for(uint32_t * p = reinterpret_cast<uint32_t *>(pStart) ; p < pEnd ; p++ )
      * p = htobe32( * p );
  else if(Endian_Scalar_Size == 2)
    for(uint16_t * p = reinterpret_cast<uint16_t *>(pStart) ; p < pEnd ; p++ )
      * p = htobe16( * p );
  w.write(pStart, TotalDataSize);
  // Swap back from network byte order
  if(Endian_Scalar_Size == 8)
    for(uint64_t * p = reinterpret_cast<uint64_t *>(pStart) ; p < pEnd ; p++ )
      * p = be64toh( * p );
  else if(Endian_Scalar_Size == 4)
    for(uint32_t * p = reinterpret_cast<uint32_t *>(pStart) ; p < pEnd ; p++ )
      * p = be32toh( * p );
  else if(Endian_Scalar_Size == 2)
    for(uint16_t * p = reinterpret_cast<uint16_t *>(pStart) ; p < pEnd ; p++ )
      * p = be16toh( * p );
  // checksum
#ifdef USE_IPP
  u32 = htobe32(GridChecksum::crc32c(tensor.data(), TotalDataSize));
#else
  u32 = htobe32(GridChecksum::crc32(tensor.data(), TotalDataSize));
#endif
  w.write(reinterpret_cast<const char *>(&u32), sizeof(u32));
}

/******************************************************************************
 Load NamedTensor binary format (NB: On-disk format is Big Endian)
 Assumes the Scalar_ objects are contiguous (no padding)
 ******************************************************************************/

template<typename Scalar_, int NumIndices_, uint16_t Endian_Scalar_Size>
void NamedTensor<Scalar_, NumIndices_, Endian_Scalar_Size>::ReadBinary(const std::string filename) {
  LOG(Message) << "Reading NamedTensor from \"" << filename << "\"" << std::endl;
  std::ifstream r(filename, std::ios::binary);
  // Enforce assumption that the scalar is composed of fundamental elements of size Endian_Scalar_Size
  assert((Endian_Scalar_Size == 1 || Endian_Scalar_Size == 2 || Endian_Scalar_Size == 4 || Endian_Scalar_Size == 8 )
         && "NamedTensor error: Endian_Scalar_Size should be 1, 2, 4 or 8");
  assert((sizeof(Scalar_) % Endian_Scalar_Size) == 0 && "NamedTensor error: Scalar_ is not composed of Endian_Scalar_Size" );
  // Size of the data in bytes
  const uint32_t Scalar_Size{sizeof(Scalar_)};
  Index NumElements{tensor.size()};
  std::streamsize TotalDataSize{static_cast<std::streamsize>(NumElements * Scalar_Size)};
  uint64_t u64;
  r.read(reinterpret_cast<char *>(&u64), sizeof(u64));
  assert( TotalDataSize == 0 || TotalDataSize == be64toh( u64 ) && "NamedTensor error: Size of the data in bytes" );
  // Size of a Scalar_
  uint32_t u32;
  r.read(reinterpret_cast<char *>(&u32), sizeof(u32));
  assert( Scalar_Size == be32toh( u32 ) && "NamedTensor error: sizeof(Scalar_)");
  // Endian_Scalar_Size
  uint16_t u16;
  r.read(reinterpret_cast<char *>(&u16), sizeof(u16));
  assert( Endian_Scalar_Size == be16toh( u16 ) && "NamedTensor error: Scalar_Unit_size");
  // number of dimensions which aren't 1
  uint16_t NumFileDimensions;
  r.read(reinterpret_cast<char *>(&NumFileDimensions), sizeof(NumFileDimensions));
  NumFileDimensions = be16toh( NumFileDimensions );
  /*for( auto dim : tensor.dimensions() )
    if( dim == 1 )
      u16++;*/
  assert( ( TotalDataSize == 0 && this->NumIndices >= NumFileDimensions || this->NumIndices == NumFileDimensions )
         && "NamedTensor error: number of dimensions which aren't 1" );
  if( TotalDataSize == 0 ) {
    // Read each dimension, using names to skip past dimensions == 1
    std::array<Index,NumIndices_> NewDimensions;
    for( Index &i : NewDimensions ) i = 1;
    int d = 0;
    for( int FileDimension = 0; FileDimension < NumFileDimensions; FileDimension++ ) {
      // read dimension
      uint16_t thisDim;
      r.read(reinterpret_cast<char *>(&thisDim), sizeof(thisDim));
      // read dimension name
      r.read(reinterpret_cast<char *>(&u16), sizeof(u16));
      size_t l = be16toh( u16 );
      std::string s( l, '?' );
      r.read(&s[0], l);
      // skip forward to matching name
      while( IndexNames[d].size() > 0 && s != IndexNames[d] )
        assert(++d < NumIndices && "NamedTensor error: dimension name" );
      if( IndexNames[d].size() == 0 )
        IndexNames[d] = s;
      NewDimensions[d++] = be16toh( thisDim );
    }
    tensor.resize(NewDimensions);
    NumElements = 1;
    for( Index i : NewDimensions ) NumElements *= i;
    TotalDataSize = NumElements * Scalar_Size;
  } else {
    // dimensions together with names
    const auto & TensorDims{tensor.dimensions()};
    for( int d = 0; d < NumIndices_; d++ ) {
      // size of dimension
      r.read(reinterpret_cast<char *>(&u16), sizeof(u16));
      u16 = be16toh( u16 );
      assert( TensorDims[d] == u16 && "size of dimension" );
      // length of dimension name
      r.read(reinterpret_cast<char *>(&u16), sizeof(u16));
      size_t l = be16toh( u16 );
      assert( l == IndexNames[d].size() && "NamedTensor error: length of dimension name" );
      // dimension name
      std::string s( l, '?' );
      r.read(&s[0], l);
      assert( s == IndexNames[d] && "NamedTensor error: dimension name" );
    }
  }
  // Actual data
  char * const pStart{reinterpret_cast<char *>(tensor.data())};
  void * const pEnd{pStart + TotalDataSize};
  r.read(pStart,TotalDataSize);
  // Swap back from network byte order
  if(Endian_Scalar_Size == 8)
    for(uint64_t * p = reinterpret_cast<uint64_t *>(pStart) ; p < pEnd ; p++ )
      * p = be64toh( * p );
  else if(Endian_Scalar_Size == 4)
    for(uint32_t * p = reinterpret_cast<uint32_t *>(pStart) ; p < pEnd ; p++ )
      * p = be32toh( * p );
  else if(Endian_Scalar_Size == 2)
    for(uint16_t * p = reinterpret_cast<uint16_t *>(pStart) ; p < pEnd ; p++ )
      * p = be16toh( * p );
  // checksum
  r.read(reinterpret_cast<char *>(&u32), sizeof(u32));
  u32 = be32toh( u32 );
#ifdef USE_IPP
  u32 -= GridChecksum::crc32c(tensor.data(), TotalDataSize);
#else
  u32 -= GridChecksum::crc32(tensor.data(), TotalDataSize);
#endif
  assert( u32 == 0 && "NamedTensor error: PerambTensor checksum invalid");
}

/******************************************************************************
 Write NamedTensor
 ******************************************************************************/

template<typename Scalar_, int NumIndices_, uint16_t Endian_Scalar_Size>
template<typename Writer>
void NamedTensor<Scalar_, NumIndices_, Endian_Scalar_Size>::write(Writer &w, const char * pszTag)const{
  if( pszTag == nullptr )
    pszTag = "NamedTensor";
  LOG(Message) << "Writing NamedTensor to tag  " << pszTag << std::endl;
  write(w, pszTag, *this);
}

template<typename Scalar_, int NumIndices_, uint16_t Endian_Scalar_Size>
void NamedTensor<Scalar_, NumIndices_, Endian_Scalar_Size>::write(const char * filename, const char * pszTag)const{
  std::string sFileName{filename};
  sFileName.append( MDistil::FileExtension );
  LOG(Message) << "Writing NamedTensor to file " << sFileName << std::endl;
  MDistil::Default_Writer w( sFileName );
  write( w, pszTag );
}

/******************************************************************************
 Validate named tensor index names
 ******************************************************************************/

template<typename Scalar_, int NumIndices_, uint16_t Endian_Scalar_Size>
bool NamedTensor<Scalar_, NumIndices_, Endian_Scalar_Size>::ValidateIndexNames( int iNumNames, const std::string * MatchNames ) const {
  bool bSame{ iNumNames == NumIndices_ && IndexNames.size() == NumIndices_ };
  for( int i = 0; bSame && i < NumIndices_; i++ ) {
    // Case insensitive name compare
    const int iMatchLen{ static_cast<int>( MatchNames[i].size() ) };
    const int iNewLen  { static_cast<int>( IndexNames[i].size() ) };
    bSame = ( iMatchLen == iNewLen );
    for( int j = 0; bSame && j < iMatchLen; j++ ) {
      wchar_t c1 = MatchNames[i][j];
      if( c1 >= 'a' && c1 <= 'z' )
        c1 -= 'a' - 'A';
      wchar_t c2 = IndexNames[i][j];
      if( c2 >= 'a' && c1 <= 'z' )
        c2 -= 'a' - 'A';
      bSame = ( c1 == c2 );
    }
  }
  return bSame;
}

/******************************************************************************
 Read NamedTensor
 ******************************************************************************/

template<typename Scalar_, int NumIndices_, uint16_t Endian_Scalar_Size>
template<typename Reader>
void NamedTensor<Scalar_, NumIndices_, Endian_Scalar_Size>::read(Reader &r, const char * pszTag) {
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

template<typename Scalar_, int NumIndices_, uint16_t Endian_Scalar_Size>
void NamedTensor<Scalar_, NumIndices_, Endian_Scalar_Size>::read(const char * filename, const char * pszTag) {
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
  std::vector<int> latt_size   = gridHD->_gdimensions;
  latt_size[nd-1] = 1;
  std::vector<int> simd_layout = GridDefaultSimd(nd-1, vComplex::Nsimd());
  simd_layout.push_back( 1 );
  std::vector<int> mpi_layout  = gridHD->_processors;
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
  auto grid = evec[0]._grid;
  std::vector<int> siteFirst(grid->Nd(),0);
  peekSite(cv0, evec[0], siteFirst);
  auto & cplx0 = cv0()()(0);
  if( std::imag(cplx0) == 0 )
    std::cout << GridLogMessage << "RotateEigen() : Site 0 : " << cplx0 << " => already meets phase convention" << std::endl;
  else {
    const auto cplx0_mag{std::abs(cplx0)};
    const auto phase{std::conj(cplx0 / cplx0_mag)};
    std::cout << GridLogMessage << "RotateEigen() : Site 0 : |" << cplx0 << "|=" << cplx0_mag << " => phase=" << (std::arg(phase) / 3.14159265) << " pi" << std::endl;
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
