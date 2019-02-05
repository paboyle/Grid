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

/******************************************************************************
 Needed to make sure envCreate() (see Hadrons) works with specialisations
 with more than one parameter, eg obj<T1 COMMA T2>
 I imagine this exists already?
 ******************************************************************************/

#ifndef COMMA
#define COMMA ,
#endif

/******************************************************************************
 A consistent set of cross-platform methods for big endian <-> host byte ordering
 I imagine this exists already?
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
 Turns out I don't actually need this when running inside hadrons
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
 
 -Grad^2 (Peardon, 2009, pg 2, equation 3)
 Field      Type of field the operator will be applied to
 GaugeField Gauge field the operator will smear using
 
 TODO CANDIDATE for integration into laplacian operator
 should just require adding number of dimensions to act on to constructor,
 where the default=all dimensions, but we could specify 3 spatial dimensions
 
 *************************************************************************************/

template<typename Field, typename GaugeField=LatticeGaugeFieldD>
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
  
  LinOpPeardonNablaHerm(OperatorFunction<Field> & poly,LinearOperatorBase<Field>& linop) : _poly(poly), _Linop(linop) {
  }
  
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

typedef Grid::Hadrons::EigenPack<LatticeColourVector>       DistilEP;
typedef std::vector<std::vector<std::vector<SpinVector> > > DistilNoises;

/******************************************************************************
 Make a lower dimensional grid
 ******************************************************************************/

inline GridCartesian * MakeLowerDimGrid( GridCartesian * gridHD )
{
  //LOG(Message) << "MakeLowerDimGrid() begin" << std::endl;
  int nd{static_cast<int>(gridHD->_ndimension)};
  std::vector<int> latt_size   = gridHD->_gdimensions;
  latt_size[nd-1] = 1;

  std::vector<int> simd_layout = GridDefaultSimd(nd-1, vComplex::Nsimd());
  simd_layout.push_back( 1 );

  std::vector<int> mpi_layout  = gridHD->_processors;
  mpi_layout[nd-1] = 1;
  GridCartesian * gridLD = new GridCartesian(latt_size,simd_layout,mpi_layout,*gridHD);
  //LOG(Message) << "MakeLowerDimGrid() end" << std::endl;
  return gridLD;
}

/******************************************************************************
 Perambulator object
 This is an Eigen::Tensor of type Scalar_ and rank NumIndices_ (row-major order)
 They can be persisted to disk, with the on-disk format being big endian.
 Scalar_ objects are assumed to be composite objects of size Endian_Scalar_Size.
 (Disable big-endian by setting Endian_Scalar_Size=1)
 IndexNames contains one name for each index, and IndexNames are validated on load.
 (NB: Indices of dimension 1 are not saved, and not validated on load)
 ******************************************************************************/

template<typename Scalar_, int NumIndices_, uint16_t Endian_Scalar_Size = sizeof(Scalar_)>
class NamedTensor : public Eigen::Tensor<Scalar_, NumIndices_, Eigen::RowMajor>
{
public:
  typedef Eigen::Tensor<Scalar_, NumIndices_, Eigen::RowMajor> ET;
  std::array<std::string,NumIndices_> IndexNames;
public:
  template<typename... IndexTypes>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE NamedTensor(std::array<std::string,NumIndices_> &IndexNames_, Eigen::Index firstDimension, IndexTypes... otherDimensions)
  : IndexNames{IndexNames_}, ET(firstDimension, otherDimensions...)
  {
    // The number of dimensions used to construct a tensor must be equal to the rank of the tensor.
    assert(sizeof...(otherDimensions) + 1 == NumIndices_
           && "NamedTensor error: dimensions in constructor != tensor rank");
  }

  // Share data for timeslices we calculated with other nodes
  inline void SliceShare( GridCartesian * gridLowDim, GridCartesian * gridHighDim ) {
    Grid::SliceShare( gridLowDim, gridHighDim, this->data(), (int) (this->size() * sizeof(Scalar_)));
  }

  // load and save - not virtual - probably all changes
  inline void load(const std::string filename);
  inline void save(const std::string filename) const;
  inline void ReadBinary(const std::string filename);
  inline void WriteBinary(const std::string filename);
};

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
  const auto NumElements{this->size()};
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
  for( auto dim : this->dimensions() )
    if( dim == 1 )
      u16--;
  u16 = htobe16( u16 );
  w.write(reinterpret_cast<const char *>(&u16), sizeof(u16));
  // dimensions together with names
  int d = 0;
  for( auto dim : this->dimensions() ) {
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
  char * const pStart{reinterpret_cast<char *>(this->data())};
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
  u32 = htobe32(GridChecksum::crc32c(this->data(), TotalDataSize));
#else
  u32 = htobe32(GridChecksum::crc32(this->data(), TotalDataSize));
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
  const auto NumElements{this->size()};
  const std::streamsize TotalDataSize{static_cast<std::streamsize>(NumElements * Scalar_Size)};
  uint64_t u64;
  r.read(reinterpret_cast<char *>(&u64), sizeof(u64));
  assert( TotalDataSize == be64toh( u64 ) && "NamedTensor error: Size of the data in bytes" );
  // Size of a Scalar_
  uint32_t u32;
  r.read(reinterpret_cast<char *>(&u32), sizeof(u32));
  assert( Scalar_Size == be32toh( u32 ) && "NamedTensor error: sizeof(Scalar_)");
  // Endian_Scalar_Size
  uint16_t u16;
  r.read(reinterpret_cast<char *>(&u16), sizeof(u16));
  assert( Endian_Scalar_Size == be16toh( u16 ) && "NamedTensor error: Scalar_Unit_size");
  // number of dimensions which aren't 1
  r.read(reinterpret_cast<char *>(&u16), sizeof(u16));
  u16 = be16toh( u16 );
  for( auto dim : this->dimensions() )
    if( dim == 1 )
      u16++;
  assert( this->NumIndices == u16 && "NamedTensor error: number of dimensions which aren't 1" );
  // dimensions together with names
  int d = 0;
  for( auto dim : this->dimensions() ) {
    if( dim != 1 ) {
      // size of dimension
      r.read(reinterpret_cast<char *>(&u16), sizeof(u16));
      assert( dim == be16toh( u16 ) && "size of dimension" );
      // length of dimension name
      r.read(reinterpret_cast<char *>(&u16), sizeof(u16));
      size_t l = be16toh( u16 );
      assert( l == IndexNames[d].size() && "NamedTensor error: length of dimension name" );
      // dimension name
      std::string s( l, '?' );
      r.read(&s[0], l);
      assert( s == IndexNames[d] && "NamedTensor error: dimension name" );
    }
    d++;
  }
  // Actual data
  char * const pStart{reinterpret_cast<char *>(this->data())};
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
  u32 -= GridChecksum::crc32c(this->data(), TotalDataSize);
#else
  u32 -= GridChecksum::crc32(this->data(), TotalDataSize);
#endif
  assert( u32 == 0 && "NamedTensor error: Perambulator checksum invalid");
}

/******************************************************************************
 Save NamedTensor Hdf5 format
 ******************************************************************************/

template<typename Scalar_, int NumIndices_, uint16_t Endian_Scalar_Size>
void NamedTensor<Scalar_, NumIndices_, Endian_Scalar_Size>::save(const std::string filename) const {
  LOG(Message) << "Writing NamedTensor to \"" << filename << "\"" << std::endl;
#ifndef HAVE_HDF5
  LOG(Message) << "Error: I/O for NamedTensor requires HDF5" << std::endl;
#else
  Hdf5Writer w(filename);
  //w << this->NumIndices << this->dimensions() << this->IndexNames;
#endif
}

/******************************************************************************
 Load NamedTensor Hdf5 format
 ******************************************************************************/

template<typename Scalar_, int NumIndices_, uint16_t Endian_Scalar_Size>
void NamedTensor<Scalar_, NumIndices_, Endian_Scalar_Size>::load(const std::string filename) {
  LOG(Message) << "Reading NamedTensor from \"" << filename << "\"" << std::endl;
#ifndef HAVE_HDF5
  LOG(Message) << "Error: I/O for NamedTensor requires HDF5" << std::endl;
#else
  Hdf5Reader r(filename);
  typename ET::Dimensions d;
  std::array<std::string,NumIndices_> n;
  //r >> this->NumIndices >> d >> n;
  //this->IndexNames = n;
#endif
}

/******************************************************************************
 Perambulator object
 ******************************************************************************/

template<typename Scalar_, int NumIndices_, uint16_t Endian_Scalar_Size = sizeof(Scalar_)>
using Perambulator = NamedTensor<Scalar_, NumIndices_, Endian_Scalar_Size>;

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
