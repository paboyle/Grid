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
END_MODULE_NAMESPACE // Grid

/******************************************************************************
 Common elements for distillation
 ******************************************************************************/

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 Make a lower dimensional grid
 ******************************************************************************/

inline GridCartesian * MakeLowerDimGrid( GridCartesian * gridHD )
{
  int nd{static_cast<int>(gridHD->_ndimension)};
  std::vector<int> latt_size   = gridHD->_fdimensions;
  latt_size[nd-1] = 1;

  std::vector<int> simd_layout = GridDefaultSimd(nd-1, vComplex::Nsimd());
  simd_layout.push_back( 1 );

  std::vector<int> mpi_layout  = gridHD->_processors;
  mpi_layout[nd-1] = 1;
  return new GridCartesian(latt_size,simd_layout,mpi_layout,*gridHD);
}

/******************************************************************************
 Perambulator object
 ******************************************************************************/

template<typename LatticeObj>
class Perambulator : Serializable{
  // TODO: The next line makes friends across all combinations
  //     (not much of a problem given all public anyway ...)
  //      FYI, the bug here was that I forgot that the friend is templated
  template<typename T> friend std::ostream & operator<<(std::ostream &os, const Perambulator<T>& p);
protected:
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS( Perambulator,
                                  std::string,             ID,          // Allows owner to specialise
                                  std::string,             Provenance,  // For info only
                                  std::vector<int>,        dimensions,
                                  std::vector<LatticeObj>, perambulator,
                                  // Following items are redundant, but useful
                                  int,     nd,          // Number of dimensions
                                  size_t,  NumElements); // Number of elements
protected:
  // Constructor common code
  inline void ConstructCommon(const int * Dimensions) {
    assert(nd > 0);
    dimensions.resize(nd);
    NumElements = 1;
    for(int i = 0 ; i < nd ; i++) {
      assert(Dimensions[i] > 0);
      NumElements *= (size_t) Dimensions[i];
      dimensions[i] = Dimensions[i];
    }
    //const LatticeObj perambulatorDefault;
    perambulator.resize(NumElements);//,perambulatorDefault);
  }
public:
  // Constructor with dimensions passed as std::vector<int>
  inline Perambulator(const std::vector<int> & Dimensions)
  : nd {(int) Dimensions.size()} {
    ConstructCommon( &Dimensions[0] ); }
  
  // Constructor with dimensions passed as std::vector<int>
  inline Perambulator(const std::vector<int> & Dimensions, const std::string sID)
  : nd {(int) Dimensions.size()}, ID(sID) {
    ConstructCommon( &Dimensions[0] ); }
  
  // Constructor with dimensions passed as std::vector<int>
  inline Perambulator(const std::vector<int> & Dimensions, const std::string sID, const std::string sProvenance)
  : nd {(int) Dimensions.size()}, ID(sID), Provenance(sProvenance) {
    ConstructCommon( &Dimensions[0] ); }
  
  // Constructor with dimensions passed as individual parameters
  // FYI: The caller is free to ignore the names and use the indices however they see fit
  inline Perambulator(int NumNoise, int NumEvec=1, int NumTime=1, int NumSpin=1, int I_k=1, int I_t=1, int I_s=1) {
    int Dimensions[]={NumNoise,NumEvec,NumTime,NumSpin,I_k,I_t,I_s};
    nd = sizeof(Dimensions)/sizeof(Dimensions[0]);
    while( nd > 1 && Dimensions[nd-1] == 1 )
      nd--;
    ConstructCommon( Dimensions );
  }

  inline Perambulator(const std::string sID, int NumNoise, int NumEvec=1, int NumTime=1, int NumSpin=1, int I_k=1, int I_t=1, int I_s=1) : ID{sID} {
    int Dimensions[]={NumNoise,NumEvec,NumTime,NumSpin,I_k,I_t,I_s};
    nd = sizeof(Dimensions)/sizeof(Dimensions[0]);
    while( nd > 1 && Dimensions[nd-1] == 1 )
      nd--;
    ConstructCommon( Dimensions );
  }
  
  inline Perambulator(const std::string sID, const std::string sProvenance, int NumNoise, int NumEvec=1, int NumTime=1, int NumSpin=1, int I_k=1, int I_t=1, int I_s=1) : ID{sID}, Provenance{sProvenance} {
    int Dimensions[]={NumNoise,NumEvec,NumTime,NumSpin,I_k,I_t,I_s};
    nd = sizeof(Dimensions)/sizeof(Dimensions[0]);
    while( nd > 1 && Dimensions[nd-1] == 1 )
      nd--;
    ConstructCommon( Dimensions );
  }

  inline LatticeObj & operator()(size_t count, const int * Coord) {
    assert( count == nd );
    assert( Coord );
    size_t idx = 0;
    // C memory order (???)
    for( int d = 0 ; d < nd ; d++ ) {
      assert( Coord[d] < dimensions[d] );
      idx *= (size_t) dimensions[d];
      idx += (size_t) Coord[d];
    }
    return perambulator[idx];
  }
  
  inline LatticeObj & operator()(const std::vector<int> Coord) {
    return operator()(Coord.size(), &Coord[0]);
  }
  
  inline LatticeObj & operator()(int idxNoise, int idxEvec=0, int idxTime=0, int idxSpin=0, int I_k=0, int I_t=0, int I_s=0) {
    int MyIndex[]={idxNoise,idxEvec,idxTime,idxSpin,I_k,I_t,I_s};
    int i = sizeof(MyIndex)/sizeof(MyIndex[0]);
    assert( i >= nd );
    while( i > nd )
      assert(MyIndex[--i] == 0);
    return operator()(i, MyIndex);
  }

  // Share data for timeslices we calculated with other nodes
  inline void SliceShare( GridCartesian * gridLowDim, GridCartesian * gridHighDim ) {
    Grid::SliceShare( gridLowDim, gridHighDim, &perambulator[0],
                                               (int) perambulator.size() * sizeof(perambulator[0]) );
  }
  
  /*************************************************************************************
   
   Write/Read perambulator to/from disk
   
   Temporary version - keep the code running until such time as correct format written
   
   *************************************************************************************/
  
  inline void WriteTemporary(const std::string filename) const
  {
    std::cout << GridLogMessage << "Writing perambulator ID \"" << ID << "\" to " << filename << std::endl;
    BinaryWriter myPhDThesis( filename + ".tmp" );
    write( myPhDThesis, "Perambulator", *this );
  }
  
  inline bool ReadTemporary(const std::string filename)
  {
    std::string _filename{filename};
    _filename.append( ".tmp" );
    bool bReturnValue = false;
    std::fstream f;
    f.open(_filename,std::ios_base::in);
    if( !f.is_open() )
      std::cout << GridLogMessage << "Cached perambulator file " << _filename << " does not exist" << std::endl;
    else {
      f.close();
      bReturnValue = true;
      auto MyID{ID};
      std::cout << GridLogMessage << "Reading perambulator ID \"" << ID << "\" from " << _filename << std::endl;
      BinaryReader reader( _filename );
      read( reader, "Perambulator", *this );
      std::cout << GridLogMessage << "Perambulator ID read from " << _filename << " was \"" << ID << "\"" << std::endl;
      assert(MyID == ID);
    }
    return bReturnValue;
  }
  
  /*************************************************************************************
   
   Write perambulator to disk
   TODO 1. Ensure precision on disk can be specified independently of in memory
   2. Validate format with Peter, Antonin et al
   3. This object "works" for small lattii (lattices),
   BUT, the final block is written as XML, and while write is fast, read is painfully slow
   i.e. on a 24^3 x 64 lattice, I abandoned the perambulator read after 1 hour on Tesseract
   
   *************************************************************************************/
  
  inline void WritePoorly(LatticeGaugeField &field, const std::string filename)
  {
    std::cout << GridLogMessage << "Writing perambulator ID \"" << ID << "\" to " << filename << std::endl;
    assert(nd>=2); // Really should be a little bigger
    GridBase * gridHighDim = field._grid;
    //ScidacWriterPerambulator binWriter(gridHighDim->IsBoss());
    ScidacWriter binWriter(gridHighDim->IsBoss());
    //makeFileDir(filename, gridHighDim); // Assume this makes directory ... but why pass it the grid?
    binWriter.open(filename);
    // Write the header
    {
      XmlWriter xmlWriter("", "perambulatorPar");
      xmlWriter.pushXmlString("<ID>" + ID + "</ID>");
      xmlWriter.pushXmlString("<Provenance>" + Provenance + "</Provenance>");
      // TODO add all the perambulator parameters here
      binWriter.writeLimeObject(1, 1, xmlWriter, "parameters", SCIDAC_FILE_XML);
      std::cout << GridLogMessage << "Perambulator header written" << std::endl;
    }
    // Now write the local portion of the Perambulator
    {
      //binWriter.writeScidacPerambulatorRecord(field, *this);
      //std::cout << GridLogMessage << "Perambulator body written" << std::endl;
    }
    {
      ////////////////////////////////////////
      // fill the Grid header
      ////////////////////////////////////////
      FieldMetaData header;
      scidacRecord  _scidacRecord;
      scidacFile    _scidacFile;
      ScidacMetaData(field,header,_scidacRecord,_scidacFile);
      
      //////////////////////////////////////////////
      // Fill the Lime file record by record
      //////////////////////////////////////////////
      constexpr auto precision = std::numeric_limits<Real>::digits10;
      binWriter.writeLimeObject(1,0,header ,std::string("FieldMetaData"),std::string(GRID_FORMAT)); // Open message
      binWriter.writeLimeObject(0, 0, _scidacRecord, _scidacRecord.SerialisableClassName(),
                                std::string(SCIDAC_PRIVATE_RECORD_XML));
      binWriter.writeLimeObject(0,1,*this,this->SerialisableClassName(),std::string("Perambulator"),precision);
    }
  }
  
  /*************************************************************************************
   
   Read perambulator from disk
   TODO 1. Ensure precision on disk can be specified independently of in memory
   2. Validate format with Peter
   3. Abandoning for now because of synchronisation during write.
   Object small enough to send to root and write from there
   
   *************************************************************************************/
  
  struct PerambHeader{
    std::string ID, Provenance;
  };
  
  inline bool ReadPoorly(LatticeGaugeField &field, const std::string filename)
  {
    assert(nd>=2); // Really should be a little bigger
    bool bReturnValue = false;
    std::fstream f;
    f.open(filename,std::ios_base::in);
    if( !f.is_open() )
      std::cout << GridLogMessage << "Cached perambulator file " << filename << " does not exist" << std::endl;
    else {
      f.close();
      ScidacReader binReader;
      binReader.open(filename);
      PerambHeader header;
      // Read the header
      {
        std::string recordXml;
        std::cout << GridLogMessage << "Reading perambulator header from " << filename << std::endl;
        binReader.readLimeObject(recordXml, SCIDAC_FILE_XML);
        XmlReader xmlReader(recordXml, true, "perambulatorPar");
        xmlReader.push("perambulatorPar");
        //xmlReader.readCurrentSubtree(header.ID);
        xmlReader.readDefault("ID",header.ID);
        std::cout << GridLogMessage << "Perambulator ID=" << header.ID << std::endl;
        //xmlReader.nextElement();
        //xmlReader.readCurrentSubtree(header.Provenance);
        xmlReader.readDefault("Provenance",header.Provenance);
        std::cout << GridLogMessage << "Perambulator Provenance=" << header.Provenance << std::endl;
        assert( header.ID == ID );
        bReturnValue = true;
      }
      // Now read the Perambulator
      {
        ////////////////////////////////////////
        // fill the Grid header
        ////////////////////////////////////////
        FieldMetaData header;
        scidacRecord  _scidacRecord;
        
        //////////////////////////////////////////////
        // Fill the Lime file record by record
        //////////////////////////////////////////////
        binReader.readLimeObject(header ,std::string("FieldMetaData"),std::string(GRID_FORMAT)); // Open message
        binReader.readLimeObject( _scidacRecord, _scidacRecord.SerialisableClassName(),
                                 std::string(SCIDAC_PRIVATE_RECORD_XML));
        binReader.readLimeObject(*this,this->SerialisableClassName(),std::string("Perambulator"));
      }
      // TODO Add validation that the field matches what we read in
    }
    return bReturnValue;
  }
};

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_Distil_hpp_
