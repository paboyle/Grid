/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/parallelIO/IldgIO.h

Copyright (C) 2015

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef GRID_ILDG_IO_H
#define GRID_ILDG_IO_H

#ifdef HAVE_LIME
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>

#include <pwd.h>
#include <sys/utsname.h>
#include <unistd.h>

//Lime is a must have for this functionality
extern "C" {  // for linkage
#include "lime.h"
}

namespace Grid {
namespace QCD {

 template<class word> inline std::string ScidacWordMnemonic(void){ return std::string("unknown"); }
 template<> inline std::string ScidacWordMnemonic<double>  (void){ return std::string("D"); }
 template<> inline std::string ScidacWordMnemonic<float>   (void){ return std::string("F"); }
 template<> inline std::string ScidacWordMnemonic< int32_t>(void){ return std::string("I32_t"); }
 template<> inline std::string ScidacWordMnemonic<uint32_t>(void){ return std::string("U32_t"); }
 template<> inline std::string ScidacWordMnemonic< int64_t>(void){ return std::string("I64_t"); }
 template<> inline std::string ScidacWordMnemonic<uint64_t>(void){ return std::string("U64_t"); }

 template<class vobj> std::string ScidacRecordTypeString(int &colors, int &spins, int & typesize,int &datacount) { 

   typedef typename getPrecision<vobj>::real_scalar_type stype;

   int _ColourN       = indexRank<ColourIndex,vobj>();
   int _ColourScalar  =  isScalar<ColourIndex,vobj>();
   int _ColourVector  =  isVector<ColourIndex,vobj>();
   int _ColourMatrix  =  isMatrix<ColourIndex,vobj>();

   int _SpinN       = indexRank<SpinIndex,vobj>();
   int _SpinScalar  =  isScalar<SpinIndex,vobj>();
   int _SpinVector  =  isVector<SpinIndex,vobj>();
   int _SpinMatrix  =  isMatrix<SpinIndex,vobj>();

   int _LorentzN       = indexRank<LorentzIndex,vobj>();
   int _LorentzScalar  =  isScalar<LorentzIndex,vobj>();
   int _LorentzVector  =  isVector<LorentzIndex,vobj>();
   int _LorentzMatrix  =  isMatrix<LorentzIndex,vobj>();

   std::stringstream stream;

   stream << "GRID_";
   stream << ScidacWordMnemonic<stype>();

   //   std::cout << " Lorentz N/S/V/M : " << _LorentzN<<" "<<_LorentzScalar<<"/"<<_LorentzVector<<"/"<<_LorentzMatrix<<std::endl;
   //   std::cout << " Spin    N/S/V/M : " << _SpinN   <<" "<<_SpinScalar   <<"/"<<_SpinVector   <<"/"<<_SpinMatrix<<std::endl;
   //   std::cout << " Colour  N/S/V/M : " << _ColourN <<" "<<_ColourScalar <<"/"<<_ColourVector <<"/"<<_ColourMatrix<<std::endl;

   if ( _LorentzVector )   stream << "_LorentzVector"<<_LorentzN;
   if ( _LorentzMatrix )   stream << "_LorentzMatrix"<<_LorentzN;

   if ( _SpinVector )   stream << "_SpinVector"<<_SpinN;
   if ( _SpinMatrix )   stream << "_SpinMatrix"<<_SpinN;

   if ( _ColourVector )   stream << "_ColourVector"<<_ColourN;
   if ( _ColourMatrix )   stream << "_ColourMatrix"<<_ColourN;

   if ( _ColourScalar && _LorentzScalar && _SpinScalar )   stream << "_Complex";


   typesize = sizeof(typename vobj::scalar_type);

   if ( _ColourMatrix ) typesize*= _ColourN*_ColourN;
   else                 typesize*= _ColourN;

   if ( _SpinMatrix )   typesize*= _SpinN*_SpinN;
   else                 typesize*= _SpinN;

   colors    = _ColourN;
   spins     = _SpinN;
   datacount = _LorentzN;

   return stream.str();
 }
 
 template<class vobj> std::string ScidacRecordTypeString(Lattice<vobj> & lat,int &colors, int &spins, int & typesize,int &datacount) { 
   return ScidacRecordTypeString<vobj>(colors,spins,typesize,datacount);
 };

 template<class vobj> void ScidacMetaData(Lattice<vobj> & field,
					  FieldMetaData &header,
					  scidacRecord & _scidacRecord,
					  scidacFile   & _scidacFile) 
 {
   typedef typename getPrecision<vobj>::real_scalar_type stype;

   /////////////////////////////////////
   // Pull Grid's metadata
   /////////////////////////////////////
   PrepareMetaData(field,header);

   /////////////////////////////////////
   // Scidac Private File structure
   /////////////////////////////////////
   _scidacFile              = scidacFile(field._grid);

   /////////////////////////////////////
   // Scidac Private Record structure
   /////////////////////////////////////
   scidacRecord sr;
   sr.datatype   = ScidacRecordTypeString(field,sr.colors,sr.spins,sr.typesize,sr.datacount);
   sr.date       = header.creation_date;
   sr.precision  = ScidacWordMnemonic<stype>();
   sr.recordtype = GRID_IO_FIELD;

   _scidacRecord = sr;

   std::cout << GridLogMessage << "Build SciDAC datatype " <<sr.datatype<<std::endl;
 }
 
 ///////////////////////////////////////////////////////
 // Scidac checksum
 ///////////////////////////////////////////////////////
 static int scidacChecksumVerify(scidacChecksum &scidacChecksum_,uint32_t scidac_csuma,uint32_t scidac_csumb)
 {
   uint32_t scidac_checksuma = stoull(scidacChecksum_.suma,0,16);
   uint32_t scidac_checksumb = stoull(scidacChecksum_.sumb,0,16);
   if ( scidac_csuma !=scidac_checksuma) return 0;
   if ( scidac_csumb !=scidac_checksumb) return 0;
    return 1;
 }

////////////////////////////////////////////////////////////////////////////////////
// Lime, ILDG and Scidac I/O classes
////////////////////////////////////////////////////////////////////////////////////
class LimeIO : public BinaryIO {
 public:

   ///////////////////////////////////////////////////
   // FIXME: format for RNG? Now just binary out instead
   // FIXME: Make interface able to write multiple records
   // FIXME: Split into LimeReader and LimeWriter
   ///////////////////////////////////////////////////
   /*
   FILE * File;
   LimeWriter LimeW;
   LimeReader LimeR;
   template<class serialisable_object>
   int readObject(serialisable_object &object,std::string object_name,std::string record_name)

  int createLimeRecordHeader(std::string message, int MB, int ME, size_t PayloadSize);
  template<class serialisable_object>
  int writeObject(int MB,int ME,serialisable_object &object,std::string object_name,std::string record_name)
  template<class vobj>
  int writeLimeLatticeBinaryObject(Lattice<vobj> &field,std::string filename,std::string record_name)
   */
  ///////////////////////////////////////////////////////
  // Lime utility functions
  ///////////////////////////////////////////////////////

  static int createLimeRecordHeader(std::string message, int MB, int ME, size_t PayloadSize, LimeWriter* L)
  {
    LimeRecordHeader *h;
    h = limeCreateHeader(MB, ME, const_cast<char *>(message.c_str()), PayloadSize);
    assert(limeWriteRecordHeader(h, L) >= 0);
    limeDestroyHeader(h);
    return LIME_SUCCESS;
  }

  ////////////////////////////////////////////
  // Write a generic serialisable object
  ////////////////////////////////////////////
  template<class serialisable_object>
  static void writeLimeObject(int MB,int ME,serialisable_object &object,std::string object_name,std::string record_name, LimeWriter *LimeW)
  {
    std::string xmlstring;
    {
      XmlWriter WR("","");
      write(WR,object_name,object);
      xmlstring = WR.XmlString();
    }
    uint64_t nbytes = xmlstring.size();
    LimeRecordHeader *h = limeCreateHeader(MB, ME,(char *)record_name.c_str(), nbytes);
    int err=limeWriteRecordHeader(h, LimeW); assert(err>=0);
    err=limeWriteRecordData(&xmlstring[0], &nbytes, LimeW); assert(err>=0);
    err=limeWriterCloseRecord(LimeW);  assert(err>=0);
    limeDestroyHeader(h);
  }
  ////////////////////////////////////////////
  // Read a generic serialisable object
  ////////////////////////////////////////////
  template<class serialisable_object>
  static void readLimeObject(serialisable_object &object,std::string object_name,std::string record_name, LimeReader *LimeR)
  {
    std::string xmlstring;
    // should this be a do while; can we miss a first record??
    while ( limeReaderNextRecord(LimeR) == LIME_SUCCESS ) { 

      uint64_t nbytes = limeReaderBytes(LimeR);//size of this record (configuration)

      if ( strncmp(limeReaderType(LimeR), record_name.c_str(),strlen(record_name.c_str()) )  ) {
	std::vector<char> xmlc(nbytes+1,'\0');
	limeReaderReadData((void *)&xmlc[0], &nbytes, LimeR);    
	XmlReader RD(&xmlc[0],"");
	read(RD,object_name,object);
	return;
      }

    }  
    assert(0);
  }

  ////////////////////////////////////////////
  // Read a generic lattice field and verify checksum
  ////////////////////////////////////////////
  template<class vobj>
  static void readLimeLatticeBinaryObject(Lattice<vobj> &field,std::string filename,std::string record_name,FILE *File, LimeReader *LimeR)
  {
    typedef typename vobj::scalar_object sobj;
    scidacChecksum scidacChecksum_;
    uint32_t nersc_csum,scidac_csuma,scidac_csumb;

    std::string format = getFormatString<vobj>();

    while ( limeReaderNextRecord(LimeR) == LIME_SUCCESS ) { 

      std::cout << GridLogMessage << limeReaderType(LimeR) <<std::endl;
	
      if ( strncmp(limeReaderType(LimeR), record_name.c_str(),strlen(record_name.c_str()) )  ) {


	off_t offset= ftell(File);
	BinarySimpleMunger<sobj,sobj> munge;
	BinaryIO::readLatticeObject< sobj, sobj >(field, filename, munge, offset, format,nersc_csum,scidac_csuma,scidac_csumb);

	/////////////////////////////////////////////
	// Insist checksum is next record
	/////////////////////////////////////////////
	readLimeObject(scidacChecksum_,std::string("scidacChecksum"),record_name,LimeR);

	/////////////////////////////////////////////
	// Verify checksums
	/////////////////////////////////////////////
	scidacChecksumVerify(scidacChecksum_,scidac_csuma,scidac_csumb);
	return;
      }
    }
  }

  ////////////////////////////////////////////
  // Write a generic lattice field and csum
  ////////////////////////////////////////////
  template<class vobj>
  static void writeLimeLatticeBinaryObject(Lattice<vobj> &field,std::string filename,std::string record_name,FILE *File, LimeWriter *LimeW)
  {

    ////////////////////////////////////////////
    // Create record header
    ////////////////////////////////////////////
    typedef typename vobj::scalar_object sobj;
    int err;
    uint32_t nersc_csum,scidac_csuma,scidac_csumb;
    uint64_t PayloadSize = sizeof(sobj) * field._grid->_gsites;
    createLimeRecordHeader(record_name, 0, 0, PayloadSize, LimeW);

    ////////////////////////////////////////////////////////////////////
    // NB: FILE and iostream are jointly writing disjoint sequences in the
    // the same file through different file handles (integer units).
    // 
    // These are both buffered, so why I think this code is right is as follows.
    //
    // i)  write record header to FILE *File, telegraphing the size. 
    // ii) ftell reads the offset from FILE *File .
    // iii) iostream / MPI Open independently seek this offset. Write sequence direct to disk.
    //      Closes iostream and flushes.
    // iv) fseek on FILE * to end of this disjoint section.
    //  v) Continue writing scidac record.
    ////////////////////////////////////////////////////////////////////
    off_t offset = ftell(File);
    std::string format = getFormatString<vobj>();
    BinarySimpleMunger<sobj,sobj> munge;
    BinaryIO::writeLatticeObject<vobj,sobj>(field, filename, munge, offset, format,nersc_csum,scidac_csuma,scidac_csumb);
    err=limeWriterCloseRecord(LimeW);  assert(err>=0);
    ////////////////////////////////////////
    // Write checksum element, propagaing forward from the BinaryIO
    // Always pair a checksum with a binary object, and close message
    ////////////////////////////////////////
    scidacChecksum checksum;
    std::stringstream streama; streama << std::hex << scidac_csuma;
    std::stringstream streamb; streamb << std::hex << scidac_csumb;
    checksum.suma= streama.str();
    checksum.sumb= streamb.str();
    std::cout << GridLogMessage<<" writing scidac checksums "<<std::hex<<scidac_csuma<<"/"<<scidac_csumb<<std::dec<<std::endl;
    writeLimeObject(0,1,checksum,std::string("scidacChecksum"    ),std::string(SCIDAC_CHECKSUM),LimeW);
  }
  // Could end the LIME base class here
};

class ScidacIO : public LimeIO {
 public:
   /*
    LimeWriter *LimeW;
    LimeReader *LimeR;
    FILE *File;
  template<class userFile>
  int open(std::string filename,GridBase *grid,userFile &_userFile,int volfmt) {
 
  }
  void close(void) {
 
  }
  template<class vobj,class userRecord>
  int writeScidacField(Lattice<vobj> &field,userRecord &_userRecord,int volfmt) 
  template<class vobj,class userRecord>
  int  readScidacField(Lattice<vobj> &field,userRecord &_userRecord,int volfmt) 
   */
  ////////////////////////////////////////////////
  // Write generic lattice field in scidac format
  ////////////////////////////////////////////////
  template <class vobj,class userFile, class userRecord>
  static void writeScidacField(std::string filename,Lattice<vobj> &field,userFile _userFile,userRecord _userRecord) 
  {
    typedef typename vobj::scalar_object sobj;
    uint64_t nbytes;
    GridBase * grid = field._grid;

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
    FILE *File = fopen(filename.c_str(), "w");
    LimeWriter *LimeW = limeCreateWriter(File);
    assert(LimeW != NULL );

    writeLimeObject(1,0,header ,std::string("FieldMetaData"),std::string(GRID_FORMAT),LimeW); // Open message 
    writeLimeObject(0,0,_userFile,_userFile.SerialisableClassName(),std::string(SCIDAC_FILE_XML),LimeW);
    writeLimeObject(0,0,_scidacFile,_scidacFile.SerialisableClassName(),std::string(SCIDAC_PRIVATE_FILE_XML),LimeW);
    writeLimeObject(0,0,_userRecord,_userRecord.SerialisableClassName(),std::string(SCIDAC_RECORD_XML),LimeW);
    writeLimeObject(0,0,_scidacRecord,_scidacRecord.SerialisableClassName(),std::string(SCIDAC_PRIVATE_RECORD_XML),LimeW);
    writeLimeLatticeBinaryObject(field,filename,std::string(ILDG_BINARY_DATA),File,LimeW);      // Closes message with checksum

    limeDestroyWriter(LimeW);
    fclose(File);
  }
};

class IldgIO : public ScidacIO {
 public:

  ///////////////////////////////////
  // A little helper
  ///////////////////////////////////
  static void writeLimeIldgLFN(std::string &LFN,LimeWriter *LimeW)
  {
    uint64_t PayloadSize = LFN.size();
    int err;
    createLimeRecordHeader(ILDG_DATA_LFN, 0 , 0, PayloadSize, LimeW);
    err=limeWriteRecordData(const_cast<char*>(LFN.c_str()), &PayloadSize, LimeW); assert(err>=0);
    err=limeWriterCloseRecord(LimeW); assert(err>=0);
  }

  ////////////////////////////////////////////////////////////////
  // Special ILDG operations ; gauge configs only.
  // Don't require scidac records EXCEPT checksum
  // Use Grid MetaData object if present.
  ////////////////////////////////////////////////////////////////
  template <class vsimd>
  static void writeConfiguration(std::string filename,Lattice<iLorentzColourMatrix<vsimd> > &Umu) 
  {
    GridBase * grid = Umu._grid;
    typedef Lattice<iLorentzColourMatrix<vsimd> > GaugeField;
    typedef iLorentzColourMatrix<vsimd> vobj;
    typedef typename vobj::scalar_object sobj;

    uint64_t nbytes;

    ////////////////////////////////////////
    // fill the Grid header
    ////////////////////////////////////////
    FieldMetaData header;
    scidacRecord  _scidacRecord;
    scidacFile    _scidacFile;

    ScidacMetaData(Umu,header,_scidacRecord,_scidacFile);

    std::string format = header.floating_point;

    assert ( (format == std::string("IEEE32BIG"))  
           ||(format == std::string("IEEE64BIG")) );

    //////////////////////////////////////////////////////
    // Fill ILDG header data struct
    //////////////////////////////////////////////////////
    ildgFormat ildgfmt ;
    ildgfmt.field     = std::string("su3gauge");

    if ( format == std::string("IEEE32BIG") ) { 
      ildgfmt.precision = 32;
    } else { 
      ildgfmt.precision = 64;
    }
    ildgfmt.version = 1.0;
    ildgfmt.lx = header.dimension[0];
    ildgfmt.ly = header.dimension[1];
    ildgfmt.lz = header.dimension[2];
    ildgfmt.lt = header.dimension[3];
    assert(header.nd==4);
    assert(header.nd==header.dimension.size());

    //////////////////////////////////////////////////////////////////////////////
    // Fill the USQCD info field
    //////////////////////////////////////////////////////////////////////////////
    usqcdInfo info;
    info.version=1.0;
    info.plaq   = header.plaquette;
    info.linktr = header.link_trace;

    std::cout << GridLogMessage << " Writing config; IldgIO "<<std::endl;
    //////////////////////////////////////////////
    // Fill the Lime file record by record
    //////////////////////////////////////////////

    FILE *File = fopen(filename.c_str(), "w");
    LimeWriter *LimeW = limeCreateWriter(File); assert(LimeW != NULL);
    writeLimeObject(1,0,header ,std::string("FieldMetaData"),std::string(GRID_FORMAT),LimeW); // Open message 
    writeLimeObject(0,0,info,info.SerialisableClassName(),std::string(SCIDAC_FILE_XML),LimeW);
    writeLimeObject(0,0,_scidacFile,_scidacFile.SerialisableClassName(),std::string(SCIDAC_PRIVATE_FILE_XML),LimeW);
    writeLimeObject(0,0,info,info.SerialisableClassName(),std::string(SCIDAC_RECORD_XML),LimeW);
    writeLimeObject(0,0,_scidacRecord,_scidacRecord.SerialisableClassName(),std::string(SCIDAC_PRIVATE_RECORD_XML),LimeW);
    writeLimeObject(0,0,ildgfmt,std::string("ildgFormat")   ,std::string(ILDG_FORMAT),LimeW); // rec
    writeLimeIldgLFN(header.ildg_lfn, LimeW);                                                 // rec
    writeLimeLatticeBinaryObject(Umu,filename,std::string(ILDG_BINARY_DATA),File,LimeW);      // Closes message with checksum
    limeDestroyWriter(LimeW);
    fclose(File);
  }

  ////////////////////////////////////////////////////////////////
  // Read either Grid/SciDAC/ILDG configuration
  // Don't require scidac records EXCEPT checksum
  // Use Grid MetaData object if present.
  // Else use ILDG MetaData object if present.
  // Else use SciDAC MetaData object if present.
  ////////////////////////////////////////////////////////////////
  template <class vsimd>
  static void readConfiguration(std::string filename,Lattice<iLorentzColourMatrix<vsimd> > &Umu, FieldMetaData &FieldMetaData_) {

    typedef Lattice<iLorentzColourMatrix<vsimd> > GaugeField;
    typedef typename GaugeField::vector_object  vobj;
    typedef typename vobj::scalar_object sobj;

    typedef LorentzColourMatrixF fobj;
    typedef LorentzColourMatrixD dobj;

    GridBase *grid = Umu._grid;

    std::vector<int> dims = Umu._grid->FullDimensions();

    assert(dims.size()==4);

    FILE *File = fopen(filename.c_str(), "r");
    LimeReader *LimeR = limeCreateReader(File);

    // Metadata holders
    ildgFormat     ildgFormat_    ;
    std::string    ildgLFN_       ;
    scidacChecksum scidacChecksum_; 
    usqcdInfo      usqcdInfo_     ;

    // track what we read from file
    int found_ildgFormat    =0;
    int found_ildgLFN       =0;
    int found_scidacChecksum=0;
    int found_usqcdInfo     =0;
    int found_ildgBinary =0;
    int found_FieldMetaData =0;

    uint32_t nersc_csum;
    uint32_t scidac_csuma;
    uint32_t scidac_csumb;

    // Binary format
    std::string format;

    //////////////////////////////////////////////////////////////////////////
    // Loop over all records
    // -- Order is poorly guaranteed except ILDG header preceeds binary section.
    // -- Run like an event loop.
    // -- Impose trust hierarchy. Grid takes precedence & look for ILDG, and failing
    //    that Scidac. 
    // -- Insist on Scidac checksum record.
    //////////////////////////////////////////////////////////////////////////

    while ( limeReaderNextRecord(LimeR) == LIME_SUCCESS ) { 

      uint64_t nbytes = limeReaderBytes(LimeR);//size of this record (configuration)
      
      //////////////////////////////////////////////////////////////////
      // If not BINARY_DATA read a string and parse
      //////////////////////////////////////////////////////////////////
      if ( strncmp(limeReaderType(LimeR), ILDG_BINARY_DATA,strlen(ILDG_BINARY_DATA) )  ) {
	
	// Copy out the string
	std::vector<char> xmlc(nbytes+1,'\0');
	limeReaderReadData((void *)&xmlc[0], &nbytes, LimeR);    
	std::cout << GridLogMessage<< "Non binary record :" <<limeReaderType(LimeR) <<std::endl; //<<"\n"<<(&xmlc[0])<<std::endl;

	//////////////////////////////////
	// ILDG format record
	if ( !strncmp(limeReaderType(LimeR), ILDG_FORMAT,strlen(ILDG_FORMAT)) ) { 

	  XmlReader RD(&xmlc[0],"");
	  read(RD,"ildgFormat",ildgFormat_);

	  if ( ildgFormat_.precision == 64 ) format = std::string("IEEE64BIG");
	  if ( ildgFormat_.precision == 32 ) format = std::string("IEEE32BIG");

	  assert( ildgFormat_.lx == dims[0]);
	  assert( ildgFormat_.ly == dims[1]);
	  assert( ildgFormat_.lz == dims[2]);
	  assert( ildgFormat_.lt == dims[3]);

	  found_ildgFormat = 1;
	}

	if ( !strncmp(limeReaderType(LimeR), ILDG_DATA_LFN,strlen(ILDG_DATA_LFN)) ) {
	  FieldMetaData_.ildg_lfn = std::string(&xmlc[0]);
	  found_ildgLFN = 1;
	}

	if ( !strncmp(limeReaderType(LimeR), GRID_FORMAT,strlen(ILDG_FORMAT)) ) { 

	  XmlReader RD(&xmlc[0],"");
	  read(RD,"FieldMetaData",FieldMetaData_);

	  format = FieldMetaData_.floating_point;

	  assert(FieldMetaData_.dimension[0] == dims[0]);
	  assert(FieldMetaData_.dimension[1] == dims[1]);
	  assert(FieldMetaData_.dimension[2] == dims[2]);
	  assert(FieldMetaData_.dimension[3] == dims[3]);

	  found_FieldMetaData = 1;
	}

	if ( !strncmp(limeReaderType(LimeR), SCIDAC_RECORD_XML,strlen(SCIDAC_RECORD_XML)) ) { 
	  XmlReader RD(&xmlc[0],"");
	  read(RD,"usqcdInfo",usqcdInfo_);
	  found_usqcdInfo = 1;
	}

	if ( !strncmp(limeReaderType(LimeR), SCIDAC_CHECKSUM,strlen(SCIDAC_CHECKSUM)) ) { 
	  XmlReader RD(&xmlc[0],"");
	  read(RD,"scidacChecksum",scidacChecksum_);
	  found_scidacChecksum = 1;
	}

      } else {  
	/////////////////////////////////
	// Binary data
	/////////////////////////////////
	std::cout << GridLogMessage << "ILDG Binary record found : "  ILDG_BINARY_DATA << std::endl;
	off_t offset= ftell(File);

	if ( format == std::string("IEEE64BIG") ) {
	  GaugeSimpleMunger<dobj, sobj> munge;
	  BinaryIO::readLatticeObject< vobj, dobj >(Umu, filename, munge, offset, format,nersc_csum,scidac_csuma,scidac_csumb);
	} else { 
	  GaugeSimpleMunger<fobj, sobj> munge;
	  BinaryIO::readLatticeObject< vobj, fobj >(Umu, filename, munge, offset, format,nersc_csum,scidac_csuma,scidac_csumb);
	}

	found_ildgBinary = 1;
      }

    }

    //////////////////////////////////////////////////////
    // Minimally must find binary segment and checksum
    // Since this is an ILDG reader require ILDG format
    //////////////////////////////////////////////////////
    assert(found_ildgBinary);
    assert(found_ildgFormat);
    assert(found_scidacChecksum);

    // Must find something with the lattice dimensions
    assert(found_FieldMetaData||found_ildgFormat);

    if ( found_FieldMetaData ) {

      std::cout << GridLogMessage<<"Grid MetaData was record found: configuration was probably written by Grid ! Yay ! "<<std::endl;

    } else { 

      assert(found_ildgFormat);
      assert ( ildgFormat_.field == std::string("su3gauge") );

      ///////////////////////////////////////////////////////////////////////////////////////
      // Populate our Grid metadata as best we can
      ///////////////////////////////////////////////////////////////////////////////////////

      std::ostringstream vers; vers << ildgFormat_.version;
      FieldMetaData_.hdr_version = vers.str();
      FieldMetaData_.data_type = std::string("4D_SU3_GAUGE_3X3");

      FieldMetaData_.nd=4;
      FieldMetaData_.dimension.resize(4);

      FieldMetaData_.dimension[0] = ildgFormat_.lx ;
      FieldMetaData_.dimension[1] = ildgFormat_.ly ;
      FieldMetaData_.dimension[2] = ildgFormat_.lz ;
      FieldMetaData_.dimension[3] = ildgFormat_.lt ;

      if ( found_usqcdInfo ) { 
	FieldMetaData_.plaquette = usqcdInfo_.plaq;
	FieldMetaData_.link_trace= usqcdInfo_.linktr;
	std::cout << GridLogMessage <<"This configuration was probably written by USQCD "<<std::endl;
	std::cout << GridLogMessage <<"USQCD xml record Plaquette : "<<FieldMetaData_.plaquette<<std::endl;
	std::cout << GridLogMessage <<"USQCD xml record LinkTrace : "<<FieldMetaData_.link_trace<<std::endl;
      } else { 
	FieldMetaData_.plaquette = 0.0;
	FieldMetaData_.link_trace= 0.0;
	std::cout << GridLogWarning << "This configuration is unsafe with no plaquette records that can verify it !!! "<<std::endl;
      }
    }

    ////////////////////////////////////////////////////////////
    // Really really want to mandate a scidac checksum
    ////////////////////////////////////////////////////////////
    if ( found_scidacChecksum ) {
      FieldMetaData_.scidac_checksuma = stoull(scidacChecksum_.suma,0,16);
      FieldMetaData_.scidac_checksumb = stoull(scidacChecksum_.sumb,0,16);
      scidacChecksumVerify(scidacChecksum_,scidac_csuma,scidac_csumb);
      assert( scidac_csuma ==FieldMetaData_.scidac_checksuma);
      assert( scidac_csumb ==FieldMetaData_.scidac_checksumb);
      std::cout << GridLogMessage<<"SciDAC checksums match " << std::endl;
    } else { 
      std::cout << GridLogWarning<<"SciDAC checksums not found. This is unsafe. " << std::endl;
      assert(0); // Can I insist always checksum ?
    }

    if ( found_FieldMetaData || found_usqcdInfo ) {
      FieldMetaData checker;
      GaugeStatistics(Umu,checker);
      assert(fabs(checker.plaquette  - FieldMetaData_.plaquette )<1.0e-5);
      assert(fabs(checker.link_trace - FieldMetaData_.link_trace)<1.0e-5);
      std::cout << GridLogMessage<<"Plaquette and link trace match " << std::endl;
    }
  }
 };

}}

//HAVE_LIME
#endif

#endif
