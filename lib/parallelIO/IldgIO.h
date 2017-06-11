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

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>

#include <pwd.h>
#include <sys/utsname.h>
#include <unistd.h>

#ifdef HAVE_LIME

extern "C" {  // for linkage
#include "lime.h"
}


// Unused SCIDAC records names
// SCIDAC_PRIVATE_FILE_XML   "scidac-private-file-xml"
// SCIDAC_SITELIST           "scidac-sitelist"
// SCIDAC_FILE_XML           "scidac-file-xml"
// SCIDAC_RIVATE_RECORD_XML "scidac-private-record-xml"
// SCIDAC_RECORD_XML         "scidac-record-xml"
// SCIDAC_BINARY_DATA        "scidac-binary-data"
//
// Scidac checksum: CRC32 every site, xor reduce some hash of this.
// https://github.com/usqcd-software/qio/blob/master/lib/dml/DML_utils.c

namespace Grid {
namespace QCD {

class IldgIO : public BinaryIO {
 public:

  static int createHeader(std::string message, int MB, int ME, size_t PayloadSize, LimeWriter* L)
  {
    LimeRecordHeader *h;
    h = limeCreateHeader(MB, ME, const_cast<char *>(message.c_str()), PayloadSize);
    assert(limeWriteRecordHeader(h, L) >= 0);
    limeDestroyHeader(h);
    return LIME_SUCCESS;
  }

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
    assert(limeWriteRecordHeader(h, LimeW)>=0);
    assert(limeWriteRecordData(&xmlstring[0], &nbytes, LimeW)>=0);
    limeWriterCloseRecord(LimeW);
    limeDestroyHeader(h);
  }

  static unsigned int writeHeader(FieldMetaData &header, LimeWriter *LimeW) {

    uint64_t nbytes;

    ildgFormat ildgfmt ;
    usqcdInfo info;

    //////////////////////////////////////////////////////
    // Fill ILDG header data struct
    //////////////////////////////////////////////////////
    ildgfmt.field     = std::string("su3gauge");
    ildgfmt.precision = 64;
    ildgfmt.version = 1.0;
    ildgfmt.lx = header.dimension[0];
    ildgfmt.ly = header.dimension[1];
    ildgfmt.lz = header.dimension[2];
    ildgfmt.lt = header.dimension[3];
    assert(header.nd==4);
    assert(header.nd==header.dimension.size());

    info.version=1.0;
    info.plaq   = header.plaquette;
    info.linktr = header.link_trace;

    // Following scidac file downloaded from NERSC under MILC
    // Begin message, keep open on successive records
    //Message 1
    // Type:           scidac-private-file-xml <scidacFile><version>1.1</version><spacetime>4</spacetime><dims>16 16 16 48 </dims><volfmt>0</volfmt></scidacFile>
    // Type:           scidac-file-xml         <title>MILC ILDG archival gauge configuration</title>
    //Message 2
    // Type:           scidac-private-record-xml <scidacRecord><version>1.0</version><date>Thu May 11 00:11:33 2006 UTC</date><globaldata>0</globaldata>
    //                    <datatype>QDP_F3_ColorMatrix</datatype><precision>F</precision><colors>3</colors><typesize>72</typesize><datacount>4</datacount></scidacRecord>
    // Type:           scidac-record-xml 
    // Type:           ildg-format
    // Type:           ildg-data-lfn
    // Type:           ildg-binary-data
    // Type:           scidac-checksum

    writeLimeObject(1,0,header ,std::string("FieldMetaData"),std::string(GRID_FORMAT),LimeW);
    writeLimeObject(0,0,info   ,std::string("usqcdInfo"    ),std::string(USQCD_INFO ),LimeW);
    writeLimeObject(0,0,ildgfmt,std::string("ildgFormat")   ,std::string(ILDG_FORMAT),LimeW);
    // LFN is not a serializable object
    {
      std::string LFN = header.ildg_lfn; 
      uint64_t PayloadSize = LFN.size();
      createHeader(ILDG_DATA_LFN, 0 , 0, PayloadSize, LimeW);
      limeWriteRecordData(const_cast<char*>(LFN.c_str()), &PayloadSize, LimeW);
      limeWriterCloseRecord(LimeW);
    }
    return 0;
  }

  template <class vsimd>
  static void writeConfiguration(std::string filename,Lattice<iLorentzColourMatrix<vsimd> > &Umu, std::string format) {

    FILE *File = fopen(filename.c_str(), "w");
    LimeWriter *LimeW = limeCreateWriter(File);

    typedef Lattice<iLorentzColourMatrix<vsimd> > GaugeField;
    typedef iLorentzColourMatrix<vsimd> vobj;
    typedef typename vobj::scalar_object sobj;
    typedef LorentzColourMatrixD fobj;

    GridBase * grid = Umu._grid;

    ////////////////////////////////////////
    // fill the headers
    ////////////////////////////////////////
    FieldMetaData header;

    GridMetaData(grid,header); 
    GaugeStatistics<GaugeField>(Umu,header);
    MachineCharacteristics(header);

    assert( (format=="IEEE64BIG") || (format=="IEEE32BIG"));
    header.floating_point = format;
    header.checksum = 0x0; // unused in ILDG
    writeHeader(header,LimeW);

    ////////////////////////////////////////
    // Write data record header
    ////////////////////////////////////////
    uint64_t PayloadSize = sizeof(fobj) * Umu._grid->_gsites;
    createHeader(ILDG_BINARY_DATA, 0, 0, PayloadSize, LimeW);
    
    off_t offset = ftell(File);
    uint32_t nersc_csum,scidac_csuma,scidac_csumb;
    GaugeSimpleMunger<sobj, fobj> munge;
    BinaryIO::writeLatticeObject<vobj, fobj >(Umu, filename, munge, offset, header.floating_point,
					      nersc_csum,scidac_csuma,scidac_csumb);
    limeWriterCloseRecord(LimeW);

    ////////////////////////////////////////
    // Write checksum element, propagaing forward from the BinaryIO
    ////////////////////////////////////////
    scidacChecksum checksum;
    checksum.suma= scidac_csuma;
    checksum.sumb= scidac_csumb;
    //    std::cout << " writing scidac checksums "<<std::hex<<scidac_csuma<<"/"<<scidac_csumb<<std::dec<<std::endl;
    writeLimeObject(0,1,checksum,std::string("scidacChecksum"    ),std::string(SCIDAC_CHECKSUM),LimeW);

    fclose(File);
  }

  template <class vsimd>
  static void readConfiguration(std::string filename,Lattice<iLorentzColourMatrix<vsimd> > &Umu, FieldMetaData &FieldMetaData_) {

    typedef Lattice<iLorentzColourMatrix<vsimd> > GaugeField;
    typedef LorentzColourMatrixD sobjd;
    typedef LorentzColourMatrixF sobjf;
    typedef iLorentzColourMatrix<vsimd> itype;
    typedef LorentzColourMatrix sobj;

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

	  //	  std::cout << "This is an ILDG format record : "<<format<<std::endl;

	  assert( ildgFormat_.lx == dims[0]);
	  assert( ildgFormat_.ly == dims[1]);
	  assert( ildgFormat_.lz == dims[2]);
	  assert( ildgFormat_.lt == dims[3]);

	  found_ildgFormat = 1;
	}

	if ( !strncmp(limeReaderType(LimeR), ILDG_DATA_LFN,strlen(ILDG_DATA_LFN)) ) {
	  FieldMetaData_.ildg_lfn = std::string(&xmlc[0]);
	  //	  std::cout << "ILDG logical file name "<< FieldMetaData_.ildg_lfn << std::endl;
	  found_ildgLFN = 1;
	}

	if ( !strncmp(limeReaderType(LimeR), GRID_FORMAT,strlen(ILDG_FORMAT)) ) { 

	  XmlReader RD(&xmlc[0],"");
	  read(RD,"FieldMetaData",FieldMetaData_);

	  //	  std::cout << "Grid header found : format is "<<FieldMetaData_.floating_point<<std::endl;

	  format = FieldMetaData_.floating_point;

	  assert(FieldMetaData_.dimension[0] == dims[0]);
	  assert(FieldMetaData_.dimension[1] == dims[1]);
	  assert(FieldMetaData_.dimension[2] == dims[2]);
	  assert(FieldMetaData_.dimension[3] == dims[3]);

	  found_FieldMetaData = 1;
	}

	if ( !strncmp(limeReaderType(LimeR), USQCD_INFO,strlen(USQCD_INFO)) ) { 
	  XmlReader RD(&xmlc[0],"");
	  read(RD,USQCD_INFO,usqcdInfo_);
	  //	  std::cout << "USQCD info record found " <<std::endl;
	  found_usqcdInfo = 1;
	}

	if ( !strncmp(limeReaderType(LimeR), SCIDAC_CHECKSUM,strlen(SCIDAC_CHECKSUM)) ) { 
	  XmlReader RD(&xmlc[0],"");
	  read(RD,"scidacChecksum",scidacChecksum_);
	  FieldMetaData_.scidac_checksuma = scidacChecksum_.suma;
	  FieldMetaData_.scidac_checksumb = scidacChecksum_.sumb;
	  //std::cout << " Read Out "<<scidacChecksum_.version<<"/"<< scidacChecksum_.suma<<"/"<<scidacChecksum_.sumb<<std::endl;
	  found_scidacChecksum = 1;
	}

      } else {  
	/////////////////////////////////
	// Binary data
	/////////////////////////////////
	std::cout << GridLogMessage << ILDG_BINARY_DATA << std::endl;
	off_t offset= ftell(File);
	GaugeSimpleMunger<sobjd, sobj> munge;
	BinaryIO::readLatticeObject< itype, sobjd >(Umu, filename, munge, offset, format,
						    nersc_csum,scidac_csuma,scidac_csumb);
	found_ildgBinary = 1;
      }

    }

    //////////////////////////////////////////////////////
    // Minimally must find binary segment and checksum
    //////////////////////////////////////////////////////
    assert(found_ildgBinary);
    assert(found_scidacChecksum);

    // Must find something with the lattice dimensions
    assert(found_FieldMetaData||found_ildgFormat);

    if ( found_FieldMetaData ) {

      std::cout << GridLogMessage<<"a Grid MetaData was record found: configuration was probably written by Grid ! Yay ! "<<std::endl;
      //      std::cout << "Read Grid Plaqette  "<<FieldMetaData_.plaquette<<std::endl;
      //      std::cout << "Read Grid LinkTrace "<<FieldMetaData_.link_trace<<std::endl;

    } else { 

      assert(found_ildgFormat);
      assert ( ildgFormat_.field == std::string("su3gauge") );

      ///////////////////////////////////////////////////////////////////////////////////////
      // Populate our Grid metadata as best we can
      ///////////////////////////////////////////////////////////////////////////////////////

      std::ostringstream vers; vers << ildgFormat_.version;
      FieldMetaData_.hdr_version = vers.str();
      FieldMetaData_.data_type = std::string("4D_SU3_GAUGE_3X3");

      assert(FieldMetaData_.nd==4);
      assert(FieldMetaData_.dimension.size()==4);

      FieldMetaData_.dimension[0] = ildgFormat_.lx ;
      FieldMetaData_.dimension[1] = ildgFormat_.ly ;
      FieldMetaData_.dimension[2] = ildgFormat_.lz ;
      FieldMetaData_.dimension[3] = ildgFormat_.lt ;

      if ( found_usqcdInfo ) { 
	FieldMetaData_.plaquette = usqcdInfo_.plaq;
	FieldMetaData_.link_trace= usqcdInfo_.linktr;
	//	std::cout << "This configuration was probably written by USQCD and not Grid "<<std::endl;
	//	std::cout << "Read USQCD Plaquette  "<<FieldMetaData_.plaquette<<std::endl;
	//	std::cout << "Read USQCD LinkTrace  "<<FieldMetaData_.link_trace<<std::endl;
      } else { 
	FieldMetaData_.plaquette = 0.0;
	FieldMetaData_.link_trace= 0.0;
	std::cout << "Uhoh... This configuration is unsafe and contains no recognised checksum or physics records that can verify it !!! "<<std::endl;
      }
    }

    if ( found_scidacChecksum ) {
      assert( scidac_csuma ==FieldMetaData_.scidac_checksuma);
      assert( scidac_csumb ==FieldMetaData_.scidac_checksumb);
      std::cout << GridLogMessage<<"SciDAC checksums match " << std::endl;
    }

    if ( found_FieldMetaData || found_usqcdInfo ) {
      FieldMetaData checker;
      GaugeStatistics<GaugeField>(Umu,checker);
      assert(fabs(checker.plaquette  - FieldMetaData_.plaquette )<1.0e-5);
      assert(fabs(checker.link_trace - FieldMetaData_.link_trace)<1.0e-5);
      std::cout << GridLogMessage<<"Plaquette and link trace match " << std::endl;
    }
  }

  // format for RNG? Now just binary out
};
}
}

//HAVE_LIME
#endif

#endif
