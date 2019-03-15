/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/parallelIO/NerscIO.h

    Copyright (C) 2015

    Author: Matt Spraggs <matthew.spraggs@gmail.com>
    Author: Peter Boyle <paboyle@ph.ed.ac.uk>
    Author: paboyle <paboyle@ph.ed.ac.uk>

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
#ifndef GRID_NERSC_IO_H
#define GRID_NERSC_IO_H

namespace Grid {
  namespace QCD {

    using namespace Grid;

    ////////////////////////////////////////////////////////////////////////////////
    // Write and read from fstream; comput header offset for payload
    ////////////////////////////////////////////////////////////////////////////////
    class NerscIO : public BinaryIO { 
    public:

      static inline void truncate(std::string file){
	std::ofstream fout(file,std::ios::out);
      }
  
      static inline unsigned int writeHeader(FieldMetaData &field,std::string file)
      {
      std::ofstream fout(file,std::ios::out|std::ios::in);
      fout.seekp(0,std::ios::beg);
      dump_meta_data(field, fout);
      field.data_start = fout.tellp();
      return field.data_start;
    }

      // for the header-reader
      static inline int readHeader(std::string file,GridBase *grid,  FieldMetaData &field)
      {
      uint64_t offset=0;
      std::map<std::string,std::string> header;
      std::string line;

      //////////////////////////////////////////////////
      // read the header
      //////////////////////////////////////////////////
      std::ifstream fin(file);

      getline(fin,line); // read one line and insist is 

      removeWhitespace(line);
      std::cout << GridLogMessage << "* " << line << std::endl;

      assert(line==std::string("BEGIN_HEADER"));

      do {
      getline(fin,line); // read one line
      std::cout << GridLogMessage << "* "<<line<< std::endl;
      int eq = line.find("=");
      if(eq >0) {
      std::string key=line.substr(0,eq);
      std::string val=line.substr(eq+1);
      removeWhitespace(key);
      removeWhitespace(val);
      
      header[key] = val;
    }
    } while( line.find("END_HEADER") == std::string::npos );

      field.data_start = fin.tellg();

      //////////////////////////////////////////////////
      // chomp the values
      //////////////////////////////////////////////////
      field.hdr_version    = header["HDR_VERSION"];
      field.data_type      = header["DATATYPE"];
      field.storage_format = header["STORAGE_FORMAT"];
  
      field.dimension[0] = std::stol(header["DIMENSION_1"]);
      field.dimension[1] = std::stol(header["DIMENSION_2"]);
      field.dimension[2] = std::stol(header["DIMENSION_3"]);
      field.dimension[3] = std::stol(header["DIMENSION_4"]);

      assert(grid->_ndimension == 4);
      for(int d=0;d<4;d++){
      assert(grid->_fdimensions[d]==field.dimension[d]);
    }

      field.link_trace = std::stod(header["LINK_TRACE"]);
      field.plaquette  = std::stod(header["PLAQUETTE"]);

      field.boundary[0] = header["BOUNDARY_1"];
      field.boundary[1] = header["BOUNDARY_2"];
      field.boundary[2] = header["BOUNDARY_3"];
      field.boundary[3] = header["BOUNDARY_4"];

      field.checksum = std::stoul(header["CHECKSUM"],0,16);
      field.ensemble_id      = header["ENSEMBLE_ID"];
      field.ensemble_label   = header["ENSEMBLE_LABEL"];
      field.sequence_number  = std::stol(header["SEQUENCE_NUMBER"]);
      field.creator          = header["CREATOR"];
      field.creator_hardware = header["CREATOR_HARDWARE"];
      field.creation_date    = header["CREATION_DATE"];
      field.archive_date     = header["ARCHIVE_DATE"];
      field.floating_point   = header["FLOATING_POINT"];

      return field.data_start;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Now the meat: the object readers
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template<class vsimd>
    static inline void readConfiguration(Lattice<iLorentzColourMatrix<vsimd> > &Umu,
					 FieldMetaData& header,
					 std::string file)
    {
      typedef Lattice<iLorentzColourMatrix<vsimd> > GaugeField;

      GridBase *grid = Umu._grid;
      uint64_t offset = readHeader(file,Umu._grid,header);

      FieldMetaData clone(header);

      std::string format(header.floating_point);

      int ieee32big = (format == std::string("IEEE32BIG"));
      int ieee32    = (format == std::string("IEEE32"));
      int ieee64big = (format == std::string("IEEE64BIG"));
      int ieee64    = (format == std::string("IEEE64"));

      uint32_t nersc_csum,scidac_csuma,scidac_csumb;
      // depending on datatype, set up munger;
      // munger is a function of <floating point, Real, data_type>
      if ( header.data_type == std::string("4D_SU3_GAUGE") ) {
	if ( ieee32 || ieee32big ) {
	  BinaryIO::readLatticeObject<iLorentzColourMatrix<vsimd>, LorentzColour2x3F> 
	    (Umu,file,Gauge3x2munger<LorentzColour2x3F,LorentzColourMatrix>(), offset,format,
	     nersc_csum,scidac_csuma,scidac_csumb);
	}
	if ( ieee64 || ieee64big ) {
	  BinaryIO::readLatticeObject<iLorentzColourMatrix<vsimd>, LorentzColour2x3D> 
	    (Umu,file,Gauge3x2munger<LorentzColour2x3D,LorentzColourMatrix>(),offset,format,
	     nersc_csum,scidac_csuma,scidac_csumb);
	}
      } else if ( header.data_type == std::string("4D_SU3_GAUGE_3x3") ) {
	if ( ieee32 || ieee32big ) {
	  BinaryIO::readLatticeObject<iLorentzColourMatrix<vsimd>,LorentzColourMatrixF>
	    (Umu,file,GaugeSimpleMunger<LorentzColourMatrixF,LorentzColourMatrix>(),offset,format,
	     nersc_csum,scidac_csuma,scidac_csumb);
	}
	if ( ieee64 || ieee64big ) {
	  BinaryIO::readLatticeObject<iLorentzColourMatrix<vsimd>,LorentzColourMatrixD>
	    (Umu,file,GaugeSimpleMunger<LorentzColourMatrixD,LorentzColourMatrix>(),offset,format,
	     nersc_csum,scidac_csuma,scidac_csumb);
	}
      } else {
	assert(0);
      }

      GaugeStatistics(Umu,clone);

      std::cout<<GridLogMessage <<"NERSC Configuration "<<file<<" checksum "<<std::hex<<nersc_csum<< std::dec
	       <<" header   "<<std::hex<<header.checksum<<std::dec <<std::endl;
      std::cout<<GridLogMessage <<"NERSC Configuration "<<file<<" plaquette "<<clone.plaquette
	       <<" header    "<<header.plaquette<<std::endl;
      std::cout<<GridLogMessage <<"NERSC Configuration "<<file<<" link_trace "<<clone.link_trace
	       <<" header    "<<header.link_trace<<std::endl;

      if ( fabs(clone.plaquette -header.plaquette ) >=  1.0e-5 ) { 
	std::cout << " Plaquette mismatch "<<std::endl;
	std::cout << Umu[0]<<std::endl;
	std::cout << Umu[1]<<std::endl;
      }
      if ( nersc_csum != header.checksum ) { 
	std::cerr << " checksum mismatch " << std::endl;
	std::cerr << " plaqs " << clone.plaquette << " " << header.plaquette << std::endl;
	std::cerr << " trace " << clone.link_trace<< " " << header.link_trace<< std::endl;
	std::cerr << " nersc_csum  " <<std::hex<< nersc_csum << " " << header.checksum<< std::dec<< std::endl;
	exit(0);
      }
      assert(fabs(clone.plaquette -header.plaquette ) < 1.0e-5 );
      assert(fabs(clone.link_trace-header.link_trace) < 1.0e-6 );
      assert(nersc_csum == header.checksum );
      
      std::cout<<GridLogMessage <<"NERSC Configuration "<<file<< " and plaquette, link trace, and checksum agree"<<std::endl;
    }

      template<class vsimd>
      static inline void writeConfiguration(Lattice<iLorentzColourMatrix<vsimd> > &Umu,
					    std::string file, 
					    int two_row,
					    int bits32)
      {
	typedef Lattice<iLorentzColourMatrix<vsimd> > GaugeField;

	typedef iLorentzColourMatrix<vsimd> vobj;
	typedef typename vobj::scalar_object sobj;

	FieldMetaData header;
	///////////////////////////////////////////
	// Following should become arguments
	///////////////////////////////////////////
	header.sequence_number = 1;
	header.ensemble_id     = "UKQCD";
	header.ensemble_label  = "DWF";

	typedef LorentzColourMatrixD fobj3D;
	typedef LorentzColour2x3D    fobj2D;
  
	GridBase *grid = Umu._grid;

	GridMetaData(grid,header);
	assert(header.nd==4);
	GaugeStatistics(Umu,header);
	MachineCharacteristics(header);

	uint64_t offset;

	// Sod it -- always write 3x3 double
	header.floating_point = std::string("IEEE64BIG");
	header.data_type      = std::string("4D_SU3_GAUGE_3x3");
	GaugeSimpleUnmunger<fobj3D,sobj> munge;
	if ( grid->IsBoss() ) { 
	  truncate(file);
	  offset = writeHeader(header,file);
	}
	grid->Broadcast(0,(void *)&offset,sizeof(offset));

	uint32_t nersc_csum,scidac_csuma,scidac_csumb;
	BinaryIO::writeLatticeObject<vobj,fobj3D>(Umu,file,munge,offset,header.floating_point,
								  nersc_csum,scidac_csuma,scidac_csumb);
	header.checksum = nersc_csum;
	if ( grid->IsBoss() ) { 
	  writeHeader(header,file);
	}

	std::cout<<GridLogMessage <<"Written NERSC Configuration on "<< file << " checksum "
		 <<std::hex<<header.checksum
		 <<std::dec<<" plaq "<< header.plaquette <<std::endl;

      }
      ///////////////////////////////
      // RNG state
      ///////////////////////////////
      static inline void writeRNGState(GridSerialRNG &serial,GridParallelRNG &parallel,std::string file)
      {
	typedef typename GridParallelRNG::RngStateType RngStateType;

	// Following should become arguments
	FieldMetaData header;
	header.sequence_number = 1;
	header.ensemble_id     = "UKQCD";
	header.ensemble_label  = "DWF";

	GridBase *grid = parallel._grid;

	GridMetaData(grid,header);
	assert(header.nd==4);
	header.link_trace=0.0;
	header.plaquette=0.0;
	MachineCharacteristics(header);

	uint64_t offset;
  
#ifdef RNG_RANLUX
	header.floating_point = std::string("UINT64");
	header.data_type      = std::string("RANLUX48");
#endif
#ifdef RNG_MT19937
	header.floating_point = std::string("UINT32");
	header.data_type      = std::string("MT19937");
#endif
#ifdef RNG_SITMO
	header.floating_point = std::string("UINT64");
	header.data_type      = std::string("SITMO");
#endif

	if ( grid->IsBoss() ) { 
	  truncate(file);
	  offset = writeHeader(header,file);
	}
	grid->Broadcast(0,(void *)&offset,sizeof(offset));
	
	uint32_t nersc_csum,scidac_csuma,scidac_csumb;
	BinaryIO::writeRNG(serial,parallel,file,offset,nersc_csum,scidac_csuma,scidac_csumb);
	header.checksum = nersc_csum;
	if ( grid->IsBoss() ) { 
	  offset = writeHeader(header,file);
	}

	std::cout<<GridLogMessage 
		 <<"Written NERSC RNG STATE "<<file<< " checksum "
		 <<std::hex<<header.checksum
		 <<std::dec<<std::endl;

      }
    
      static inline void readRNGState(GridSerialRNG &serial,GridParallelRNG & parallel,FieldMetaData& header,std::string file)
      {
	typedef typename GridParallelRNG::RngStateType RngStateType;

	GridBase *grid = parallel._grid;

	uint64_t offset = readHeader(file,grid,header);

	FieldMetaData clone(header);

	std::string format(header.floating_point);
	std::string data_type(header.data_type);

#ifdef RNG_RANLUX
	assert(format == std::string("UINT64"));
	assert(data_type == std::string("RANLUX48"));
#endif
#ifdef RNG_MT19937
	assert(format == std::string("UINT32"));
	assert(data_type == std::string("MT19937"));
#endif
#ifdef RNG_SITMO
	assert(format == std::string("UINT64"));
	assert(data_type == std::string("SITMO"));
#endif

	// depending on datatype, set up munger;
	// munger is a function of <floating point, Real, data_type>
	uint32_t nersc_csum,scidac_csuma,scidac_csumb;
	BinaryIO::readRNG(serial,parallel,file,offset,nersc_csum,scidac_csuma,scidac_csumb);

	if ( nersc_csum != header.checksum ) { 
	  std::cerr << "checksum mismatch "<<std::hex<< nersc_csum <<" "<<header.checksum<<std::dec<<std::endl;
	  exit(0);
	}
	assert(nersc_csum == header.checksum );

	std::cout<<GridLogMessage <<"Read NERSC RNG file "<<file<< " format "<< data_type <<std::endl;
      }

    };

  }}
#endif
