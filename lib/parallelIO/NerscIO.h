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

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>

#include <unistd.h>
#include <sys/utsname.h>
#include <pwd.h>

namespace Grid {
namespace QCD {

using namespace Grid;

////////////////////////////////////////////////////////////////////////////////
// Some data types for intermediate storage
////////////////////////////////////////////////////////////////////////////////
  template<typename vtype> using iLorentzColour2x3 = iVector<iVector<iVector<vtype, Nc>, 2>, 4 >;

  typedef iLorentzColour2x3<Complex>  LorentzColour2x3;
  typedef iLorentzColour2x3<ComplexF> LorentzColour2x3F;
  typedef iLorentzColour2x3<ComplexD> LorentzColour2x3D;

////////////////////////////////////////////////////////////////////////////////
// header specification/interpretation
////////////////////////////////////////////////////////////////////////////////
class NerscField {
 public:
    // header strings (not in order)
    int dimension[4];
    std::string boundary[4]; 
    int data_start;
    std::string hdr_version;
    std::string storage_format;
    // Checks on data
    double link_trace;
    double plaquette;
    uint32_t checksum;
    unsigned int sequence_number;
    std::string data_type;
    std::string ensemble_id ;
    std::string ensemble_label ;
    std::string creator ;
    std::string creator_hardware ;
    std::string creation_date ;
    std::string archive_date ;
    std::string floating_point;
};

//////////////////////////////////////////////////////////////////////
// Bit and Physical Checksumming and QA of data
//////////////////////////////////////////////////////////////////////

inline void NerscGrid(GridBase *grid,NerscField &header)
{
  assert(grid->_ndimension==4);
  for(int d=0;d<4;d++) {
    header.dimension[d] = grid->_fdimensions[d];
  }
  for(int d=0;d<4;d++) {
    header.boundary[d] = std::string("PERIODIC");
  }
}
template<class GaugeField>
inline void NerscStatistics(GaugeField & data,NerscField &header)
{
  // How to convert data precision etc...
  header.link_trace=Grid::QCD::WilsonLoops<PeriodicGimplR>::linkTrace(data);
  header.plaquette =Grid::QCD::WilsonLoops<PeriodicGimplR>::avgPlaquette(data);
}

inline void NerscMachineCharacteristics(NerscField &header)
{
  // Who
  struct passwd *pw = getpwuid (getuid());
  if (pw) header.creator = std::string(pw->pw_name); 

  // When
  std::time_t t = std::time(nullptr);
  std::tm tm = *std::localtime(&t);
  std::ostringstream oss; 
  //  oss << std::put_time(&tm, "%c %Z");
  header.creation_date = oss.str();
  header.archive_date  = header.creation_date;

  // What
  struct utsname name;  uname(&name);
  header.creator_hardware = std::string(name.nodename)+"-";
  header.creator_hardware+= std::string(name.machine)+"-";
  header.creator_hardware+= std::string(name.sysname)+"-";
  header.creator_hardware+= std::string(name.release);

}
//////////////////////////////////////////////////////////////////////
// Utilities ; these are QCD aware
//////////////////////////////////////////////////////////////////////
    inline void NerscChecksum(uint32_t *buf,uint32_t buf_size_bytes,uint32_t &csum)
    {
      BinaryIO::Uint32Checksum(buf,buf_size_bytes,csum);
    }
    inline void reconstruct3(LorentzColourMatrix & cm)
    {
      const int x=0;
      const int y=1;
      const int z=2;
      for(int mu=0;mu<4;mu++){
	cm(mu)()(2,x) = adj(cm(mu)()(0,y)*cm(mu)()(1,z)-cm(mu)()(0,z)*cm(mu)()(1,y)); //x= yz-zy
	cm(mu)()(2,y) = adj(cm(mu)()(0,z)*cm(mu)()(1,x)-cm(mu)()(0,x)*cm(mu)()(1,z)); //y= zx-xz
	cm(mu)()(2,z) = adj(cm(mu)()(0,x)*cm(mu)()(1,y)-cm(mu)()(0,y)*cm(mu)()(1,x)); //z= xy-yx
      }
    }

    template<class fobj,class sobj>
    struct NerscSimpleMunger{

      void operator() (fobj &in,sobj &out,uint32_t &csum){

      for(int mu=0;mu<4;mu++){
      for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	out(mu)()(i,j) = in(mu)()(i,j);
      }}}
      NerscChecksum((uint32_t *)&in,sizeof(in),csum); 
      };
    };

    template<class fobj,class sobj>
    struct NerscSimpleUnmunger{
      void operator() (sobj &in,fobj &out,uint32_t &csum){
	for(int mu=0;mu<Nd;mu++){
	for(int i=0;i<Nc;i++){
	for(int j=0;j<Nc;j++){
	  out(mu)()(i,j) = in(mu)()(i,j);
	}}}
	NerscChecksum((uint32_t *)&out,sizeof(out),csum); 
      };
    };
 
    template<class fobj,class sobj>
    struct Nersc3x2munger{
      void operator() (fobj &in,sobj &out,uint32_t &csum){
     
	NerscChecksum((uint32_t *)&in,sizeof(in),csum); 

	for(int mu=0;mu<4;mu++){
	  for(int i=0;i<2;i++){
	    for(int j=0;j<3;j++){
	      out(mu)()(i,j) = in(mu)(i)(j);
	    }}
	}
	reconstruct3(out);
      }
    };

    template<class fobj,class sobj>
    struct Nersc3x2unmunger{

      void operator() (sobj &in,fobj &out,uint32_t &csum){


	for(int mu=0;mu<4;mu++){
	  for(int i=0;i<2;i++){
	    for(int j=0;j<3;j++){
	      out(mu)(i)(j) = in(mu)()(i,j);
	    }}
	}

	NerscChecksum((uint32_t *)&out,sizeof(out),csum); 

      }
    };


////////////////////////////////////////////////////////////////////////////////
// Write and read from fstream; comput header offset for payload
////////////////////////////////////////////////////////////////////////////////
class NerscIO : public BinaryIO { 
 public:

  static inline void truncate(std::string file){
    std::ofstream fout(file,std::ios::out);
  }
  static inline unsigned int writeHeader(NerscField &field,std::string file)
  {
    std::ofstream fout(file,std::ios::out|std::ios::in);
  
    fout.seekp(0,std::ios::beg);
    fout << "BEGIN_HEADER"      << std::endl;
    fout << "HDR_VERSION = "    << field.hdr_version    << std::endl;
    fout << "DATATYPE = "       << field.data_type      << std::endl;
    fout << "STORAGE_FORMAT = " << field.storage_format << std::endl;

    for(int i=0;i<4;i++){
      fout << "DIMENSION_" << i+1 << " = " << field.dimension[i] << std::endl ;
    }
    // just to keep the space and write it later
    fout << "LINK_TRACE = " << std::setprecision(10) << field.link_trace << std::endl;
    fout << "PLAQUETTE  = " << std::setprecision(10) << field.plaquette  << std::endl;
    for(int i=0;i<4;i++){
      fout << "BOUNDARY_"<<i+1<<" = " << field.boundary[i] << std::endl;
    }

    fout << "CHECKSUM = "<< std::hex << std::setw(10) << field.checksum << std::dec<<std::endl;

    fout << "ENSEMBLE_ID = "     << field.ensemble_id      << std::endl;
    fout << "ENSEMBLE_LABEL = "  << field.ensemble_label   << std::endl;
    fout << "SEQUENCE_NUMBER = " << field.sequence_number  << std::endl;
    fout << "CREATOR = "         << field.creator          << std::endl;
    fout << "CREATOR_HARDWARE = "<< field.creator_hardware << std::endl;
    fout << "CREATION_DATE = "   << field.creation_date    << std::endl;
    fout << "ARCHIVE_DATE = "    << field.archive_date     << std::endl;
    fout << "FLOATING_POINT = "  << field.floating_point   << std::endl;
    fout << "END_HEADER"         << std::endl;
    field.data_start = fout.tellp();
    return field.data_start;
}

// for the header-reader
static inline int readHeader(std::string file,GridBase *grid,  NerscField &field)
{
  int offset=0;
  std::map<std::string,std::string> header;
  std::string line;

  //////////////////////////////////////////////////
  // read the header
  //////////////////////////////////////////////////
  std::ifstream fin(file);

  getline(fin,line); // read one line and insist is 

  removeWhitespace(line);
  assert(line==std::string("BEGIN_HEADER"));

  do {
    getline(fin,line); // read one line
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
static inline void readConfiguration(Lattice<iLorentzColourMatrix<vsimd> > &Umu,NerscField& header,std::string file)
{
  typedef Lattice<iLorentzColourMatrix<vsimd> > GaugeField;

  GridBase *grid = Umu._grid;
  int offset = readHeader(file,Umu._grid,header);

  NerscField clone(header);

  std::string format(header.floating_point);

  int ieee32big = (format == std::string("IEEE32BIG"));
  int ieee32    = (format == std::string("IEEE32"));
  int ieee64big = (format == std::string("IEEE64BIG"));
  int ieee64    = (format == std::string("IEEE64"));

  uint32_t csum;
  // depending on datatype, set up munger;
  // munger is a function of <floating point, Real, data_type>
  if ( header.data_type == std::string("4D_SU3_GAUGE") ) {
    if ( ieee32 || ieee32big ) {
      //      csum=BinaryIO::readObjectSerial<iLorentzColourMatrix<vsimd>, LorentzColour2x3F> 
      csum=BinaryIO::readObjectParallel<iLorentzColourMatrix<vsimd>, LorentzColour2x3F> 
	(Umu,file,Nersc3x2munger<LorentzColour2x3F,LorentzColourMatrix>(), offset,format);
    }
    if ( ieee64 || ieee64big ) {
      //      csum=BinaryIO::readObjectSerial<iLorentzColourMatrix<vsimd>, LorentzColour2x3D> 
      csum=BinaryIO::readObjectParallel<iLorentzColourMatrix<vsimd>, LorentzColour2x3D> 
	(Umu,file,Nersc3x2munger<LorentzColour2x3D,LorentzColourMatrix>(),offset,format);
    }
  } else if ( header.data_type == std::string("4D_SU3_GAUGE_3X3") ) {
    if ( ieee32 || ieee32big ) {
      //      csum=BinaryIO::readObjectSerial<iLorentzColourMatrix<vsimd>,LorentzColourMatrixF>
      csum=BinaryIO::readObjectParallel<iLorentzColourMatrix<vsimd>,LorentzColourMatrixF>
	(Umu,file,NerscSimpleMunger<LorentzColourMatrixF,LorentzColourMatrix>(),offset,format);
    }
    if ( ieee64 || ieee64big ) {
      //      csum=BinaryIO::readObjectSerial<iLorentzColourMatrix<vsimd>,LorentzColourMatrixD>
      csum=BinaryIO::readObjectParallel<iLorentzColourMatrix<vsimd>,LorentzColourMatrixD>
	(Umu,file,NerscSimpleMunger<LorentzColourMatrixD,LorentzColourMatrix>(),offset,format);
    }
  } else {
    assert(0);
  }

  NerscStatistics<GaugeField>(Umu,clone);

  assert(fabs(clone.plaquette -header.plaquette ) < 1.0e-5 );
  assert(fabs(clone.link_trace-header.link_trace) < 1.0e-6 );
  assert(csum == header.checksum );

  std::cout<<GridLogMessage <<"Read NERSC Configuration "<<file<< " and plaquette, link trace, and checksum agree"<<std::endl;
}

template<class vsimd>
static inline void writeConfiguration(Lattice<iLorentzColourMatrix<vsimd> > &Umu,std::string file, int two_row,int bits32)
{
  typedef Lattice<iLorentzColourMatrix<vsimd> > GaugeField;

  typedef iLorentzColourMatrix<vsimd> vobj;
  typedef typename vobj::scalar_object sobj;

  // Following should become arguments
  NerscField header;
  header.sequence_number = 1;
  header.ensemble_id     = "UKQCD";
  header.ensemble_label  = "DWF";

  typedef LorentzColourMatrixD fobj3D;
  typedef LorentzColour2x3D    fobj2D;
  typedef LorentzColourMatrixF fobj3f;
  typedef LorentzColour2x3F    fobj2f;

  GridBase *grid = Umu._grid;

  NerscGrid(grid,header);
  NerscStatistics<GaugeField>(Umu,header);
  NerscMachineCharacteristics(header);

  uint32_t csum;
  int offset;
  
  truncate(file);

  if ( two_row ) { 

    header.floating_point = std::string("IEEE64BIG");
    header.data_type      = std::string("4D_SU3_GAUGE");
    Nersc3x2unmunger<fobj2D,sobj> munge;
    BinaryIO::Uint32Checksum<vobj,fobj2D>(Umu, munge,header.checksum);
    offset = writeHeader(header,file);
    csum=BinaryIO::writeObjectSerial<vobj,fobj2D>(Umu,file,munge,offset,header.floating_point);

    std::string file1 = file+"para";
    int offset1 = writeHeader(header,file1);
    int csum1=BinaryIO::writeObjectParallel<vobj,fobj2D>(Umu,file1,munge,offset,header.floating_point);

    
    std::cout << GridLogMessage << " TESTING PARALLEL WRITE offsets " << offset1 << " "<< offset << std::endl;
    std::cout << GridLogMessage << " TESTING PARALLEL WRITE csums   " << csum1 << " "<<std::hex<< csum << std::dec<< std::endl;

    assert(offset1==offset);  
    assert(csum1==csum);  

  } else { 
    header.floating_point = std::string("IEEE64BIG");
    header.data_type      = std::string("4D_SU3_GAUGE_3X3");
    NerscSimpleUnmunger<fobj3D,sobj> munge;
    BinaryIO::Uint32Checksum<vobj,fobj3D>(Umu, munge,header.checksum);
    offset = writeHeader(header,file);
    csum=BinaryIO::writeObjectSerial<vobj,fobj3D>(Umu,file,munge,offset,header.floating_point);
  }

  std::cout<<GridLogMessage <<"Written NERSC Configuration "<<file<< " checksum "<<std::hex<<csum<< std::dec<<" plaq "<< header.plaquette <<std::endl;

 }


    ///////////////////////////////
    // RNG state
    ///////////////////////////////
static inline void writeRNGState(GridSerialRNG &serial,GridParallelRNG &parallel,std::string file)
{
  typedef typename GridParallelRNG::RngStateType RngStateType;

  // Following should become arguments
  NerscField header;
  header.sequence_number = 1;
  header.ensemble_id     = "UKQCD";
  header.ensemble_label  = "DWF";

  GridBase *grid = parallel._grid;

  NerscGrid(grid,header);
  header.link_trace=0.0;
  header.plaquette=0.0;
  NerscMachineCharacteristics(header);

  uint32_t csum;
  int offset;
  
#ifdef RNG_RANLUX
    header.floating_point = std::string("UINT64");
    header.data_type      = std::string("RANLUX48");
#else
    header.floating_point = std::string("UINT32");
    header.data_type      = std::string("MT19937");
#endif

  truncate(file);
  offset = writeHeader(header,file);
  csum=BinaryIO::writeRNGSerial(serial,parallel,file,offset);
  header.checksum = csum;
  offset = writeHeader(header,file);

  std::cout<<GridLogMessage <<"Written NERSC RNG STATE "<<file<< " checksum "<<std::hex<<csum<<std::dec<<std::endl;

 }
    
static inline void readRNGState(GridSerialRNG &serial,GridParallelRNG & parallel,NerscField& header,std::string file)
{
  typedef typename GridParallelRNG::RngStateType RngStateType;

  GridBase *grid = parallel._grid;

  int offset = readHeader(file,grid,header);

  NerscField clone(header);

  std::string format(header.floating_point);
  std::string data_type(header.data_type);

#ifdef RNG_RANLUX
  assert(format == std::string("UINT64"));
  assert(data_type == std::string("RANLUX48"));
#else
  assert(format == std::string("UINT32"));
  assert(data_type == std::string("MT19937"));
#endif

  // depending on datatype, set up munger;
  // munger is a function of <floating point, Real, data_type>
  uint32_t csum=BinaryIO::readRNGSerial(serial,parallel,file,offset);

  assert(csum == header.checksum );

  std::cout<<GridLogMessage <<"Read NERSC RNG file "<<file<< " format "<< data_type <<std::endl;
}

};


}}
#endif
