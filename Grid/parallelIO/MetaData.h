/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/parallelIO/NerscIO.h

    Copyright (C) 2015


    Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <unistd.h>
#include <sys/utsname.h>
#include <pwd.h>

namespace Grid {

  ///////////////////////////////////////////////////////
  // Precision mapping
  ///////////////////////////////////////////////////////
  template<class vobj> static std::string getFormatString (void)
  {
    std::string format;
    typedef typename getPrecision<vobj>::real_scalar_type stype;
    if ( sizeof(stype) == sizeof(float) ) {
      format = std::string("IEEE32BIG");
    }
    if ( sizeof(stype) == sizeof(double) ) {
      format = std::string("IEEE64BIG");
    }
    return format;
  }
  ////////////////////////////////////////////////////////////////////////////////
  // header specification/interpretation
  ////////////////////////////////////////////////////////////////////////////////
    class FieldNormMetaData : Serializable {
    public:
      GRID_SERIALIZABLE_CLASS_MEMBERS(FieldNormMetaData, double, norm2);
    };
    class FieldMetaData : Serializable {
    public:

      GRID_SERIALIZABLE_CLASS_MEMBERS(FieldMetaData,
				      int, nd,
				      std::vector<int>, dimension,
				      std::vector<std::string>, boundary,
				      int, data_start,
				      std::string, hdr_version,
				      std::string, storage_format,
				      double, link_trace,
				      double, plaquette,
				      uint32_t, checksum,
				      uint32_t, scidac_checksuma,
				      uint32_t, scidac_checksumb,
				      unsigned int, sequence_number,
				      std::string, data_type,
				      std::string, ensemble_id,
				      std::string, ensemble_label,
				      std::string, ildg_lfn,
				      std::string, creator,
				      std::string, creator_hardware,
				      std::string, creation_date,
				      std::string, archive_date,
				      std::string, floating_point);
      // WARNING: non-initialised values might lead to twisted parallel IO
      // issues, std::string are fine because they initliase to size 0
      // as per C++ standard.
      FieldMetaData(void) 
      : nd(4), dimension(4,0), boundary(4, ""), data_start(0),
      link_trace(0.), plaquette(0.), checksum(0),
      scidac_checksuma(0), scidac_checksumb(0), sequence_number(0)
      {}
    };

  namespace QCD {

    using namespace Grid;


    //////////////////////////////////////////////////////////////////////
    // Bit and Physical Checksumming and QA of data
    //////////////////////////////////////////////////////////////////////
    inline void GridMetaData(GridBase *grid,FieldMetaData &header)
    {
      int nd = grid->_ndimension;
      header.nd = nd;
      header.dimension.resize(nd);
      header.boundary.resize(nd);
      header.data_start = 0;
      for(int d=0;d<nd;d++) {
	header.dimension[d] = grid->_fdimensions[d];
      }
      for(int d=0;d<nd;d++) {
	header.boundary[d] = std::string("PERIODIC");
      }
    }

    inline void MachineCharacteristics(FieldMetaData &header)
    {
      // Who
      struct passwd *pw = getpwuid (getuid());
      if (pw) header.creator = std::string(pw->pw_name); 

      // When
      std::time_t t = std::time(nullptr);
      std::tm tm_ = *std::localtime(&t);
      std::ostringstream oss; 
      //      oss << std::put_time(&tm_, "%c %Z");
      header.creation_date = oss.str();
      header.archive_date  = header.creation_date;

      // What
      struct utsname name;  uname(&name);
      header.creator_hardware = std::string(name.nodename)+"-";
      header.creator_hardware+= std::string(name.machine)+"-";
      header.creator_hardware+= std::string(name.sysname)+"-";
      header.creator_hardware+= std::string(name.release);
    }

#define dump_meta_data(field, s)					\
      s << "BEGIN_HEADER"      << std::endl;				\
      s << "HDR_VERSION = "    << field.hdr_version    << std::endl;	\
      s << "DATATYPE = "       << field.data_type      << std::endl;	\
      s << "STORAGE_FORMAT = " << field.storage_format << std::endl;	\
      for(int i=0;i<4;i++){						\
	s << "DIMENSION_" << i+1 << " = " << field.dimension[i] << std::endl ; \
      }									\
      s << "LINK_TRACE = " << std::setprecision(10) << field.link_trace << std::endl; \
      s << "PLAQUETTE  = " << std::setprecision(10) << field.plaquette  << std::endl; \
      for(int i=0;i<4;i++){						\
	s << "BOUNDARY_"<<i+1<<" = " << field.boundary[i] << std::endl;	\
      }									\
									\
      s << "CHECKSUM = "<< std::hex << std::setw(10) << field.checksum << std::dec<<std::endl; \
      s << "SCIDAC_CHECKSUMA = "<< std::hex << std::setw(10) << field.scidac_checksuma << std::dec<<std::endl; \
      s << "SCIDAC_CHECKSUMB = "<< std::hex << std::setw(10) << field.scidac_checksumb << std::dec<<std::endl; \
      s << "ENSEMBLE_ID = "     << field.ensemble_id      << std::endl;	\
      s << "ENSEMBLE_LABEL = "  << field.ensemble_label   << std::endl;	\
      s << "SEQUENCE_NUMBER = " << field.sequence_number  << std::endl;	\
      s << "CREATOR = "         << field.creator          << std::endl;	\
      s << "CREATOR_HARDWARE = "<< field.creator_hardware << std::endl;	\
      s << "CREATION_DATE = "   << field.creation_date    << std::endl;	\
      s << "ARCHIVE_DATE = "    << field.archive_date     << std::endl;	\
      s << "FLOATING_POINT = "  << field.floating_point   << std::endl;	\
      s << "END_HEADER"         << std::endl;

template<class vobj> inline void PrepareMetaData(Lattice<vobj> & field, FieldMetaData &header)
{
  GridBase *grid = field._grid;
  std::string format = getFormatString<vobj>();
   header.floating_point = format;
   header.checksum = 0x0; // Nersc checksum unused in ILDG, Scidac
   GridMetaData(grid,header); 
   MachineCharacteristics(header);
 }
 inline void GaugeStatistics(Lattice<vLorentzColourMatrixF> & data,FieldMetaData &header)
 {
   // How to convert data precision etc...
   header.link_trace=Grid::QCD::WilsonLoops<PeriodicGimplF>::linkTrace(data);
   header.plaquette =Grid::QCD::WilsonLoops<PeriodicGimplF>::avgPlaquette(data);
 }
 inline void GaugeStatistics(Lattice<vLorentzColourMatrixD> & data,FieldMetaData &header)
 {
   // How to convert data precision etc...
   header.link_trace=Grid::QCD::WilsonLoops<PeriodicGimplD>::linkTrace(data);
   header.plaquette =Grid::QCD::WilsonLoops<PeriodicGimplD>::avgPlaquette(data);
 }
 template<> inline void PrepareMetaData<vLorentzColourMatrixF>(Lattice<vLorentzColourMatrixF> & field, FieldMetaData &header)
 {
   
   GridBase *grid = field._grid;
   std::string format = getFormatString<vLorentzColourMatrixF>();
   header.floating_point = format;
   header.checksum = 0x0; // Nersc checksum unused in ILDG, Scidac
   GridMetaData(grid,header); 
   GaugeStatistics(field,header);
   MachineCharacteristics(header);
 }
 template<> inline void PrepareMetaData<vLorentzColourMatrixD>(Lattice<vLorentzColourMatrixD> & field, FieldMetaData &header)
 {
   GridBase *grid = field._grid;
   std::string format = getFormatString<vLorentzColourMatrixD>();
   header.floating_point = format;
   header.checksum = 0x0; // Nersc checksum unused in ILDG, Scidac
   GridMetaData(grid,header); 
   GaugeStatistics(field,header);
   MachineCharacteristics(header);
 }

    //////////////////////////////////////////////////////////////////////
    // Utilities ; these are QCD aware
    //////////////////////////////////////////////////////////////////////
    inline void reconstruct3(LorentzColourMatrix & cm)
    {
      const int x=0;
      const int y=1;
      const int z=2;
      for(int mu=0;mu<Nd;mu++){
	cm(mu)()(2,x) = adj(cm(mu)()(0,y)*cm(mu)()(1,z)-cm(mu)()(0,z)*cm(mu)()(1,y)); //x= yz-zy
	cm(mu)()(2,y) = adj(cm(mu)()(0,z)*cm(mu)()(1,x)-cm(mu)()(0,x)*cm(mu)()(1,z)); //y= zx-xz
	cm(mu)()(2,z) = adj(cm(mu)()(0,x)*cm(mu)()(1,y)-cm(mu)()(0,y)*cm(mu)()(1,x)); //z= xy-yx
      }
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Some data types for intermediate storage
    ////////////////////////////////////////////////////////////////////////////////
    template<typename vtype> using iLorentzColour2x3 = iVector<iVector<iVector<vtype, Nc>, 2>, Nd >;

    typedef iLorentzColour2x3<Complex>  LorentzColour2x3;
    typedef iLorentzColour2x3<ComplexF> LorentzColour2x3F;
    typedef iLorentzColour2x3<ComplexD> LorentzColour2x3D;

/////////////////////////////////////////////////////////////////////////////////
// Simple classes for precision conversion
/////////////////////////////////////////////////////////////////////////////////
template <class fobj, class sobj>
struct BinarySimpleUnmunger {
  typedef typename getPrecision<fobj>::real_scalar_type fobj_stype;
  typedef typename getPrecision<sobj>::real_scalar_type sobj_stype;
  
  void operator()(sobj &in, fobj &out) {
    // take word by word and transform accoding to the status
    fobj_stype *out_buffer = (fobj_stype *)&out;
    sobj_stype *in_buffer = (sobj_stype *)&in;
    size_t fobj_words = sizeof(out) / sizeof(fobj_stype);
    size_t sobj_words = sizeof(in) / sizeof(sobj_stype);
    assert(fobj_words == sobj_words);
    
    for (unsigned int word = 0; word < sobj_words; word++)
      out_buffer[word] = in_buffer[word];  // type conversion on the fly
    
  }
};

template <class fobj, class sobj>
struct BinarySimpleMunger {
  typedef typename getPrecision<fobj>::real_scalar_type fobj_stype;
  typedef typename getPrecision<sobj>::real_scalar_type sobj_stype;

  void operator()(fobj &in, sobj &out) {
    // take word by word and transform accoding to the status
    fobj_stype *in_buffer = (fobj_stype *)&in;
    sobj_stype *out_buffer = (sobj_stype *)&out;
    size_t fobj_words = sizeof(in) / sizeof(fobj_stype);
    size_t sobj_words = sizeof(out) / sizeof(sobj_stype);
    assert(fobj_words == sobj_words);
    
    for (unsigned int word = 0; word < sobj_words; word++)
      out_buffer[word] = in_buffer[word];  // type conversion on the fly
    
  }
};


    template<class fobj,class sobj>
    struct GaugeSimpleMunger{
      void operator()(fobj &in, sobj &out) {
        for (int mu = 0; mu < Nd; mu++) {
          for (int i = 0; i < Nc; i++) {
          for (int j = 0; j < Nc; j++) {
	    out(mu)()(i, j) = in(mu)()(i, j);
	  }}
        }
      };
    };

    template <class fobj, class sobj>
    struct GaugeSimpleUnmunger {

      void operator()(sobj &in, fobj &out) {
        for (int mu = 0; mu < Nd; mu++) {
          for (int i = 0; i < Nc; i++) {
          for (int j = 0; j < Nc; j++) {
	    out(mu)()(i, j) = in(mu)()(i, j);
	  }}
        }
      };
    };

    template<class fobj,class sobj>
    struct Gauge3x2munger{
      void operator() (fobj &in,sobj &out){
	for(int mu=0;mu<Nd;mu++){
	  for(int i=0;i<2;i++){
	  for(int j=0;j<3;j++){
	    out(mu)()(i,j) = in(mu)(i)(j);
	  }}
	}
	reconstruct3(out);
      }
    };

    template<class fobj,class sobj>
    struct Gauge3x2unmunger{
      void operator() (sobj &in,fobj &out){
	for(int mu=0;mu<Nd;mu++){
	  for(int i=0;i<2;i++){
	  for(int j=0;j<3;j++){
	    out(mu)(i)(j) = in(mu)()(i,j);
	  }}
	}
      }
    };
  }


}
