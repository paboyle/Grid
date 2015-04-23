#ifndef GRID_NERSC_IO_H
#define GRID_NERSC_IO_H

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>

#ifdef HAVE_ENDIAN_H
#include <endian.h>
#include <arpa/inet.h>
#define ntohll be64toh
#else 
#include <arpa/inet.h>
#endif

namespace Grid {

  using namespace QCD;

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


////////////////////////////////////////////////////////////////////////////////
// Write and read from fstream; comput header offset for payload
////////////////////////////////////////////////////////////////////////////////
inline unsigned int writeNerscHeader(NerscField &field,std::string file)
{
  std::ofstream fout(file,std::ios::out);
  
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
  fout << "CHECKSUM = "<< std::hex << std::setw(16) << 0 << field.checksum << std::endl;

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

// A little helper
inline void removeWhitespace(std::string &key)
{
  key.erase(std::remove_if(key.begin(), key.end(), ::isspace),key.end());
}
// for the header-reader
inline int readNerscHeader(std::string file,GridBase *grid,  NerscField &field)
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

  field.hdr_version    = header[std::string("HDR_VERSION")];
  field.data_type      = header[std::string("DATATYPE")];
  field.storage_format = header[std::string("STORAGE_FORMAT")];
  
  field.dimension[0] = std::stol(header[std::string("DIMENSION_1")]);
  field.dimension[1] = std::stol(header[std::string("DIMENSION_2")]);
  field.dimension[2] = std::stol(header[std::string("DIMENSION_3")]);
  field.dimension[3] = std::stol(header[std::string("DIMENSION_4")]);

  assert(grid->_ndimension == 4);
  for(int d=0;d<4;d++){
    assert(grid->_fdimensions[d]==field.dimension[d]);
  }

  field.link_trace = std::stod(header[std::string("LINK_TRACE")]);
  field.plaquette  = std::stod(header[std::string("PLAQUETTE")]);

  field.boundary[0] = header[std::string("BOUNDARY_1")];
  field.boundary[1] = header[std::string("BOUNDARY_2")];
  field.boundary[2] = header[std::string("BOUNDARY_3")];
  field.boundary[3] = header[std::string("BOUNDARY_4")];

  field.checksum = std::stoul(header[std::string("CHECKSUM")],0,16);
  field.ensemble_id      = header[std::string("ENSEMBLE_ID")];
  field.ensemble_label   = header[std::string("ENSEMBLE_LABEL")];
  field.sequence_number  = std::stol(header[std::string("SEQUENCE_NUMBER")]);
  field.creator          = header[std::string("CREATOR")];
  field.creator_hardware = header[std::string("CREATOR_HARDWARE")];
  field.creation_date    = header[std::string("CREATION_DATE")];
  field.archive_date     = header[std::string("ARCHIVE_DATE")];
  field.floating_point   = header[std::string("FLOATING_POINT")];

  return field.data_start;
}


//////////////////////////////////////////////////////////////////////
// Utilities
//////////////////////////////////////////////////////////////////////
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


 void inline be32toh_v(void *file_object,uint32_t bytes)
 {
   uint32_t * f = (uint32_t *)file_object;
   for(int i=0;i*sizeof(uint32_t)<bytes;i++){  
     f[i] = ntohl(f[i]);
   }
 }
 void inline le32toh_v(void *file_object,uint32_t bytes)
 {
   uint32_t *fp = (uint32_t *)file_object;

   uint32_t f;

   for(int i=0;i*sizeof(uint32_t)<bytes;i++){  
     f = fp[i];
     // got network order and the network to host
     f = ((f&0xFF)<<24) | ((f&0xFF00)<<8) | ((f&0xFF0000)>>8) | ((f&0xFF000000UL)>>24) ; 
     fp[i] = ntohl(f);
   }
 }
 void inline be64toh_v(void *file_object,uint32_t bytes)
 {
   uint64_t * f = (uint64_t *)file_object;
   for(int i=0;i*sizeof(uint64_t)<bytes;i++){  
     f[i] = ntohll(f[i]);
   }
 }
 void inline le64toh_v(void *file_object,uint32_t bytes)
 {
   uint64_t *fp = (uint64_t *)file_object;
   uint64_t f,g;

   for(int i=0;i*sizeof(uint64_t)<bytes;i++){  
     f = fp[i];
     // got network order and the network to host
     g = ((f&0xFF)<<24) | ((f&0xFF00)<<8) | ((f&0xFF0000)>>8) | ((f&0xFF000000UL)>>24) ; 
     g = g << 32;
     f = f >> 32;
     g|= ((f&0xFF)<<24) | ((f&0xFF00)<<8) | ((f&0xFF0000)>>8) | ((f&0xFF000000UL)>>24) ; 
     fp[i] = ntohl(g);
   }
 }

inline void NerscChecksum(uint32_t *buf,uint32_t buf_size,uint32_t &csum)
{
  for(int i=0;i*sizeof(uint32_t)<buf_size;i++){
    csum=csum+buf[i];
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
      for(int mu=0;mu<4;mu++){
      for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
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

    NerscChecksum((uint32_t *)&out,sizeof(out),csum); 

    for(int mu=0;mu<4;mu++){
      for(int i=0;i<2;i++){
      for(int j=0;j<3;j++){
	out(mu)(i)(j) = in(mu)()(i,j);
      }}
    }
  }
};

////////////////////////////////////////////////////////////////////////////
// Template wizardry to map types to strings for NERSC in an extensible way
////////////////////////////////////////////////////////////////////////////
 template<class vobj> struct NerscDataType {
   static void DataType     (std::string &str) { str = std::string("4D_BINARY_UNKNOWN"); };
   static void FloatingPoint(std::string &str) { str = std::string("IEEE64BIG"); };
 };

 template<> struct NerscDataType<iColourMatrix<ComplexD> > {
   static void DataType     (std::string &str) { str = std::string("4D_SU3_GAUGE_3X3"); };
   static void FloatingPoint(std::string &str) { str = std::string("IEEE64BIG");};
 };

 template<> struct NerscDataType<iColourMatrix<ComplexF> > {
   static void DataType     (std::string &str) { str = std::string("4D_SU3_GAUGE_3X3"); };
   static void FloatingPoint(std::string &str) { str = std::string("IEEE32BIG");};
 };

//////////////////////////////////////////////////////////////////////
// Bit and Physical Checksumming and QA of data
//////////////////////////////////////////////////////////////////////
/*
template<class vobj> inline uint32_t NerscChecksum(Lattice<vobj> & data)
{
  uint32_t sum;
  for(int ss=0;ss<data._grid->Osites();ss++){
    uint32_t *iptr = (uint32_t *)& data._odata[0] ;
    for(int i=0;i<sizeof(vobj);i+=sizeof(uint32_t)){
      sum=sum+iptr[i];
    }
  }
  data._grid->globalSum(sum);
  return sum;
}
*/
template<class vobj> inline void NerscPhysicalCharacteristics(Lattice<vobj> & data,NerscField &header)
{
  header.data_type      = NerscDataType<vobj>::DataType;
  header.floating_point = NerscDataType<vobj>::FloatingPoint;
  return;
}

 template<> inline void NerscPhysicalCharacteristics(LatticeGaugeField & data,NerscField &header)
{
  NerscDataType<decltype(data._odata[0])>::DataType(header.data_type);
  NerscDataType<decltype(data._odata[0])>::FloatingPoint(header.floating_point);
  header.link_trace=1.0;
  header.plaquette =1.0;
}

template<class vobj> inline void NerscStatisics(Lattice<vobj> & data,NerscField &header)
{
  assert(data._grid->_ndimension==4);

  for(int d=0;d<4;d++)
    header.dimension[d] = data._grid->_fdimensions[d];

  // compute checksum and any physical properties contained for this type
  //  header.checksum = NerscChecksum(data);

  NerscPhysicalCharacteristics(data,header);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Now the meat: the object readers
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class vobj,class sobj,class fobj,class munger>
inline void readNerscObject(Lattice<vobj> &Umu,std::string file,munger munge,int offset,std::string &format)
{
  GridBase *grid = Umu._grid;

  int ieee32big = (format == std::string("IEEE32BIG"));
  int ieee32    = (format == std::string("IEEE32"));
  int ieee64big = (format == std::string("IEEE64BIG"));
  int ieee64    = (format == std::string("IEEE64"));

  // Find the location of each site and send to primary node
  //  for(int site=0; site < Layout::vol(); ++site){
  //     multi1d<int> coord = crtesn(site, Layout::lattSize());
  //     for(int dd=0; dd<Nd; dd++){        /* dir */
  //        cfg_in.readArray(su3_buffer, float_size, mat_size);
  //
  // Above from Chroma; defines loop order now that NERSC doc no longer
  // available (how short sighted is that?)
  {
    std::ifstream fin(file,std::ios::binary|std::ios::in);
    fin.seekg(offset);

    Umu = zero;
    uint32_t csum=0;
    fobj file_object;
    sobj munged;
    
    for(int t=0;t<grid->_fdimensions[3];t++){
    for(int z=0;z<grid->_fdimensions[2];z++){
    for(int y=0;y<grid->_fdimensions[1];y++){
    for(int x=0;x<grid->_fdimensions[0];x++){

      std::vector<int> site({x,y,z,t});

      if ( grid->IsBoss() ) {
	fin.read((char *)&file_object,sizeof(file_object));

	if(ieee32big) be32toh_v((void *)&file_object,sizeof(file_object));
	if(ieee32)    le32toh_v((void *)&file_object,sizeof(file_object));
	if(ieee64big) be64toh_v((void *)&file_object,sizeof(file_object));
	if(ieee64)    le64toh_v((void *)&file_object,sizeof(file_object));

	munge(file_object,munged,csum);
      }
      // The boss who read the file has their value poked
      pokeSite(munged,Umu,site);
    }}}}
  }
}

template<class vobj,class sobj,class fobj,class munger>
inline void writeNerscObject(Lattice<vobj> &Umu,std::string file,munger munge,int offset,
                       int sequence,double lt,double pl)
{
  GridBase *grid = Umu._grid;
  NerscField header;
  
  //////////////////////////////////////////////////
  // First write the header; this is in wrong place
  //////////////////////////////////////////////////
  assert(grid->_ndimension == 4);
  for(int d=0;d<4;d++){
    header.dimension[d]=grid->_fdimensions[d];
    header.boundary [d]=std::string("PERIODIC");; 
  }
  header.hdr_version=std::string("WHATDAHECK");
  //  header.storage_format=storage_format<vobj>::string; // use template specialisation
  //  header.data_type=data_type<vobj>::string;
  header.storage_format=std::string("debug");
  header.data_type     =std::string("debug");

  //FIXME; use template specialisation to fill these out
  header.link_trace   =lt;
  header.plaquette    =pl;
  header.checksum     =0;

  //
  header.sequence_number =sequence;
  header.ensemble_id     =std::string("UKQCD");
  header.ensemble_label  =std::string("UKQCD");
  header.creator         =std::string("Tadahito");
  header.creator_hardware=std::string("BlueGene/Q");
  header.creation_date   =std::string("AnnoDomini");
  header.archive_date    =std::string("AnnoDomini");
  header.floating_point  =std::string("IEEE64BIG");
  //  header.data_start=;
  //  unsigned int checksum;

  //////////////////////////////////////////////////
  // Now write the body
  //////////////////////////////////////////////////
  {
    std::ofstream fout(file,std::ios::binary|std::ios::in);
    fout.seekp(offset);

    Umu = zero;
    uint32_t csum=0;
    fobj file_object;
    sobj unmunged;
    for(int t=0;t<grid->_fdimensions[3];t++){
    for(int z=0;z<grid->_fdimensions[2];z++){
    for(int y=0;y<grid->_fdimensions[1];y++){
    for(int x=0;x<grid->_fdimensions[0];x++){
      std::vector<int> site({x,y,z,t});
      peekSite(unmunged,Umu,site);
      munge(unmunged,file_object,csum);
      // broadcast & insert
      fout.write((char *)&file_object,sizeof(file_object));
    }}}}
  }
}



inline void readNerscConfiguration(LatticeGaugeField &Umu,NerscField& header,std::string file)
{
  GridBase *grid = Umu._grid;

  int offset = readNerscHeader(file,Umu._grid,header);

  std::string format(header.floating_point);

  int ieee32big = (format == std::string("IEEE32BIG"));
  int ieee32    = (format == std::string("IEEE32"));
  int ieee64big = (format == std::string("IEEE64BIG"));
  int ieee64    = (format == std::string("IEEE64"));

  // depending on datatype, set up munger;
  // munger is a function of <floating point, Real, data_type>
  if ( header.data_type == std::string("4D_SU3_GAUGE") ) {
    if ( ieee32 || ieee32big ) {
      readNerscObject<vLorentzColourMatrix, LorentzColourMatrix, LorentzColour2x3F> 
	(Umu,file,
	 Nersc3x2munger<LorentzColour2x3F,LorentzColourMatrix>(),
	 offset,format);
    }
    if ( ieee64 || ieee64big ) {
      readNerscObject<vLorentzColourMatrix, LorentzColourMatrix, LorentzColour2x3D> 
	(Umu,file,
	 Nersc3x2munger<LorentzColour2x3D,LorentzColourMatrix>(),
	 offset,format);
    }
  } else if ( header.data_type == std::string("4D_SU3_GAUGE_3X3") ) {
    if ( ieee32 || ieee32big ) {
      readNerscObject<vLorentzColourMatrix,LorentzColourMatrix,LorentzColourMatrixF>
	(Umu,file,NerscSimpleMunger<LorentzColourMatrixF,LorentzColourMatrix>(),offset,format);
    }
    if ( ieee64 || ieee64big ) {
      readNerscObject<vLorentzColourMatrix,LorentzColourMatrix,LorentzColourMatrixD>
	(Umu,file,NerscSimpleMunger<LorentzColourMatrixD,LorentzColourMatrix>(),offset,format);
    }
  } else {
    assert(0);
  }

}

template<class vobj>
inline void writeNerscConfiguration(Lattice<vobj> &Umu,NerscField &header,std::string file)
{
  GridBase &grid = Umu._grid;
  
  NerscStatisics(Umu,header);

  int offset = writeNerscHeader(header,file);

  writeNerscObject(Umu,NerscSimpleMunger<vobj,vobj>(),offset);
}


}
#endif
