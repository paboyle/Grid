    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/parallelIO/BinaryIO.h

    Copyright (C) 2015

    Author: Peter Boyle <paboyle@ph.ed.ac.uk>
    Author: Guido Cossu<guido.cossu@ed.ac.uk>

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
#ifndef GRID_BINARY_IO_H
#define GRID_BINARY_IO_H


#include "IldgIOtypes.h"

#ifdef HAVE_ENDIAN_H
#include <endian.h>
#endif
#include <arpa/inet.h>
#include <algorithm>


inline uint32_t byte_reverse32(uint32_t f) { 
      f = ((f&0xFF)<<24) | ((f&0xFF00)<<8) | ((f&0xFF0000)>>8) | ((f&0xFF000000UL)>>24) ; 
      return f;
}
inline uint64_t byte_reverse64(uint64_t f) { 
  uint64_t g;
  g = ((f&0xFF)<<24) | ((f&0xFF00)<<8) | ((f&0xFF0000)>>8) | ((f&0xFF000000UL)>>24) ; 
  g = g << 32;
  f = f >> 32;
  g|= ((f&0xFF)<<24) | ((f&0xFF00)<<8) | ((f&0xFF0000)>>8) | ((f&0xFF000000UL)>>24) ; 
  return g;
}

#if BYTE_ORDER == BIG_ENDIAN 
inline uint64_t Grid_ntohll(uint64_t A) { return A; }
#else
inline uint64_t Grid_ntohll(uint64_t A) { 
  return byte_reverse64(A);
}
#endif

namespace Grid { 

  // A little helper
  inline void removeWhitespace(std::string &key)
  {
    key.erase(std::remove_if(key.begin(), key.end(), ::isspace),key.end());
  }

class BinaryIO {

 public:


  // Network is big endian
  static inline void htobe32_v(void *file_object,uint32_t bytes){ be32toh_v(file_object,bytes);} 
  static inline void htobe64_v(void *file_object,uint32_t bytes){ be64toh_v(file_object,bytes);} 
  static inline void htole32_v(void *file_object,uint32_t bytes){ le32toh_v(file_object,bytes);} 
  static inline void htole64_v(void *file_object,uint32_t bytes){ le64toh_v(file_object,bytes);} 

  static inline void be32toh_v(void *file_object,uint32_t bytes)
  {
    uint32_t * f = (uint32_t *)file_object;
    for(int i=0;i*sizeof(uint32_t)<bytes;i++){  
      f[i] = ntohl(f[i]);
    }
  }

  // LE must Swap and switch to host
  static inline void le32toh_v(void *file_object,uint32_t bytes)
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

  // BE is same as network
  static inline void be64toh_v(void *file_object,uint32_t bytes)
  {
    uint64_t * f = (uint64_t *)file_object;
    for(int i=0;i*sizeof(uint64_t)<bytes;i++){  
      f[i] = Grid_ntohll(f[i]);
    }
  }
  
  // LE must swap and switch;
  static inline void le64toh_v(void *file_object,uint32_t bytes)
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
      fp[i] = Grid_ntohll(g);
    }
  }

  template<class vobj,class fobj,class munger> static inline void Uint32Checksum(Lattice<vobj> &lat,munger munge,uint32_t &csum)
  {
    typedef typename vobj::scalar_object sobj;
    GridBase *grid = lat._grid ;
    std::cout <<GridLogMessage<< "Uint32Checksum "<<norm2(lat)<<std::endl;
    sobj siteObj;
    fobj fileObj;

    csum = 0;
    std::vector<int> lcoor;
    for(int l=0;l<grid->lSites();l++){
      Lexicographic::CoorFromIndex(lcoor,l,grid->_ldimensions);
      peekLocalSite(siteObj,lat,lcoor);
      munge(siteObj,fileObj,csum);
    }
    grid->GlobalSum(csum);
  }
    
  static inline void Uint32Checksum(uint32_t *buf,uint32_t buf_size_bytes,uint32_t &csum)
  {
    for(int i=0;i*sizeof(uint32_t)<buf_size_bytes;i++){
      csum=csum+buf[i];
    }
  }

  // Simple classes for precision conversion
  template <class fobj, class sobj>
  struct BinarySimpleUnmunger {
    typedef typename getPrecision<fobj>::real_scalar_type fobj_stype;
    typedef typename getPrecision<sobj>::real_scalar_type sobj_stype;

    void operator()(sobj &in, fobj &out, uint32_t &csum) {
      // take word by word and transform accoding to the status
      fobj_stype *out_buffer = (fobj_stype *)&out;
      sobj_stype *in_buffer = (sobj_stype *)&in;
      size_t fobj_words = sizeof(out) / sizeof(fobj_stype);
      size_t sobj_words = sizeof(in) / sizeof(sobj_stype);
      assert(fobj_words == sobj_words);

      for (unsigned int word = 0; word < sobj_words; word++)
        out_buffer[word] = in_buffer[word];  // type conversion on the fly

      BinaryIO::Uint32Checksum((uint32_t *)&out, sizeof(out), csum);
    }
  };

  template <class fobj, class sobj>
  struct BinarySimpleMunger {
    typedef typename getPrecision<fobj>::real_scalar_type fobj_stype;
    typedef typename getPrecision<sobj>::real_scalar_type sobj_stype;

    void operator()(fobj &in, sobj &out, uint32_t &csum) {
      // take word by word and transform accoding to the status
      fobj_stype *in_buffer = (fobj_stype *)&in;
      sobj_stype *out_buffer = (sobj_stype *)&out;
      size_t fobj_words = sizeof(in) / sizeof(fobj_stype);
      size_t sobj_words = sizeof(out) / sizeof(sobj_stype);
      assert(fobj_words == sobj_words);

      for (unsigned int word = 0; word < sobj_words; word++)
        out_buffer[word] = in_buffer[word];  // type conversion on the fly

      BinaryIO::Uint32Checksum((uint32_t *)&in, sizeof(in), csum);
    }
  };

  template<class vobj,class fobj,class munger>
  static inline uint32_t readObjectSerial(Lattice<vobj> &Umu,std::string file,munger munge,int offset,const std::string &format)
  {
    typedef typename vobj::scalar_object sobj;

    GridBase *grid = Umu._grid;

    std::cout<< GridLogMessage<< "Serial read I/O "<< file<< std::endl;
    GridStopWatch timer; timer.Start();

    int ieee32big = (format == std::string("IEEE32BIG"));
    int ieee32    = (format == std::string("IEEE32"));
    int ieee64big = (format == std::string("IEEE64BIG"));
    int ieee64    = (format == std::string("IEEE64"));

    // Find the location of each site and send to primary node
    // Take loop order from Chroma; defines loop order now that NERSC doc no longer
    // available (how short sighted is that?)
    std::ifstream fin(file,std::ios::binary|std::ios::in);
    fin.seekg(offset);

    Umu = zero;
    uint32_t csum=0;
    uint64_t bytes=0;

    int lx = grid->_fdimensions[0];
    std::vector<fobj> file_object(lx);
    std::vector<sobj> munged(lx);
    for(int t=0;t<grid->_fdimensions[3];t++){
    for(int z=0;z<grid->_fdimensions[2];z++){
    for(int y=0;y<grid->_fdimensions[1];y++){
    {
      bytes += sizeof(fobj)*lx;
      if (grid->IsBoss()) {
        fin.read((char *)&file_object[0], sizeof(fobj)*lx); assert( fin.fail()==0);
	for(int x=0;x<lx;x++){
	  if (ieee32big) be32toh_v((void *)&file_object[x], sizeof(fobj));
	  if (ieee32)    le32toh_v((void *)&file_object[x], sizeof(fobj));
	  if (ieee64big) be64toh_v((void *)&file_object[x], sizeof(fobj));
	  if (ieee64)    le64toh_v((void *)&file_object[x], sizeof(fobj));
	  munge(file_object[x], munged[x], csum);
	}
      }
      for(int x=0;x<lx;x++){
	std::vector<int> site({x,y,z,t});
	// The boss who read the file has their value poked
	pokeSite(munged[x],Umu,site);
      }
    }}}}
    timer.Stop();
    std::cout<<GridLogPerformance<<"readObjectSerial: read "<< bytes <<" bytes in "<<timer.Elapsed() <<" "
	     << (double)bytes/ (double)timer.useconds() <<" MB/s "  <<std::endl;

    grid->Broadcast(0,(void *)&csum,sizeof(csum));
    return csum;
  }

  template<class vobj,class fobj,class munger> 
  static inline uint32_t writeObjectSerial(Lattice<vobj> &Umu,std::string file,munger munge,int offset,
					   const std::string & format)
  {
        typedef typename vobj::scalar_object sobj;

    GridBase *grid = Umu._grid;

    int ieee32big = (format == std::string("IEEE32BIG"));
    int ieee32    = (format == std::string("IEEE32"));
    int ieee64big = (format == std::string("IEEE64BIG"));
    int ieee64    = (format == std::string("IEEE64"));

    //////////////////////////////////////////////////
    // Serialise through node zero
    //////////////////////////////////////////////////
    std::cout<< GridLogMessage<< "Serial write I/O "<< file<<std::endl;
    GridStopWatch timer; timer.Start();
    
    std::ofstream fout;
    if ( grid->IsBoss() ) {
      fout.open(file,std::ios::binary|std::ios::out|std::ios::in);
      fout.seekp(offset);
    }
    uint64_t bytes=0;
    uint32_t csum=0;
    int lx = grid->_fdimensions[0];
    std::vector<fobj> file_object(lx);
    std::vector<sobj> unmunged(lx);
    for(int t=0;t<grid->_fdimensions[3];t++){
    for(int z=0;z<grid->_fdimensions[2];z++){
    for(int y=0;y<grid->_fdimensions[1];y++){
    {

      std::vector<int> site({0,y,z,t});
      // peek & write
      for(int x=0;x<lx;x++){
	site[0]=x;
	peekSite(unmunged[x],Umu,site);
      }
      
      if ( grid->IsBoss() ) {
	for(int x=0;x<lx;x++){
	  munge(unmunged[x],file_object[x],csum);
	  if(ieee32big) htobe32_v((void *)&file_object[x],sizeof(fobj));
	  if(ieee32)    htole32_v((void *)&file_object[x],sizeof(fobj));
	  if(ieee64big) htobe64_v((void *)&file_object[x],sizeof(fobj));
	  if(ieee64)    htole64_v((void *)&file_object[x],sizeof(fobj));
	}
	fout.write((char *)&file_object[0],sizeof(fobj)*lx);assert( fout.fail()==0);
	bytes+=sizeof(fobj)*lx;
      }
    }}}}

    timer.Stop();
    std::cout<<GridLogPerformance<<"writeObjectSerial: wrote "<< bytes <<" bytes in "<<timer.Elapsed() <<" "
	     << (double)bytes/timer.useconds() <<" MB/s "  <<std::endl;
    
    grid->Broadcast(0,(void *)&csum,sizeof(csum));
    return csum;
  }
  
  static inline uint32_t writeRNGSerial(GridSerialRNG &serial,GridParallelRNG &parallel,std::string file,int offset)
  {
    typedef typename GridSerialRNG::RngStateType RngStateType;
    const int RngStateCount = GridSerialRNG::RngStateCount;

    GridBase *grid = parallel._grid;
    int gsites = grid->_gsites;

    GridStopWatch timer; timer.Start();
    //////////////////////////////////////////////////
    // Serialise through node zero
    //////////////////////////////////////////////////
    std::ofstream fout;
    if (grid->IsBoss()) {
      fout.open(file, std::ios::binary | std::ios::out);
      if (!fout.is_open()) {
        std::cout << GridLogMessage << "writeRNGSerial: Error opening file " << file << std::endl;
        exit(0);// write better error handling
      }
      fout.seekp(offset);
    }

    std::cout << GridLogMessage << "Serial RNG write I/O on file " << file << std::endl;
    uint32_t csum = 0;
    std::vector<RngStateType> saved(RngStateCount);
    int bytes = sizeof(RngStateType) * saved.size();
    std::cout << GridLogDebug << "RngStateCount: " << RngStateCount << std::endl;
    std::cout << GridLogDebug << "Type has " << bytes << " bytes" << std::endl;
    std::vector<int> gcoor;

    for(int gidx=0;gidx<gsites;gidx++){

      int rank,o_idx,i_idx;
      grid->GlobalIndexToGlobalCoor(gidx,gcoor);
      grid->GlobalCoorToRankIndex(rank,o_idx,i_idx,gcoor);
      int l_idx=parallel.generator_idx(o_idx,i_idx);

      if( rank == grid->ThisRank() ){
	//  std::cout << "rank" << rank<<" Getting state for index "<<l_idx<<std::endl;
	parallel.GetState(saved,l_idx);
      }
      grid->Broadcast(rank, (void *)&saved[0], bytes);

      if ( grid->IsBoss() ) {
	Uint32Checksum((uint32_t *)&saved[0],bytes,csum);
	fout.write((char *)&saved[0],bytes);assert( fout.fail()==0);
      }
      
    }

    if ( grid->IsBoss() ) {
      serial.GetState(saved,0);
      Uint32Checksum((uint32_t *)&saved[0],bytes,csum);
      fout.write((char *)&saved[0],bytes);assert( fout.fail()==0);
    }

    grid->Broadcast(0, (void *)&csum, sizeof(csum));

    if (grid->IsBoss()) 
      fout.close();

    timer.Stop();

    std::cout << GridLogMessage << "RNG file checksum " << std::hex << csum << std::dec << std::endl;
    std::cout << GridLogMessage << "RNG state saved in " << timer.Elapsed() << std::endl;
    return csum;
  }


  static inline uint32_t readRNGSerial(GridSerialRNG &serial,GridParallelRNG &parallel,std::string file,int offset)
  {
    typedef typename GridSerialRNG::RngStateType RngStateType;
    const int RngStateCount = GridSerialRNG::RngStateCount;

    GridBase *grid = parallel._grid;
    int gsites = grid->_gsites;

    //////////////////////////////////////////////////
    // Serialise through node zero
    //////////////////////////////////////////////////
    std::cout<< GridLogMessage<< "Serial RNG read I/O of file "<<file<<std::endl;

    std::ifstream fin;
    if (grid->IsBoss()) {
      fin.open(file, std::ios::binary | std::ios::in);
      if (!fin.is_open()) {
        std::cout << GridLogMessage << "readRNGSerial: Error opening file " << file << std::endl;
        exit(0);// write better error handling
      }
      fin.seekg(offset);
    }

    
    uint32_t csum=0;
    std::vector<RngStateType> saved(RngStateCount);
    int bytes = sizeof(RngStateType)*saved.size();
    std::cout << GridLogDebug << "RngStateCount: " << RngStateCount << std::endl;
    std::cout << GridLogDebug << "Type has " << bytes << " bytes" << std::endl;
    std::vector<int> gcoor;

    std::cout << GridLogDebug << "gsites: " << gsites << " loop" << std::endl;
    for(int gidx=0;gidx<gsites;gidx++){

      int rank,o_idx,i_idx;
      grid->GlobalIndexToGlobalCoor(gidx,gcoor);
      grid->GlobalCoorToRankIndex(rank,o_idx,i_idx,gcoor);
      int l_idx=parallel.generator_idx(o_idx,i_idx);

      if ( grid->IsBoss() ) {
	fin.read((char *)&saved[0],bytes);assert( fin.fail()==0);
	Uint32Checksum((uint32_t *)&saved[0],bytes,csum);
      }
      
      grid->Broadcast(0,(void *)&saved[0],bytes);

      if( rank == grid->ThisRank() ){
        parallel.SetState(saved,l_idx);
      }
    }

    if ( grid->IsBoss() ) {
      fin.read((char *)&saved[0],bytes);assert( fin.fail()==0);
      serial.SetState(saved,0);
      Uint32Checksum((uint32_t *)&saved[0],bytes,csum);
    }

    std::cout << GridLogMessage << "RNG file checksum " << std::hex << csum << std::dec << std::endl;

    grid->Broadcast(0,(void *)&csum,sizeof(csum));

    return csum;
  }


  template <class vobj, class fobj, class munger>
  static inline uint32_t readObjectParallel(Lattice<vobj> &Umu,
                                            std::string file, 
                                            munger munge,
                                            int offset,
                                            const std::string &format,
                                            ILDGtype ILDG = ILDGtype()) {
    typedef typename vobj::scalar_object sobj;

    GridBase *grid = Umu._grid;

    int ieee32big = (format == std::string("IEEE32BIG"));
    int ieee32    = (format == std::string("IEEE32"));
    int ieee64big = (format == std::string("IEEE64BIG"));
    int ieee64    = (format == std::string("IEEE64"));


    // Ideally one reader/writer per xy plane and read these contiguously
    // with comms from nominated I/O nodes.
    std::ifstream fin;

    int nd = grid->_ndimension;
    std::vector<int> parallel(nd,1); parallel[0] = 0;
    std::vector<int> ioproc  (nd);
    std::vector<int> start(nd);
    std::vector<int> range(nd);

    for(int d=0;d<nd;d++){
      assert(grid->CheckerBoarded(d) == 0);
    }

    uint64_t slice_vol = 1;

    int IOnode = 1;
    int gstrip = grid->_gdimensions[0];
    int lstrip = grid->_ldimensions[0];

    int chunk ;
    if ( nd==1) chunk = gstrip;
    else        chunk = gstrip*grid->_ldimensions[1];

    for(int d=0;d<grid->_ndimension;d++) {
      
      if (parallel[d]) {
	range[d] = grid->_ldimensions[d];
	start[d] = grid->_processor_coor[d]*range[d];
	ioproc[d]= grid->_processor_coor[d];
      } else {
	range[d] = grid->_gdimensions[d];
	start[d] = 0;
	ioproc[d]= 0;
	
	if ( grid->_processor_coor[d] != 0 ) IOnode = 0;
      }
      slice_vol = slice_vol * range[d];
    }

    {
      uint32_t tmp = IOnode;
      grid->GlobalSum(tmp);
      std::cout<< std::dec ;
      std::cout<< GridLogMessage<< "Parallel read I/O from "<< file << " with " <<tmp<< " IOnodes for subslice ";
      for(int d=0;d<grid->_ndimension;d++){
	std::cout<< range[d];
	if( d< grid->_ndimension-1 ) 
	  std::cout<< " x ";
      }
      std::cout << std::endl;
      std::cout<< GridLogMessage<< "Parallel I/O local  strip size is "<< lstrip <<std::endl;
      std::cout<< GridLogMessage<< "Parallel I/O global strip size is "<< gstrip <<std::endl;
      std::cout<< GridLogMessage<< "Parallel I/O chunk size is "<< chunk  <<std::endl;
    }

    GridStopWatch timer; timer.Start();
    uint64_t bytes=0;

    int myrank = grid->ThisRank();
    int iorank = grid->RankFromProcessorCoor(ioproc);

    if (!ILDG.is_ILDG) {
      if ( IOnode ) { 
	fin.open(file,std::ios::binary|std::ios::in);
      }
    }

    //////////////////////////////////////////////////////////
    // Find the location of each site and send to primary node
    // Take loop order from Chroma; defines loop order now that NERSC doc no longer
    // available (how short sighted is that?)
    //////////////////////////////////////////////////////////
    Umu = zero;
    static uint32_t csum; csum=0;//static for SHMEM

    std::vector<fobj> fileObj(chunk); // FIXME
    std::vector<sobj> siteObj(chunk); // Use comm allocator to place in symmetric region for SHMEM

    // need to implement these loops in Nd independent way with a lexico conversion
    for(int tlex=0;tlex<slice_vol;tlex+=chunk){

      std::vector<int> tsite(nd); // temporary mixed up site
      std::vector<int> gsite(nd);
      std::vector<int> lsite(nd);

      Lexicographic::CoorFromIndex(tsite,tlex,range);

      for(int d=0;d<nd;d++){
	lsite[d] = tsite[d]%grid->_ldimensions[d];  // local site
	gsite[d] = tsite[d]+start[d];               // global site
      }

      ///////////////////////////////////////////
      // Get the global lexico base of the chunk
      ///////////////////////////////////////////
      int rank, o_idx,i_idx, g_idx;
      grid->GlobalCoorToRankIndex(rank,o_idx,i_idx,gsite);
      grid->GlobalCoorToGlobalIndex(gsite,g_idx);

      ////////////////////////////////
      // iorank reads from the seek
      ////////////////////////////////
      if (myrank == iorank) {

      	if (ILDG.is_ILDG){
#ifdef HAVE_LIME
	  // use C-LIME to populate the record
          uint64_t sizeFO = sizeof(fobj)*chunk;
          limeReaderSeek(ILDG.LR, g_idx*sizeFO, SEEK_SET);
          int status = limeReaderReadData((void *)&fileObj[0], &sizeFO, ILDG.LR);
#endif
        } else{
          fin.seekg(offset+g_idx*sizeof(fobj));
          fin.read((char *)&fileObj[0],sizeof(fobj)*chunk);
        }
        bytes+=sizeof(fobj)*chunk;

        if(ieee32big) be32toh_v((void *)&fileObj[0],sizeof(fobj)*chunk);
        if(ieee32)    le32toh_v((void *)&fileObj[0],sizeof(fobj)*chunk);
        if(ieee64big) be64toh_v((void *)&fileObj[0],sizeof(fobj)*chunk);
        if(ieee64)    le64toh_v((void *)&fileObj[0],sizeof(fobj)*chunk);

	for(int c=0;c<chunk;c++) munge(fileObj[c],siteObj[c],csum);

      } 
     
      // Possibly do transport through pt2pt 
      for(int cc=0;cc<chunk;cc+=lstrip){

	/////////////////////////////////
	// Get the rank of owner of strip
	/////////////////////////////////
	Lexicographic::CoorFromIndex(tsite,tlex+cc,range);

	for(int d=0;d<nd;d++){
	  lsite[d] = tsite[d]%grid->_ldimensions[d];  // local site
	  gsite[d] = tsite[d]+start[d];               // global site
	}
	grid->GlobalCoorToRankIndex(rank,o_idx,i_idx,gsite);

	if ( rank != iorank ) { 
	  if ( (myrank == rank) || (myrank==iorank) ) {
	    grid->SendRecvPacket((void *)&siteObj[cc],(void *)&siteObj[cc],iorank,rank,sizeof(sobj)*lstrip);
	  }
	}
	// Poke at destination
	if ( myrank == rank ) {
	  for(int x=0;x<lstrip;x++){
	    lsite[0]=x;
	    pokeLocalSite(siteObj[cc+x],Umu,lsite);
	  }
	}
	grid->Barrier(); // necessary?
      }
    }

    grid->GlobalSum(csum);
    grid->GlobalSum(bytes);
    grid->Barrier();

    timer.Stop();
    std::cout<<GridLogPerformance<<"readObjectParallel: read "<< bytes <<" bytes in "<<timer.Elapsed() <<" "
	     <<(double)bytes/timer.useconds() <<" MB/s "  <<std::endl;
    return csum;
  }

  //////////////////////////////////////////////////////////
  // Parallel writer
  //////////////////////////////////////////////////////////
  template <class vobj, class fobj, class munger>
  static inline uint32_t writeObjectParallel(Lattice<vobj> &Umu,
                                             std::string file, munger munge,
                                             int offset,
                                             const std::string &format,
                                             ILDGtype ILDG = ILDGtype()) {
    typedef typename vobj::scalar_object sobj;
    GridBase *grid = Umu._grid;

    int ieee32big = (format == std::string("IEEE32BIG"));
    int ieee32 = (format == std::string("IEEE32"));
    int ieee64big = (format == std::string("IEEE64BIG"));
    int ieee64 = (format == std::string("IEEE64"));

    if (!(ieee32big || ieee32 || ieee64big || ieee64)) {
      std::cout << GridLogError << "Unrecognized file format " << format << std::endl;
      std::cout << GridLogError << "Allowed: IEEE32BIG | IEEE32 | IEEE64BIG | IEEE64" << std::endl;
      exit(0);
    }

    int nd = grid->_ndimension;
    for (int d = 0; d < nd; d++) {
      assert(grid->CheckerBoarded(d) == 0);
    }

    std::vector<int> parallel(nd, 1);
    std::vector<int> ioproc(nd);
    std::vector<int> start(nd);
    std::vector<int> range(nd);

    uint64_t slice_vol = 1;

    int IOnode = 1;

    for (int d = 0; d < grid->_ndimension; d++) {
      if (d != grid->_ndimension - 1) parallel[d] = 0;

      if (parallel[d]) {
	range[d] = grid->_ldimensions[d];
	start[d] = grid->_processor_coor[d]*range[d];
	ioproc[d]= grid->_processor_coor[d];
      } else {
	range[d] = grid->_gdimensions[d];
	start[d] = 0;
	ioproc[d]= 0;

	if ( grid->_processor_coor[d] != 0 ) IOnode = 0;
      }

      slice_vol = slice_vol * range[d];
    }

    {
      uint32_t tmp = IOnode;
      grid->GlobalSum(tmp);
      std::cout<< GridLogMessage<< "Parallel write I/O from "<< file
	       << " with " <<tmp<< " IOnodes for subslice ";
      for(int d=0;d<grid->_ndimension;d++){
	std::cout<< range[d];
	if( d< grid->_ndimension-1 ) 
	  std::cout<< " x ";
      }
      std::cout << std::endl;
    }
    
    GridStopWatch timer;
    timer.Start();
    uint64_t bytes=0;

    int myrank = grid->ThisRank();
    int iorank = grid->RankFromProcessorCoor(ioproc);

    // Take into account block size of parallel file systems want about
    // 4-16MB chunks.
    // Ideally one reader/writer per xy plane and read these contiguously
    // with comms from nominated I/O nodes.
    std::ofstream fout;
    if (!ILDG.is_ILDG)
    	if (IOnode){
    		fout.open(file, std::ios::binary | std::ios::in | std::ios::out);
    		if (!fout.is_open()) {
    			std::cout << GridLogMessage << "writeObjectParallel: Error opening file " << file
    			<< std::endl;
    			exit(0);
    		}
    	}


    //////////////////////////////////////////////////////////
    // Find the location of each site and send to primary node
    // Take loop order from Chroma; defines loop order now that NERSC doc no
    // longer
    // available (how short sighted is that?)
    //////////////////////////////////////////////////////////

    uint32_t csum = 0;
    fobj fileObj;
    static sobj siteObj;  // static for SHMEM target; otherwise dynamic allocate
                          // with AlignedAllocator

    // should aggregate a whole chunk and then write.
    // need to implement these loops in Nd independent way with a lexico
    // conversion
    for (int tlex = 0; tlex < slice_vol; tlex++) {

      std::vector<int> tsite(nd);  // temporary mixed up site
      std::vector<int> gsite(nd);
      std::vector<int> lsite(nd);

      Lexicographic::CoorFromIndex(tsite, tlex, range);

      for(int d = 0;d < nd; d++){
	lsite[d] = tsite[d] % grid->_ldimensions[d];  // local site
	gsite[d] = tsite[d] + start[d];               // global site
      }

      /////////////////////////
      // Get the rank of owner of data
      /////////////////////////
      int rank, o_idx, i_idx, g_idx;
      grid->GlobalCoorToRankIndex(rank, o_idx, i_idx, gsite);
      grid->GlobalCoorToGlobalIndex(gsite, g_idx);

      ////////////////////////////////
      // iorank writes from the seek
      ////////////////////////////////

      // Owner of data peeks it
      peekLocalSite(siteObj, Umu, lsite);

      // Pair of nodes may need to do pt2pt send
      if ( rank != iorank ) { // comms is necessary
	if ( (myrank == rank) || (myrank==iorank) ) { // and we have to do it
	  // Send to IOrank 
	  grid->SendRecvPacket((void *)&siteObj,(void *)&siteObj,rank,iorank,sizeof(siteObj));
	}
      }

      grid->Barrier();  // necessary?

      if (myrank == iorank) {
        munge(siteObj, fileObj, csum);

        if (ieee32big) htobe32_v((void *)&fileObj, sizeof(fileObj));
        if (ieee32) htole32_v((void *)&fileObj, sizeof(fileObj));
        if (ieee64big) htobe64_v((void *)&fileObj, sizeof(fileObj));
        if (ieee64) htole64_v((void *)&fileObj, sizeof(fileObj));


        if (ILDG.is_ILDG) {
          #ifdef HAVE_LIME
          uint64_t sizeFO = sizeof(fileObj);
 					limeWriterSeek(ILDG.LW, g_idx*sizeFO, SEEK_SET);
          int status = limeWriteRecordData((void *)&fileObj, &sizeFO, ILDG.LW);
          #endif
        } 

        else {
          fout.seekp(offset + g_idx * sizeof(fileObj));
          fout.write((char *)&fileObj, sizeof(fileObj));assert( fout.fail()==0);
        }
        bytes += sizeof(fileObj);
      }
    }
    
    grid->GlobalSum(csum);
    grid->GlobalSum(bytes);
    
    timer.Stop();
    std::cout << GridLogPerformance << "writeObjectParallel: wrote " << bytes
              << " bytes in " << timer.Elapsed() << " "
              << (double)bytes / timer.useconds() << " MB/s " << std::endl;



     grid->Barrier();  // necessary?
     if (IOnode) 
      fout.close();


    return csum;
  }
};
}

#endif
