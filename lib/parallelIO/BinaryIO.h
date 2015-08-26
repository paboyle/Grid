#ifndef GRID_BINARY_IO_H
#define GRID_BINARY_IO_H

#ifdef HAVE_ENDIAN_H
#include <endian.h>
#endif


#include <arpa/inet.h>

// 64bit endian swap is a portability pain
#ifndef __has_builtin         // Optional of course.
#define __has_builtin(x) 0  // Compatibility with non-clang compilers.
#endif

#if HAVE_DECL_BE64TOH 
#undef Grid_ntohll
#define Grid_ntohll be64toh
#endif

#if HAVE_DECL_NTOHLL
#undef  Grid_ntohll
#define Grid_ntohll ntohll
#endif

#ifndef Grid_ntohll

#if BYTE_ORDER == BIG_ENDIAN 

#define Grid_ntohll(A) (A)

#else 

#if __has_builtin(__builtin_bswap64)
#define Grid_ntohll(A) __builtin_bswap64(A)
#else
#error
#endif

#endif

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

    sobj siteObj;
    fobj fileObj;

    csum = 0;
    std::vector<int> lcoor;
    for(int l=0;l<grid->lSites();l++){
      grid->CoorFromIndex(lcoor,l,grid->_ldimensions);
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
  
  template<class vobj,class fobj,class munger>
  static inline uint32_t readObjectSerial(Lattice<vobj> &Umu,std::string file,munger munge,int offset,const std::string &format)
  {
    typedef typename vobj::scalar_object sobj;

    GridBase *grid = Umu._grid;

    std::cout<< GridLogMessage<< "Serial read I/O "<< file<< std::endl;

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
    return csum;
  }

  template<class vobj,class fobj,class munger> 
  static inline uint32_t writeObjectSerial(Lattice<vobj> &Umu,std::string file,munger munge,int offset,const std::string & format)
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

    std::ofstream fout;
    if ( grid->IsBoss() ) {
      fout.open(file,std::ios::binary|std::ios::out|std::ios::in);
      fout.seekp(offset);
    }
    
    uint32_t csum=0;
    fobj file_object;
    sobj unmunged;
    for(int t=0;t<grid->_fdimensions[3];t++){
    for(int z=0;z<grid->_fdimensions[2];z++){
    for(int y=0;y<grid->_fdimensions[1];y++){
    for(int x=0;x<grid->_fdimensions[0];x++){

      std::vector<int> site({x,y,z,t});
      // peek & write
      peekSite(unmunged,Umu,site);

      munge(unmunged,file_object,csum);

      
      if ( grid->IsBoss() ) {
	
	if(ieee32big) htobe32_v((void *)&file_object,sizeof(file_object));
	if(ieee32)    htole32_v((void *)&file_object,sizeof(file_object));
	if(ieee64big) htobe64_v((void *)&file_object,sizeof(file_object));
	if(ieee64)    htole64_v((void *)&file_object,sizeof(file_object));
	
	fout.write((char *)&file_object,sizeof(file_object));
      }
    }}}}

    return csum;
  }

  template<class vobj,class fobj,class munger>
  static inline uint32_t readObjectParallel(Lattice<vobj> &Umu,std::string file,munger munge,int offset,const std::string &format)
  {
    typedef typename vobj::scalar_object sobj;

    GridBase *grid = Umu._grid;

    int ieee32big = (format == std::string("IEEE32BIG"));
    int ieee32    = (format == std::string("IEEE32"));
    int ieee64big = (format == std::string("IEEE64BIG"));
    int ieee64    = (format == std::string("IEEE64"));


    // Take into account block size of parallel file systems want about
    // 4-16MB chunks.
    // Ideally one reader/writer per xy plane and read these contiguously
    // with comms from nominated I/O nodes.
    std::ifstream fin;

    int nd = grid->_ndimension;
    std::vector<int> parallel(nd,1);
    std::vector<int> ioproc  (nd);
    std::vector<int> start(nd);
    std::vector<int> range(nd);

    for(int d=0;d<nd;d++){
      assert(grid->CheckerBoarded(d) == 0);
    }

    uint64_t slice_vol = 1;

    int IOnode = 1;
    for(int d=0;d<grid->_ndimension;d++) {

      if ( d==0 ) parallel[d] = 0;
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
      std::cout<< GridLogMessage<< "Parallel read I/O to "<< file << " with " <<tmp<< " IOnodes for subslice ";
      for(int d=0;d<grid->_ndimension;d++){
	std::cout<< range[d];
	if( d< grid->_ndimension-1 ) 
	  std::cout<< " x ";
      }
      std::cout << std::endl;
    }

    int myrank = grid->ThisRank();
    int iorank = grid->RankFromProcessorCoor(ioproc);

    if ( IOnode ) { 
      fin.open(file,std::ios::binary|std::ios::in);
    }

    //////////////////////////////////////////////////////////
    // Find the location of each site and send to primary node
    // Take loop order from Chroma; defines loop order now that NERSC doc no longer
    // available (how short sighted is that?)
    //////////////////////////////////////////////////////////
    Umu = zero;
    uint32_t csum=0;
    fobj fileObj;
    sobj siteObj;

      // need to implement these loops in Nd independent way with a lexico conversion
    for(int tlex=0;tlex<slice_vol;tlex++){
	
      std::vector<int> tsite(nd); // temporary mixed up site
      std::vector<int> gsite(nd);
      std::vector<int> lsite(nd);
      std::vector<int> iosite(nd);

      grid->CoorFromIndex(tsite,tlex,range);

      for(int d=0;d<nd;d++){
	lsite[d] = tsite[d]%grid->_ldimensions[d];  // local site
	gsite[d] = tsite[d]+start[d];               // global site
      }

      /////////////////////////
      // Get the rank of owner of data
      /////////////////////////
	int rank, o_idx,i_idx, g_idx;
      grid->GlobalCoorToRankIndex(rank,o_idx,i_idx,gsite);
      grid->GlobalCoorToGlobalIndex(gsite,g_idx);
      
      ////////////////////////////////
      // iorank reads from the seek
      ////////////////////////////////
      if (myrank == iorank) {
	
	fin.seekg(offset+g_idx*sizeof(fileObj));
	fin.read((char *)&fileObj,sizeof(fileObj));
	
	if(ieee32big) be32toh_v((void *)&fileObj,sizeof(fileObj));
	if(ieee32)    le32toh_v((void *)&fileObj,sizeof(fileObj));
	if(ieee64big) be64toh_v((void *)&fileObj,sizeof(fileObj));
	if(ieee64)    le64toh_v((void *)&fileObj,sizeof(fileObj));
	
	munge(fileObj,siteObj,csum);
	
	if ( rank != myrank ) {
	  grid->SendTo((void *)&siteObj,rank,sizeof(siteObj));
	} else { 
	  pokeLocalSite(siteObj,Umu,lsite);
	}
	 
      } else { 
	if ( myrank == rank ) {
	  grid->RecvFrom((void *)&siteObj,iorank,sizeof(siteObj));
	  pokeLocalSite(siteObj,Umu,lsite);
	} 
      }
      grid->Barrier(); // necessary?
    }

    grid->GlobalSum(csum);
    
    return csum;
  }

  //////////////////////////////////////////////////////////
  // Parallel writer
  //////////////////////////////////////////////////////////
  template<class vobj,class fobj,class munger>
  static inline uint32_t writeObjectParallel(Lattice<vobj> &Umu,std::string file,munger munge,int offset,const std::string & format)
  {
    typedef typename vobj::scalar_object sobj;
    GridBase *grid = Umu._grid;

    int ieee32big = (format == std::string("IEEE32BIG"));
    int ieee32    = (format == std::string("IEEE32"));
    int ieee64big = (format == std::string("IEEE64BIG"));
    int ieee64    = (format == std::string("IEEE64"));

    int nd = grid->_ndimension;
    for(int d=0;d<nd;d++){
      assert(grid->CheckerBoarded(d) == 0);
    }

    std::vector<int> parallel(nd,1);
    std::vector<int> ioproc  (nd);
    std::vector<int> start(nd);
    std::vector<int> range(nd);

    uint64_t slice_vol = 1;

    int IOnode = 1;

    for(int d=0;d<grid->_ndimension;d++) {

      if ( d==0 ) parallel[d] = 0;

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
      std::cout<< GridLogMessage<< "Parallel write I/O from "<< file << " with " <<tmp<< " IOnodes for subslice ";
      for(int d=0;d<grid->_ndimension;d++){
	std::cout<< range[d];
	if( d< grid->_ndimension-1 ) 
	  std::cout<< " x ";
      }
      std::cout << std::endl;
    }

    int myrank = grid->ThisRank();
    int iorank = grid->RankFromProcessorCoor(ioproc);

    // Take into account block size of parallel file systems want about
    // 4-16MB chunks.
    // Ideally one reader/writer per xy plane and read these contiguously
    // with comms from nominated I/O nodes.
    std::ofstream fout;
    if ( IOnode ) fout.open(file,std::ios::binary|std::ios::in|std::ios::out);

    //////////////////////////////////////////////////////////
    // Find the location of each site and send to primary node
    // Take loop order from Chroma; defines loop order now that NERSC doc no longer
    // available (how short sighted is that?)
    //////////////////////////////////////////////////////////

    uint32_t csum=0;
    fobj fileObj;
    sobj siteObj;


      // need to implement these loops in Nd independent way with a lexico conversion
    for(int tlex=0;tlex<slice_vol;tlex++){
	
      std::vector<int> tsite(nd); // temporary mixed up site
      std::vector<int> gsite(nd);
      std::vector<int> lsite(nd);
      std::vector<int> iosite(nd);

      grid->CoorFromIndex(tsite,tlex,range);

      for(int d=0;d<nd;d++){
	lsite[d] = tsite[d]%grid->_ldimensions[d];  // local site
	gsite[d] = tsite[d]+start[d];               // global site
      }


      /////////////////////////
      // Get the rank of owner of data
      /////////////////////////
      int rank, o_idx,i_idx, g_idx;
      grid->GlobalCoorToRankIndex(rank,o_idx,i_idx,gsite);
      grid->GlobalCoorToGlobalIndex(gsite,g_idx);

      ////////////////////////////////
      // iorank writes from the seek
      ////////////////////////////////
      if (myrank == iorank) {

	if ( rank != myrank ) {
	  grid->RecvFrom((void *)&siteObj,rank,sizeof(siteObj));
	} else { 
	  peekLocalSite(siteObj,Umu,lsite);
	}
	
	munge(siteObj,fileObj,csum);

	if(ieee32big) htobe32_v((void *)&fileObj,sizeof(fileObj));
	if(ieee32)    htole32_v((void *)&fileObj,sizeof(fileObj));
	if(ieee64big) htobe64_v((void *)&fileObj,sizeof(fileObj));
	if(ieee64)    htole64_v((void *)&fileObj,sizeof(fileObj));
	
	fout.seekp(offset+g_idx*sizeof(fileObj));
	fout.write((char *)&fileObj,sizeof(fileObj));

      } else { 
	if ( myrank == rank ) {
	  peekLocalSite(siteObj,Umu,lsite);
	  grid->SendTo((void *)&siteObj,iorank,sizeof(siteObj));
	} 
      }
      grid->Barrier(); // necessary// or every 16 packets to rate throttle??
    }

    grid->GlobalSum(csum);

    return csum;
  }

};

}

#endif
