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

#if defined(GRID_COMMS_MPI) || defined(GRID_COMMS_MPI3) || defined(GRID_COMMS_MPIT) 
#define USE_MPI_IO
#else
#undef  USE_MPI_IO
#endif

#ifdef HAVE_ENDIAN_H
#include <endian.h>
#endif

#include <arpa/inet.h>
#include <algorithm>

namespace Grid { 


/////////////////////////////////////////////////////////////////////////////////
// Byte reversal garbage
/////////////////////////////////////////////////////////////////////////////////
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

// A little helper
inline void removeWhitespace(std::string &key)
{
  key.erase(std::remove_if(key.begin(), key.end(), ::isspace),key.end());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Static class holding the parallel IO code
// Could just use a namespace
///////////////////////////////////////////////////////////////////////////////////////////////////
class BinaryIO {
 public:

  /////////////////////////////////////////////////////////////////////////////
  // more byte manipulation helpers
  /////////////////////////////////////////////////////////////////////////////

  template<class vobj> static inline void Uint32Checksum(Lattice<vobj> &lat,uint32_t &nersc_csum)
  {
    typedef typename vobj::scalar_object sobj;

    GridBase *grid = lat._grid;
    int lsites = grid->lSites();

    std::vector<sobj> scalardata(lsites); 
    unvectorizeToLexOrdArray(scalardata,lat);    

    NerscChecksum(grid,scalardata,nersc_csum);
  }

  template <class fobj>
  static inline void NerscChecksum(GridBase *grid, std::vector<fobj> &fbuf, uint32_t &nersc_csum)
  {
    const uint64_t size32 = sizeof(fobj) / sizeof(uint32_t);

    uint64_t lsites = grid->lSites();
    if (fbuf.size() == 1)
    {
      lsites = 1;
    }

    #pragma omp parallel
    {
      uint32_t nersc_csum_thr = 0;

      #pragma omp for
      for (uint64_t local_site = 0; local_site < lsites; local_site++)
      {
        uint32_t *site_buf = (uint32_t *)&fbuf[local_site];
        for (uint64_t j = 0; j < size32; j++)
        {
          nersc_csum_thr = nersc_csum_thr + site_buf[j];
        }
      }

      #pragma omp critical
      {
        nersc_csum += nersc_csum_thr;
      }
    }
  }

  template<class fobj> static inline void ScidacChecksum(GridBase *grid,std::vector<fobj> &fbuf,uint32_t &scidac_csuma,uint32_t &scidac_csumb)
  {
    const uint64_t size32 = sizeof(fobj)/sizeof(uint32_t);


    int nd = grid->_ndimension;

    uint64_t lsites              =grid->lSites();
    if (fbuf.size()==1) {
      lsites=1;
    }
    std::vector<int> local_vol   =grid->LocalDimensions();
    std::vector<int> local_start =grid->LocalStarts();
    std::vector<int> global_vol  =grid->FullDimensions();

#pragma omp parallel
    { 
      std::vector<int> coor(nd);
      uint32_t scidac_csuma_thr=0;
      uint32_t scidac_csumb_thr=0;
      uint32_t site_crc=0;

#pragma omp for
      for(uint64_t local_site=0;local_site<lsites;local_site++){

	uint32_t * site_buf = (uint32_t *)&fbuf[local_site];

	/* 
	 * Scidac csum  is rather more heavyweight
	 */
	int global_site;

	Lexicographic::CoorFromIndex(coor,local_site,local_vol);

	for(int d=0;d<nd;d++) {
	  coor[d] = coor[d]+local_start[d];
	}

	Lexicographic::IndexFromCoor(coor,global_site,global_vol);

	uint32_t gsite29   = global_site%29;
	uint32_t gsite31   = global_site%31;
	
	site_crc = crc32(0,(unsigned char *)site_buf,sizeof(fobj));
	//	std::cout << "Site "<<local_site << " crc "<<std::hex<<site_crc<<std::dec<<std::endl;
	//	std::cout << "Site "<<local_site << std::hex<<site_buf[0] <<site_buf[1]<<std::dec <<std::endl;
	scidac_csuma_thr ^= site_crc<<gsite29 | site_crc>>(32-gsite29);
	scidac_csumb_thr ^= site_crc<<gsite31 | site_crc>>(32-gsite31);
      }

#pragma omp critical
      {
	scidac_csuma^= scidac_csuma_thr;
	scidac_csumb^= scidac_csumb_thr;
      }
    }
  }

  // Network is big endian
  static inline void htobe32_v(void *file_object,uint32_t bytes){ be32toh_v(file_object,bytes);} 
  static inline void htobe64_v(void *file_object,uint32_t bytes){ be64toh_v(file_object,bytes);} 
  static inline void htole32_v(void *file_object,uint32_t bytes){ le32toh_v(file_object,bytes);} 
  static inline void htole64_v(void *file_object,uint32_t bytes){ le64toh_v(file_object,bytes);} 

  static inline void be32toh_v(void *file_object,uint64_t bytes)
  {
    uint32_t * f = (uint32_t *)file_object;
    uint64_t count = bytes/sizeof(uint32_t);
    parallel_for(uint64_t i=0;i<count;i++){  
      f[i] = ntohl(f[i]);
    }
  }
  // LE must Swap and switch to host
  static inline void le32toh_v(void *file_object,uint64_t bytes)
  {
    uint32_t *fp = (uint32_t *)file_object;
    uint32_t f;

    uint64_t count = bytes/sizeof(uint32_t);
    parallel_for(uint64_t i=0;i<count;i++){  
      f = fp[i];
      // got network order and the network to host
      f = ((f&0xFF)<<24) | ((f&0xFF00)<<8) | ((f&0xFF0000)>>8) | ((f&0xFF000000UL)>>24) ; 
      fp[i] = ntohl(f);
    }
  }

  // BE is same as network
  static inline void be64toh_v(void *file_object,uint64_t bytes)
  {
    uint64_t * f = (uint64_t *)file_object;
    uint64_t count = bytes/sizeof(uint64_t);
    parallel_for(uint64_t i=0;i<count;i++){  
      f[i] = Grid_ntohll(f[i]);
    }
  }
  
  // LE must swap and switch;
  static inline void le64toh_v(void *file_object,uint64_t bytes)
  {
    uint64_t *fp = (uint64_t *)file_object;
    uint64_t f,g;
    
    uint64_t count = bytes/sizeof(uint64_t);
    parallel_for(uint64_t i=0;i<count;i++){  
      f = fp[i];
      // got network order and the network to host
      g = ((f&0xFF)<<24) | ((f&0xFF00)<<8) | ((f&0xFF0000)>>8) | ((f&0xFF000000UL)>>24) ; 
      g = g << 32;
      f = f >> 32;
      g|= ((f&0xFF)<<24) | ((f&0xFF00)<<8) | ((f&0xFF0000)>>8) | ((f&0xFF000000UL)>>24) ; 
      fp[i] = Grid_ntohll(g);
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  // Real action:
  // Read or Write distributed lexico array of ANY object to a specific location in file 
  //////////////////////////////////////////////////////////////////////////////////////

  static const int BINARYIO_MASTER_APPEND = 0x10;
  static const int BINARYIO_UNORDERED     = 0x08;
  static const int BINARYIO_LEXICOGRAPHIC = 0x04;
  static const int BINARYIO_READ          = 0x02;
  static const int BINARYIO_WRITE         = 0x01;

  template<class word,class fobj>
  static inline void IOobject(word w,
			      GridBase *grid,
			      std::vector<fobj> &iodata,
			      std::string file,
			      Integer offset,
			      const std::string &format, int control,
			      uint32_t &nersc_csum,
			      uint32_t &scidac_csuma,
			      uint32_t &scidac_csumb)
  {
    grid->Barrier();
    GridStopWatch timer; 
    GridStopWatch bstimer;
    
    nersc_csum=0;
    scidac_csuma=0;
    scidac_csumb=0;

    int ndim                 = grid->Dimensions();
    int nrank                = grid->ProcessorCount();
    int myrank               = grid->ThisRank();

    std::vector<int>  psizes = grid->ProcessorGrid(); 
    std::vector<int>  pcoor  = grid->ThisProcessorCoor();
    std::vector<int> gLattice= grid->GlobalDimensions();
    std::vector<int> lLattice= grid->LocalDimensions();

    std::vector<int> lStart(ndim);
    std::vector<int> gStart(ndim);

    // Flatten the file
    uint64_t lsites = grid->lSites();
    if ( control & BINARYIO_MASTER_APPEND )  {
      assert(iodata.size()==1);
    } else {
      assert(lsites==iodata.size());
    }
    for(int d=0;d<ndim;d++){
      gStart[d] = lLattice[d]*pcoor[d];
      lStart[d] = 0;
    }

#ifdef USE_MPI_IO
    std::vector<int> distribs(ndim,MPI_DISTRIBUTE_BLOCK);
    std::vector<int> dargs   (ndim,MPI_DISTRIBUTE_DFLT_DARG);
    MPI_Datatype mpiObject;
    MPI_Datatype fileArray;
    MPI_Datatype localArray;
    MPI_Datatype mpiword;
    MPI_Offset disp = offset;
    MPI_File fh ;
    MPI_Status status;
    int numword;

    if ( sizeof( word ) == sizeof(float ) ) {
      numword = sizeof(fobj)/sizeof(float);
      mpiword = MPI_FLOAT;
    } else {
      numword = sizeof(fobj)/sizeof(double);
      mpiword = MPI_DOUBLE;
    }

    //////////////////////////////////////////////////////////////////////////////
    // Sobj in MPI phrasing
    //////////////////////////////////////////////////////////////////////////////
    int ierr;
    ierr = MPI_Type_contiguous(numword,mpiword,&mpiObject);    assert(ierr==0);
    ierr = MPI_Type_commit(&mpiObject);

    //////////////////////////////////////////////////////////////////////////////
    // File global array data type
    //////////////////////////////////////////////////////////////////////////////
    ierr=MPI_Type_create_subarray(ndim,&gLattice[0],&lLattice[0],&gStart[0],MPI_ORDER_FORTRAN, mpiObject,&fileArray);    assert(ierr==0);
    ierr=MPI_Type_commit(&fileArray);    assert(ierr==0);

    //////////////////////////////////////////////////////////////////////////////
    // local lattice array
    //////////////////////////////////////////////////////////////////////////////
    ierr=MPI_Type_create_subarray(ndim,&lLattice[0],&lLattice[0],&lStart[0],MPI_ORDER_FORTRAN, mpiObject,&localArray);    assert(ierr==0);
    ierr=MPI_Type_commit(&localArray);    assert(ierr==0);
#endif

    //////////////////////////////////////////////////////////////////////////////
    // Byte order
    //////////////////////////////////////////////////////////////////////////////
    int ieee32big = (format == std::string("IEEE32BIG"));
    int ieee32    = (format == std::string("IEEE32"));
    int ieee64big = (format == std::string("IEEE64BIG"));
    int ieee64    = (format == std::string("IEEE64"));

    //////////////////////////////////////////////////////////////////////////////
    // Do the I/O
    //////////////////////////////////////////////////////////////////////////////
    if ( control & BINARYIO_READ ) { 

      timer.Start();

      if ( (control & BINARYIO_LEXICOGRAPHIC) && (nrank > 1) ) {
#ifdef USE_MPI_IO
	std::cout<< GridLogMessage<<"IOobject: MPI read I/O "<< file<< std::endl;
	ierr=MPI_File_open(grid->communicator,(char *) file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);    assert(ierr==0);
	ierr=MPI_File_set_view(fh, disp, mpiObject, fileArray, "native", MPI_INFO_NULL);    assert(ierr==0);
	ierr=MPI_File_read_all(fh, &iodata[0], 1, localArray, &status);    assert(ierr==0);
	MPI_File_close(&fh);
	MPI_Type_free(&fileArray);
	MPI_Type_free(&localArray);
#else 
	assert(0);
#endif
      } else {
	std::cout << GridLogMessage <<"IOobject: C++ read I/O " << file << " : "
                  << iodata.size() * sizeof(fobj) << " bytes" << std::endl;
        std::ifstream fin;
        fin.open(file, std::ios::binary | std::ios::in);
        if (control & BINARYIO_MASTER_APPEND)
        {
          fin.seekg(-sizeof(fobj), fin.end);
        }
        else
        {
          fin.seekg(offset + myrank * lsites * sizeof(fobj));
        }
        fin.read((char *)&iodata[0], iodata.size() * sizeof(fobj));
        assert(fin.fail() == 0);
        fin.close();
      }
      timer.Stop();

      grid->Barrier();

      bstimer.Start();
      ScidacChecksum(grid,iodata,scidac_csuma,scidac_csumb);
      if (ieee32big) be32toh_v((void *)&iodata[0], sizeof(fobj)*iodata.size());
      if (ieee32)    le32toh_v((void *)&iodata[0], sizeof(fobj)*iodata.size());
      if (ieee64big) be64toh_v((void *)&iodata[0], sizeof(fobj)*iodata.size());
      if (ieee64)    le64toh_v((void *)&iodata[0], sizeof(fobj)*iodata.size());
      NerscChecksum(grid,iodata,nersc_csum);
      bstimer.Stop();
    }
    
    if ( control & BINARYIO_WRITE ) { 

      bstimer.Start();
      NerscChecksum(grid,iodata,nersc_csum);
      if (ieee32big) htobe32_v((void *)&iodata[0], sizeof(fobj)*iodata.size());
      if (ieee32)    htole32_v((void *)&iodata[0], sizeof(fobj)*iodata.size());
      if (ieee64big) htobe64_v((void *)&iodata[0], sizeof(fobj)*iodata.size());
      if (ieee64)    htole64_v((void *)&iodata[0], sizeof(fobj)*iodata.size());
      ScidacChecksum(grid,iodata,scidac_csuma,scidac_csumb);
      bstimer.Stop();

      grid->Barrier();

      timer.Start();
      if ( (control & BINARYIO_LEXICOGRAPHIC) && (nrank > 1) ) {
#ifdef USE_MPI_IO
        std::cout << GridLogMessage <<"IOobject: MPI write I/O " << file << std::endl;
        ierr = MPI_File_open(grid->communicator, (char *)file.c_str(), MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	//        std::cout << GridLogMessage << "Checking for errors" << std::endl;
        if (ierr != MPI_SUCCESS)
        {
          char error_string[BUFSIZ];
          int length_of_error_string, error_class;

          MPI_Error_class(ierr, &error_class);
          MPI_Error_string(error_class, error_string, &length_of_error_string);
          fprintf(stderr, "%3d: %s\n", myrank, error_string);
          MPI_Error_string(ierr, error_string, &length_of_error_string);
          fprintf(stderr, "%3d: %s\n", myrank, error_string);
          MPI_Abort(MPI_COMM_WORLD, 1); //assert(ierr == 0);
        }

        std::cout << GridLogDebug << "MPI read I/O set view " << file << std::endl;
        ierr = MPI_File_set_view(fh, disp, mpiObject, fileArray, "native", MPI_INFO_NULL);
        assert(ierr == 0);

        std::cout << GridLogDebug << "MPI read I/O write all " << file << std::endl;
        ierr = MPI_File_write_all(fh, &iodata[0], 1, localArray, &status);
        assert(ierr == 0);

        MPI_File_close(&fh);
        MPI_Type_free(&fileArray);
        MPI_Type_free(&localArray);
#else 
	assert(0);
#endif
      } else { 

        std::cout << GridLogMessage << "IOobject: C++ write I/O " << file << " : "
                  << iodata.size() * sizeof(fobj) << " bytes" << std::endl;
        
	std::ofstream fout; 
	fout.exceptions ( std::fstream::failbit | std::fstream::badbit );
	try {
	  fout.open(file,std::ios::binary|std::ios::out|std::ios::in);
	} catch (const std::fstream::failure& exc) {
	  std::cout << GridLogError << "Error in opening the file " << file << " for output" <<std::endl;
	  std::cout << GridLogError << "Exception description: " << exc.what() << std::endl;
	  std::cout << GridLogError << "Probable cause: wrong path, inaccessible location "<< std::endl;
#ifdef USE_MPI_IO
	  MPI_Abort(MPI_COMM_WORLD,1);
#else
	  exit(1);
#endif
	}
	
	if ( control & BINARYIO_MASTER_APPEND )  {
	  try {
	    fout.seekp(0,fout.end);
	  } catch (const std::fstream::failure& exc) {
	    std::cout << "Exception in seeking file end " << file << std::endl;
	  }
	} else {
	  try { 
	    fout.seekp(offset+myrank*lsites*sizeof(fobj));
	  } catch (const std::fstream::failure& exc) {
	    std::cout << "Exception in seeking file " << file <<" offset "<< offset << std::endl;
	  }
	}

	try {
	  fout.write((char *)&iodata[0],iodata.size()*sizeof(fobj));//assert( fout.fail()==0);
	}
	catch (const std::fstream::failure& exc) {
	  std::cout << "Exception in writing file " << file << std::endl;
	  std::cout << GridLogError << "Exception description: "<< exc.what() << std::endl;
#ifdef USE_MPI_IO
	  MPI_Abort(MPI_COMM_WORLD,1);
#else
	  exit(1);
#endif
	}
	fout.close();
      }
      timer.Stop();
    }
    
    std::cout<<GridLogMessage<<"IOobject: ";
    if ( control & BINARYIO_READ) std::cout << " read  ";
    else                          std::cout << " write ";
    uint64_t bytes = sizeof(fobj)*iodata.size()*nrank;
    std::cout<< bytes <<" bytes in "<<timer.Elapsed() <<" "
	     << (double)bytes/ (double)timer.useconds() <<" MB/s "<<std::endl;

    std::cout<<GridLogMessage<<"IOobject: endian and checksum overhead "<<bstimer.Elapsed()  <<std::endl;

    //////////////////////////////////////////////////////////////////////////////
    // Safety check
    //////////////////////////////////////////////////////////////////////////////
    // if the data size is 1 we do not want to sum over the MPI ranks
    if (iodata.size() != 1){
      grid->Barrier();
      grid->GlobalSum(nersc_csum);
      grid->GlobalXOR(scidac_csuma);
      grid->GlobalXOR(scidac_csumb);
      grid->Barrier();
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // Read a Lattice of object
  //////////////////////////////////////////////////////////////////////////////////////
  template<class vobj,class fobj,class munger>
  static inline void readLatticeObject(Lattice<vobj> &Umu,
				       std::string file,
				       munger munge,
				       Integer offset,
				       const std::string &format,
				       uint32_t &nersc_csum,
				       uint32_t &scidac_csuma,
				       uint32_t &scidac_csumb)
  {
    typedef typename vobj::scalar_object sobj;
    typedef typename vobj::Realified::scalar_type word;    word w=0;

    GridBase *grid = Umu._grid;
    int lsites = grid->lSites();

    std::vector<sobj> scalardata(lsites); 
    std::vector<fobj>     iodata(lsites); // Munge, checksum, byte order in here
    
    IOobject(w,grid,iodata,file,offset,format,BINARYIO_READ|BINARYIO_LEXICOGRAPHIC,
	     nersc_csum,scidac_csuma,scidac_csumb);

    GridStopWatch timer; 
    timer.Start();

    parallel_for(int x=0;x<lsites;x++) munge(iodata[x], scalardata[x]);

    vectorizeFromLexOrdArray(scalardata,Umu);    
    grid->Barrier();

    timer.Stop();
    std::cout<<GridLogMessage<<"readLatticeObject: vectorize overhead "<<timer.Elapsed()  <<std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////
  // Write a Lattice of object
  //////////////////////////////////////////////////////////////////////////////////////
  template<class vobj,class fobj,class munger>
    static inline void writeLatticeObject(Lattice<vobj> &Umu,
					  std::string file,
					  munger munge,
					  Integer offset,
					  const std::string &format,
					  uint32_t &nersc_csum,
					  uint32_t &scidac_csuma,
					  uint32_t &scidac_csumb)
  {
    typedef typename vobj::scalar_object sobj;
    typedef typename vobj::Realified::scalar_type word;    word w=0;
    GridBase *grid = Umu._grid;
    int lsites = grid->lSites();

    std::vector<sobj> scalardata(lsites); 
    std::vector<fobj>     iodata(lsites); // Munge, checksum, byte order in here

    //////////////////////////////////////////////////////////////////////////////
    // Munge [ .e.g 3rd row recon ]
    //////////////////////////////////////////////////////////////////////////////
    GridStopWatch timer; timer.Start();
    unvectorizeToLexOrdArray(scalardata,Umu);    

    parallel_for(int x=0;x<lsites;x++) munge(scalardata[x],iodata[x]);

    grid->Barrier();
    timer.Stop();

    IOobject(w,grid,iodata,file,offset,format,BINARYIO_WRITE|BINARYIO_LEXICOGRAPHIC,
	     nersc_csum,scidac_csuma,scidac_csumb);

    std::cout<<GridLogMessage<<"writeLatticeObject: unvectorize overhead "<<timer.Elapsed()  <<std::endl;
  }
  
  /////////////////////////////////////////////////////////////////////////////
  // Read a RNG;  use IOobject and lexico map to an array of state 
  //////////////////////////////////////////////////////////////////////////////////////
  static inline void readRNG(GridSerialRNG &serial,
			     GridParallelRNG &parallel,
			     std::string file,
			     Integer offset,
			     uint32_t &nersc_csum,
			     uint32_t &scidac_csuma,
			     uint32_t &scidac_csumb)
  {
    typedef typename GridSerialRNG::RngStateType RngStateType;
    const int RngStateCount = GridSerialRNG::RngStateCount;
    typedef std::array<RngStateType,RngStateCount> RNGstate;
    typedef RngStateType word;    word w=0;

    std::string format = "IEEE32BIG";

    GridBase *grid = parallel._grid;
    int gsites = grid->gSites();
    int lsites = grid->lSites();

    uint32_t nersc_csum_tmp   = 0;
    uint32_t scidac_csuma_tmp = 0;
    uint32_t scidac_csumb_tmp = 0;

    GridStopWatch timer;

    std::cout << GridLogMessage << "RNG read I/O on file " << file << std::endl;

    std::vector<RNGstate> iodata(lsites);
    IOobject(w,grid,iodata,file,offset,format,BINARYIO_READ|BINARYIO_LEXICOGRAPHIC,
	     nersc_csum,scidac_csuma,scidac_csumb);

    timer.Start();
    parallel_for(int lidx=0;lidx<lsites;lidx++){
      std::vector<RngStateType> tmp(RngStateCount);
      std::copy(iodata[lidx].begin(),iodata[lidx].end(),tmp.begin());
      parallel.SetState(tmp,lidx);
    }
    timer.Stop();

    iodata.resize(1);
    IOobject(w,grid,iodata,file,offset,format,BINARYIO_READ|BINARYIO_MASTER_APPEND,
	     nersc_csum_tmp,scidac_csuma_tmp,scidac_csumb_tmp);

    {
      std::vector<RngStateType> tmp(RngStateCount);
      std::copy(iodata[0].begin(),iodata[0].end(),tmp.begin());
      serial.SetState(tmp,0);
    }

    nersc_csum   = nersc_csum   + nersc_csum_tmp;
    scidac_csuma = scidac_csuma ^ scidac_csuma_tmp;
    scidac_csumb = scidac_csumb ^ scidac_csumb_tmp;

    std::cout << GridLogMessage << "RNG file nersc_checksum   " << std::hex << nersc_csum << std::dec << std::endl;
    std::cout << GridLogMessage << "RNG file scidac_checksuma " << std::hex << scidac_csuma << std::dec << std::endl;
    std::cout << GridLogMessage << "RNG file scidac_checksumb " << std::hex << scidac_csumb << std::dec << std::endl;

    std::cout << GridLogMessage << "RNG state overhead " << timer.Elapsed() << std::endl;
  }
  /////////////////////////////////////////////////////////////////////////////
  // Write a RNG; lexico map to an array of state and use IOobject
  //////////////////////////////////////////////////////////////////////////////////////
  static inline void writeRNG(GridSerialRNG &serial,
			      GridParallelRNG &parallel,
			      std::string file,
			      Integer offset,
			      uint32_t &nersc_csum,
			      uint32_t &scidac_csuma,
			      uint32_t &scidac_csumb)
  {
    typedef typename GridSerialRNG::RngStateType RngStateType;
    typedef RngStateType word; word w=0;
    const int RngStateCount = GridSerialRNG::RngStateCount;
    typedef std::array<RngStateType,RngStateCount> RNGstate;

    GridBase *grid = parallel._grid;
    int gsites = grid->gSites();
    int lsites = grid->lSites();

    uint32_t nersc_csum_tmp;
    uint32_t scidac_csuma_tmp;
    uint32_t scidac_csumb_tmp;

    GridStopWatch timer;
    std::string format = "IEEE32BIG";

    std::cout << GridLogMessage << "RNG write I/O on file " << file << std::endl;

    timer.Start();
    std::vector<RNGstate> iodata(lsites);
    parallel_for(int lidx=0;lidx<lsites;lidx++){
      std::vector<RngStateType> tmp(RngStateCount);
      parallel.GetState(tmp,lidx);
      std::copy(tmp.begin(),tmp.end(),iodata[lidx].begin());
    }
    timer.Stop();

    IOobject(w,grid,iodata,file,offset,format,BINARYIO_WRITE|BINARYIO_LEXICOGRAPHIC,
	     nersc_csum,scidac_csuma,scidac_csumb);

    iodata.resize(1);
    {
      std::vector<RngStateType> tmp(RngStateCount);
      serial.GetState(tmp,0);
      std::copy(tmp.begin(),tmp.end(),iodata[0].begin());
    }
    IOobject(w,grid,iodata,file,offset,format,BINARYIO_WRITE|BINARYIO_MASTER_APPEND,
	     nersc_csum_tmp,scidac_csuma_tmp,scidac_csumb_tmp);

    nersc_csum   = nersc_csum   + nersc_csum_tmp;
    scidac_csuma = scidac_csuma ^ scidac_csuma_tmp;
    scidac_csumb = scidac_csumb ^ scidac_csumb_tmp;
    
    std::cout << GridLogMessage << "RNG file checksum " << std::hex << nersc_csum    << std::dec << std::endl;
    std::cout << GridLogMessage << "RNG file checksuma " << std::hex << scidac_csuma << std::dec << std::endl;
    std::cout << GridLogMessage << "RNG file checksumb " << std::hex << scidac_csumb << std::dec << std::endl;
    std::cout << GridLogMessage << "RNG state overhead " << timer.Elapsed() << std::endl;
  }
};
}
#endif
