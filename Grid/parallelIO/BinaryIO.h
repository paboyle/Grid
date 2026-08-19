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
#pragma once

#if defined(GRID_COMMS_MPI) || defined(GRID_COMMS_MPI3) || defined(GRID_COMMS_MPIT) 
#define USE_MPI_IO
#else
#undef  USE_MPI_IO
#endif

#ifdef HAVE_ENDIAN_H
#include <endian.h>
#endif

#include <arpa/inet.h>
#include <sys/stat.h>
#include <algorithm>

NAMESPACE_BEGIN(Grid);

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
  struct IoPerf
  {
    uint64_t size{0},time{0};
    double   mbytesPerSecond{0.};
  };

  static IoPerf lastPerf;
  static int latticeWriteMaxRetry;
  static uint64_t aggregateTargetBytes;

  /////////////////////////////////////////////////////////////////////////////
  // more byte manipulation helpers
  /////////////////////////////////////////////////////////////////////////////

  template<class vobj> static inline void Uint32Checksum(Lattice<vobj> &lat,uint32_t &nersc_csum)
  {
    typedef typename vobj::scalar_object sobj;

    GridBase *grid = lat.Grid();
    uint64_t lsites = grid->lSites();

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

    thread_region
    {
      uint32_t nersc_csum_thr = 0;

      thread_for_in_region( local_site, lsites, 
      {
        uint32_t *site_buf = (uint32_t *)&fbuf[local_site];
        for (uint64_t j = 0; j < size32; j++)
        {
          nersc_csum_thr = nersc_csum_thr + site_buf[j];
        }
      });

      thread_critical
      {
        nersc_csum += nersc_csum_thr;
      }
    }
  }

  template<class fobj> static inline void ScidacChecksum(GridBase *grid,std::vector<fobj> &fbuf,uint32_t &scidac_csuma,uint32_t &scidac_csumb)
  {
    int nd = grid->_ndimension;

    uint64_t lsites              =grid->lSites();
    if (fbuf.size()==1) {
      lsites=1;
    }
    Coordinate local_vol   =grid->LocalDimensions();
    Coordinate local_start =grid->LocalStarts();
    Coordinate global_vol  =grid->FullDimensions();

    thread_region
    { 
      Coordinate coor(nd);
      uint32_t scidac_csuma_thr=0;
      uint32_t scidac_csumb_thr=0;
      uint32_t site_crc=0;

      thread_for_in_region( local_site, lsites, 
      {

	uint32_t * site_buf = (uint32_t *)&fbuf[local_site];

	/* 
	 * Scidac csum  is rather more heavyweight
	 * FIXME -- 128^3 x 256 x 16 will overflow.
	 */
	
	int64_t global_site;

	Lexicographic::CoorFromIndex(coor,local_site,local_vol);

	for(int d=0;d<nd;d++) {
	  coor[d] = coor[d]+local_start[d];
	}

	Lexicographic::IndexFromCoor(coor,global_site,global_vol);

	uint64_t gsite29   = global_site%29;
	uint64_t gsite31   = global_site%31;
	
	site_crc = crc32(0,(unsigned char *)site_buf,sizeof(fobj));
	//	std::cout << "Site "<<local_site << " crc "<<std::hex<<site_crc<<std::dec<<std::endl;
	//	std::cout << "Site "<<local_site << std::hex<<site_buf[0] <<site_buf[1]<<std::dec <<std::endl;
	scidac_csuma_thr ^= site_crc<<gsite29 | site_crc>>(32-gsite29);
	scidac_csumb_thr ^= site_crc<<gsite31 | site_crc>>(32-gsite31);
      });

      thread_critical
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
    thread_for( i, count, {  
      f[i] = ntohl(f[i]);
    });
  }
  // LE must Swap and switch to host
  static inline void le32toh_v(void *file_object,uint64_t bytes)
  {
    uint32_t *fp = (uint32_t *)file_object;

    uint64_t count = bytes/sizeof(uint32_t);
    thread_for(i,count,{
      uint32_t f;
      f = fp[i];
      // got network order and the network to host
      f = ((f&0xFF)<<24) | ((f&0xFF00)<<8) | ((f&0xFF0000)>>8) | ((f&0xFF000000UL)>>24) ; 
      fp[i] = ntohl(f);
    });
  }

  // BE is same as network
  static inline void be64toh_v(void *file_object,uint64_t bytes)
  {
    uint64_t * f = (uint64_t *)file_object;
    uint64_t count = bytes/sizeof(uint64_t);
    thread_for( i, count, {
      f[i] = Grid_ntohll(f[i]);
    });
  }
  
  // LE must swap and switch;
  static inline void le64toh_v(void *file_object,uint64_t bytes)
  {
    uint64_t *fp = (uint64_t *)file_object;
    uint64_t count = bytes/sizeof(uint64_t);
    thread_for( i, count, {
      uint64_t f,g;
      f = fp[i];
      // got network order and the network to host
      g = ((f&0xFF)<<24) | ((f&0xFF00)<<8) | ((f&0xFF0000)>>8) | ((f&0xFF000000UL)>>24) ; 
      g = g << 32;
      f = f >> 32;
      g|= ((f&0xFF)<<24) | ((f&0xFF00)<<8) | ((f&0xFF0000)>>8) | ((f&0xFF000000UL)>>24) ; 
      fp[i] = Grid_ntohll(g);
    });
  }
  /////////////////////////////////////////////////////////////////////////////
  // Real action:
  // Read or Write distributed lexico array of ANY object to a specific location in file 
  //////////////////////////////////////////////////////////////////////////////////////

  static const int BINARYIO_AGGREGATE     = 0x20;
  static const int BINARYIO_MASTER_APPEND = 0x10;
  static const int BINARYIO_UNORDERED     = 0x08;
  static const int BINARYIO_LEXICOGRAPHIC = 0x04;
  static const int BINARYIO_READ          = 0x02;
  static const int BINARYIO_WRITE         = 0x01;

  // Single point of control for the aggregated path. Setting
  // GRID_BINARYIO_NOAGGREGATE falls back to plain lexicographic I/O.
  static int DefaultControl(void)
  {
    static int ctrl = getenv("GRID_BINARYIO_NOAGGREGATE")
                    ? BINARYIO_LEXICOGRAPHIC
                    : BINARYIO_LEXICOGRAPHIC|BINARYIO_AGGREGATE;
    return ctrl;
  }

#ifdef USE_MPI_IO
  /////////////////////////////////////////////////////////////////////////////
  // Aggregation: self controlled transposition onto an I/O friendly layout.
  //
  // Under BINARYIO_LEXICOGRAPHIC the subarray file view handed to MPI-IO has
  // contiguous runs of only lLattice[0]*sizeof(fobj) bytes -- a few KB for
  // typical local volumes. Rather than rely on collective buffering to repair
  // that, redistribute the payload ourselves so every rank owns a contiguous
  // range of the global lexicographic site ordering, then issue large plain
  // contiguous writes.
  //
  // "Un-splitting" the nunsplit fastest dimensions means the row of ranks
  // sharing the remaining process coordinates collectively owns whole global
  // hyperplanes. All data movement is then confined to that row communicator.
  // Every rank still owns exactly lSites() sites afterwards, so the exchange is
  // a pure permutation and needs no divisibility condition on the process grid.
  /////////////////////////////////////////////////////////////////////////////
  struct AggregationPlan {
    int      nunsplit{0};  // number of fastest dimensions un-split
    int      rowsize{0};   // ranks in the aggregation (row) communicator
    int      rowrank{0};   // our logical (lexicographic) index within the row
    uint64_t lsites{0};    // sites per rank -- invariant under the permutation
    uint64_t chunk{0};     // sites in one globally contiguous run owned by the row
    std::unique_ptr<CartesianCommunicator> rowcomm;
    // counts and displacements are indexed by rank within rowcomm
    std::vector<int>      sendcounts, senddispls, recvcounts, recvdispls;
    std::vector<uint64_t> scatter;     // recv slot -> slot in the aggregated buffer
    std::vector<uint64_t> extentGsite; // global lex site index of extent start
    std::vector<uint64_t> extentLocal; // offset of extent within aggregated buffer
    std::vector<uint64_t> extentSites; // sites in this extent
  };

  static inline void BuildAggregationPlan(GridBase *grid,uint64_t fobjSize,AggregationPlan &p)
  {
    int ndim           = grid->Dimensions();
    Coordinate psizes  = grid->ProcessorGrid();
    Coordinate pcoor   = grid->ThisProcessorCoor();
    Coordinate gLattice= grid->GlobalDimensions();
    Coordinate lLattice= grid->LocalDimensions();
    Coordinate lstart  = grid->LocalStarts();

    uint64_t lsites = grid->lSites();
    p.lsites = lsites;

    //////////////////////////////////////////////////////////////////////////
    // Un-splitting dims 0..k-1 gives the row a contiguous run of
    //     chunk(k) = prod_{d<k} gLattice[d] * lLattice[k]
    // sites, and each rank writes extents of min(chunk,lsites). Take the
    // smallest k that reaches the target so we disturb as few dimensions --
    // and move as little data -- as possible.
    //////////////////////////////////////////////////////////////////////////
    int k = ndim-1;
    for(int trial=1; trial<ndim; trial++){
      uint64_t chunk = lLattice[trial];
      for(int d=0; d<trial; d++) chunk *= gLattice[d];
      if ( std::min(chunk,lsites)*fobjSize >= aggregateTargetBytes ) { k = trial; break; }
    }
    p.nunsplit = k;

    //////////////////////////////////////////////////////////////////////////
    // The box the row collectively owns, expressed in global coordinates.
    // Restricting the global lexicographic order to this box preserves the
    // ordering, so the row index below is monotone in the global index.
    //////////////////////////////////////////////////////////////////////////
    Coordinate B(ndim), S(ndim);
    for(int d=0; d<ndim; d++){
      if ( d<k ) { B[d] = gLattice[d]; S[d] = 0;         }
      else       { B[d] = lLattice[d]; S[d] = lstart[d]; }
    }

    uint64_t chunk = lLattice[k];
    for(int d=0; d<k; d++) chunk *= gLattice[d];
    p.chunk = chunk;

    //////////////////////////////////////////////////////////////////////////
    // Row communicator: the ranks sharing the process coordinates of the slow
    // (still split) dimensions. This is the sub-division the Cartesian
    // communicator already performs for AllToAll(dim,...), widened from one
    // dimension to the k fastest.
    //////////////////////////////////////////////////////////////////////////
    Coordinate row(ndim,1);
    for(int d=0; d<k; d++) row[d] = psizes[d];
    int srank;
    p.rowcomm.reset(new CartesianCommunicator(row,*grid,srank));
    p.rowsize = p.rowcomm->ProcessorCount();

    //////////////////////////////////////////////////////////////////////////
    // Our logical index in the row is the forward lexicographic index of the
    // un-split process coordinates, so that increasing logical index means
    // increasing global lexicographic position in the file. The communicator
    // numbers its own ranks by the reversed (MPI) convention, so build the map
    // between the two rather than assuming either.
    //////////////////////////////////////////////////////////////////////////
    int64_t logical=0, lstride=1;
    for(int d=0; d<k; d++){ logical += pcoor[d]*lstride; lstride *= psizes[d]; }
    GRID_ASSERT(lstride == (int64_t)p.rowsize);
    p.rowrank = (int)logical;

    std::vector<uint64_t> commOf(p.rowsize,0);
    commOf[p.rowrank] = (uint64_t)p.rowcomm->ThisRank();
    p.rowcomm->GlobalSumVector(&commOf[0],p.rowsize);

    uint64_t mystart = (uint64_t)p.rowrank * lsites;
    uint64_t myend   = mystart + lsites;

    Coordinate lcoor(ndim), bcoor(ndim), gcoor(ndim);

    //////////////////////////////////////////////////////////////////////////
    // Send side. Walking our local sites in local lexicographic order walks
    // the row index monotonically, so the send buffer is iodata untouched and
    // we need only the per destination counts.
    //////////////////////////////////////////////////////////////////////////
    std::vector<int> sendLogical(p.rowsize,0);
    for(uint64_t L=0; L<lsites; L++){
      Lexicographic::CoorFromIndex(lcoor,L,lLattice);
      for(int d=0; d<ndim; d++) bcoor[d] = (d<k) ? (lstart[d]+lcoor[d]) : lcoor[d];
      int64_t ri; Lexicographic::IndexFromCoor(bcoor,ri,B);
      sendLogical[ ri/(int64_t)lsites ]++;
    }
    p.sendcounts.assign(p.rowsize,0);
    p.senddispls.assign(p.rowsize,0);
    { int64_t disp=0;
      for(int d=0; d<p.rowsize; d++){          // send buffer is in logical order
        int c = (int)commOf[d];
        p.sendcounts[c] = sendLogical[d];
        p.senddispls[c] = (int)disp;
        disp += sendLogical[d];
      }
    }

    //////////////////////////////////////////////////////////////////////////
    // Receive side. For each slot of our aggregated range work out which rank
    // of the row owns it. Within one source the slots arrive in increasing row
    // index order, which is the order the source sends them in.
    //////////////////////////////////////////////////////////////////////////
    std::vector<int> recvLogical(p.rowsize,0), recvDisplLogical(p.rowsize,0);
    std::vector<int> source(lsites);
    for(uint64_t pos=0; pos<lsites; pos++){
      Lexicographic::CoorFromIndex(bcoor,(int64_t)(mystart+pos),B);
      int64_t j=0, jstride=1;
      for(int d=0; d<k; d++){ j += (bcoor[d]/lLattice[d])*jstride; jstride *= psizes[d]; }
      source[pos] = (int)j;
      recvLogical[j]++;
    }
    p.recvcounts.assign(p.rowsize,0);
    p.recvdispls.assign(p.rowsize,0);
    { int64_t disp=0;
      for(int s=0; s<p.rowsize; s++){         // recv buffer is in logical order
        int c = (int)commOf[s];
        recvDisplLogical[s] = (int)disp;
        p.recvcounts[c] = recvLogical[s];
        p.recvdispls[c] = (int)disp;
        disp += recvLogical[s];
      }
    }
    p.scatter.resize(lsites);
    {
      std::vector<int> fill(p.rowsize,0);
      for(uint64_t pos=0; pos<lsites; pos++){
        int j = source[pos];
        p.scatter[ recvDisplLogical[j] + fill[j]++ ] = pos;
      }
    }

    //////////////////////////////////////////////////////////////////////////
    // The two sides are derived independently; make them check each other.
    //////////////////////////////////////////////////////////////////////////
    {
      std::vector<uint64_t> sendc(p.rowsize),recvc(p.rowsize);
      for(int c=0;c<p.rowsize;c++) sendc[c]=(uint64_t)p.sendcounts[c];
      p.rowcomm->AllToAll(&sendc[0],&recvc[0],1,sizeof(uint64_t));
      for(int c=0;c<p.rowsize;c++) GRID_ASSERT((int)recvc[c]==p.recvcounts[c]);
    }

    //////////////////////////////////////////////////////////////////////////
    // Decompose our range into globally contiguous file extents.
    //////////////////////////////////////////////////////////////////////////
    for(uint64_t c = mystart/chunk; c <= (myend-1)/chunk; c++){
      uint64_t lo = std::max(mystart, c*chunk);
      uint64_t hi = std::min(myend, (c+1)*chunk);
      Lexicographic::CoorFromIndex(bcoor,(int64_t)(c*chunk),B);
      for(int d=0;d<ndim;d++) gcoor[d] = (d<k) ? bcoor[d] : bcoor[d]+S[d];
      int64_t gbase; Lexicographic::IndexFromCoor(gcoor,gbase,gLattice);
      p.extentGsite.push_back( (uint64_t)gbase + (lo - c*chunk) );
      p.extentLocal.push_back( lo - mystart );
      p.extentSites.push_back( hi - lo );
    }
  }

  static inline void ReportAggregationPlan(GridBase *grid,const AggregationPlan &p,uint64_t fobjSize,const char *what)
  {
    if ( !grid->IsBoss() ) return;
    std::cout << GridLogMessage << "IOobject: aggregate " << what
              << " un-splitting " << p.nunsplit << " fastest dimensions, row of "
              << p.rowsize << " ranks" << std::endl;
    std::cout << GridLogMessage << "IOobject: aggregate " << p.extentSites.size()
              << " extent(s)/rank, first " << p.extentSites[0]*fobjSize/1024./1024. << " MB"
              << " (target " << aggregateTargetBytes/1024./1024. << " MB)" << std::endl;
    std::cout << GridLogMessage << "IOobject: aggregate buffer overhead "
              << p.lsites*fobjSize/1024./1024. << " MB/rank" << std::endl;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Stage timings.  The interesting quantity is the slowest rank, since every
  // stage is followed sooner or later by a synchronisation, so reduce with
  // GlobalMax rather than reporting whatever the boss happened to see.
  ////////////////////////////////////////////////////////////////////////////
  static inline void ReportStages(GridBase *grid,const char *what,
                                  const std::vector<const char *> &names,
                                  std::vector<RealD> &useconds)
  {
    GRID_ASSERT(names.size()==useconds.size());
    for(uint64_t i=0;i<useconds.size();i++) grid->GlobalMax(useconds[i]);
    if ( grid->IsBoss() ) {
      std::cout << GridLogMessage << "IOobject: aggregate " << what << " stages (max over ranks, s):";
      for(uint64_t i=0;i<names.size();i++)
        std::cout << " " << names[i] << " " << useconds[i]/1.0e6;
      std::cout << std::endl;
    }
  }

  template<class fobj>
  static inline void AggregateExchange(GridBase *grid,AggregationPlan &p,std::vector<fobj> &iodata,
                                       std::vector<fobj> &aggregated,int forward)
  {
    uint64_t lsites = p.lsites;
    GridStopWatch talloc,tperm,tcomm;

    talloc.Start();
    std::vector<fobj> tmp(lsites);
    talloc.Stop();

    if ( forward ) { // iodata (local order) -> aggregated (lexicographic order)
      tcomm.Start();
      p.rowcomm->AllToAllV(&iodata[0],p.sendcounts,p.senddispls,
                           &tmp[0],   p.recvcounts,p.recvdispls,sizeof(fobj));
      tcomm.Stop();
      tperm.Start();
      thread_for(s,lsites,{ aggregated[p.scatter[s]] = tmp[s]; });
      tperm.Stop();
    } else {         // aggregated -> iodata, the exact mirror
      tperm.Start();
      thread_for(s,lsites,{ tmp[s] = aggregated[p.scatter[s]]; });
      tperm.Stop();
      tcomm.Start();
      p.rowcomm->AllToAllV(&tmp[0],   p.recvcounts,p.recvdispls,
                           &iodata[0],p.sendcounts,p.senddispls,sizeof(fobj));
      tcomm.Stop();
    }

    std::vector<RealD> us = { (RealD)talloc.useconds(), (RealD)tperm.useconds(), (RealD)tcomm.useconds() };
    ReportStages(grid,forward?"exchange (write)":"exchange (read)",
                 {"alloc","permute","alltoallv"},us);
  }

  template<class fobj>
  static inline void AggregateWrite(GridBase *grid,AggregationPlan &p,std::vector<fobj> &aggregated,
                                    std::string file,uint64_t offset)
  {
    //////////////////////////////////////////////////////////////////////////
    // All ranks write concurrently into a shared file, so the file must exist
    // before any of them open it for update, but it does NOT have to be the
    // right length first: the extents tile the record exactly, so writing them
    // extends a short file to precisely offset+payload.
    //
    // Records are created in sequence, so this payload ends the file: the
    // length must end up precisely offset+payload. Anything beyond is left
    // over from whatever the file previously held and must not survive -- a
    // shorter new record written over a longer old one would otherwise leave
    // a trailing fragment of the previous contents masquerading as data.
    // That is the only case needing a truncate, so stat first and truncate
    // afterwards only when the size actually came out wrong. Measured on
    // Frontier, an unconditional truncate up front cost 0.22 to 5.4 s per
    // record -- 15 to 25% of a 19 GB write and 100% of a small one -- while
    // create, open and close together cost a few milliseconds. It is per
    // record, so multi record files do not amortise it away.
    //
    // ::truncate is used because the C++ standard library cannot express this.
    // std::filebuf has no length operation at all; ios::trunc only truncates to
    // zero; seeking past the end and writing a byte can grow a file but never
    // shrink one; and there is no portable way to recover a descriptor from a
    // stream in order to call ftruncate. C++17 does finally offer
    // std::filesystem::resize_file, but that would be Grid's first <filesystem>
    // dependency and needs -lstdc++fs on the older toolchains still in use.
    //////////////////////////////////////////////////////////////////////////
    GridStopWatch tcreate,ttrunc,tbar,topen,twrite,tclose,tskew;
    uint64_t need = offset + (uint64_t)grid->_gsites*sizeof(fobj);

    tcreate.Start();
    if ( grid->IsBoss() ) {
      // opening for update needs the file to exist; create one only if not
      std::fstream probe(file,std::ios::binary|std::ios::out|std::ios::in);
      if ( !probe.is_open() ) {
        std::ofstream create(file,std::ios::binary|std::ios::out);
        create.close();
      }
    }
    tcreate.Stop();

    tbar.Start();
    grid->Barrier();
    tbar.Stop();

    std::ofstream fout;
    fout.exceptions( std::fstream::failbit | std::fstream::badbit );
    try {
      topen.Start();
      fout.open(file,std::ios::binary|std::ios::out|std::ios::in);
      topen.Stop();
      twrite.Start();
      for(uint64_t e=0;e<p.extentSites.size();e++){
        fout.seekp(offset + p.extentGsite[e]*sizeof(fobj));
        fout.write((char *)&aggregated[p.extentLocal[e]],p.extentSites[e]*sizeof(fobj));
      }
      twrite.Stop();
      tclose.Start();
      fout.close();   // flushes the stream buffer; does not force writeback
      tclose.Stop();
    } catch (const std::fstream::failure& exc) {
      std::cout << GridLogError << "Error in aggregate write to " << file << std::endl;
      std::cout << GridLogError << "Exception description: " << exc.what() << std::endl;
      GridAbort();
    }

    ////////////////////////////////////////////////////////////////////////
    // Timed apart from the truncate that follows it.  seek+write above is the
    // slowest rank; this barrier is what the fastest rank then waits, so the
    // pair separates the write cost from the spread across ranks. Folding it
    // into the truncate makes a millisecond stat look like a second.
    ////////////////////////////////////////////////////////////////////////
    tskew.Start();
    grid->Barrier();                 // every extent must be on its way first
    tskew.Stop();

    ttrunc.Start();
    if ( grid->IsBoss() ) {
      struct stat sb;
      int ierr = ::stat(file.c_str(),&sb);
      GRID_ASSERT(ierr==0);
      if ( (uint64_t)sb.st_size != need ) {   // only when a longer record preceded us
        ierr = ::truncate(file.c_str(),(off_t)need);
        GRID_ASSERT(ierr==0);
      }
    }
    grid->Barrier();
    ttrunc.Stop();

    std::vector<RealD> us = { (RealD)tcreate.useconds(), (RealD)tbar.useconds(),
                              (RealD)topen.useconds(),   (RealD)twrite.useconds(),
                              (RealD)tclose.useconds(),  (RealD)tskew.useconds(),
                              (RealD)ttrunc.useconds() };
    ReportStages(grid,"write",{"create","barrier","open","seek+write","close","skew","stat+truncate"},us);
  }

  template<class fobj>
  static inline void AggregateRead(GridBase *grid,AggregationPlan &p,std::vector<fobj> &aggregated,
                                   std::string file,uint64_t offset)
  {
    GridStopWatch topen,tread,tclose;
    std::ifstream fin;
    topen.Start();
    fin.open(file,std::ios::binary|std::ios::in);
    topen.Stop();
    tread.Start();
    for(uint64_t e=0;e<p.extentSites.size();e++){
      fin.seekg(offset + p.extentGsite[e]*sizeof(fobj));
      fin.read((char *)&aggregated[p.extentLocal[e]],p.extentSites[e]*sizeof(fobj));
      GRID_ASSERT(fin.fail()==0);
    }
    tread.Stop();
    tclose.Start();
    fin.close();
    tclose.Stop();

    std::vector<RealD> us = { (RealD)topen.useconds(), (RealD)tread.useconds(), (RealD)tclose.useconds() };
    ReportStages(grid,"read",{"open","seek+read","close"},us);
  }
#endif

  template<class word,class fobj>
  static inline void IOobject(word w,
			      GridBase *grid,
			      std::vector<fobj> &iodata,
			      std::string file,
			      uint64_t& offset,
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

    Coordinate  psizes = grid->ProcessorGrid(); 
    Coordinate  pcoor  = grid->ThisProcessorCoor();
    Coordinate gLattice= grid->GlobalDimensions();
    Coordinate lLattice= grid->LocalDimensions();

    Coordinate lStart(ndim);
    Coordinate gStart(ndim);

    // Flatten the file
    uint64_t lsites = grid->lSites();
    if ( control & BINARYIO_MASTER_APPEND )  {
      GRID_ASSERT(iodata.size()==1);
    } else {
      GRID_ASSERT(lsites==iodata.size());
    }
    for(int d=0;d<ndim;d++){
      gStart[d] = lLattice[d]*pcoor[d];
      lStart[d] = 0;
    }

    //////////////////////////////////////////////////////////////////////////////
    // Aggregate the lexicographic layout onto contiguous per rank extents
    // ourselves rather than leaving it to MPI-IO collective buffering
    //////////////////////////////////////////////////////////////////////////////
    int aggregate = (control & BINARYIO_AGGREGATE)
                 && (control & BINARYIO_LEXICOGRAPHIC)
                 && !(control & BINARYIO_MASTER_APPEND)
                 && (nrank > 1);
#ifndef USE_MPI_IO
    GRID_ASSERT(aggregate==0); // BINARYIO_AGGREGATE requires MPI
#endif

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
    ierr = MPI_Type_contiguous(numword,mpiword,&mpiObject);    GRID_ASSERT(ierr==0);
    ierr = MPI_Type_commit(&mpiObject);

    // The subarray view is what aggregation exists to avoid; do not build it
    if ( !aggregate ) {
    //////////////////////////////////////////////////////////////////////////////
    // File global array data type
    //////////////////////////////////////////////////////////////////////////////
    ierr=MPI_Type_create_subarray(ndim,&gLattice[0],&lLattice[0],&gStart[0],MPI_ORDER_FORTRAN, mpiObject,&fileArray);    GRID_ASSERT(ierr==0);
    ierr=MPI_Type_commit(&fileArray);    GRID_ASSERT(ierr==0);

    //////////////////////////////////////////////////////////////////////////////
    // local lattice array
    //////////////////////////////////////////////////////////////////////////////
    ierr=MPI_Type_create_subarray(ndim,&lLattice[0],&lLattice[0],&lStart[0],MPI_ORDER_FORTRAN, mpiObject,&localArray);    GRID_ASSERT(ierr==0);
    ierr=MPI_Type_commit(&localArray);    GRID_ASSERT(ierr==0);
    }
#endif

    //////////////////////////////////////////////////////////////////////////////
    // Byte order
    //////////////////////////////////////////////////////////////////////////////
    int ieee32big = (format == std::string("IEEE32BIG"));
    int ieee32    = (format == std::string("IEEE32"));
    int ieee64big = (format == std::string("IEEE64BIG"));
    int ieee64    = (format == std::string("IEEE64") || format == std::string("IEEE64LITTLE"));
    GRID_ASSERT(ieee64||ieee32|ieee64big||ieee32big);
    GRID_ASSERT((ieee64+ieee32+ieee64big+ieee32big)==1);
    //////////////////////////////////////////////////////////////////////////////
    // Do the I/O
    //////////////////////////////////////////////////////////////////////////////
    if ( control & BINARYIO_READ ) { 

      timer.Start();

      if ( aggregate ) {
#ifdef USE_MPI_IO
        std::cout<< GridLogMessage<<"IOobject: aggregate read I/O "<< file<< std::endl;
        AggregationPlan plan;
        BuildAggregationPlan(grid,sizeof(fobj),plan);
        ReportAggregationPlan(grid,plan,sizeof(fobj),"read");
        std::vector<fobj> aggregated(lsites);
        AggregateRead(grid,plan,aggregated,file,offset);
        AggregateExchange(grid,plan,iodata,aggregated,0);
#else
        GRID_ASSERT(0);
#endif
      } else if ( (control & BINARYIO_LEXICOGRAPHIC) && (nrank > 1) ) {
#ifdef USE_MPI_IO
	std::cout<< GridLogMessage<<"IOobject: MPI read I/O "<< file<< std::endl;
	ierr=MPI_File_open(grid->communicator,(char *) file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);    GRID_ASSERT(ierr==0);
	ierr=MPI_File_set_view(fh, disp, mpiObject, fileArray, "native", MPI_INFO_NULL);    GRID_ASSERT(ierr==0);
	ierr=MPI_File_read_all(fh, &iodata[0], 1, localArray, &status);    GRID_ASSERT(ierr==0);
	MPI_File_close(&fh);
	MPI_Type_free(&fileArray);
	MPI_Type_free(&localArray);
#else 
	assert(0);
#endif
      } else {
	std::cout << GridLogMessage <<"IOobject: C++ read I/O " << file << " : "
                  << iodata.size() * sizeof(fobj) << " bytes and offset " << offset << std::endl;
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
        GRID_ASSERT(fin.fail() == 0);
        fin.close();
      }
      
      grid->Barrier();

	  timer.Stop();

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
      if ( aggregate ) {
#ifdef USE_MPI_IO
        std::cout << GridLogMessage <<"IOobject: aggregate write I/O " << file << std::endl;
        AggregationPlan plan;
        BuildAggregationPlan(grid,sizeof(fobj),plan);
        ReportAggregationPlan(grid,plan,sizeof(fobj),"write");
        std::vector<fobj> aggregated(lsites);
        AggregateExchange(grid,plan,iodata,aggregated,1);
        AggregateWrite(grid,plan,aggregated,file,offset);
        ////////////////////////////////////////////////////////////////////////
        // Not every rank ends at the end of the payload, so the position can
        // not be recovered from a file handle. Callers (Lime record chaining)
        // rely on this being the first byte past the record.
        ////////////////////////////////////////////////////////////////////////
        offset = offset + (uint64_t)grid->_gsites*sizeof(fobj);
#else
        GRID_ASSERT(0);
#endif
      } else if ( (control & BINARYIO_LEXICOGRAPHIC) && (nrank > 1) ) {
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

        std::cout << GridLogDebug << "MPI write I/O set view " << file << std::endl;
        ierr = MPI_File_set_view(fh, disp, mpiObject, fileArray, "native", MPI_INFO_NULL);
        GRID_ASSERT(ierr == 0);

        std::cout << GridLogDebug << "MPI write I/O write all " << file << std::endl;
        ierr = MPI_File_write_all(fh, &iodata[0], 1, localArray, &status);
        GRID_ASSERT(ierr == 0);

        MPI_Offset os;
        MPI_File_get_position(fh, &os);
        MPI_File_get_byte_offset(fh, os, &disp);
        offset = disp;


        MPI_File_close(&fh);
        MPI_Type_free(&fileArray);
        MPI_Type_free(&localArray);
#else 
	assert(0);
#endif
      } else { 

        std::cout << GridLogMessage << "IOobject: C++ write I/O " << file << " : "
                  << iodata.size() * sizeof(fobj) << " bytes and offset " << offset << std::endl;
        
	std::ofstream fout; 
	fout.exceptions ( std::fstream::failbit | std::fstream::badbit );

        ////////////////////////////////////////////////////////////////////
        // Grid's model is that the boss rank performs the metadata
        // operations and every other rank only seeks and writes into a file
        // that already exists. Opening with ios::out on all ranks broke that:
        // it is O_TRUNC, so a rank opening late truncated the file back to
        // zero after an earlier rank had written its segment, leaving a hole
        // in its place. The barriers around this block are outside it and do
        // not order the opens against the writes. Let the boss create and
        // empty the file, then everyone opens for update only. Same resulting
        // length, one metadata operation instead of one per rank, no race.
        ////////////////////////////////////////////////////////////////////
        if ( !offset && grid->IsBoss() ) { // offset zero: this record starts the file
          std::ofstream create(file,std::ios::binary|std::ios::out);
          create.close();
        }
        grid->Barrier();

	try {
	  fout.open(file,std::ios::binary|std::ios::out|std::ios::in);
	} catch (const std::fstream::failure& exc) {
	  std::cout << GridLogError << "Error in opening the file " << file << " for output" <<std::endl;
	  std::cout << GridLogError << "Exception description: " << exc.what() << std::endl;
	  //	  std::cout << GridLogError << "Probable cause: wrong path, inaccessible location "<< std::endl;
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
  offset  = fout.tellp();
	fout.close();
      }
      grid->Barrier();
      timer.Stop();
    }
    
    lastPerf.size            = sizeof(fobj)*iodata.size()*nrank;
    lastPerf.time            = timer.useconds();
    lastPerf.mbytesPerSecond = lastPerf.size/1024./1024./(lastPerf.time/1.0e6);
    std::cout<<GridLogMessage<<"IOobject: ";
    if ( control & BINARYIO_READ) std::cout << " read  ";
    else                          std::cout << " write ";
    uint64_t bytes = sizeof(fobj)*iodata.size()*nrank;
    std::cout<< lastPerf.size <<" bytes in "<< timer.Elapsed() <<" "
	     << lastPerf.mbytesPerSecond <<" MB/s "<<std::endl;

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
				       uint64_t offset,
				       const std::string &format,
				       uint32_t &nersc_csum,
				       uint32_t &scidac_csuma,
				       uint32_t &scidac_csumb,
				       int control=DefaultControl()
				       )
  {
    typedef typename vobj::scalar_object sobj;
    typedef typename vobj::Realified::scalar_type word;    word w=0;

    GridBase *grid = Umu.Grid();
    uint64_t lsites = grid->lSites();

    std::vector<sobj> scalardata(lsites); 
    std::vector<fobj>     iodata(lsites); // Munge, checksum, byte order in here
    
    IOobject(w,grid,iodata,file,offset,format,BINARYIO_READ|control,
	     nersc_csum,scidac_csuma,scidac_csumb);

    GridStopWatch timer; 
    timer.Start();

    thread_for(x,lsites, { munge(iodata[x], scalardata[x]); });

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
					  uint64_t offset,
					  const std::string &format,
					  uint32_t &nersc_csum,
					  uint32_t &scidac_csuma,
					  uint32_t &scidac_csumb,
					  int control=DefaultControl())
  {
    typedef typename vobj::scalar_object sobj;
    typedef typename vobj::Realified::scalar_type word;    word w=0;
    GridBase *grid = Umu.Grid();
    uint64_t lsites = grid->lSites(), offsetCopy = offset;
    int attemptsLeft = std::max(0, BinaryIO::latticeWriteMaxRetry);
    bool checkWrite = (BinaryIO::latticeWriteMaxRetry >= 0);

    std::vector<sobj> scalardata(lsites); 
    std::vector<fobj>     iodata(lsites); // Munge, checksum, byte order in here

    //////////////////////////////////////////////////////////////////////////////
    // Munge [ .e.g 3rd row recon ]
    //////////////////////////////////////////////////////////////////////////////
    GridStopWatch timer; timer.Start();
    unvectorizeToLexOrdArray(scalardata,Umu);    

    thread_for(x, lsites, { munge(scalardata[x],iodata[x]); });

    grid->Barrier();
    timer.Stop();
    while (attemptsLeft >= 0)
    {
      grid->Barrier();
      IOobject(w,grid,iodata,file,offset,format,BINARYIO_WRITE|control,
	             nersc_csum,scidac_csuma,scidac_csumb);
      if (checkWrite)
      {
        std::vector<fobj> ckiodata(lsites);
        uint32_t          cknersc_csum, ckscidac_csuma, ckscidac_csumb;
        uint64_t          ckoffset = offsetCopy;

        std::cout << GridLogMessage << "writeLatticeObject: read back object" << std::endl;
        grid->Barrier();
        IOobject(w,grid,ckiodata,file,ckoffset,format,BINARYIO_READ|control,
	               cknersc_csum,ckscidac_csuma,ckscidac_csumb);
        if ((cknersc_csum != nersc_csum) or (ckscidac_csuma != scidac_csuma) or (ckscidac_csumb != scidac_csumb))
        {
          std::cout << GridLogMessage << "writeLatticeObject: read test checksum failure, re-writing (" << attemptsLeft << " attempt(s) remaining)" << std::endl;
          offset = offsetCopy;
          thread_for(x,lsites, { munge(scalardata[x],iodata[x]); });
        }
        else
        {
          std::cout << GridLogMessage << "writeLatticeObject: read test checksum correct" << std::endl;
          break;
        }
      }
      attemptsLeft--;
    }
    

    std::cout<<GridLogMessage<<"writeLatticeObject: unvectorize overhead "<<timer.Elapsed()  <<std::endl;
  }
  
  /////////////////////////////////////////////////////////////////////////////
  // Read a RNG;  use IOobject and lexico map to an array of state 
  //////////////////////////////////////////////////////////////////////////////////////
  static inline void readRNG(GridSerialRNG &serial_rng,
			     GridParallelRNG &parallel_rng,
			     std::string file,
			     uint64_t offset,
			     uint32_t &nersc_csum,
			     uint32_t &scidac_csuma,
			     uint32_t &scidac_csumb)
  {
    typedef typename GridSerialRNG::RngStateType RngStateType;
    const int RngStateCount = GridSerialRNG::RngStateCount;
    typedef std::array<RngStateType,RngStateCount> RNGstate;
    typedef RngStateType word;    word w=0;

    std::string format = "IEEE32BIG";

    GridBase *grid = parallel_rng.Grid();
    uint64_t gsites = grid->gSites();
    uint64_t lsites = grid->lSites();

    uint32_t nersc_csum_tmp   = 0;
    uint32_t scidac_csuma_tmp = 0;
    uint32_t scidac_csumb_tmp = 0;

    GridStopWatch timer;

    std::cout << GridLogMessage << "RNG read I/O on file " << file << std::endl;

    std::vector<RNGstate> iodata(lsites);
    IOobject(w,grid,iodata,file,offset,format,BINARYIO_READ|DefaultControl(),
	     nersc_csum,scidac_csuma,scidac_csumb);

    timer.Start();
    thread_for(lidx,lsites,{  // FIX ME, suboptimal implementation
      std::vector<RngStateType> tmp(RngStateCount);
      std::copy(iodata[lidx].begin(),iodata[lidx].end(),tmp.begin());
      Coordinate lcoor;
      grid->LocalIndexToLocalCoor(lidx, lcoor);
      int o_idx=grid->oIndex(lcoor);
      int i_idx=grid->iIndex(lcoor);
      int gidx=parallel_rng.generator_idx(o_idx,i_idx);
      parallel_rng.SetState(tmp,gidx);
      });
    timer.Stop();

    iodata.resize(1);
    IOobject(w,grid,iodata,file,offset,format,BINARYIO_READ|BINARYIO_MASTER_APPEND,
	     nersc_csum_tmp,scidac_csuma_tmp,scidac_csumb_tmp);

    {
      std::vector<RngStateType> tmp(RngStateCount);
      std::copy(iodata[0].begin(),iodata[0].end(),tmp.begin());
      serial_rng.SetState(tmp,0);
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
  static inline void writeRNG(GridSerialRNG &serial_rng,
			      GridParallelRNG &parallel_rng,
			      std::string file,
			      uint64_t offset,
			      uint32_t &nersc_csum,
			      uint32_t &scidac_csuma,
			      uint32_t &scidac_csumb)
  {
    typedef typename GridSerialRNG::RngStateType RngStateType;
    typedef RngStateType word; word w=0;
    const int RngStateCount = GridSerialRNG::RngStateCount;
    typedef std::array<RngStateType,RngStateCount> RNGstate;

    GridBase *grid = parallel_rng.Grid();
    uint64_t gsites = grid->gSites();
    uint64_t lsites = grid->lSites();

    uint32_t nersc_csum_tmp;
    uint32_t scidac_csuma_tmp;
    uint32_t scidac_csumb_tmp;

    GridStopWatch timer;
    std::string format = "IEEE32BIG";

    std::cout << GridLogMessage << "RNG write I/O on file " << file << std::endl;

    timer.Start();
    std::vector<RNGstate> iodata(lsites);
    thread_for(lidx,lsites,{
      std::vector<RngStateType> tmp(RngStateCount);
      Coordinate lcoor;
      grid->LocalIndexToLocalCoor(lidx, lcoor);
      int o_idx=grid->oIndex(lcoor);
      int i_idx=grid->iIndex(lcoor);
      int gidx=parallel_rng.generator_idx(o_idx,i_idx);
      parallel_rng.GetState(tmp,gidx);
      std::copy(tmp.begin(),tmp.end(),iodata[lidx].begin());
    });
    timer.Stop();

    IOobject(w,grid,iodata,file,offset,format,BINARYIO_WRITE|DefaultControl(),
	     nersc_csum,scidac_csuma,scidac_csumb);
    iodata.resize(1);
    {
      std::vector<RngStateType> tmp(RngStateCount);
      serial_rng.GetState(tmp,0);
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

NAMESPACE_END(Grid);
