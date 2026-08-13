/**************************************************************************
 * MPI only parallel file I/O reproducer -- no Grid, no accelerator.
 *
 * Writes an Nd lattice of NWORD-double site objects to a single shared file
 * in *global lexicographic* order (dimension 0 fastest), which is the order
 * that makes a file independent of the process decomposition that wrote it.
 * Three paths produce that file, and are timed and validated against each
 * other:
 *
 *   raw        each rank writes one disjoint contiguous segment with pwrite.
 *              NOT lexicographic -- the file depends on the decomposition.
 *              This is the zero-transposition control: the bandwidth the
 *              filesystem and the clients can sustain with no reordering.
 *
 *   mpiio      MPI_Type_create_subarray + MPI_File_set_view +
 *              MPI_File_write_all / MPI_File_read_all.  The canonical
 *              two-phase collective.  This is the path under test.
 *
 *   aggregate  user level two-phase: MPI_Alltoallv within a row
 *              subcommunicator to give each rank a few large contiguous
 *              file extents, then independent pwrite / pread of those.
 *              Produces a byte identical file to mpiio.
 *
 * Validation.  Site content is a deterministic function of the *global*
 * lexicographic site index, so a byte that ends up in the wrong place is
 * detectable, not merely a checksum difference.  Three independent checks:
 *
 *   1. crc32 per site, combined across sites with rotate-and-xor keyed by
 *      the global index (so the combination is order independent) and
 *      reduced with MPI_BXOR.  Proves the right sites with the right
 *      contents are present.  Does NOT prove they are in the right place.
 *   2. On read, every rank recomputes the expected content of each global
 *      site it believes it owns and memcmps.  Proves placement, in
 *      parallel, at any volume.  Reports the first offending global site.
 *   3. Cross validation: write with mpiio, read with aggregate, and vice
 *      versa.  Both paths must agree on the same file, so a shared
 *      misunderstanding of the layout cannot hide.
 *   4. Optionally (--serial-crc) rank 0 streams the whole file and crc32s
 *      it.  Slow, incontrovertible, and the two lexicographic files must
 *      give the same value.
 *
 * Build:   mpicxx -O2 -std=c++11 io_mpi.cc -o io_mpi
 * Run:     mpirun -n 32 ./io_mpi --grid 32.32.64.128 --mpi 4.4.2.1
 **************************************************************************/

#include <algorithm>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>

#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>

#include <mpi.h>

/**************************************************************
 * Globals
 **************************************************************/
static MPI_Comm WorldComm;
static int WorldRank, WorldSize;
static int ShmRank,   ShmSize;

static int      Nd       = 4;
static int64_t  NWORD    = 72;   /* doubles per site; 72 -> 576 B, a LorentzColourMatrixD */
static int64_t  fobjSize = 72*sizeof(double);

#define CHECK(ierr) do {                                                \
    int _e = (ierr);                                                    \
    if ( _e != MPI_SUCCESS ) {                                          \
      char _s[MPI_MAX_ERROR_STRING]; int _l;                            \
      MPI_Error_string(_e,_s,&_l);                                      \
      printf("rank %d: MPI error at %s:%d : %s\n",WorldRank,__FILE__,__LINE__,_s); \
      fflush(stdout); MPI_Abort(WorldComm,1);                           \
    } } while(0)

#define SYSCHECK(cond,what) do {                                        \
    if ( !(cond) ) {                                                    \
      printf("rank %d: %s failed at %s:%d : %s\n",WorldRank,what,__FILE__,__LINE__,strerror(errno)); \
      fflush(stdout); MPI_Abort(WorldComm,1);                           \
    } } while(0)

static inline double usecond(void) {
  struct timeval tv; gettimeofday(&tv,NULL);
  return 1.0e6*tv.tv_sec + 1.0*tv.tv_usec;
}

/**************************************************************
 * crc32 (IEEE 802.3, reflected), table built once.  Deliberately
 * the slow obvious implementation: this is a correctness tool and
 * we would rather it be readable than fast.
 **************************************************************/
static uint32_t crc_table[256];
static void crc32_init(void)
{
  for(uint32_t n=0;n<256;n++){
    uint32_t c = n;
    for(int k=0;k<8;k++) c = (c&1) ? (0xEDB88320u ^ (c>>1)) : (c>>1);
    crc_table[n] = c;
  }
}
static inline uint32_t crc32_buf(const void *buf,size_t len,uint32_t crc)
{
  const unsigned char *p = (const unsigned char *)buf;
  crc = crc ^ 0xFFFFFFFFu;
  for(size_t i=0;i<len;i++) crc = crc_table[(crc ^ p[i]) & 0xFF] ^ (crc>>8);
  return crc ^ 0xFFFFFFFFu;
}
static inline uint32_t rotl32(uint32_t x,int n){ n &= 31; return n ? ((x<<n)|(x>>(32-n))) : x; }

/**************************************************************
 * Lexicographic index <-> coordinate, dimension 0 fastest.
 * This is MPI_ORDER_FORTRAN, and it is the file order.
 **************************************************************/
static inline void CoorFromIndex(int64_t *coor,int64_t idx,const int64_t *dims,int nd)
{
  for(int d=0;d<nd;d++){ coor[d] = idx % dims[d]; idx /= dims[d]; }
}
static inline int64_t IndexFromCoor(const int64_t *coor,const int64_t *dims,int nd)
{
  int64_t idx=0, stride=1;
  for(int d=0;d<nd;d++){ idx += coor[d]*stride; stride *= dims[d]; }
  return idx;
}

/**************************************************************
 * Site content as a pure function of the GLOBAL site index, so a
 * misplaced site is detectable and identifiable, not just a
 * checksum difference.  splitmix64 mantissa -> exactly
 * representable doubles, so byte comparison is well defined.
 **************************************************************/
static inline void FillSite(double *s,uint64_t g)
{
  for(int64_t w=0;w<NWORD;w++){
    uint64_t x = g*0x9E3779B97F4A7C15ull + (uint64_t)w*0xBF58476D1CE4E5B9ull;
    x ^= x>>30; x *= 0xBF58476D1CE4E5B9ull;
    x ^= x>>27; x *= 0x94D049BB133111EBull;
    x ^= x>>31;
    s[w] = (double)(int64_t)(x>>11) * (1.0/9007199254740992.0);
  }
}

/**************************************************************
 * Layout.  Everything the three paths need to agree on.
 **************************************************************/
struct Layout {
  std::vector<int64_t> gLattice, lLattice, psizes, pcoor, gStart;
  int64_t lsites, gsites;
  uint64_t rawOffset;      /* this rank's byte offset for the raw control */
};

static void BuildLayout(Layout &L,const std::vector<int64_t> &g,const std::vector<int64_t> &p)
{
  L.gLattice = g; L.psizes = p;
  L.lLattice.resize(Nd); L.pcoor.resize(Nd); L.gStart.resize(Nd);

  /* process coordinate: dimension 0 fastest, matching the lattice order.
     Note we do NOT use MPI_Cart_create -- it ranks with the last dimension
     fastest, and mixing the two conventions is the classic way to get a
     subtly transposed file. */
  int64_t r = WorldRank;
  for(int d=0;d<Nd;d++){ L.pcoor[d] = r % p[d]; r /= p[d]; }

  L.lsites = 1; L.gsites = 1;
  for(int d=0;d<Nd;d++){
    if ( g[d] % p[d] ) {
      if(!WorldRank) printf("global dimension %d (%lld) not divisible by %lld\n",
                            d,(long long)g[d],(long long)p[d]);
      MPI_Abort(WorldComm,1);
    }
    L.lLattice[d] = g[d]/p[d];
    L.gStart[d]   = L.pcoor[d]*L.lLattice[d];
    L.lsites *= L.lLattice[d];
    L.gsites *= g[d];
  }
  L.rawOffset = (uint64_t)WorldRank * (uint64_t)L.lsites * (uint64_t)fobjSize;
}

/**************************************************************
 * Aggregation plan.
 *
 * Un-splitting the k fastest dimensions of the process grid gives the
 * resulting row of ranks a globally contiguous run of
 *      chunk(k) = prod_{d<k} gLattice[d] * lLattice[k]
 * sites.  Choose the smallest k that reaches the target extent, so we
 * disturb as few dimensions -- and move as little data -- as possible.
 *
 * The row's box B, restricted from the global lattice, preserves the
 * global lexicographic order.  Walking our local sites in local order
 * therefore walks the row index strictly monotonically:
 *
 *   incrementing local dim 0 moves ri by +1;
 *   a carry out of dim d moves ri by B[d] - lLattice[d] + 1 > 0, since
 *   B[d] >= lLattice[d] for every d.
 *
 * so the Alltoallv send buffer is the untouched local array and only the
 * per-destination counts are needed.  That is the whole trick.
 **************************************************************/
struct AggregationPlan {
  int      nunsplit, rowsize, rowrank;
  int64_t  lsites, chunk;
  MPI_Comm rowcomm;
  std::vector<int>      sendcounts, senddispls, recvcounts, recvdispls;  /* in sites */
  std::vector<int64_t>  scatter;
  std::vector<int64_t>  extentGsite, extentLocal, extentSites;
};

static void BuildAggregationPlan(const Layout &L,uint64_t targetBytes,AggregationPlan &p)
{
  const int64_t *g = &L.gLattice[0], *l = &L.lLattice[0], *ps = &L.psizes[0];
  int64_t lsites = L.lsites;
  p.lsites = lsites;

  int k = Nd-1;
  for(int trial=1; trial<Nd; trial++){
    int64_t chunk = l[trial];
    for(int d=0;d<trial;d++) chunk *= g[d];
    if ( (uint64_t)(std::min(chunk,lsites)*fobjSize) >= targetBytes ) { k = trial; break; }
  }
  p.nunsplit = k;

  std::vector<int64_t> B(Nd), S(Nd);
  for(int d=0;d<Nd;d++){
    if ( d<k ) { B[d] = g[d]; S[d] = 0;          }
    else       { B[d] = l[d]; S[d] = L.gStart[d]; }
  }
  int64_t chunk = l[k];
  for(int d=0;d<k;d++) chunk *= g[d];
  p.chunk = chunk;

  /* Row communicator.  Colour = the still-split (slow) process coordinates;
     key = the lexicographic index of the un-split (fast) ones, so that
     MPI_Comm_split hands back rowrank == logical index and the send/recv
     count arrays can be indexed directly.  Grid has to build an explicit
     map here because its Cartesian communicator numbers ranks the other
     way round; a standalone reproducer does not have to inherit that. */
  int colour=0, cstride=1;
  for(int d=k;d<Nd;d++){ colour += (int)(L.pcoor[d]*cstride); cstride *= (int)ps[d]; }
  int logical=0, lstride=1;
  for(int d=0;d<k;d++){ logical += (int)(L.pcoor[d]*lstride); lstride *= (int)ps[d]; }

  CHECK(MPI_Comm_split(WorldComm,colour,logical,&p.rowcomm));
  MPI_Comm_size(p.rowcomm,&p.rowsize);
  MPI_Comm_rank(p.rowcomm,&p.rowrank);
  if ( p.rowsize != lstride || p.rowrank != logical ) {
    printf("rank %d: row communicator inconsistent (%d vs %d, %d vs %d)\n",
           WorldRank,p.rowsize,lstride,p.rowrank,logical);
    fflush(stdout); MPI_Abort(WorldComm,1);
  }

  int64_t mystart = (int64_t)p.rowrank * lsites;
  int64_t myend   = mystart + lsites;

  std::vector<int64_t> lcoor(Nd), bcoor(Nd), gcoor(Nd);

  /* send side: counts only, buffer stays in local order (monotonicity above) */
  p.sendcounts.assign(p.rowsize,0);
  p.senddispls.assign(p.rowsize,0);
  for(int64_t Lx=0; Lx<lsites; Lx++){
    CoorFromIndex(&lcoor[0],Lx,l,Nd);
    for(int d=0;d<Nd;d++) bcoor[d] = (d<k) ? (L.gStart[d]+lcoor[d]) : lcoor[d];
    int64_t ri = IndexFromCoor(&bcoor[0],&B[0],Nd);
    p.sendcounts[ ri/lsites ]++;
  }
  { int64_t disp=0;
    for(int d=0;d<p.rowsize;d++){ p.senddispls[d] = (int)disp; disp += p.sendcounts[d]; } }

  /* receive side: which rank of the row owns each slot of our aggregated
     range.  Within one source the slots arrive in increasing row index
     order, which is the order that source sends them in. */
  p.recvcounts.assign(p.rowsize,0);
  p.recvdispls.assign(p.rowsize,0);
  std::vector<int> source(lsites);
  for(int64_t pos=0; pos<lsites; pos++){
    CoorFromIndex(&bcoor[0],mystart+pos,&B[0],Nd);
    int64_t j=0, jstride=1;
    for(int d=0;d<k;d++){ j += (bcoor[d]/l[d])*jstride; jstride *= ps[d]; }
    source[pos] = (int)j;
    p.recvcounts[j]++;
  }
  { int64_t disp=0;
    for(int s=0;s<p.rowsize;s++){ p.recvdispls[s] = (int)disp; disp += p.recvcounts[s]; } }

  p.scatter.resize(lsites);
  { std::vector<int> fill(p.rowsize,0);
    for(int64_t pos=0; pos<lsites; pos++){
      int j = source[pos];
      p.scatter[ p.recvdispls[j] + fill[j]++ ] = pos;
    } }

  /* the two sides were derived independently; make them check each other */
  { std::vector<int> back(p.rowsize,0);
    CHECK(MPI_Alltoall(&p.sendcounts[0],1,MPI_INT,&back[0],1,MPI_INT,p.rowcomm));
    for(int c=0;c<p.rowsize;c++){
      if ( back[c] != p.recvcounts[c] ) {
        printf("rank %d: aggregation plan inconsistent at %d: %d vs %d\n",
               WorldRank,c,back[c],p.recvcounts[c]);
        fflush(stdout); MPI_Abort(WorldComm,1);
      }
    } }

  /* decompose our aggregated range into globally contiguous file extents */
  for(int64_t c = mystart/chunk; c <= (myend-1)/chunk; c++){
    int64_t lo = std::max(mystart,  c*chunk);
    int64_t hi = std::min(myend,   (c+1)*chunk);
    CoorFromIndex(&bcoor[0],c*chunk,&B[0],Nd);
    for(int d=0;d<Nd;d++) gcoor[d] = (d<k) ? bcoor[d] : bcoor[d]+S[d];
    int64_t gbase = IndexFromCoor(&gcoor[0],g,Nd);
    p.extentGsite.push_back(gbase + (lo - c*chunk));
    p.extentLocal.push_back(lo - mystart);
    p.extentSites.push_back(hi - lo);
  }
}

static void ReportPlan(const AggregationPlan &p,uint64_t targetBytes)
{
  /* Whether the counts differ between destinations decides which branch of
     MPI_Alltoallv is exercised.  Say so, rather than leaving the reader to
     infer it from the geometry: a sweep that never leaves the uniform case
     has not tested the interesting one. */
  int lo = p.sendcounts[0], hi = p.sendcounts[0];
  for(int c=0;c<p.rowsize;c++){ lo = std::min(lo,p.sendcounts[c]); hi = std::max(hi,p.sendcounts[c]); }
  int uniform = (lo==hi), alluniform;
  CHECK(MPI_Allreduce(&uniform,&alluniform,1,MPI_INT,MPI_MIN,WorldComm));

  if ( WorldRank ) return;
  printf("= aggregate: un-splitting %d fastest dimension(s), row of %d ranks, counts %s\n",
         p.nunsplit,p.rowsize,alluniform?"UNIFORM":"NON-UNIFORM");
  printf("= aggregate: %d extent(s)/rank, first %.1f MB (target %.1f MB)\n",
         (int)p.extentSites.size(),p.extentSites[0]*fobjSize/1048576.0,targetBytes/1048576.0);
  printf("= aggregate: buffer overhead %.1f MB/rank\n",p.lsites*fobjSize/1048576.0);
  fflush(stdout);
}

/* forward: iodata (local order) -> aggregated (lexicographic order) */
static void AggregateExchange(AggregationPlan &p,double *iodata,double *aggregated,
                              MPI_Datatype siteType,int forward,double *t_comm,double *t_perm)
{
  int64_t lsites = p.lsites;
  std::vector<double> tmp(lsites*NWORD);
  double t0,t1,t2;

  if ( forward ) {
    t0 = usecond();
    CHECK(MPI_Alltoallv(iodata, &p.sendcounts[0],&p.senddispls[0],siteType,
                        &tmp[0],&p.recvcounts[0],&p.recvdispls[0],siteType,p.rowcomm));
    t1 = usecond();
    for(int64_t s=0;s<lsites;s++)
      memcpy(aggregated + p.scatter[s]*NWORD, &tmp[s*NWORD], fobjSize);
    t2 = usecond();
    *t_comm += t1-t0; *t_perm += t2-t1;
  } else {
    t0 = usecond();
    for(int64_t s=0;s<lsites;s++)
      memcpy(&tmp[s*NWORD], aggregated + p.scatter[s]*NWORD, fobjSize);
    t1 = usecond();
    CHECK(MPI_Alltoallv(&tmp[0],&p.recvcounts[0],&p.recvdispls[0],siteType,
                        iodata, &p.sendcounts[0],&p.senddispls[0],siteType,p.rowcomm));
    t2 = usecond();
    *t_perm += t1-t0; *t_comm += t2-t1;
  }
}

/**************************************************************
 * Checksums.
 *
 * Order independent site checksum: crc32 of each site rotated by the
 * global index and xor-combined, then MPI_BXOR reduced.  This is the
 * SciDAC csuma/csumb scheme and it proves the right site contents are
 * present; it says nothing about where in the file they ended up.
 * Placement is proved separately, by ValidateLocal below.
 **************************************************************/
struct Checksum { uint32_t a,b; };

static Checksum SiteChecksum(const double *data,int64_t nsite,
                             const int64_t *gindex)   /* gindex[i] = global index of site i */
{
  Checksum c; c.a = 0; c.b = 0;
  for(int64_t i=0;i<nsite;i++){
    uint32_t crc = crc32_buf(data + i*NWORD,fobjSize,0);
    uint64_t g   = (uint64_t)gindex[i];
    c.a ^= rotl32(crc,(int)(g%29));
    c.b ^= rotl32(crc,(int)(g%31));
  }
  uint32_t r[2] = { c.a, c.b }, s[2];
  CHECK(MPI_Allreduce(r,s,2,MPI_UINT32_T,MPI_BXOR,WorldComm));
  c.a = s[0]; c.b = s[1];
  return c;
}

/* every rank recomputes what it should have read and compares.  Parallel,
   affordable at any volume, and it proves placement because each rank
   checks the bytes it found at the file positions it believes it owns. */
static int64_t ValidateLocal(const double *data,int64_t nsite,const int64_t *gindex)
{
  std::vector<double> want(NWORD);
  for(int64_t i=0;i<nsite;i++){
    FillSite(&want[0],(uint64_t)gindex[i]);
    if ( memcmp(&want[0],data + i*NWORD,fobjSize) ) return gindex[i];
  }
  return -1;
}

/* rank 0 streams the file and crc32s it.  Slow; the incontrovertible check. */
static uint32_t SerialFileCrc(const char *file,uint64_t offset,uint64_t bytes)
{
  uint32_t crc = 0;
  if ( !WorldRank ) {
    int fd = open(file,O_RDONLY);
    SYSCHECK(fd>=0,"open for serial crc");
    const size_t bufsz = 8*1024*1024;
    std::vector<char> buf(bufsz);
    uint64_t done=0;
    while ( done < bytes ) {
      size_t want = (size_t)std::min((uint64_t)bufsz,bytes-done);
      ssize_t got = pread(fd,&buf[0],want,(off_t)(offset+done));
      SYSCHECK(got==(ssize_t)want,"pread for serial crc");
      crc = crc32_buf(&buf[0],want,crc);
      done += want;
    }
    close(fd);
  }
  CHECK(MPI_Bcast(&crc,1,MPI_UINT32_T,0,WorldComm));
  return crc;
}

/**************************************************************
 * The three write / read paths.  Each returns seconds, measured as the
 * maximum over ranks of the interval spanning open..close, with a
 * barrier before the clock starts.
 **************************************************************/
static int    do_fsync   = 0;
static int    do_memsub  = 0;
static int    use_stdio  = 0;
static int    reuse_plan = 0;
static uint64_t file_offset = 0;
static MPI_Info io_hints = MPI_INFO_NULL;
static int    hints_reported = 0;

/* Grid's bracket, exactly: barrier, start, work, barrier, stop, and quote
   the boss rank's own stopwatch.  The trailing barrier makes that the
   slowest rank plus the barrier itself, which is what BinaryIO.h reports;
   an allreduced max of per-rank spans is a slightly different -- and
   slightly smaller -- number, and the difference would show up as the two
   tools disagreeing for no reason to do with I/O. */
static double Elapsed(double t0)
{
  CHECK(MPI_Barrier(WorldComm));
  double dt = (usecond()-t0)/1.0e6;
  CHECK(MPI_Bcast(&dt,1,MPI_DOUBLE,0,WorldComm));
  return dt;
}

/**************************************************************
 * File sink / source.
 *
 * Grid writes through std::ofstream (buffered) and this reproducer
 * defaults to pwrite (unbuffered).  --stdio selects the former, so that
 * if the two tools disagree on the POSIX paths the buffering can be
 * ruled in or out without editing anything.
 **************************************************************/
struct Sink   { int fd; std::ofstream *os; };
struct Source { int fd; std::ifstream *is; };

/* boss creates, everyone else only seeks and writes into a file that
   already exists.  No rank but 0 ever uses O_CREAT/O_TRUNC or ios::out:
   an unsynchronised truncation from a late opener silently erases an
   earlier rank's segment. */
static Sink SinkOpen(const char *file)
{
  Sink s; s.fd = -1; s.os = NULL;
  if ( !WorldRank ) {
    int fd = open(file,O_WRONLY|O_CREAT,0644);
    SYSCHECK(fd>=0,"create");
    close(fd);
  }
  CHECK(MPI_Barrier(WorldComm));
  if ( use_stdio ) {
    s.os = new std::ofstream(file,std::ios::binary|std::ios::out|std::ios::in);
    SYSCHECK(s.os->is_open(),"ofstream open for write");
  } else {
    s.fd = open(file,O_WRONLY);
    SYSCHECK(s.fd>=0,"open for write");
  }
  return s;
}
static void SinkWrite(Sink &s,const char *buf,uint64_t bytes,uint64_t off)
{
  if ( s.os ) {
    s.os->seekp((std::streamoff)off);
    s.os->write(buf,(std::streamsize)bytes);
    SYSCHECK(!s.os->fail(),"ofstream write");
    return;
  }
  uint64_t done=0;
  while ( done < bytes ) {
    ssize_t w = pwrite(s.fd,buf+done,bytes-done,(off_t)(off+done));
    SYSCHECK(w>0,"pwrite");
    done += (uint64_t)w;
  }
}
static void SinkClose(Sink &s,const char *file,uint64_t need)
{
  if ( s.os ) {
    s.os->flush();
    if ( do_fsync ) { /* no portable fd from ofstream; reopen to sync */
      int fd = open(file,O_WRONLY); if ( fd>=0 ) { fsync(fd); close(fd); }
    }
    s.os->close(); delete s.os; s.os = NULL;
  } else {
    if ( do_fsync ) SYSCHECK(fsync(s.fd)==0,"fsync");
    SYSCHECK(close(s.fd)==0,"close");
    s.fd = -1;
  }
  CHECK(MPI_Barrier(WorldComm));          /* every extent must be on its way */
  if ( !WorldRank ) {
    struct stat sb;
    SYSCHECK(stat(file,&sb)==0,"stat");
    if ( (uint64_t)sb.st_size != need )   /* only when a longer record preceded us */
      SYSCHECK(truncate(file,(off_t)need)==0,"truncate");
  }
  CHECK(MPI_Barrier(WorldComm));
}

static Source SourceOpen(const char *file)
{
  Source s; s.fd = -1; s.is = NULL;
  if ( use_stdio ) {
    s.is = new std::ifstream(file,std::ios::binary|std::ios::in);
    SYSCHECK(s.is->is_open(),"ifstream open for read");
  } else {
    s.fd = open(file,O_RDONLY);
    SYSCHECK(s.fd>=0,"open for read");
  }
  return s;
}
static void SourceRead(Source &s,char *buf,uint64_t bytes,uint64_t off)
{
  if ( s.is ) {
    s.is->seekg((std::streamoff)off);
    s.is->read(buf,(std::streamsize)bytes);
    SYSCHECK(!s.is->fail(),"ifstream read");
    return;
  }
  uint64_t done=0;
  while ( done < bytes ) {
    ssize_t r = pread(s.fd,buf+done,bytes-done,(off_t)(off+done));
    SYSCHECK(r>0,"pread");
    done += (uint64_t)r;
  }
}
static void SourceClose(Source &s)
{
  if ( s.is ) { s.is->close(); delete s.is; s.is = NULL; }
  else        { close(s.fd); s.fd = -1; }
}

/*-------------------- raw: one contiguous segment per rank ----------------*/
static double WriteRaw(const Layout &L,const double *data,const char *file)
{
  uint64_t bytes = (uint64_t)L.lsites*fobjSize;
  CHECK(MPI_Barrier(WorldComm));
  double t0 = usecond();
  Sink s = SinkOpen(file);
  SinkWrite(s,(const char *)data,bytes,file_offset+L.rawOffset);
  SinkClose(s,file,file_offset+(uint64_t)L.gsites*fobjSize);
  return Elapsed(t0);
}
static double ReadRaw(const Layout &L,double *data,const char *file)
{
  uint64_t bytes = (uint64_t)L.lsites*fobjSize;
  CHECK(MPI_Barrier(WorldComm));
  double t0 = usecond();
  Source s = SourceOpen(file);
  SourceRead(s,(char *)data,bytes,file_offset+L.rawOffset);
  SourceClose(s);
  return Elapsed(t0);
}

/*-------------------- mpiio: subarray view + collective -------------------*/
static void ReportHints(MPI_File fh)
{
  if ( WorldRank || hints_reported ) return;
  hints_reported = 1;
  MPI_Info info;
  if ( MPI_File_get_info(fh,&info) != MPI_SUCCESS ) return;
  int nkeys; MPI_Info_get_nkeys(info,&nkeys);
  printf("= ROMIO hints actually in force (%d):\n",nkeys);
  for(int i=0;i<nkeys;i++){
    char key[MPI_MAX_INFO_KEY], val[MPI_MAX_INFO_VAL]; int flag;
    MPI_Info_get_nthkey(info,i,key);
    MPI_Info_get(info,key,MPI_MAX_INFO_VAL-1,val,&flag);
    if ( flag ) printf("=   %-32s %s\n",key,val);
  }
  MPI_Info_free(&info);
  fflush(stdout);
}

static void MPIIOTypes(const Layout &L,MPI_Datatype siteType,
                       MPI_Datatype *fileArray,MPI_Datatype *memType)
{
  std::vector<int> g(Nd),l(Nd),s(Nd),z(Nd,0);
  for(int d=0;d<Nd;d++){ g[d]=(int)L.gLattice[d]; l[d]=(int)L.lLattice[d]; s[d]=(int)L.gStart[d]; }

  CHECK(MPI_Type_create_subarray(Nd,&g[0],&l[0],&s[0],MPI_ORDER_FORTRAN,siteType,fileArray));
  CHECK(MPI_Type_commit(fileArray));

  if ( do_memsub ) {   /* the degenerate subarray Grid builds, for exactness */
    CHECK(MPI_Type_create_subarray(Nd,&l[0],&l[0],&z[0],MPI_ORDER_FORTRAN,siteType,memType));
    CHECK(MPI_Type_commit(memType));
  } else {
    *memType = MPI_DATATYPE_NULL;
  }
}

static double MPIIO(const Layout &L,double *data,const char *file,
                    MPI_Datatype siteType,int write)
{
  MPI_Datatype fileArray,memType;
  MPIIOTypes(L,siteType,&fileArray,&memType);

  MPI_File fh; MPI_Status st;
  int amode = write ? (MPI_MODE_RDWR|MPI_MODE_CREATE) : MPI_MODE_RDONLY;

  CHECK(MPI_Barrier(WorldComm));
  double t0 = usecond();
  CHECK(MPI_File_open(WorldComm,(char *)file,amode,io_hints,&fh));
  CHECK(MPI_File_set_view(fh,(MPI_Offset)file_offset,siteType,fileArray,
                          (char *)"native",io_hints));
  if ( write ) {
    if ( do_memsub ) CHECK(MPI_File_write_all(fh,data,1,memType,&st));
    else             CHECK(MPI_File_write_all(fh,data,(int)L.lsites,siteType,&st));
    if ( do_fsync )  CHECK(MPI_File_sync(fh));
    /* Grid recovers the next record's displacement this way, and it is
       inside its timed region, so it is inside ours. */
    MPI_Offset os,disp;
    CHECK(MPI_File_get_position(fh,&os));
    CHECK(MPI_File_get_byte_offset(fh,os,&disp));
  } else {
    if ( do_memsub ) CHECK(MPI_File_read_all(fh,data,1,memType,&st));
    else             CHECK(MPI_File_read_all(fh,data,(int)L.lsites,siteType,&st));
  }
  ReportHints(fh);
  CHECK(MPI_File_close(&fh));
  double dt = Elapsed(t0);

  MPI_Type_free(&fileArray);
  if ( do_memsub ) MPI_Type_free(&memType);
  return dt;
}

/*-------------------- aggregate: Alltoallv + big pwrites ------------------*/
static double AggregateIO(const Layout &L,AggregationPlan &hoisted,uint64_t target,
                          double *data,const char *file,
                          MPI_Datatype siteType,int write,
                          double *t_comm,double *t_perm,double *t_plan)
{
  uint64_t need = file_offset + (uint64_t)L.gsites*fobjSize;

  CHECK(MPI_Barrier(WorldComm));
  double t0 = usecond();

  /* Grid rebuilds the plan on every IOobject call, inside its timed region,
     so by default so do we -- otherwise the two tools are not measuring the
     same thing.  --reuse-plan hoists it out and reports what it cost, which
     is the way to show the plan is not where the time goes. */
  AggregationPlan local;
  AggregationPlan &p = reuse_plan ? hoisted : local;
  if ( !reuse_plan ) {
    double tp0 = usecond();
    BuildAggregationPlan(L,target,local);
    *t_plan += usecond()-tp0;
  }
  std::vector<double> aggregated(L.lsites*NWORD);

  if ( write ) {
    AggregateExchange(p,data,&aggregated[0],siteType,1,t_comm,t_perm);
    Sink s = SinkOpen(file);
    for(size_t e=0;e<p.extentSites.size();e++)
      SinkWrite(s,(const char *)(&aggregated[0] + p.extentLocal[e]*NWORD),
                (uint64_t)p.extentSites[e]*fobjSize,
                file_offset + (uint64_t)p.extentGsite[e]*fobjSize);
    SinkClose(s,file,need);
  } else {
    Source s = SourceOpen(file);
    for(size_t e=0;e<p.extentSites.size();e++)
      SourceRead(s,(char *)(&aggregated[0] + p.extentLocal[e]*NWORD),
                 (uint64_t)p.extentSites[e]*fobjSize,
                 file_offset + (uint64_t)p.extentGsite[e]*fobjSize);
    SourceClose(s);
    AggregateExchange(p,data,&aggregated[0],siteType,0,t_comm,t_perm);
  }
  double dt = Elapsed(t0);
  if ( !reuse_plan ) MPI_Comm_free(&local.rowcomm);
  return dt;
}

/**************************************************************
 * Drop whatever of the file this client is caching.  Advisory only,
 * and Lustre is free to ignore it, but without it every read after a
 * write in the same job is a page cache measurement.
 **************************************************************/
static void DropCache(const char *file)
{
  int fd = open(file,O_RDONLY);
  if ( fd>=0 ) {
#ifdef POSIX_FADV_DONTNEED
    posix_fadvise(fd,0,0,POSIX_FADV_DONTNEED);
#endif
    close(fd);
  }
  CHECK(MPI_Barrier(WorldComm));
}

/**************************************************************
 * Command line
 **************************************************************/
static std::string CmdPayload(char **b,char **e,const std::string &opt)
{
  char **itr = std::find(b,e,opt);
  if ( itr != e && ++itr != e ) return std::string(*itr);
  return std::string("");
}
static bool CmdExists(char **b,char **e,const std::string &opt){ return std::find(b,e,opt) != e; }
static void CmdIntVector(const std::string &str,std::vector<int64_t> &vec)
{
  vec.resize(0);
  std::stringstream ss(str);
  long long i;
  while ( ss >> i ) { vec.push_back((int64_t)i); if ( ispunct(ss.peek()) ) ss.ignore(); }
}
static void ParseHints(const std::string &str)
{
  if ( str.empty() ) return;
  CHECK(MPI_Info_create(&io_hints));
  std::stringstream ss(str);
  std::string kv;
  while ( std::getline(ss,kv,',') ) {
    size_t eq = kv.find('=');
    if ( eq == std::string::npos ) continue;
    std::string k = kv.substr(0,eq), v = kv.substr(eq+1);
    CHECK(MPI_Info_set(io_hints,(char *)k.c_str(),(char *)v.c_str()));
    if ( !WorldRank ) printf("= hint %s = %s\n",k.c_str(),v.c_str());
  }
}

/**************************************************************
 * Reporting
 **************************************************************/
struct Samples { std::vector<double> s; };

/* MiB/s, dividing by 1024*1024 -- which is what BinaryIO.h computes for
   lastPerf.mbytesPerSecond and prints as "MB/s".  Quoting true decimal MB/s
   here would make this tool read 4.86% faster than Grid on identical work,
   and that discrepancy would be blamed on something real. */
static void Record(Samples &S,double bytes,double secs){ S.s.push_back(bytes/1024.0/1024.0/secs); }

static void ReportSamples(const char *what,const char *path,const Samples &S)
{
  if ( WorldRank || S.s.empty() ) return;
  double best=0, sum=0;
  for(size_t i=0;i<S.s.size();i++){ best = std::max(best,S.s[i]); sum += S.s[i]; }
  printf("PERF %-5s %-9s best %10.1f MiB/s mean %10.1f MiB/s first %10.1f MiB/s samples",
         what,path,best,sum/S.s.size(),S.s[0]);
  for(size_t i=0;i<S.s.size();i++) printf(" %.1f",S.s[i]);
  printf("\n"); fflush(stdout);
}

/**************************************************************
 * Main
 **************************************************************/
int main(int argc,char **argv)
{
  MPI_Init(&argc,&argv);
  WorldComm = MPI_COMM_WORLD;
  MPI_Comm_rank(WorldComm,&WorldRank);
  MPI_Comm_size(WorldComm,&WorldSize);
  MPI_Comm shm;
  MPI_Comm_split_type(WorldComm,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,&shm);
  MPI_Comm_rank(shm,&ShmRank); MPI_Comm_size(shm,&ShmSize);
  MPI_Comm_free(&shm);

  crc32_init();

  std::vector<int64_t> g,p;
  if ( CmdExists(argv,argv+argc,"--grid") ) CmdIntVector(CmdPayload(argv,argv+argc,"--grid"),g);
  if ( CmdExists(argv,argv+argc,"--mpi")  ) CmdIntVector(CmdPayload(argv,argv+argc,"--mpi"), p);
  if ( g.size()==0 || p.size()==0 || g.size()!=p.size() ) {
    if ( !WorldRank ) {
      printf("usage: io_mpi --grid n1.n2.n3.n4 --mpi p1.p2.p3.p4 [options]\n");
      printf("  --words N          doubles per site        (default 72, = 576 B)\n");
      printf("  --target BYTES     aggregate extent target (default 4194304)\n");
      printf("  --reps N           timed repetitions       (default 3)\n");
      printf("  --offset BYTES     record displacement in the file (default 0)\n");
      printf("  --hints k=v,k=v    MPI_Info passed to open and set_view\n");
      printf("  --fsync            include fsync / MPI_File_sync in the timed region\n");
      printf("  --drop-cache       posix_fadvise(DONTNEED) between write and read\n");
      printf("  --mem-subarray     use a degenerate memory subarray, as Grid does\n");
      printf("  --stdio            POSIX paths via std::ofstream/ifstream, as Grid does\n");
      printf("  --reuse-plan       build the aggregation plan once, outside the timed region\n");
      printf("  --serial-crc       rank 0 streams and crc32s each file (slow)\n");
      printf("  --no-validate      skip correctness, timing only\n");
    }
    MPI_Finalize(); return 0;
  }
  Nd = (int)g.size();

  int64_t reps = 3;
  uint64_t target = 4*1024*1024;
  int validate = !CmdExists(argv,argv+argc,"--no-validate");
  int serialcrc =  CmdExists(argv,argv+argc,"--serial-crc");
  int dropcache =  CmdExists(argv,argv+argc,"--drop-cache");
  do_fsync      =  CmdExists(argv,argv+argc,"--fsync");
  do_memsub     =  CmdExists(argv,argv+argc,"--mem-subarray");
  use_stdio     =  CmdExists(argv,argv+argc,"--stdio");
  reuse_plan    =  CmdExists(argv,argv+argc,"--reuse-plan");
  if ( CmdExists(argv,argv+argc,"--words")  ) NWORD  = atoll(CmdPayload(argv,argv+argc,"--words").c_str());
  if ( CmdExists(argv,argv+argc,"--target") ) target = strtoull(CmdPayload(argv,argv+argc,"--target").c_str(),NULL,0);
  if ( CmdExists(argv,argv+argc,"--reps")   ) reps   = atoll(CmdPayload(argv,argv+argc,"--reps").c_str());
  if ( CmdExists(argv,argv+argc,"--offset") ) file_offset = strtoull(CmdPayload(argv,argv+argc,"--offset").c_str(),NULL,0);
  fobjSize = NWORD*sizeof(double);

  Layout L;
  BuildLayout(L,g,p);

  { int64_t prod=1; for(int d=0;d<Nd;d++) prod *= p[d];
    if ( prod != WorldSize ) {
      if(!WorldRank) printf("--mpi product %lld != %d ranks\n",(long long)prod,WorldSize);
      MPI_Finalize(); return 1;
    } }

  ParseHints(CmdPayload(argv,argv+argc,"--hints"));

  MPI_Datatype siteType;
  CHECK(MPI_Type_contiguous((int)fobjSize,MPI_BYTE,&siteType));
  CHECK(MPI_Type_commit(&siteType));

  AggregationPlan plan;
  BuildAggregationPlan(L,target,plan);

  uint64_t payload = (uint64_t)L.gsites*fobjSize;

  if ( !WorldRank ) {
    printf("===========================================================\n");
    printf("= io_mpi : %d ranks, %d ranks/node, %d nodes\n",
           WorldSize,ShmSize,WorldSize/ShmSize);
    printf("= global "); for(int d=0;d<Nd;d++) printf("%lld ",(long long)g[d]);
    printf("  mpi ");    for(int d=0;d<Nd;d++) printf("%lld ",(long long)p[d]);
    printf("  local ");  for(int d=0;d<Nd;d++) printf("%lld ",(long long)L.lLattice[d]);
    printf("\n");
    printf("= site %lld B, %lld sites/rank = %.1f MB/rank, record %.2f GB\n",
           (long long)fobjSize,(long long)L.lsites,
           L.lsites*fobjSize/1048576.0,payload/1073741824.0);
    printf("= file view: %lld contiguous run(s)/rank of %lld B\n",
           (long long)(L.lsites/L.lLattice[0]),(long long)(L.lLattice[0]*fobjSize));
    printf("= offset %llu, fsync %s, mem-subarray %s, stdio %s, plan %s, reps %lld\n",
           (unsigned long long)file_offset,do_fsync?"on":"off",
           do_memsub?"on":"off",use_stdio?"on":"off",
           reuse_plan?"hoisted":"timed (as Grid)",(long long)reps);
    printf("===========================================================\n");
    fflush(stdout);
  }
  ReportPlan(plan,target);

  /* global index of each local site, in local lexicographic order.  The
     three paths all present data in this order, so one table serves all. */
  std::vector<int64_t> gindex(L.lsites);
  { std::vector<int64_t> lc(Nd), gc(Nd);
    for(int64_t i=0;i<L.lsites;i++){
      CoorFromIndex(&lc[0],i,&L.lLattice[0],Nd);
      for(int d=0;d<Nd;d++) gc[d] = L.gStart[d]+lc[d];
      gindex[i] = IndexFromCoor(&gc[0],&L.gLattice[0],Nd);
    } }

  std::vector<double> data(L.lsites*NWORD), back(L.lsites*NWORD);
  for(int64_t i=0;i<L.lsites;i++) FillSite(&data[i*NWORD],(uint64_t)gindex[i]);

  Checksum ref = SiteChecksum(&data[0],L.lsites,&gindex[0]);
  if ( !WorldRank ) printf("= reference site checksum %08x %08x\n",ref.a,ref.b);

  const char *fraw = "io_raw.bin", *fmpi = "io_mpiio.bin", *fagg = "io_agg.bin";
  int failures = 0;

  /* Start from nothing.  MPI-IO never shortens a file, so a leftover file
     from a larger configuration would survive underneath this record and
     the length check below would fire on the stale tail rather than on
     anything this run did.  (That the tail survives at all is the reason
     the aggregate path truncates; it is not what we are testing here.) */
  if ( !WorldRank ) { unlink(fraw); unlink(fmpi); unlink(fagg); }
  CHECK(MPI_Barrier(WorldComm));

  /*=================== correctness ===================*/
  if ( validate ) {
    if ( !WorldRank ) { printf("\n= CORRECTNESS\n"); fflush(stdout); }

    WriteRaw(L,&data[0],fraw);
    MPIIO(L,&data[0],fmpi,siteType,1);
    { double tc=0,tp=0,tb=0; AggregateIO(L,plan,target,&data[0],fagg,siteType,1,&tc,&tp,&tb); }

    /* the two lexicographic files must be identical, and both must read
       back correctly through BOTH lexicographic paths.  Cross reading is
       the point: a shared misunderstanding of the layout cannot survive
       mpiio-write/aggregate-read agreeing with aggregate-write/mpiio-read. */
    struct { const char *file; const char *wpath; } lex[2] = { {fmpi,"mpiio"}, {fagg,"aggregate"} };
    for(int w=0;w<2;w++){
      for(int r=0;r<2;r++){
        memset(&back[0],0,back.size()*sizeof(double));
        double tc=0,tp=0,tb=0;
        if ( r==0 ) MPIIO(L,&back[0],lex[w].file,siteType,0);
        else        AggregateIO(L,plan,target,&back[0],lex[w].file,siteType,0,&tc,&tp,&tb);

        Checksum c  = SiteChecksum(&back[0],L.lsites,&gindex[0]);
        int64_t bad = ValidateLocal(&back[0],L.lsites,&gindex[0]);
        int64_t worst;
        CHECK(MPI_Allreduce(&bad,&worst,1,MPI_INT64_T,MPI_MAX,WorldComm));

        int ok = (c.a==ref.a) && (c.b==ref.b) && (worst<0);
        if ( !ok ) failures++;
        if ( !WorldRank ) {
          printf("=   write %-9s read %-9s : checksum %08x %08x %s, placement ",
                 lex[w].wpath, r==0?"mpiio":"aggregate", c.a,c.b,
                 (c.a==ref.a&&c.b==ref.b)?"OK":"MISMATCH");
          if ( worst<0 ) printf("OK\n");
          else           printf("WRONG, first bad global site %lld\n",(long long)worst);
        }
      }
    }

    /* raw is a different layout, so it validates only against itself */
    { memset(&back[0],0,back.size()*sizeof(double));
      ReadRaw(L,&back[0],fraw);
      Checksum c  = SiteChecksum(&back[0],L.lsites,&gindex[0]);
      int64_t bad = ValidateLocal(&back[0],L.lsites,&gindex[0]);
      int64_t worst; CHECK(MPI_Allreduce(&bad,&worst,1,MPI_INT64_T,MPI_MAX,WorldComm));
      int ok = (c.a==ref.a) && (c.b==ref.b) && (worst<0);
      if ( !ok ) failures++;
      if ( !WorldRank ) printf("=   write raw       read raw       : %s (control, non-lexicographic)\n",
                               ok?"OK":"FAILED"); }

    /* file lengths must be exactly offset + payload */
    if ( !WorldRank ) {
      const char *f[3] = { fraw,fmpi,fagg };
      for(int i=0;i<3;i++){
        struct stat sb;
        if ( stat(f[i],&sb)==0 && (uint64_t)sb.st_size != file_offset+payload ) {
          printf("=   LENGTH WRONG %s : %llu, expected %llu\n",
                 f[i],(unsigned long long)sb.st_size,
                 (unsigned long long)(file_offset+payload));
          failures++;
        }
      }
    }

    if ( serialcrc ) {
      uint32_t cm = SerialFileCrc(fmpi,file_offset,payload);
      uint32_t ca = SerialFileCrc(fagg,file_offset,payload);
      if ( !WorldRank ) {
        printf("=   serial whole-file crc32: mpiio %08x  aggregate %08x  %s\n",
               cm,ca,(cm==ca)?"IDENTICAL":"DIFFER");
        fflush(stdout);
      }
      if ( cm != ca ) failures++;
    }

    int f; CHECK(MPI_Allreduce(&failures,&f,1,MPI_INT,MPI_MAX,WorldComm)); failures = f;
    if ( !WorldRank ) printf("= correctness: %s\n",failures?"FAILED":"all checks passed");
  }

  /*=================== performance ===================*/
  if ( reps > 0 ) {
    if ( !WorldRank ) { printf("\n= PERFORMANCE (%lld reps)\n",(long long)reps); fflush(stdout); }

    Samples wraw,wmpi,wagg,rraw,rmpi,ragg;
    double bytes = (double)payload;
    double tc=0,tp=0,tb=0;

    for(int64_t n=0;n<reps;n++){
      Record(wraw,bytes,WriteRaw(L,&data[0],fraw));
      Record(wmpi,bytes,MPIIO(L,&data[0],fmpi,siteType,1));
      Record(wagg,bytes,AggregateIO(L,plan,target,&data[0],fagg,siteType,1,&tc,&tp,&tb));
    }
    if ( dropcache ) { DropCache(fraw); DropCache(fmpi); DropCache(fagg); }
    for(int64_t n=0;n<reps;n++){
      Record(rraw,bytes,ReadRaw(L,&back[0],fraw));
      Record(rmpi,bytes,MPIIO(L,&back[0],fmpi,siteType,0));
      Record(ragg,bytes,AggregateIO(L,plan,target,&back[0],fagg,siteType,0,&tc,&tp,&tb));
    }

    ReportSamples("write","raw",      wraw);
    ReportSamples("write","mpiio",    wmpi);
    ReportSamples("write","aggregate",wagg);
    ReportSamples("read", "raw",      rraw);
    ReportSamples("read", "mpiio",    rmpi);
    ReportSamples("read", "aggregate",ragg);

    if ( !WorldRank )
      printf("PERF alltoallv %.3f s, permute %.3f s, plan %.3f s over all %lld aggregate calls "
             "(rank 0; the redistribution cost, inside the timed region)\n",
             tc/1.0e6,tp/1.0e6,tb/1.0e6,(long long)(2*reps));
  }

  if ( !WorldRank ) { printf("= done\n"); fflush(stdout); }

  MPI_Type_free(&siteType);
  MPI_Comm_free(&plan.rowcomm);
  if ( io_hints != MPI_INFO_NULL ) MPI_Info_free(&io_hints);
  MPI_Finalize();
  return failures ? 1 : 0;
}
