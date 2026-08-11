    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/IO/Test_aggregate_io.cc

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

// Correctness and performance test for BINARYIO_AGGREGATE.
//
// Correctness, per aggregateTargetBytes:
//   1. write via the MPI-IO lexicographic path      -> ref.bin
//   2. write via the aggregate path                 -> agg.bin
//   3. the two files must be byte identical         <- proves the layout matches
//   4. read agg.bin back through the aggregate path <- proves the mirror inverts
//   5. write and read the non-lexicographic path    -> raw.bin.  Its layout is
//      different by construction (each rank owns one contiguous segment in
//      rank order) so it cannot be compared byte for byte, but the NERSC and
//      SciDAC checksums are computed from the global site index and are
//      therefore layout independent: they must match the other two paths.
//   6. a record written over a longer pre-existing file must leave the file at
//      exactly offset+payload, with no trailing fragment of the old contents
//
// Performance: three paths, both directions, timed with the client page cache
// dropped before every read so that a read back reports filesystem bandwidth
// rather than memory bandwidth.  The non-lexicographic path is the zero
// overhead reference: no transposition, no layout independence, one disjoint
// contiguous segment per rank, which is the arrangement that reaches full
// filesystem bandwidth on a leadership machine.  It is the upper bound the
// other two are trying to approach.
//
// Options:
//   --aggregate-target <bytes>   sweep this one target only (default: sweep
//                                1, 1024, 64K, 4M)
//   --io-reps <n>                repetitions in the performance section
//                                (default 3; 0 disables it)
//   --io-no-correctness          skip the correctness section, which reads the
//                                whole file on one rank and is not affordable
//                                at very large volume
//   --io-read-only               time reads only, of files left in place by an
//                                earlier job.  Reading back what this job just
//                                wrote measures the client page cache; a fresh
//                                allocation pointed at the same directory is
//                                the only way to get a cold read without root.
//
// The exchange is only meaningfully exercised when the fast dimensions are
// split across ranks; --mpi 1.1.X.Y leaves the rows of size one and the test
// then passes vacuously.  Non-uniform AllToAllV counts additionally need an
// odd process factor in a fast dimension and a small local volume.

#include <Grid/Grid.h>
#include <fcntl.h>

using namespace Grid;

/////////////////////////////////////////////////////////////////////////////
// Compare in chunks.  Slurping both files into memory is fine for a few MB
// and fatal for the multi-GB records this test is meant to reach.
/////////////////////////////////////////////////////////////////////////////
static bool FilesIdentical(std::string a,std::string b)
{
  std::ifstream fa(a,std::ios::binary), fb(b,std::ios::binary);
  if ( !fa.good() || !fb.good() ) {
    std::cout<<GridLogMessage<<"  could not open "<<a<<" and/or "<<b<<std::endl;
    return false;
  }
  fa.seekg(0,std::ios::end); fb.seekg(0,std::ios::end);
  uint64_t sa = (uint64_t)fa.tellg(), sb = (uint64_t)fb.tellg();
  if ( sa != sb ) {
    std::cout<<GridLogMessage<<"  size mismatch "<<sa<<" vs "<<sb<<std::endl;
    return false;
  }
  fa.seekg(0,std::ios::beg); fb.seekg(0,std::ios::beg);

  const uint64_t chunk = 8*1024*1024;
  std::vector<char> va(chunk), vb(chunk);
  uint64_t done=0;
  while ( done < sa ) {
    uint64_t n = std::min(chunk,sa-done);
    fa.read(&va[0],n);
    fb.read(&vb[0],n);
    for(uint64_t i=0;i<n;i++){
      if ( va[i]!=vb[i] ) {
        std::cout<<GridLogMessage<<"  first differing byte at "<<done+i<<" of "<<sa<<std::endl;
        return false;
      }
    }
    done += n;
  }
  return true;
}

/////////////////////////////////////////////////////////////////////////////
// Reading back a file we have just written measures the client page cache,
// not the filesystem: the earlier runs of this test reported 8 GB/s on reads
// and ~1 GB/s on writes for the same data.  POSIX_FADV_DONTNEED asks the
// kernel to drop the cached pages for the file.  It is advisory and every
// rank must do it, since each client caches independently, so treat this as
// best effort rather than a guarantee of a cold read.
/////////////////////////////////////////////////////////////////////////////
static void DropCache(GridBase *grid,std::string file)
{
  grid->Barrier();
  int fd = ::open(file.c_str(),O_RDONLY);
  if ( fd >= 0 ) {
#ifdef POSIX_FADV_DONTNEED
    ::posix_fadvise(fd,0,0,POSIX_FADV_DONTNEED);
#endif
    ::close(fd);
  }
  grid->Barrier();
}

static uint64_t OptionU64(int argc,char **argv,const char *opt,uint64_t def)
{
  if ( GridCmdOptionExists(argv,argv+argc,opt) ) {
    std::string arg = GridCmdOptionPayload(argv,argv+argc,opt);
    return (uint64_t)std::stoull(arg);
  }
  return def;
}

int main(int argc,char **argv)
{
  Grid_init(&argc,&argv);

  Coordinate latt = GridDefaultLatt();
  Coordinate simd = GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate mpi  = GridDefaultMpi();
  GridCartesian grid(latt,simd,mpi);

  typedef vLorentzColourMatrixD vobj;
  typedef LorentzColourMatrixD  sobj;

  GridParallelRNG pRNG(&grid);
  pRNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));
  LatticeGaugeFieldD Umu(&grid);
  random(pRNG,Umu);

  BinarySimpleMunger<sobj,sobj> munge;
  const std::string format("IEEE64BIG");

  const int lex = BinaryIO::BINARYIO_LEXICOGRAPHIC;
  const int agg = BinaryIO::BINARYIO_LEXICOGRAPHIC|BinaryIO::BINARYIO_AGGREGATE;
  const int raw = 0;   // no BINARYIO_LEXICOGRAPHIC: contiguous segment per rank

  uint64_t payload = (uint64_t)grid._gsites*sizeof(sobj);

  std::vector<uint64_t> targets = {1, 1024, 64*1024, 4*1024*1024};
  if ( GridCmdOptionExists(argv,argv+argc,"--aggregate-target") ) {
    targets.clear();
    targets.push_back(OptionU64(argc,argv,"--aggregate-target",4*1024*1024));
  }
  uint64_t reps        = OptionU64(argc,argv,"--io-reps",3);
  bool     correctness = !GridCmdOptionExists(argv,argv+argc,"--io-no-correctness");
  // Read only: time reads of files left by an earlier job.  The only way to
  // get a cold client cache without root is to read on an allocation that did
  // not write the data, so run one job to write and a second, pointed at the
  // same directory, with this flag.
  bool     readonly    = GridCmdOptionExists(argv,argv+argc,"--io-read-only");
  if ( readonly ) correctness = false;

  std::cout<<GridLogMessage<<"Record payload "<<payload<<" bytes = "
           <<payload/1024./1024.<<" MB, "
           <<payload/(RealD)grid.ProcessorCount()/1024./1024.<<" MB/rank"<<std::endl;

  int failures=0;

  //////////////////////////////////////////////////////////////////////////
  // Correctness
  //////////////////////////////////////////////////////////////////////////
  if ( correctness ) for(auto target : targets){

    std::cout<<GridLogMessage<<"=== correctness, aggregateTargetBytes = "<<target<<" ==="<<std::endl;

    uint32_t n1,a1,b1, n2,a2,b2, n3,a3,b3;
    uint64_t off;

    // Start from a clean slate.  The aggregate path sets the file length to
    // exactly offset+payload; the MPI-IO path (MPI_MODE_CREATE) leaves any
    // pre-existing tail in place.  Comparing stale files would therefore
    // report a size mismatch that says nothing about the payload.
    if ( grid.IsBoss() ) { ::unlink("ref.bin"); ::unlink("agg.bin"); ::unlink("raw.bin"); }
    grid.Barrier();

    off=0;
    BinaryIO::writeLatticeObject<vobj,sobj>(Umu,"ref.bin",munge,off,format,n1,a1,b1,lex);

    BinaryIO::aggregateTargetBytes = target;
    off=0;
    BinaryIO::writeLatticeObject<vobj,sobj>(Umu,"agg.bin",munge,off,format,n2,a2,b2,agg);

    grid.Barrier();

    if ( grid.IsBoss() ) {
      if ( !FilesIdentical("ref.bin","agg.bin") ) {
        std::cout<<GridLogError<<"  FAIL: aggregate file differs from lexicographic file"<<std::endl;
        failures++;
      } else {
        std::cout<<GridLogMessage<<"  files byte identical"<<std::endl;
      }
    }
    if ( (n1!=n2)||(a1!=a2)||(b1!=b2) ) {
      std::cout<<GridLogError<<"  FAIL: checksum mismatch between paths"<<std::endl;
      failures++;
    }
    // writeLatticeObject takes offset by value, so the out-parameter that
    // IOobject sets never reaches us here and cannot be checked directly.
    // The observable equivalent is the file length: both paths must leave the
    // record ending at exactly offset+payload.
    if ( grid.IsBoss() ) {
      for(auto f : {std::string("ref.bin"),std::string("agg.bin")}){
        std::ifstream fs(f,std::ios::binary|std::ios::ate);
        uint64_t sz = (uint64_t)fs.tellg();
        if ( sz != payload ) {
          std::cout<<GridLogError<<"  FAIL: "<<f<<" is "<<sz<<" bytes, expected "<<payload<<std::endl;
          failures++;
        }
      }
    }

    LatticeGaugeFieldD Uchk(&grid);
    DropCache(&grid,"agg.bin");
    off=0;
    BinaryIO::readLatticeObject<vobj,sobj>(Uchk,"agg.bin",munge,off,format,n3,a3,b3,agg);
    if ( (n3!=n1)||(a3!=a1)||(b3!=b1) ) {
      std::cout<<GridLogError<<"  FAIL: read back checksum mismatch"<<std::endl;
      failures++;
    }
    Uchk = Uchk - Umu;
    RealD residual = norm2(Uchk);
    std::cout<<GridLogMessage<<"  read back residual "<<residual<<std::endl;
    if ( residual != 0.0 ) {
      std::cout<<GridLogError<<"  FAIL: read back does not reproduce the field"<<std::endl;
      failures++;
    }

    ////////////////////////////////////////////////////////////////////////
    // Non-lexicographic.  Different file layout by construction, so compare
    // by checksum and by round trip rather than by bytes.
    ////////////////////////////////////////////////////////////////////////
    uint32_t n4,a4,b4, n5,a5,b5;
    off=0;
    BinaryIO::writeLatticeObject<vobj,sobj>(Umu,"raw.bin",munge,off,format,n4,a4,b4,raw);
    grid.Barrier();
    if ( (n4!=n1)||(a4!=a1)||(b4!=b1) ) {
      std::cout<<GridLogError<<"  FAIL: raw path checksum differs; the NERSC and"
               <<" SciDAC checksums are layout independent and must agree"<<std::endl;
      failures++;
    }
    DropCache(&grid,"raw.bin");
    off=0;
    BinaryIO::readLatticeObject<vobj,sobj>(Uchk,"raw.bin",munge,off,format,n5,a5,b5,raw);
    if ( (n5!=n1)||(a5!=a1)||(b5!=b1) ) {
      std::cout<<GridLogError<<"  FAIL: raw read back checksum mismatch"<<std::endl;
      failures++;
    }
    Uchk = Uchk - Umu;
    RealD rawresidual = norm2(Uchk);
    std::cout<<GridLogMessage<<"  raw read back residual "<<rawresidual<<std::endl;
    if ( rawresidual != 0.0 ) {
      std::cout<<GridLogError<<"  FAIL: raw read back does not reproduce the field"<<std::endl;
      failures++;
    }
  }

  //////////////////////////////////////////////////////////////////////////
  // Truncation.  offset!=0 is the case that matters: it is what ILDG and
  // NERSC use, and it is the branch that used to only ever grow the file.
  // The oversized starting file is made by extending a sparse one rather
  // than writing padding from a single rank, which does not scale.
  //////////////////////////////////////////////////////////////////////////
  if ( correctness ) {
    BinaryIO::aggregateTargetBytes = 4*1024*1024;
    for(uint64_t testOffset : {(uint64_t)0, (uint64_t)1024}){

      uint64_t expect = testOffset + payload;

      if ( grid.IsBoss() ) {
        { std::ofstream create("trunc.bin",std::ios::binary|std::ios::out); create.close(); }
        int ierr = ::truncate("trunc.bin",(off_t)(expect+65536));
        GRID_ASSERT(ierr==0);
      }
      grid.Barrier();

      uint32_t n,a,b;
      uint64_t off = testOffset;
      BinaryIO::writeLatticeObject<vobj,sobj>(Umu,"trunc.bin",munge,off,format,n,a,b,agg);
      grid.Barrier();

      if ( grid.IsBoss() ) {
        std::ifstream f("trunc.bin",std::ios::binary|std::ios::ate);
        uint64_t sz = (uint64_t)f.tellg();
        f.close();
        if ( sz != expect ) {
          std::cout<<GridLogError<<"  FAIL: offset "<<testOffset<<" left file "<<sz
                   <<" bytes, expected "<<expect<<std::endl;
          failures++;
        } else {
          std::cout<<GridLogMessage<<"  truncation ok at offset "<<testOffset
                   <<": file is exactly "<<sz<<" bytes"<<std::endl;
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////
  // Performance.  Four numbers per repetition: write and read, old path and
  // new.  Reads are preceded by a cache drop; writes are not, so a write
  // number is "time to hand the data to the client cache and close", the
  // same convention for both paths.
  //////////////////////////////////////////////////////////////////////////
  if ( reps ) {

    LatticeGaugeFieldD Uio(&grid);
    uint32_t n,a,b;

    // The NERSC and SciDAC checksums are computed from the global site index,
    // so all three layouts must produce the same values.  This costs nothing
    // and is the only correctness check available at a volume where the byte
    // for byte comparison (single rank, whole file) is unaffordable.
    uint32_t cn[6],ca[6],cb[6];
    auto agreeing = [&](const char *what,int lo,int hi){
      for(int i=lo+1;i<=hi;i++){
        if ( (cn[i]!=cn[lo])||(ca[i]!=ca[lo])||(cb[i]!=cb[lo]) ) {
          std::cout<<GridLogError<<"  FAIL: "<<what<<" checksums disagree between paths"<<std::endl;
          return false;
        }
      }
      return true;
    };

    for(auto target : targets){

      BinaryIO::aggregateTargetBytes = target;
      std::cout<<GridLogMessage<<"=== performance, aggregateTargetBytes = "<<target<<" ==="<<std::endl;

      std::vector<RealD> wref,wagg,wraw,rref,ragg,rraw;

      for(uint64_t rep=0;rep<reps;rep++){

        uint64_t off;

        // Unlink only before the first repetition.  Lustre metadata cost is
        // per file, not per byte, so rep 0 reports "create the file and write
        // it" and the later reps report the steady state of overwriting an
        // existing file -- which is what a multi record file does for every
        // record after the first, and what production actually looks like.
        if ( (rep==0) && !readonly && grid.IsBoss() ) {
          ::unlink("ref.bin"); ::unlink("agg.bin"); ::unlink("raw.bin");
        }
        grid.Barrier();

        if ( !readonly ) {

        off=0;
        BinaryIO::writeLatticeObject<vobj,sobj>(Umu,"ref.bin",munge,off,format,n,a,b,lex);
        wref.push_back(BinaryIO::lastPerf.mbytesPerSecond);
        cn[0]=n; ca[0]=a; cb[0]=b;

        off=0;
        BinaryIO::writeLatticeObject<vobj,sobj>(Umu,"agg.bin",munge,off,format,n,a,b,agg);
        wagg.push_back(BinaryIO::lastPerf.mbytesPerSecond);
        cn[1]=n; ca[1]=a; cb[1]=b;

        off=0;
        BinaryIO::writeLatticeObject<vobj,sobj>(Umu,"raw.bin",munge,off,format,n,a,b,raw);
        wraw.push_back(BinaryIO::lastPerf.mbytesPerSecond);
        cn[2]=n; ca[2]=a; cb[2]=b;
        if ( !agreeing("write",0,2) ) failures++;

        }  // !readonly

        DropCache(&grid,"ref.bin");
        off=0;
        BinaryIO::readLatticeObject<vobj,sobj>(Uio,"ref.bin",munge,off,format,n,a,b,lex);
        rref.push_back(BinaryIO::lastPerf.mbytesPerSecond);
        cn[3]=n; ca[3]=a; cb[3]=b;

        DropCache(&grid,"agg.bin");
        off=0;
        BinaryIO::readLatticeObject<vobj,sobj>(Uio,"agg.bin",munge,off,format,n,a,b,agg);
        ragg.push_back(BinaryIO::lastPerf.mbytesPerSecond);
        cn[4]=n; ca[4]=a; cb[4]=b;

        DropCache(&grid,"raw.bin");
        off=0;
        BinaryIO::readLatticeObject<vobj,sobj>(Uio,"raw.bin",munge,off,format,n,a,b,raw);
        rraw.push_back(BinaryIO::lastPerf.mbytesPerSecond);
        cn[5]=n; ca[5]=a; cb[5]=b;
        if ( !agreeing("read back",readonly?3:0,5) ) failures++;
      }

      if ( grid.IsBoss() ) {
        auto report = [&](const char *name,std::vector<RealD> &v){
          if ( v.empty() ) return;
          RealD best=0, sum=0;
          for(auto x : v){ if(x>best) best=x; sum+=x; }
          // First sample includes file creation, later ones do not; quote both
          // rather than a mean that mixes the two.
          std::cout<<GridLogMessage<<"  PERF target="<<target<<" "<<name
                   <<" best "<<best<<" MB/s, mean "<<sum/v.size()
                   <<" MB/s, first(cold create) "<<v[0]<<" MB/s, samples";
          for(auto x : v) std::cout<<" "<<x;
          std::cout<<std::endl;
        };
        report("write raw       ",wraw);   // zero overhead reference
        report("write MPI-IO    ",wref);
        report("write aggregate ",wagg);
        report("read  raw       ",rraw);
        report("read  MPI-IO    ",rref);
        report("read  aggregate ",ragg);

        // Fraction of the zero overhead reference that each layout preserving
        // path achieves.  This is the number the whole exercise is about.
        auto best = [](std::vector<RealD> &v){ RealD m=0; for(auto x:v) if(x>m) m=x; return m; };
        if ( !wraw.empty() && best(wraw) > 0 ) {
          std::cout<<GridLogMessage<<"  PERF target="<<target
                   <<" write fraction of raw: MPI-IO "<<best(wref)/best(wraw)
                   <<"  aggregate "<<best(wagg)/best(wraw)<<std::endl;
        }
        if ( !rraw.empty() && best(rraw) > 0 ) {
          std::cout<<GridLogMessage<<"  PERF target="<<target
                   <<" read  fraction of raw: MPI-IO "<<best(rref)/best(rraw)
                   <<"  aggregate "<<best(ragg)/best(rraw)<<std::endl;
        }
      }
    }
  }

  if ( grid.IsBoss() ) {
    if ( failures ) std::cout<<GridLogError  <<failures<<" FAILURE(S)"<<std::endl;
    else            std::cout<<GridLogMessage<<"ALL AGGREGATE IO TESTS PASSED"<<std::endl;
  }

  Grid_finalize();
  return failures!=0;
}
