/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/algorithms/multigrid/DenseCoarseMatrix.h

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

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

#include <Grid/algorithms/blas/BatchedBlas.h>
#include <Grid/algorithms/blas/BatchedInverse.h>
#include <Grid/algorithms/multigrid/RecursiveSchurInverse.h>
#include <Grid/algorithms/multigrid/BlockCyclicSchurInverse.h>
#include <Grid/algorithms/multigrid/BlockCyclicRedistribute.h>

#include <unordered_map>

NAMESPACE_BEGIN(Grid);

//////////////////////////////////////////////////////////////////////////////////////
// DenseCoarseMatrix: a coarsened operator treated as a DENSE matrix -- explicit,
// row-distributed A^{-1} of a GeneralCoarsenedMatrix.  Library-grade successor of
// the example-local DistributedDenseInverse (Example_pvdagm_mrhs_3level_dense.cc,
// FROZEN as the regression baseline).
//
// What is new versus the example class:
//  - Stencil -> dense DIRECT IMPORT.  The coarse operator IS the dense matrix
//    unrolled: Dense[(s,a),(s+shift_p,b)] += A[p][s]_{a,b}.  Rows of my sites are
//    assembled from purely LOCAL _A[p] data: no operator applies, no comms -- the
//    O(N) probe assembly (93 s at N=69120) is retired.  ACCUMULATE (+=) because on
//    short axes distinct shifts wrap to the same neighbour.  An IMPORT CERTIFICATE
//    compares the dense apply against Op.M on a NON-CONSTANT vector (a constant one
//    cannot see a shift-sign error); DENSE_IMPORT_SIGN=-1 flips the convention
//    without recompiling.
//  - Split-K apply through GridBLAS.gemmBatched with EXPLICIT leading dimensions
//    (arXiv:2409.03904 fig 11): the tiny-output/huge-K GEMM Y = slab^T X becomes
//    DENSE_SPLITK chunk-GEMMs by pointer offset into the resident slab (lda = N),
//    partials reduced in one accelerator_for.  Platform-agnostic: deviceVector +
//    GridBLAS run the SAME code on HIP/CUDA/SYCL and CPU(Eigen).
//  - deviceVector everywhere in the apply path; the ONE surviving naked-HIP block
//    is the boss inversion buffer (quarantined below, documented).
//
// Setup: SLAB_FILE=<stem> loads per-rank <stem>.<rank> (header-guarded N/nrows/
// nbasis -- the interchange format shared with the frozen example; the STEM must
// encode cfg/mass/blocking/nbasis, only the header is guarded).  Absent: direct
// import -> import certificate -> chunked zero-fill+GlobalSum gather streamed to
// the boss GCD -> cgetrf_64 (ILP64) -> rows of A^{-1} via blocked identity
// cgetrs_64 + broadcast, each rank keeping the rows of its own sites -> save.
// VERIFY ||A Ainv x - x||/||x|| runs in BOTH paths (and now certifies the DEVICE
// slab + split-K path, since the single-RHS apply routes through the same core).
//
// Env: SLAB_FILE  DENSE_SPLITK (default 32, snapped to a divisor of N)
//      DENSE_DEVICE_SUM  DENSE_IMPORT_SIGN  DENSE_APPLY_PROFILE  DENSE_CC_CHECK
//      DENSE_SCHUR (0/absent: single-GCD gather-invert; 1: distributed
//      recursive Schur; 2: AUDIT -- run BOTH on the same imported A, report
//      the slab difference, keep the Schur result)  DENSE_PANEL_BYTES
//
// The DENSE_SCHUR=1 path is the RecursiveSchurInverse distributed
// factorisation: it lifts the fp32 N ~ 90k boss-HBM ceiling (the CC-grid
// 256-rank SIMD cap remains -- separate issue).  Internal only: slab layout,
// apply path, SLAB_FILE format and VERIFY are identical in every mode.
//
// Tensor-depth agnostic: site scalar objects treated as contiguous ComplexD
// (iScalar wrappers add no data), so any MG level's coarse operator imports.
//////////////////////////////////////////////////////////////////////////////////////
//
// Depends on a coarse operator only to extract its matrix elements; thereafter
// it is given a coarse vector and applies the inverse. Import() is a template
// member so any of the coarse classes will do, and the type of a
// DenseCoarseMatrix does not record which one built it.
//
template<class CComplex,int nbasis>
class DenseCoarseMatrix : public LinearFunction<Lattice<iVector<CComplex,nbasis> > > {
public:
  typedef iVector<CComplex,nbasis >           siteVector;
  typedef Lattice<siteVector>                 CoarseVector;
  typedef Lattice<iMatrix<CComplex,nbasis > > CoarseMatrix;
  typedef CoarseVector                        Field;
  using LinearFunction<Field>::operator();
  typedef typename Field::vector_object vobj;
  typedef typename vobj::scalar_object  sobj;
  typedef typename CoarseMatrix::vector_object Mvobj;
  typedef typename Mvobj::scalar_object        Msobj;

  GridBase *grid;
  int      nd;
  int64_t  N;                       // dense rank = gSites * nbasis
  int      lsites;                  // my local sites
  int64_t  nrows;                   // my rows = lsites * nbasis
  std::vector<Coordinate> myLcoor;  // local coordinate of my site ss
  std::vector<int64_t>    myGsite;  // global lex site index of my site ss
  std::vector<ComplexF>   slab;     // nrows x N row-major: A during setup, rows of A^{-1} after

  static const int64_t CHUNKROWS = 1024;  // getrs harvest block (trsm efficiency + fewer broadcasts)
  static const int MRHS_MAX = 32;

  // Apply machinery: resident slab + persistent buffers + AOT split-K pointers.
  GridBLAS BLAS;
  deviceVector<ComplexF>  dSlab;
  deviceVector<ComplexF>  dX;       // N x MRHS_MAX
  deviceVector<ComplexF>  dY;       // nrows x MRHS_MAX
  deviceVector<ComplexF>  dPartial; // NK x (nrows x MRHS_MAX)
  deviceVector<ComplexF*> aptrs;    // slab K-chunk pointers   (lda = N)
  deviceVector<ComplexF*> xptrs;    // X    K-chunk pointers   (ldb = N)
  deviceVector<ComplexF*> cptrs;    // partial buffers         (ldc = nrows)
  std::vector<ComplexF>   hX;
  std::vector<ComplexF>   hY;
  int NK;                           // split-K chunk count (divides N)
  int devSum;
  double schurAuditRel;             // DENSE_SCHUR=2: rel slab diff single-vs-schur (-1 = not run)

  DenseCoarseMatrix(GridBase *g)
    : grid(g)
  {
    GRID_ASSERT( sizeof(sobj)  == nbasis*sizeof(ComplexD) );
    GRID_ASSERT( sizeof(Msobj) == nbasis*nbasis*sizeof(ComplexD) );
    nd     = grid->_ndimension;
    N      = grid->gSites() * nbasis;
    lsites = grid->lSites();
    nrows  = (int64_t)lsites * nbasis;
    schurAuditRel = -1.0;

    std::cout << GridLogMessage << "DenseCoarseMatrix: N = " << N
              << " (" << grid->gSites() << " sites x " << nbasis << ")"
              << "  rows/rank = " << nrows
              << "  slab = " << (double)nrows*N*sizeof(ComplexF)/1024./1024. << " MB/rank"
              << std::endl;

    ////////////////////////////////////////////////////////////////////
    // Enumerate my sites: local coords and global lexicographic indices
    ////////////////////////////////////////////////////////////////////
    Coordinate ldims = grid->LocalDimensions();
    Coordinate gdims = grid->GlobalDimensions();
    myLcoor.resize(lsites);
    myGsite.resize(lsites);
    for(int ss=0; ss<lsites; ss++){
      Coordinate lcoor(nd);
      Lexicographic::CoorFromIndex(lcoor, ss, ldims);
      Coordinate gcoor(nd);
      for(int d=0; d<nd; d++) gcoor[d] = grid->_lstart[d] + lcoor[d];
      int64_t gsite;
      Lexicographic::IndexFromCoor(gcoor, gsite, gdims);
      myLcoor[ss] = lcoor;
      myGsite[ss] = gsite;
    }

    slab.resize((uint64_t)nrows * N);
  }

  ////////////////////////////////////////////////////////////////////
  // The only place a coarse operator is needed: pull its elements, invert,
  // and make the slab resident. Any class exposing Geometry() and
  // ExtractMatrix(p,A) will do -- single RHS or either multiRHS.
  ////////////////////////////////////////////////////////////////////
  template<class CoarseOp>
  void Import(CoarseOp &Op)
  {
    double t0 = usecond();
    ////////////////////////////////////////////////////////////////////
    // 0. Slab cache: SLAB_FILE=<stem> -> per-rank raw file <stem>.<rank>.
    //    SAME format as the frozen example (interchange compatible).
    ////////////////////////////////////////////////////////////////////
    bool loaded = false;
    char *sfile = getenv("SLAB_FILE");
    std::string slabfile;
    if (sfile) {
      slabfile = std::string(sfile) + "." + std::to_string(grid->ThisRank());
      FILE *f = fopen(slabfile.c_str(),"rb");
      if (f) {
        int64_t hdr[4] = {0,0,0,0};
        GRID_ASSERT( fread(hdr,sizeof(int64_t),4,f) == 4 );
        GRID_ASSERT( hdr[0] == (int64_t)0x44454E5345 );  // magic "DENSE"
        GRID_ASSERT( hdr[1] == N && hdr[2] == (int64_t)nrows && hdr[3] == (int64_t)nbasis );
        uint64_t nelem = (uint64_t)nrows * N;
        GRID_ASSERT( fread(&slab[0], sizeof(ComplexF), nelem, f) == nelem );
        fclose(f);
        loaded = true;
        std::cout << GridLogMessage << "DenseCoarseMatrix: slab loaded from "
                  << slabfile << " -- skipping import/factor/solve" << std::endl;
      } else {
        std::cout << GridLogMessage << "DenseCoarseMatrix: slab cache " << slabfile
                  << " absent -- full setup, will write it" << std::endl;
      }
    }
    if (!loaded) {
      ImportDense(Op);        // slab <- my rows of A   (LOCAL, no comms)
      ImportCertificate(Op);  // dense apply == Op.M, before inversion
      InvertDense(Op);        // slab <- my rows of A^{-1}
      double t1 = usecond();
      if (sfile) {
        FILE *f = fopen(slabfile.c_str(),"wb");
        GRID_ASSERT(f != nullptr);
        int64_t hdr[4] = { (int64_t)0x44454E5345, N, (int64_t)nrows, (int64_t)nbasis };
        GRID_ASSERT( fwrite(hdr,sizeof(int64_t),4,f) == 4 );
        uint64_t nelem = (uint64_t)nrows * N;
        GRID_ASSERT( fwrite(&slab[0], sizeof(ComplexF), nelem, f) == nelem );
        fclose(f);
        std::cout << GridLogMessage << "DenseCoarseMatrix: slab written to " << slabfile << std::endl;
      }
      std::cout << GridLogMessage << "DenseCoarseMatrix: import+invert took "
                << (t1-t0)/1.0e6 << " s" << std::endl;
    }

    ////////////////////////////////////////////////////////////////////
    // Device residency + persistent apply buffers + AOT split-K pointers
    ////////////////////////////////////////////////////////////////////
    {
      uint64_t sbytes = (uint64_t)nrows * N * sizeof(ComplexF);
      dSlab.resize((uint64_t)nrows*N);
      acceleratorCopyToDevice(&slab[0],&dSlab[0],sbytes);

      // DENSE_SPLITK: requested chunk count, snapped DOWN to a divisor of N.
      int req = getenv("DENSE_SPLITK") ? atoi(getenv("DENSE_SPLITK")) : 32;
      if (req < 1) req = 1;
      NK = 1;
      for(int j=1;j<=req;j++) if ( (N % j) == 0 ) NK = j;
      int64_t Kc = N / NK;

      dX.resize((uint64_t)N*MRHS_MAX);
      dY.resize((uint64_t)nrows*MRHS_MAX);
      dPartial.resize((uint64_t)NK*nrows*MRHS_MAX);
      hX.resize((uint64_t)N*MRHS_MAX);
      hY.resize((uint64_t)nrows*MRHS_MAX);

      aptrs.resize(NK); xptrs.resize(NK); cptrs.resize(NK);
      std::vector<ComplexF*> h(NK);
      for(int j=0;j<NK;j++) h[j] = &dSlab[0]    + (uint64_t)j*Kc;              // K-offset, lda=N
      acceleratorCopyToDevice(&h[0],&aptrs[0],NK*sizeof(ComplexF*));
      for(int j=0;j<NK;j++) h[j] = &dX[0]       + (uint64_t)j*Kc;              // K-offset, ldb=N
      acceleratorCopyToDevice(&h[0],&xptrs[0],NK*sizeof(ComplexF*));
      for(int j=0;j<NK;j++) h[j] = &dPartial[0] + (uint64_t)j*nrows*MRHS_MAX;  // compact, ldc=nrows
      acceleratorCopyToDevice(&h[0],&cptrs[0],NK*sizeof(ComplexF*));

      devSum = getenv("DENSE_DEVICE_SUM") ? atoi(getenv("DENSE_DEVICE_SUM")) : 0;
      const char *sumName[4] = {"host allreduce","DEVICE-buffer allreduce (GPU-aware MPI)",
                                "DEVICE cartesian ring allreduce (P2P)","DEVICE flat ring allreduce (P2P)"};
      GRID_ASSERT(devSum>=0 && devSum<=3);
      std::cout << GridLogMessage << "DenseCoarseMatrix: slab resident on device ("
                << sbytes/1024./1024. << " MB/rank), split-K NK=" << NK << " (Kc=" << Kc << "); "
                << sumName[devSum] << std::endl;
    }

    ////////////////////////////////////////////////////////////////////
    // VERIFY: || A (Ainv x) - x || / ||x|| through the DEVICE split-K core.
    ////////////////////////////////////////////////////////////////////
    {
      Field x(grid); Field y(grid); Field z(grid);
      x = ComplexD(1.0,0.0);
      double ta = usecond();
      (*this)(x, y);
      double tb = usecond();
      ApplyOracle(Op, y, z);
      z = z - x;
      RealD rel = std::sqrt(norm2(z)/norm2(x));
      std::cout << GridLogMessage << "DenseCoarseMatrix: VERIFY ||A Ainv x - x||/||x|| = "
                << rel << "   (one apply took " << (tb-ta)/1000.0 << " ms)" << std::endl;
      GRID_ASSERT(rel < 1.0e-2);
    }
    std::cout << GridLogMessage << "DenseCoarseMatrix: setup complete, total "
              << (usecond()-t0)/1.0e6 << " s" << std::endl;
  }

  ////////////////////////////////////////////////////////////////////
  // Apply the source operator to a D dimensional field, whichever kind it is.
  //
  // A multiRHS operator lives on the D+1 grid, so drive it with several right
  // hand sides at once: slice r carries (r+1)*in, and linearity says the
  // results must scale likewise. One apply, and unlike a single rhs check it
  // also catches rhs mixing. Cheap check, not a production path.
  ////////////////////////////////////////////////////////////////////
  template<class CoarseOp>
  void ApplyOracle(CoarseOp &Op,const Field &in, Field &out)
  {
    if ( Op.Grid() == grid ) { Op.M(in,out); return; }

    GridBase *mgrid = Op.Grid();
    GRID_ASSERT(mgrid->_ndimension == nd+1);
    int nr = mgrid->_fdimensions[0];

    Field min(mgrid), mout(mgrid);
    for(int r=0;r<nr;r++){
      Field scaled(grid);
      scaled = ComplexD(r+1.0,0.0)*in;
      InsertSliceFast(scaled,min,r,0);
    }

    Op.M(min,mout);

    ExtractSliceFast(out,mout,0,0);
    for(int r=1;r<nr;r++){
      Field sr(grid),d(grid);
      ExtractSliceFast(sr,mout,r,0);
      d = sr - ComplexD(r+1.0,0.0)*out;
      RealD rel = std::sqrt(norm2(d)/norm2(sr));
      if ( rel >= 1.0e-6 ) {
        std::cout << GridLogMessage << "DenseCoarseMatrix: oracle rhs "<<r
                  <<" inconsistent with rhs 0, rel "<<rel<<std::endl;
      }
      GRID_ASSERT( rel < 1.0e-6 );
    }
  }

  ////////////////////////////////////////////////////////////////////
  // 1. Direct stencil -> dense import of MY ROWS of A (no comms):
  //      Dense[(s,a),(wrap(s+shift_p),b)] += A[p][s]_{a,b}
  ////////////////////////////////////////////////////////////////////
  template<class CoarseOp>
  void ImportDense(CoarseOp &Op)
  {
    double t = -usecond();
    Coordinate gdims = grid->GlobalDimensions();
    int sign = getenv("DENSE_IMPORT_SIGN") ? atoi(getenv("DENSE_IMPORT_SIGN")) : 1;
    GRID_ASSERT( sign==1 || sign==-1 );

    uint64_t nelem = (uint64_t)nrows * N;
    thread_for(i, nelem, { slab[i] = ComplexF(0.0,0.0); });

    for(int p=0; p<Op.Geometry().npoint; p++){
      Coordinate shift = Op.Geometry().shifts[p];
      // _A[p] is PADDED after ExchangeCoarseLinks (end of CoarsenOperator):
      // extract the unpadded field before peeking with unpadded coordinates
      // (exactly as MultiGeneralCoarsenedMatrix::CopyMatrix does).
      CoarseMatrix Aun(grid);  Op.ExtractMatrix(p,Aun);
      if ( getenv("DENSE_IMPORT_DEBUG") ) {
        // Peek-path vs field-norm audit: sum |peekLocalSite|^2 must match norm2
        double pk = 0.0;
        autoView(Adbg, Aun, CpuRead);
        for(int ss=0; ss<lsites; ss++){
          Msobj m;
          peekLocalSite(m, Adbg, myLcoor[ss]);
          ComplexD *md = (ComplexD *)&m;
          for(int i=0; i<nbasis*nbasis; i++) pk += md[i].real()*md[i].real() + md[i].imag()*md[i].imag();
        }
        RealD gpk = pk;
        grid->GlobalSumVector(&gpk, 1);
        std::cout << GridLogMessage << "DenseCoarseMatrix: DEBUG p=" << p
                  << " norm2(_A[p]) " << norm2(Aun)
                  << " norm2(Extract) " << norm2(Aun)
                  << " sum|peek|^2 " << gpk << std::endl;
      }
      autoView(Av, Aun, CpuRead);
      thread_for(ss, lsites, {
        Coordinate ncoor(nd);
        for(int d=0; d<nd; d++){
          int64_t g = grid->_lstart[d] + myLcoor[ss][d] + sign*shift[d];
          ncoor[d] = (int)((g % gdims[d] + gdims[d]) % gdims[d]);
        }
        int64_t nsite;
        Lexicographic::IndexFromCoor(ncoor, nsite, gdims);
        Msobj m;
        peekLocalSite(m, Av, myLcoor[ss]);
        ComplexD *md = (ComplexD *)&m;
        // The operator contracts out(s,b) = sum_a A[p](s)(a,b) in(nbr,a)
        // (GeneralCoarsenedMatrix.h Mult kernel): the stored site matrix
        // acts TRANSPOSED, so element (a,b) lands at dense row (s,b),
        // column (nbr,a).  BUG LEDGER 2026-08-14: the original mapping
        // wrote (s,a),(nbr,b) -- caught by the IMPORT CERTIFICATE on its
        // FIRST fresh-import exercise (Test_schur_dense_coarse); every
        // production slab predates this path (probe-import SLAB_FILEs),
        // so no production output is suspect.
        for(int b=0; b<nbasis; b++){
          ComplexF *row = &slab[(uint64_t)(ss*nbasis+b)*N + nsite*nbasis];
          for(int a=0; a<nbasis; a++)
            row[a] += ComplexF(md[a*nbasis+b]);   // += : wrapped shifts may collide
        }
      });
    }
    t += usecond();

    // Structural diagnostic: identically-zero rows of my slab (a healthy
    // coarse operator has none; dead rows mean a rank-deficient import
    // or operator and the inverse will be NaN).
    int64_t zrows = 0;
    for(int64_t r=0; r<nrows; r++)
    {
      double mx = 0.0;
      const ComplexF *row = &slab[(uint64_t)r*N];
      for(int64_t j=0; j<N; j++)
      {
        mx = std::max(mx, (double)abs(row[j]));
      }
      if ( mx < 1.0e-30 ) zrows++;
    }
    RealD gz = (RealD)zrows;
    grid->GlobalSumVector(&gz, 1);

    std::cout << GridLogMessage << "DenseCoarseMatrix: stencil->dense import took "
              << t/1.0e6 << " s  (" << Op.Geometry().npoint << " points, local, no comms)"
              << "  zero rows " << (int64_t)gz << "/" << N << std::endl;

    // Debug: coordinate pattern of live sites (mechanism fingerprint)
    if ( (int64_t)gz > 0 )
    {
      int shown = 0;
      for(int ss=0; ss<lsites && shown<24; ss++)
      {
        double mx = 0.0;
        const ComplexF *row = &slab[(uint64_t)(ss*nbasis)*N];
        for(int64_t j=0; j<N; j++)
        {
          mx = std::max(mx, (double)abs(row[j]));
        }
        if ( mx > 1.0e-30 )
        {
          std::cout << GridLogMessage << "DenseCoarseMatrix: LIVE site ss=" << ss
                    << " lcoor " << myLcoor[ss] << std::endl;
          shown++;
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////
  // 2. IMPORT CERTIFICATE: dense rows vs Op.M on a NON-CONSTANT vector.
  //    (Constant x has x[s+d]==x[s-d]: blind to a shift-sign error.)
  ////////////////////////////////////////////////////////////////////
  template<class CoarseOp>
  void ImportCertificate(CoarseOp &Op)
  {
    Field x(grid); Field Ax(grid); Field Dx(grid);
    for(int ss=0; ss<lsites; ss++){
      sobj s;
      for(int b=0; b<nbasis; b++){
        double ph = 0.37*(double)(myGsite[ss]*nbasis+b);
        ((ComplexD *)&s)[b] = ComplexD(std::cos(ph),std::sin(0.61*ph));
      }
      pokeLocalSite(s, x, myLcoor[ss]);
    }
    // gather full x (zero-fill + exact GlobalSum), dense rows on host
    std::vector<ComplexD> xh((uint64_t)N, ComplexD(0.0,0.0));
    for(int ss=0; ss<lsites; ss++){
      sobj s;
      peekLocalSite(s, x, myLcoor[ss]);
      for(int b=0; b<nbasis; b++) xh[ myGsite[ss]*nbasis + b ] = ((ComplexD *)&s)[b];
    }
    grid->GlobalSumVector(&xh[0], (int)N);
    std::vector<ComplexD> yh(nrows);
    thread_for(r, nrows, {
      ComplexD acc(0.0,0.0);
      const ComplexF *row = &slab[(uint64_t)r * N];
      for(int64_t j=0; j<N; j++) acc += ComplexD(row[j]) * xh[j];
      yh[r] = acc;
    });
    for(int ss=0; ss<lsites; ss++){
      sobj s;
      for(int b=0; b<nbasis; b++) ((ComplexD *)&s)[b] = yh[ss*nbasis+b];
      pokeLocalSite(s, Dx, myLcoor[ss]);
    }
    ApplyOracle(Op, x, Ax);
    Field d(grid); d = Dx - Ax;
    RealD rel = std::sqrt(norm2(d)/norm2(Ax));
    std::cout << GridLogMessage << "DenseCoarseMatrix: IMPORT CERTIFICATE ||Dense x - A x||/||A x|| = "
              << rel << std::endl;
    if ( rel >= 1.0e-3 ) {
      std::cout << GridLogMessage << "DenseCoarseMatrix: IMPORT CERTIFICATE FAILED. If O(1), the "
                << "stencil shift-sign convention is opposite: rerun with DENSE_IMPORT_SIGN=-1"
                << std::endl;
    }
    GRID_ASSERT(rel < 1.0e-3);
  }

  ////////////////////////////////////////////////////////////////////
  // 3. Invert dispatcher.  slab holds my rows of A on entry, my rows
  //    of A^{-1} on exit, in every mode.
  //      DENSE_SCHUR absent/0 : single-GCD gather-invert (the oracle)
  //      DENSE_SCHUR=1        : distributed recursive Schur
  //      DENSE_SCHUR=2        : AUDIT -- both on the same A; report the
  //                             slab difference; keep the Schur result
  //                             (so VERIFY certifies the new path).
  ////////////////////////////////////////////////////////////////////
  template<class CoarseOp>
  void InvertDense(CoarseOp &Op)
  {
    char *sc  = getenv("DENSE_SCHUR");
    int  mode = sc ? atoi(sc) : 0;

    if ( mode == 0 )
    {
      InvertDenseSingle();
      return;
    }
    if ( mode == 1 )
    {
      InvertDenseSchur(Op);
      return;
    }
    GRID_ASSERT( mode == 2 );
    std::vector<ComplexF> Aimp(slab);       // imported A
    InvertDenseSingle();
    std::vector<ComplexF> ref(slab);        // Ainv, single path
    slab = Aimp;
    InvertDenseSchur(Op);                   // slab = Ainv, Schur path

    // NaN-PROOF comparison: max() masks NaN, so count non-finite
    // entries in each result explicitly.
    double  mx = 0.0;
    double  mr = 0.0;
    int64_t badschur  = 0;
    int64_t badsingle = 0;
    for(uint64_t i=0; i<(uint64_t)nrows*N; i++)
    {
      double as = abs(ComplexD(slab[i]));
      double ar = abs(ComplexD(ref[i]));
      if ( !std::isfinite(as) ) badschur++;
      if ( !std::isfinite(ar) ) badsingle++;
      if ( std::isfinite(as) && std::isfinite(ar) )
      {
        mx = std::max(mx, (double)abs(ComplexD(slab[i]) - ComplexD(ref[i])));
        mr = std::max(mr, ar);
      }
    }
    RealD gmx = mx;
    RealD gmr = mr;
    RealD gbs = (RealD)badschur;
    RealD gbr = (RealD)badsingle;
    grid->GlobalMax(gmx);
    grid->GlobalMax(gmr);
    grid->GlobalSumVector(&gbs, 1);
    grid->GlobalSumVector(&gbr, 1);
    schurAuditRel = gmx/gmr;
    std::cout << GridLogMessage << "DenseCoarseMatrix: DENSE_SCHUR=2 AUDIT "
              << "max|Ainv_schur - Ainv_single| = " << gmx
              << "  relative " << schurAuditRel
              << "  non-finite: schur " << (int64_t)gbs << " single " << (int64_t)gbr
              << "  (two fp32 roundings of the same inverse; expect ~ growth * eps32)"
              << std::endl;
    GRID_ASSERT( gbs == 0 );
    GRID_ASSERT( gbr == 0 );
  }

  ////////////////////////////////////////////////////////////////////
  // 3a. Single-GCD invert: chunked zero-fill+GlobalSum gather of A
  //    streamed to the boss GCD, cgetrf_64 (ILP64), rows of A^{-1} by
  //    blocked identity cgetrs_64 + broadcast; each rank keeps its own
  //    rows (in `slab`, overwriting A).  Proven path; the SCHUR oracle.
  ////////////////////////////////////////////////////////////////////
  void InvertDenseSingle(void)
  {
    double t1 = usecond();
    int boss = grid->IsBoss();
    std::vector<ComplexF> Afull;
#ifdef GRID_HIP
    // QUARANTINED naked HIP: the boss-only N^2 inversion buffer (34GB at
    // N=65536) must come from raw HBM; EvictAll flushes the device-copy
    // layer to make the window.  (FreePool of the allocator free-list
    // awaits the type-dispatched fix.)  Confined to setup; the apply path
    // is pure Grid primitives.
    rocblas_float_complex *dA = nullptr;
    rocblas_float_complex *dB = nullptr;
    int64_t *dIpiv = nullptr;
    uint64_t Abytes = (uint64_t)N * N * sizeof(ComplexF);
    MemoryManager::EvictAll();
    if (boss) {
      auto aerr = hipMalloc((void **)&dA, Abytes);
      if (aerr != hipSuccess) {
        std::cout << GridLogMessage << "DenseCoarseMatrix: hipMalloc of "
                  << Abytes/1024./1024./1024. << " GB FAILED -- reduce --device-mem" << std::endl;
        GRID_ASSERT(aerr == hipSuccess);
      }
      std::cout << GridLogMessage << "DenseCoarseMatrix: device inversion buffer allocated ("
                << Abytes/1024./1024./1024. << " GB)" << std::endl;
    }
#else
    if (boss) Afull.resize((uint64_t)N * N);
#endif
    {
      std::unordered_map<int64_t,int> rowmap;   // global row -> my slab row
      for(int ss=0; ss<lsites; ss++)
        for(int a=0; a<nbasis; a++)
          rowmap[ myGsite[ss]*nbasis + a ] = ss*nbasis + a;

      std::vector<ComplexF> chunk((uint64_t)CHUNKROWS * N);
      for(int64_t row0=0; row0<N; row0+=CHUNKROWS){
        int64_t nrow = std::min(CHUNKROWS, N-row0);
        uint64_t nelem = (uint64_t)nrow * N;
        for(uint64_t i=0;i<nelem;i++) chunk[i]=ComplexF(0.0,0.0);
        for(int64_t r=row0; r<row0+nrow; r++){
          auto it = rowmap.find(r);
          if (it != rowmap.end()) {
            uint64_t src = (uint64_t)(it->second) * N;
            uint64_t dst = (uint64_t)(r-row0) * N;
            for(int64_t j=0;j<N;j++) chunk[dst+j] = slab[src+j];
          }
        }
        grid->GlobalSumVector(&chunk[0], (int)nelem);
        if (boss) {
#ifdef GRID_HIP
          GRID_ASSERT( hipMemcpy((char *)dA + (uint64_t)row0*N*sizeof(ComplexF),
                                 &chunk[0], nelem*sizeof(ComplexF),
                                 hipMemcpyHostToDevice) == hipSuccess );
#else
          uint64_t dst = (uint64_t)row0 * N;
          for(uint64_t i=0;i<nelem;i++) Afull[dst+i] = chunk[i];
#endif
        }
      }
    }
    double t2 = usecond();
    std::cout << GridLogMessage << "DenseCoarseMatrix: gather to boss took "
              << (t2-t1)/1.0e6 << " s" << std::endl;

    ////////////////////////////////////////////////////////////////////
    // Factor in place on the boss (fp32, ILP64).  Row-major buffer handed
    // to column-major LAPACK => LU of A^T.
    ////////////////////////////////////////////////////////////////////
    if (boss) {
#ifdef GRID_HIP
      std::cout << GridLogMessage << "DenseCoarseMatrix: rocSOLVER cgetrf_64 (ILP64 LU) N=" << N
                << " in place on resident device buffer" << std::endl;
      rocblas_handle handle = GridBLASInverse::Handle();
      int64_t *dInfo;
      GRID_ASSERT( hipMalloc((void **)&dIpiv, N*sizeof(int64_t)) == hipSuccess );
      GRID_ASSERT( hipMalloc((void **)&dInfo, sizeof(int64_t))   == hipSuccess );
      auto st1 = rocsolver_cgetrf_64(handle, (int64_t)N, (int64_t)N, dA, (int64_t)N, dIpiv, dInfo);
      GRID_ASSERT( hipDeviceSynchronize() == hipSuccess );
      int64_t info_h = -1;
      GRID_ASSERT( hipMemcpy(&info_h, dInfo, sizeof(int64_t),
                             hipMemcpyDeviceToHost) == hipSuccess );
      std::cout << GridLogMessage << "DenseCoarseMatrix: cgetrf_64 status " << (int)st1
                << " info = " << (int)info_h << std::endl;
      GRID_ASSERT(st1 == rocblas_status_success);
      GRID_ASSERT(info_h == 0);
      GRID_ASSERT( hipFree(dInfo) == hipSuccess );
      GRID_ASSERT( hipMalloc((void **)&dB, (uint64_t)CHUNKROWS*N*sizeof(ComplexF)) == hipSuccess );
      // dA holds the LU of A^T; rows of A^{-1} are produced blockwise below via
      // cgetrs_64 on identity-column blocks: A^T X = E => X columns = rows of
      // A^{-1}, in exactly the linear layout the harvest expects.
#else
      // Eigen fallback: small local CPU tests only.
      std::cout << GridLogMessage << "DenseCoarseMatrix: Eigen fallback inversion N=" << N
                << (N > 10000 ? "  (WARNING: SLOW; use the HIP/rocSOLVER path)" : "")
                << std::endl;
      typedef Eigen::Matrix<std::complex<float>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> MatF;
      Eigen::Map<MatF> A(reinterpret_cast<std::complex<float>*>(&Afull[0]), N, N);
      MatF Ainv = A.inverse();
      A = Ainv;
#endif
    }
    double t3 = usecond();
    std::cout << GridLogMessage << "DenseCoarseMatrix: factorisation took "
              << (t3-t2)/1.0e6 << " s" << std::endl;

    ////////////////////////////////////////////////////////////////////
    // Blocked solve + broadcast: rows of A^{-1} chunk by chunk; each
    // rank keeps the rows of its own sites (ownership-aligned).
    ////////////////////////////////////////////////////////////////////
    {
      std::unordered_map<int64_t,int> rowmap;
      for(int ss=0; ss<lsites; ss++)
        for(int a=0; a<nbasis; a++)
          rowmap[ myGsite[ss]*nbasis + a ] = ss*nbasis + a;

      std::vector<ComplexF> chunk((uint64_t)CHUNKROWS * N);
      for(int64_t row0=0; row0<N; row0+=CHUNKROWS){
        int64_t nrow = std::min(CHUNKROWS, N-row0);
        uint64_t nelem = (uint64_t)nrow * N;
        if (boss) {
#ifdef GRID_HIP
          // Identity block E: column j = e_{row0+j}; solve A^T X = E so X's
          // columns are rows [row0,row0+nrow) of A^{-1}.
          for(uint64_t i=0;i<nelem;i++) chunk[i] = ComplexF(0.0,0.0);
          for(int64_t j=0;j<nrow;j++) chunk[(uint64_t)j*N + (uint64_t)(row0+j)] = ComplexF(1.0,0.0);
          GRID_ASSERT( hipMemcpy(dB, &chunk[0], nelem*sizeof(ComplexF), hipMemcpyHostToDevice) == hipSuccess );
          auto strs = rocsolver_cgetrs_64(GridBLASInverse::Handle(), rocblas_operation_none,
                                          (int64_t)N, (int64_t)nrow,
                                          dA, (int64_t)N, dIpiv, dB, (int64_t)N);
          GRID_ASSERT(strs == rocblas_status_success);
          GRID_ASSERT( hipDeviceSynchronize() == hipSuccess );
          GRID_ASSERT( hipMemcpy(&chunk[0], dB, nelem*sizeof(ComplexF), hipMemcpyDeviceToHost) == hipSuccess );
#else
          uint64_t src = (uint64_t)row0 * N;
          for(uint64_t i=0;i<nelem;i++) chunk[i] = Afull[src+i];
#endif
        }
        grid->Broadcast(0, &chunk[0], nelem*sizeof(ComplexF));
        for(int64_t r=row0; r<row0+nrow; r++){
          auto it = rowmap.find(r);
          if (it != rowmap.end()) {
            uint64_t dst = (uint64_t)(it->second) * N;
            uint64_t src = (uint64_t)(r-row0) * N;
            for(int64_t j=0;j<N;j++) slab[dst+j] = chunk[src+j];
          }
        }
      }
    }
#ifdef GRID_HIP
    if (boss) {
      if (dA)    GRID_ASSERT( hipFree(dA)    == hipSuccess );
      if (dB)    GRID_ASSERT( hipFree(dB)    == hipSuccess );
      if (dIpiv) GRID_ASSERT( hipFree(dIpiv) == hipSuccess );
    }
#endif
    double t4 = usecond();
    std::cout << GridLogMessage << "DenseCoarseMatrix: blocked getrs solve+scatter took "
              << (t4-t3)/1.0e6 << " s" << std::endl;
  }

  ////////////////////////////////////////////////////////////////////
  // 3b. Global column -> rank-major column map, computed LOCALLY.
  //    Rank-major ordering: rank q's rows/columns are the contiguous
  //    block [q*nrows, (q+1)*nrows), ordered by q's local site index
  //    (uniform local volumes make ownership arithmetic exact).
  //    MPI_Cart_rank is queried ONCE PER RANK (serial, P calls) into a
  //    lex-processor table; the per-site sweep is then pure arithmetic.
  ////////////////////////////////////////////////////////////////////
  void BuildRankMajorMap(std::vector<int64_t> &g2rm)
  {
    int P = grid->ProcessorCount();
    Coordinate pdims = grid->_processors;
    Coordinate gdims = grid->GlobalDimensions();
    Coordinate ldims = grid->LocalDimensions();

    std::vector<int> lexp2rank(P);
    for(int lp=0; lp<P; lp++)
    {
      Coordinate pcoor(nd);
      Lexicographic::CoorFromIndex(pcoor, lp, pdims);
      lexp2rank[lp] = grid->RankFromProcessorCoor(pcoor);
    }

    int64_t gsites = grid->gSites();
    g2rm.resize(N);
    thread_for(gsite, gsites, {
      Coordinate gcoor(nd);
      Coordinate pcoor(nd);
      Coordinate lcoor(nd);
      Lexicographic::CoorFromIndex(gcoor, gsite, gdims);
      for(int d=0; d<nd; d++)
      {
        pcoor[d] = gcoor[d]/ldims[d];
        lcoor[d] = gcoor[d]-pcoor[d]*ldims[d];
      }
      int64_t lexp;
      int64_t lsite;
      Lexicographic::IndexFromCoor(pcoor, lexp,  pdims);
      Lexicographic::IndexFromCoor(lcoor, lsite, ldims);
      int64_t base = (int64_t)lexp2rank[lexp]*nrows + lsite*nbasis;
      for(int b=0; b<nbasis; b++)
      {
        g2rm[(uint64_t)gsite*nbasis + b] = base + b;
      }
    });
  }

  ////////////////////////////////////////////////////////////////////
  // 3c. Direct stencil -> fp64 rank-major import of MY ROWS of A (the
  //    end-to-end fp64 path: the stencil source IS ComplexD; nothing is
  //    rounded through fp32 on the way into the inversion).  Same
  //    loop/sign/accumulate/transposed-contraction discipline as
  //    ImportDense; output is column-major rows x N with columns in
  //    rank-major order (g2rm).
  //    ALWAYS-ON CERTIFICATE: the fp64 import, rounded, must agree with
  //    the fp32 slab entry at the corresponding global column, over the
  //    WHOLE of my rows (few ulp: wrapped-shift collisions accumulate in
  //    different precision order).  NaN-proof: non-finite entries are
  //    counted explicitly since max() silently masks NaN.
  ////////////////////////////////////////////////////////////////////
  template<class CoarseOp>
  void ImportDenseFP64(CoarseOp &Op, BlockRows &S, std::vector<int64_t> &g2rm)
  {
    Coordinate gdims = grid->GlobalDimensions();
    int sign = getenv("DENSE_IMPORT_SIGN") ? atoi(getenv("DENSE_IMPORT_SIGN")) : 1;
    GRID_ASSERT( sign==1 || sign==-1 );

    std::vector<ComplexD> h((uint64_t)nrows*N, ComplexD(0.0,0.0));
    for(int p=0; p<Op.Geometry().npoint; p++)
    {
      Coordinate shift = Op.Geometry().shifts[p];
      CoarseMatrix Aun(grid);  Op.ExtractMatrix(p,Aun);
      autoView(Av, Aun, CpuRead);
      thread_for(ss, lsites, {
        Coordinate ncoor(nd);
        for(int d=0; d<nd; d++)
        {
          int64_t g = grid->_lstart[d] + myLcoor[ss][d] + sign*shift[d];
          ncoor[d] = (int)((g % gdims[d] + gdims[d]) % gdims[d]);
        }
        int64_t nsite;
        Lexicographic::IndexFromCoor(ncoor, nsite, gdims);
        Msobj m;
        peekLocalSite(m, Av, myLcoor[ss]);
        ComplexD *md = (ComplexD *)&m;
        // Transposed contraction as ImportDense: (a,b) lands at
        // row (s,b), column (nbr,a); column index in rank-major order.
        for(int a=0; a<nbasis; a++)
        {
          int64_t jj = g2rm[ nsite*nbasis + a ];
          for(int b=0; b<nbasis; b++)
          {
            h[(uint64_t)(ss*nbasis+b) + (uint64_t)jj*nrows] += md[a*nbasis+b];
          }
        }
      });
    }

    // Certificate vs the fp32 slab (slab holds A at this point)
    double  mx   = 0.0;
    int64_t nbad = 0;
    for(int64_t i=0; i<nrows; i++)
    {
      for(int64_t gcol=0; gcol<N; gcol++)
      {
        ComplexD d64 = h[(uint64_t)(i + g2rm[gcol]*nrows)];
        ComplexF f32 = slab[(uint64_t)i*N + gcol];
        double dev = abs(ComplexD(f32) - d64);
        if ( !std::isfinite(dev) ) nbad++;
        else mx = std::max(mx, dev);
      }
    }
    RealD gmx  = mx;
    RealD gbad = (RealD)nbad;
    grid->GlobalMax(gmx);
    grid->GlobalSumVector(&gbad, 1);
    std::cout << GridLogMessage << "DenseCoarseMatrix: fp64 import certificate "
              << "max|A64 - A32| = " << gmx
              << "  non-finite entries " << (int64_t)gbad << std::endl;
    GRID_ASSERT( gbad == 0 );
    GRID_ASSERT( gmx < 1.0e-5 );

    S.Resize(nrows, N);
    acceleratorCopyToDevice(&h[0], &S.data[0], (uint64_t)nrows*N*sizeof(ComplexD));
  }

  ////////////////////////////////////////////////////////////////////
  // 3d. Distributed recursive Schur invert, END-TO-END fp64 (decision
  //    2026-08-14): stencil (ComplexD) -> fp64 rank-major import ->
  //    fp64 recursion -> ONE terminal rounding into the fp32 apply
  //    slab.  Everything downstream (device residency, split-K apply,
  //    VERIFY, SLAB_FILE) is untouched.
  ////////////////////////////////////////////////////////////////////
  template<class CoarseOp>
  void InvertDenseSchur(CoarseOp &Op)
  {
    double t1 = usecond();
    int P  = grid->ProcessorCount();
    int me = grid->ThisRank();

    // Uniform local volumes => contiguous uniform ownership
    std::vector<int64_t> rowStart(P+1);
    for(int r=0; r<=P; r++)
    {
      rowStart[r] = (int64_t)r*nrows;
    }
    GRID_ASSERT( rowStart[P] == N );

    std::vector<int64_t> g2rm;
    BuildRankMajorMap(g2rm);

    // Self-certifying map: my own global rows land at my rank-major slots
    for(int ss=0; ss<lsites; ss++)
    {
      for(int a=0; a<nbasis; a++)
      {
        GRID_ASSERT( g2rm[ myGsite[ss]*nbasis + a ] == (int64_t)me*nrows + ss*nbasis + a );
      }
    }

    BlockRows S;
    ImportDenseFP64(Op, S, g2rm);

    ////////////////////////////////////////////////////////////////
    // DENSE_SCHUR2D=1 : invert via the 2D block-cyclic recursion
    // (BlockCyclicSchurInverse) instead of the 1D rank-range one.
    // The SAME imported rank-major rows S go in and come back, so the
    // import certificate above and the slab rounding / VERIFY below are
    // identical for both paths: a clean A/B on one imported operator.
    //
    // Everything in the 2D path -- redistribution, SUMMA rings, leaf --
    // is point-to-point SendToRecvFrom; no collectives at all.
    // DENSE_NB overrides the block size (default: rows-per-rank, which
    // makes the redistribution edges maximally regular).
    ////////////////////////////////////////////////////////////////
    int use2d = getenv("DENSE_SCHUR2D") ? atoi(getenv("DENSE_SCHUR2D")) : 0;
    int64_t panelBytes = getenv("DENSE_PANEL_BYTES") ? atol(getenv("DENSE_PANEL_BYTES"))
                                                     : (int64_t)1024*1024*1024;   // 1D path only
    double t2, t3;
    if ( use2d )
    {
      int Pr,Pc;
      BlockCyclicLayout::ChooseProcessGrid(P, Pr, Pc);
      int64_t nb = getenv("DENSE_NB") ? atol(getenv("DENSE_NB")) : nrows;
      GRID_ASSERT( nb >= 1 );
      std::cout << GridLogMessage << "DenseCoarseMatrix: 2D SCHUR invert, process grid "
                << Pr << " x " << Pc << "  nb " << nb
                << "  (pure P2P: redistribute + SUMMA rings + local leaves)" << std::endl;
      BlockCyclicMatrix A2(grid, N, nb, Pr, Pc);
      BlockCyclicSchurInverse RSI2;
      t2 = usecond();
      BlockCyclicRedistribute::RowsToCyclic(grid, rowStart, &S.data[0], nrows, A2);
      RSI2.Invert(A2);
      BlockCyclicRedistribute::CyclicToRows(grid, rowStart, A2, &S.data[0], nrows);
      t3 = usecond();
      RSI2.ReportTelemetry(grid);
    }
    else
    {
      RecursiveSchurInverse RSI(grid, N, rowStart, panelBytes);
      t2 = usecond();
      RSI.Invert(S);
      t3 = usecond();
      RSI.ReportTelemetry();
    }

    // The single terminal rounding: fp64 inverse -> fp32 apply slab
    // (row-major, global columns)
    {
      std::vector<ComplexD> h((uint64_t)nrows*N);
      acceleratorCopyFromDevice(&S.data[0], &h[0], (uint64_t)nrows*N*sizeof(ComplexD));
      thread_for(gcol, N, {
        int64_t jj = g2rm[gcol];
        for(int64_t i=0; i<nrows; i++)
        {
          slab[(uint64_t)i*N + gcol] = ComplexF(h[(uint64_t)(i + jj*nrows)]);
        }
      });
    }
    double t4 = usecond();
    std::cout << GridLogMessage << "DenseCoarseMatrix: SCHUR fp64 distributed invert took "
              << (t4-t1)/1.0e6 << " s (recursion " << (t3-t2)/1.0e6 << " s), panelBytes "
              << panelBytes << std::endl;
  }

  ////////////////////////////////////////////////////////////////////
  // CORE apply on packed data: hX[N x nr] (zero-filled, my sites only)
  // -> allreduce -> split-K GEMM against the resident slab -> reduce
  // partials -> hY[nrows x nr] (column major).  Platform-agnostic:
  // deviceVector + GridBLAS (Eigen fallback on CPU builds).
  // fp32 allreduce is EXACT: zero-fill assembly gives every element
  // exactly one contributing rank.
  ////////////////////////////////////////////////////////////////////
  void SlabApplyPacked(int nr, double *tprof)
  {
    GRID_TRACE("DenseSlabApply");
    GRID_ASSERT(nr <= MRHS_MAX);
    uint64_t nX = (uint64_t)N * nr;
    uint64_t nY = (uint64_t)nrows * nr;
    int64_t  Kc = N / NK;
    double t1 = usecond();
    double t2, t3;
    if (devSum) {
      { GRID_TRACE("DenseH2D");
        acceleratorCopyToDevice(&hX[0],&dX[0],nX*sizeof(ComplexF));
      }
      t2 = usecond();
      { GRID_TRACE("DenseAllreduce");
        // DENSE_DEVICE_SUM=1 : device-buffer MPI_Allreduce (Cray MPICH aborts
        //                      above ~8 MB: 12 RHS at N=138240 is 13.3 MB)
        // DENSE_DEVICE_SUM=2 : CartesianRingAllReduce, P2P only, no size cliff
        // DENSE_DEVICE_SUM=3 : flat RingAllReduce, P2P only
        if      (devSum==2) CartesianRingAllReduce(grid,(ComplexF *)&dX[0],nX);
        else if (devSum==3) RingAllReduce(grid,(ComplexF *)&dX[0],nX);
        else                grid->GlobalSumVector((ComplexF *)&dX[0], (int)nX);
      }
      t3 = usecond();
    } else {
      { GRID_TRACE("DenseAllreduce");
        grid->GlobalSumVector(&hX[0], (int)nX);
      }
      t2 = usecond();
      { GRID_TRACE("DenseH2D");
        acceleratorCopyToDevice(&hX[0],&dX[0],nX*sizeof(ComplexF));
      }
      t3 = usecond();
    }
    // Y = op(slab,T) . X : row-major slab (nrows x N) == col-major A^T
    // (N x nrows, lda=N) => transpose gives the nrows x N operator.
    // Split-K: NK chunk-GEMMs by pointer offset (AOT lists), then reduce.
    ComplexF one (1.0,0.0);
    ComplexF zero(0.0,0.0);
    { GRID_TRACE("DenseSplitKGEMM");
      BLAS.gemmBatched(GridBLAS_OP_T, GridBLAS_OP_N,
                       (int)nrows, nr, (int)Kc,
                       one,  aptrs, (int)N,
                             xptrs, (int)N,
                       zero, cptrs, (int)nrows);
      BLAS.synchronise();
      ComplexF *pp = &dPartial[0];
      ComplexF *py = &dY[0];
      uint64_t stride = (uint64_t)nrows*MRHS_MAX;
      int nk = NK;
      accelerator_for(i, nY, 1, {
        ComplexF acc(0.0,0.0);
        for(int j=0;j<nk;j++) acc += pp[(uint64_t)j*stride + i];
        py[i] = acc;
      });
    }
    double t4 = usecond();
    { GRID_TRACE("DenseD2H");
      acceleratorCopyFromDevice(&dY[0],&hY[0],nY*sizeof(ComplexF));
    }
    double t5 = usecond();
    if (tprof) {
      tprof[0] = devSum ? (t3-t2) : (t2-t1);   // allreduce
      tprof[1] = devSum ? (t2-t1) : (t3-t2);   // H2D
      tprof[2] = t4-t3;                        // gemm+reduce
      tprof[3] = t5-t4;                        // D2H
    }
  }

  ////////////////////////////////////////////////////////////////////
  // Single-RHS apply (also the VERIFY path => certifies device slab).
  ////////////////////////////////////////////////////////////////////
  virtual void operator()(const Field &src, Field &psi)
  {
    GRID_TRACE("DenseApply1");
    uint64_t nX = (uint64_t)N;
    { GRID_TRACE("DensePack");
      thread_for(i, nX, { hX[i]=ComplexF(0.0,0.0); });
      for(int ss=0; ss<lsites; ss++){
        sobj s;
        peekLocalSite(s, src, myLcoor[ss]);
        for(int b=0; b<nbasis; b++)
          hX[ myGsite[ss]*nbasis + b ] = ComplexF(((ComplexD *)&s)[b]);
      }
    }
    SlabApplyPacked(1, nullptr);
    { GRID_TRACE("DenseUnpack");
      for(int ss=0; ss<lsites; ss++){
        sobj s;
        for(int b=0; b<nbasis; b++)
          ((ComplexD *)&s)[b] = ComplexD(hY[ss*nbasis + b]);
        pokeLocalSite(s, psi, myLcoor[ss]);
      }
    }
  }

  ////////////////////////////////////////////////////////////////////
  // Defect of an applied inverse. Was a DENSE_CC_CHECK block inside
  // operator(), but that is virtual and cannot take an operator, so the
  // caller now asks for it explicitly.
  ////////////////////////////////////////////////////////////////////
  template<class CoarseOp>
  void CheckApply(CoarseOp &Op,const Field &src,const Field &psi)
  {
    Field tmp(grid);
    Op.M(psi, tmp);
    tmp = tmp - src;
    std::cout << GridLogMessage << "DenseCoarseMatrix: apply defect ||A x - b||/||b|| = "
              << std::sqrt(norm2(tmp)/norm2(src)) << std::endl;
  }

  ////////////////////////////////////////////////////////////////////
  // Batched (vector-of-fields) apply.
  ////////////////////////////////////////////////////////////////////
  void ApplyBatch(std::vector<Field> &src, std::vector<Field> &psi)
  {
    int nr = src.size();
    GRID_ASSERT(nr <= MRHS_MAX);
    double t0 = usecond();
    uint64_t nX = (uint64_t)N*nr;
    { GRID_TRACE("DensePack");
      thread_for(i, nX, { hX[i]=ComplexF(0.0,0.0); });
      for(int rr=0; rr<nr; rr++){
        for(int ss=0; ss<lsites; ss++){
          sobj s;
          peekLocalSite(s, src[rr], myLcoor[ss]);
          for(int b=0; b<nbasis; b++)
            hX[ (uint64_t)rr*N + myGsite[ss]*nbasis + b ] = ComplexF(((ComplexD *)&s)[b]);
        }
      }
    }
    SlabApplyPacked(nr, nullptr);
    { GRID_TRACE("DenseUnpack");
      for(int rr=0; rr<nr; rr++){
        for(int ss=0; ss<lsites; ss++){
          sobj s;
          for(int b=0; b<nbasis; b++)
            ((ComplexD *)&s)[b] = ComplexD(hY[(uint64_t)rr*nrows + (ss*nbasis+b)]);
          pokeLocalSite(s, psi[rr], myLcoor[ss]);
        }
      }
    }
    double t1 = usecond();
    std::cout << GridLogMessage << "DenseCoarseMatrix: batched apply " << nr << " rhs took "
              << (t1-t0)/1000.0 << " ms  (" << (t1-t0)/1000.0/nr << " ms/rhs)" << std::endl;
  }

  template<class CoarseOp>
  void CheckApplyBatch(CoarseOp &Op,std::vector<Field> &src,std::vector<Field> &psi,int nr)
  {
    Field tmp(grid);
    for(int rr=0; rr<nr; rr++){
      Op.M(psi[rr], tmp);
      tmp = tmp - src[rr];
      std::cout << GridLogMessage << "DenseCoarseMatrix: batch defect["<<rr<<"] = "
                << std::sqrt(norm2(tmp)/norm2(src[rr])) << std::endl;
    }
  }

  ////////////////////////////////////////////////////////////////////
  // 6D mrhs apply: operates DIRECTLY on the packed 6D field (rhs = dim 0).
  ////////////////////////////////////////////////////////////////////
  void ApplyBatch6D(const Field &in6, Field &out6, int nr)
  {
    GRID_ASSERT(nr <= MRHS_MAX);
    GRID_ASSERT(in6.Grid()->_ndimension == nd+1);   // {rhs, s, x,y,z,t}
    double t0 = usecond();
    Field &in = const_cast<Field &>(in6);
    uint64_t nX = (uint64_t)N * nr;
    thread_for(i, nX, { hX[i]=ComplexF(0.0,0.0); });
    { GRID_TRACE("DensePack");
      autoView(iv, in, CpuRead);
      Coordinate c6(nd+1);
      for(int ss=0; ss<lsites; ss++){
        for(int d=0; d<nd; d++) c6[d+1] = myLcoor[ss][d];
        for(int rr=0; rr<nr; rr++){
          c6[0] = rr;
          sobj s;
          peekLocalSite(s, iv, c6);
          for(int b=0; b<nbasis; b++)
            hX[(uint64_t)rr*N + myGsite[ss]*nbasis + b] = ComplexF(((ComplexD *)&s)[b]);
        }
      }
    }
    double t1 = usecond();
    double tprof[4];
    SlabApplyPacked(nr, tprof);
    double t5 = usecond();
    { GRID_TRACE("DenseUnpack");
      autoView(ov, out6, CpuWrite);
      Coordinate c6(nd+1);
      for(int ss=0; ss<lsites; ss++){
        for(int d=0; d<nd; d++) c6[d+1] = myLcoor[ss][d];
        for(int rr=0; rr<nr; rr++){
          c6[0] = rr;
          sobj s;
          for(int b=0; b<nbasis; b++)
            ((ComplexD *)&s)[b] = ComplexD(hY[(uint64_t)rr*nrows + (ss*nbasis+b)]);  // Y col-major
          pokeLocalSite(s, ov, c6);
        }
      }
    }
    double t6 = usecond();
    std::cout << GridLogMessage << "DenseCoarseMatrix: apply6D " << nr << " rhs took "
              << (t6-t0)/1000.0 << " ms" << std::endl;
    if ( getenv("DENSE_APPLY_PROFILE") ) {
      std::cout << GridLogMessage << "DenseCoarseMatrix: apply6D profile:"
                << " pack "        << (t1-t0)/1000.0
                << "  allreduce "  << tprof[0]/1000.0
                << "  H2D "        << tprof[1]/1000.0
                << "  gemm+reduce "<< tprof[2]/1000.0
                << "  D2H "        << tprof[3]/1000.0
                << "  unpack "     << (t6-t5)/1000.0
                << "  ms" << std::endl;
    }
  }
};

NAMESPACE_END(Grid);
