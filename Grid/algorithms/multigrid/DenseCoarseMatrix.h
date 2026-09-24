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
#include <Grid/algorithms/multigrid/BlockCyclicSchurInverse.h>
#include <Grid/algorithms/multigrid/BlockCyclicRedistribute.h>

#include <unordered_map>

NAMESPACE_BEGIN(Grid);

//////////////////////////////////////////////////////////////////////////////////////
// DenseCoarseMatrix: a coarsened operator treated as a DENSE matrix -- explicit,
// row-distributed A^{-1} of a GeneralCoarsenedMatrix.
//
//  - Stencil -> dense DIRECT IMPORT.  The coarse operator IS the dense matrix
//    unrolled: Dense[(s,a),(s+shift_p,b)] += A[p][s]_{a,b}.  Rows of my sites are
//    assembled from purely LOCAL _A[p] data: no operator applies, no comms.
//    ACCUMULATE (+=) because on short axes distinct shifts wrap to the same
//    neighbour.  An IMPORT CERTIFICATE compares the dense apply against Op.M on a
//    NON-CONSTANT vector (a constant one cannot see a shift-sign error).
//
//  - Inversion is END-TO-END fp64 through the 2D block-cyclic recursive Schur
//    complement (BlockCyclicSchurInverse): fp64 rank-major import ->
//    RowsToCyclic -> in-place recursion (pure point-to-point SUMMA rings and
//    local leaves; bitwise reproducible) -> CyclicToRows -> ONE terminal
//    rounding into the fp32 apply slab.  Distributed at every N and P.
//
//  - Split-K apply through GridBLAS.gemmBatched with EXPLICIT leading dimensions
//    (arXiv:2409.03904 fig 11): the tiny-output/huge-K GEMM Y = slab^T X becomes
//    SPLITK chunk-GEMMs by pointer offset into the resident slab (lda = N),
//    partials reduced in one accelerator_for.  The source vector is assembled by
//    a cartesian ring ALLGATHER of device buffers (pure P2P; ~8x fewer bytes
//    than a padded allreduce, and no collective size cliffs).  Platform-agnostic:
//    deviceVector + GridBLAS run the SAME code on HIP/CUDA/SYCL and CPU(Eigen).
//
// VERIFY ||A Ainv x - x||/||x|| certifies the DEVICE slab + split-K path at the
// end of Import, since the single-RHS apply routes through the same core.
//
// Tensor-depth agnostic: site scalar objects treated as contiguous coarse
// scalars (ComplexF or ComplexD, following the coefficient precision)
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
  // Scalar of the coarse site objects (ComplexF or ComplexD): the apply slab
  // is always fp32 and the inversion always fp64, but the SOURCE of both --
  // the coarse operator's matrix elements -- carries this precision.
  typedef typename GridTypeMapper<CComplex>::scalar_type CoarseScalar;

  GridBase *grid;
  int      nd;
  int64_t  N;                       // dense rank = gSites * nbasis
  int      lsites;                  // my local sites
  int64_t  nrows;                   // my rows = lsites * nbasis
  std::vector<Coordinate> myLcoor;  // local coordinate of my site ss
  std::vector<int64_t>    myGsite;  // global lex site index of my site ss
  std::vector<ComplexF>   slab;     // nrows x N row-major: A during setup, rows of A^{-1} after

  static const int MRHS_MAX = 32;
  static const int SPLITK   = 32;   // requested split-K chunk count, snapped DOWN to a divisor of N

  // Apply machinery: resident slab + persistent buffers + AOT split-K pointers.
  GridBLAS BLAS;
  deviceVector<ComplexF>  dSlab;
  deviceVector<ComplexF>  dX;       // N x MRHS_MAX
  deviceVector<ComplexF>  dY;       // nrows x MRHS_MAX
  deviceVector<ComplexF>  dG;       // N x MRHS_MAX lex-major staging for the allgather
  deviceVector<int>       dLex2Rank;// lex index of a process coordinate -> its rank (allgather block order -> row-block order)
  deviceVector<int64_t>   dRm2G;    // rank-major index (rank*nrows + ss*nbasis + b) -> global column (gsite*nbasis + b) of x / the slab
  int                     myLex;
  deviceVector<ComplexF>  dPartial; // NK x (nrows x MRHS_MAX)
  deviceVector<ComplexF*> aptrs;    // slab K-chunk pointers   (lda = N)
  deviceVector<ComplexF*> xptrs;    // X    K-chunk pointers   (ldb = N)
  deviceVector<ComplexF*> cptrs;    // partial buffers         (ldc = nrows)
  std::vector<ComplexF>   hX;
  std::vector<ComplexF>   hY;
  int NK;                           // split-K chunk count (divides N)

  // Resident device memory: the apply slab and its staging.  deviceVector is
  // not evictable, so this counts against the hard budget, not the cache.
  uint64_t DeviceBytes(void)
  {
    return (uint64_t)(dSlab.capacity()+dX.capacity()+dY.capacity()+dG.capacity()+dPartial.capacity())*sizeof(ComplexF)
         + (uint64_t)dLex2Rank.capacity()*sizeof(int)
         + (uint64_t)dRm2G.capacity()*sizeof(int64_t);
  }

  DenseCoarseMatrix(GridBase *g)
    : grid(g)
  {
    GRID_ASSERT( sizeof(sobj)  == nbasis*sizeof(CoarseScalar) );
    GRID_ASSERT( sizeof(Msobj) == nbasis*nbasis*sizeof(CoarseScalar) );
    nd     = grid->_ndimension;
    N      = grid->gSites() * nbasis;
    lsites = grid->lSites();
    nrows  = (int64_t)lsites * nbasis;

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
    {
      ImportDense(Op);        // slab <- my rows of A   (LOCAL, no comms)
      ImportCertificate(Op);  // dense apply == Op.M, before inversion
      InvertDense(Op);        // slab <- my rows of A^{-1}
      double t1 = usecond();
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

      // Split-K chunk count, snapped DOWN to a divisor of N.
      NK = 1;
      for(int j=1;j<=SPLITK;j++) if ( (N % j) == 0 ) NK = j;
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

      {
        dG.resize((uint64_t)N*MRHS_MAX);
        // allgather delivers blocks in lexicographic-coordinate order; the row
        // blocks of x are in RANK order.  Same table as BuildRankMajorMap.
        int P = grid->ProcessorCount();
        std::vector<int> l2r(P);
        for(int lp=0; lp<P; lp++){ Coordinate pc(nd); Lexicographic::CoorFromIndex(pc, lp, grid->_processors); l2r[lp] = grid->RankFromProcessorCoor(pc); }
        dLex2Rank.resize(P);
        acceleratorCopyToDevice(&l2r[0], &dLex2Rank[0], P*sizeof(int));
        myLex = CartesianLexIndex(grid);
        GRID_ASSERT( l2r[myLex] == grid->ThisRank() );
        // x and the slab columns are in GLOBAL-SITE order (hX[myGsite*nbasis+b]);
        // the gathered blocks are in RANK-MAJOR order (rank*nrows + ss*nbasis + b).
        // The two coincide only on one rank, so scatter through the inverse map.
        std::vector<int64_t> g2rm; BuildRankMajorMap(g2rm);
        std::vector<int64_t> rm2g(N); for(int64_t g=0; g<N; g++) rm2g[g2rm[g]] = g;
        for(int ss=0; ss<lsites; ss++) GRID_ASSERT( rm2g[(int64_t)grid->ThisRank()*nrows + (int64_t)ss*nbasis] == myGsite[ss]*nbasis );
        dRm2G.resize(N);
        acceleratorCopyToDevice(&rm2g[0], &dRm2G[0], N*sizeof(int64_t));
      }
      std::cout << GridLogMessage << "DenseCoarseMatrix: slab resident on device ("
                << sbytes/1024./1024. << " MB/rank), split-K NK=" << NK << " (Kc=" << Kc
                << "); DEVICE cartesian ring ALLGATHER (P2P)" << std::endl;
    }

    ////////////////////////////////////////////////////////////////////
    // VERIFY: || A (Ainv x) - x || / ||x|| through the DEVICE split-K core.
    ////////////////////////////////////////////////////////////////////
    {
      Field x(grid); Field y(grid); Field z(grid);
      x = CoarseScalar(1.0,0.0);
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
      scaled = CoarseScalar(r+1.0,0.0)*in;
      InsertSliceFast(scaled,min,r,0);
    }

    Op.M(min,mout);

    ExtractSliceFast(out,mout,0,0);
    for(int r=1;r<nr;r++){
      Field sr(grid),d(grid);
      ExtractSliceFast(sr,mout,r,0);
      d = sr - CoarseScalar(r+1.0,0.0)*out;
      RealD rel = std::sqrt(norm2(d)/norm2(sr));
      // Slices differ by rounding amplified by the operator's conditioning
      // when the product cancels (A applied to A^-1 x): the tolerance follows
      // the coarse precision (fp32 measured ~5e-6 here).  A genuine rhs
      // mix-up is O(1).
      const RealD otol = (sizeof(CoarseScalar)==sizeof(ComplexF)) ? 1.0e-4 : 1.0e-6;
      if ( rel >= otol ) {
        std::cout << GridLogMessage << "DenseCoarseMatrix: oracle rhs "<<r
                  <<" inconsistent with rhs 0, rel "<<rel<<std::endl;
      }
      GRID_ASSERT( rel < otol );
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

    uint64_t nelem = (uint64_t)nrows * N;
    thread_for(i, nelem, { slab[i] = ComplexF(0.0,0.0); });

    for(int p=0; p<Op.Geometry().npoint; p++){
      Coordinate shift = Op.Geometry().shifts[p];
      // _A[p] is PADDED after ExchangeCoarseLinks (end of CoarsenOperator):
      // extract the unpadded field before peeking with unpadded coordinates
      // (exactly as MultiGeneralCoarsenedMatrix::CopyMatrix does).
      CoarseMatrix Aun(grid);  Op.ExtractMatrix(p,Aun);
      autoView(Av, Aun, CpuRead);
      thread_for(ss, lsites, {
        Coordinate ncoor(nd);
        for(int d=0; d<nd; d++){
          int64_t g = grid->_lstart[d] + myLcoor[ss][d] + shift[d];
          ncoor[d] = (int)((g % gdims[d] + gdims[d]) % gdims[d]);
        }
        int64_t nsite;
        Lexicographic::IndexFromCoor(ncoor, nsite, gdims);
        Msobj m;
        peekLocalSite(m, Av, myLcoor[ss]);
        CoarseScalar *md = (CoarseScalar *)&m;
        // The operator contracts out(s,b) = sum_a A[p](s)(a,b) in(nbr,a)
        // (GeneralCoarsenedMatrix.h Mult kernel): the stored site matrix
        // acts TRANSPOSED, so element (a,b) lands at dense row (s,b),
        // column (nbr,a).
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
        ((CoarseScalar *)&s)[b] = CoarseScalar(std::cos(ph),std::sin(0.61*ph));
      }
      pokeLocalSite(s, x, myLcoor[ss]);
    }
    // gather full x (zero-fill + exact GlobalSum), dense rows on host
    std::vector<ComplexD> xh((uint64_t)N, ComplexD(0.0,0.0));
    for(int ss=0; ss<lsites; ss++){
      sobj s;
      peekLocalSite(s, x, myLcoor[ss]);
      for(int b=0; b<nbasis; b++) xh[ myGsite[ss]*nbasis + b ] = ComplexD(((CoarseScalar *)&s)[b]);
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
      for(int b=0; b<nbasis; b++) ((CoarseScalar *)&s)[b] = CoarseScalar(yh[ss*nbasis+b]);
      pokeLocalSite(s, Dx, myLcoor[ss]);
    }
    ApplyOracle(Op, x, Ax);
    Field d(grid); d = Dx - Ax;
    RealD rel = std::sqrt(norm2(d)/norm2(Ax));
    std::cout << GridLogMessage << "DenseCoarseMatrix: IMPORT CERTIFICATE ||Dense x - A x||/||A x|| = "
              << rel << std::endl;
    if ( rel >= 1.0e-3 ) {
      std::cout << GridLogMessage << "DenseCoarseMatrix: IMPORT CERTIFICATE FAILED. If O(1), the "
                << "stencil shift-sign convention of the coarse operator has changed: "
                << "the import in ImportDense/ImportDenseForInversion must change with it"
                << std::endl;
    }
    GRID_ASSERT(rel < 1.0e-3);
  }

  ////////////////////////////////////////////////////////////////////
  // 3a. Global column -> rank-major column map, computed LOCALLY.
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
  // 3b. Direct stencil -> rank-major import of MY ROWS of A in the
  //    inversion precision (DenseInverseScalar): the stencil source is
  //    read once, straight into the buffer the Schur recursion factorises,
  //    never via the fp32 apply slab.  Same
  //    loop/sign/accumulate/transposed-contraction discipline as
  //    ImportDense; output is column-major rows x N with columns in
  //    rank-major order (g2rm).
  //    ALWAYS-ON CERTIFICATE: this import, rounded, must agree with
  //    the fp32 slab entry at the corresponding global column, over the
  //    WHOLE of my rows (few ulp: wrapped-shift collisions accumulate in
  //    different precision order).  NaN-proof: non-finite entries are
  //    counted explicitly since max() silently masks NaN.
  ////////////////////////////////////////////////////////////////////
  template<class CoarseOp>
  void ImportDenseForInversion(CoarseOp &Op, BlockRows &S, std::vector<int64_t> &g2rm)
  {
    Coordinate gdims = grid->GlobalDimensions();

    std::vector<DenseInverseScalar> h((uint64_t)nrows*N, DenseInverseScalar(0.0,0.0));
    for(int p=0; p<Op.Geometry().npoint; p++)
    {
      Coordinate shift = Op.Geometry().shifts[p];
      CoarseMatrix Aun(grid);  Op.ExtractMatrix(p,Aun);
      autoView(Av, Aun, CpuRead);
      thread_for(ss, lsites, {
        Coordinate ncoor(nd);
        for(int d=0; d<nd; d++)
        {
          int64_t g = grid->_lstart[d] + myLcoor[ss][d] + shift[d];
          ncoor[d] = (int)((g % gdims[d] + gdims[d]) % gdims[d]);
        }
        int64_t nsite;
        Lexicographic::IndexFromCoor(ncoor, nsite, gdims);
        Msobj m;
        peekLocalSite(m, Av, myLcoor[ss]);
        CoarseScalar *md = (CoarseScalar *)&m;
        // Transposed contraction as ImportDense: (a,b) lands at
        // row (s,b), column (nbr,a); column index in rank-major order.
        for(int a=0; a<nbasis; a++)
        {
          int64_t jj = g2rm[ nsite*nbasis + a ];
          for(int b=0; b<nbasis; b++)
          {
            h[(uint64_t)(ss*nbasis+b) + (uint64_t)jj*nrows] += DenseInverseScalar(md[a*nbasis+b]);
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
        ComplexD d64 = ComplexD(h[(uint64_t)(i + g2rm[gcol]*nrows)]);
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
    std::cout << GridLogMessage << "DenseCoarseMatrix: inversion-source import certificate "
              << "max|A_inv - A_slab| = " << gmx
              << "  non-finite entries " << (int64_t)gbad << std::endl;
    GRID_ASSERT( gbad == 0 );
    GRID_ASSERT( gmx < 1.0e-5 );

    S.Resize(nrows, N);
    acceleratorCopyToDevice(&h[0], &S.data[0], (uint64_t)nrows*N*sizeof(DenseInverseScalar));
  }

  ////////////////////////////////////////////////////////////////////
  // 3c. The inverse: distributed recursive Schur, end to end in the
  //    inversion precision (DenseInverseScalar, a configure-time choice):
  //    stencil -> rank-major import -> recursion -> ONE terminal rounding
  //    into the fp32 apply slab.  Everything downstream (device
  //    residency, split-K apply, VERIFY) is fp32 regardless.
  ////////////////////////////////////////////////////////////////////
  template<class CoarseOp>
  void InvertDense(CoarseOp &Op)
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
    ImportDenseForInversion(Op, S, g2rm);

    ////////////////////////////////////////////////////////////////
    // The 2D block-cyclic recursion (BlockCyclicSchurInverse).
    // Everything in it -- redistribution, SUMMA rings, leaves -- is
    // point-to-point SendToRecvFrom; no collectives at all.  Block size
    // nb = rows-per-rank makes the redistribution edges maximally
    // regular.
    ////////////////////////////////////////////////////////////////
    double t2, t3;
    {
      int Pr,Pc;
      BlockCyclicLayout::ChooseProcessGrid(P, Pr, Pc);
      int64_t nb = nrows;
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

    // The single terminal rounding: fp64 inverse -> fp32 apply slab
    // (row-major, global columns)
    {
      std::vector<DenseInverseScalar> h((uint64_t)nrows*N);
      acceleratorCopyFromDevice(&S.data[0], &h[0], (uint64_t)nrows*N*sizeof(DenseInverseScalar));
      thread_for(gcol, N, {
        int64_t jj = g2rm[gcol];
        for(int64_t i=0; i<nrows; i++)
        {
          slab[(uint64_t)i*N + gcol] = ComplexF(h[(uint64_t)(i + jj*nrows)]);
        }
      });
    }
    double t4 = usecond();
    std::cout << GridLogMessage << "DenseCoarseMatrix: SCHUR "
              << (sizeof(DenseInverseScalar)==sizeof(ComplexF) ? "fp32" : "fp64")
              << " distributed invert took "
              << (t4-t1)/1.0e6 << " s (recursion " << (t3-t2)/1.0e6 << " s)" << std::endl;
  }

  ////////////////////////////////////////////////////////////////////
  // CORE apply on packed data: hX[N x nr] (zero-filled, my sites only)
  // -> ring allgather -> split-K GEMM against the resident slab ->
  // reduce partials -> hY[nrows x nr] (column major).  Platform-agnostic:
  // deviceVector + GridBLAS (Eigen fallback on CPU builds).
  // tprof (optional): per-phase microseconds {allgather, H2D, gemm+reduce,
  // D2H}, printed by the caller on GridLogPerformance.
  ////////////////////////////////////////////////////////////////////
  void SlabApplyPacked(int nr, double *tprof)
  {
    GRID_TRACE("DenseSlabApply");
    GRID_ASSERT(nr <= MRHS_MAX);
    uint64_t nY = (uint64_t)nrows * nr;
    int64_t  Kc = N / NK;
    double t1 = usecond();
    double t2, t3;
    {
      // ALLGATHER: x is not a reduction -- every rank owns the rows of x at
      // global columns myGsite[ss]*nbasis+b (scattered by site coordinate, NOT
      // a contiguous block) and needs all of it.  Only MY rows go host->device
      // (nrows x nr, ~15 KB at nr=1), packed rank-major [r][ss*nbasis+b],
      // gathered along the process grid, then scattered through rm2g into the
      // column-major dX (ld = N, global-site columns) the split-K GEMM reads.
      const uint64_t chunk = (uint64_t)nrows*nr;              // my block: [r][ss*nbasis+b]
      { GRID_TRACE("DenseH2D");
        std::vector<ComplexF> hG(chunk);
        for(int r=0;r<nr;r++)
          for(int ss=0; ss<lsites; ss++)
            memcpy(&hG[(uint64_t)r*nrows + (uint64_t)ss*nbasis], &hX[(uint64_t)r*N + (uint64_t)myGsite[ss]*nbasis], nbasis*sizeof(ComplexF));
        acceleratorCopyToDevice(&hG[0], &dG[(uint64_t)myLex*chunk], chunk*sizeof(ComplexF));   // my slot is my LEX index
      }
      t2 = usecond();
      { GRID_TRACE("DenseAllgather");
        CartesianRingAllGather(grid, (ComplexF *)&dG[0], chunk);
        // scatter lex block L=[r][i] -> dX[r*N + rm2g[rank(L)*nrows + i]]
        ComplexF *g = &dG[0]; ComplexF *x = &dX[0]; int *l2r = &dLex2Rank[0]; int64_t *rm2g = &dRm2G[0];
        const int64_t nrw = nrows; const int64_t NN = N; const int nrr = nr;
        accelerator_for(idx, (uint64_t)N*nr, 1, {
          int64_t r = idx / NN;  int64_t gi = idx - r*NN;
          int64_t L = gi / nrw;  int64_t i  = gi - L*nrw;
          x[r*NN + rm2g[(int64_t)l2r[L]*nrw + i]] = g[L*(nrw*nrr) + r*nrw + i];
        });
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
      tprof[0] = t3-t2;                        // allgather
      tprof[1] = t2-t1;                        // H2D
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
          hX[ myGsite[ss]*nbasis + b ] = ComplexF(((CoarseScalar *)&s)[b]);
      }
    }
    SlabApplyPacked(1, nullptr);
    { GRID_TRACE("DenseUnpack");
      for(int ss=0; ss<lsites; ss++){
        sobj s;
        for(int b=0; b<nbasis; b++)
          ((CoarseScalar *)&s)[b] = CoarseScalar(hY[ss*nbasis + b]);
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
            hX[ (uint64_t)rr*N + myGsite[ss]*nbasis + b ] = ComplexF(((CoarseScalar *)&s)[b]);
        }
      }
    }
    SlabApplyPacked(nr, nullptr);
    { GRID_TRACE("DenseUnpack");
      for(int rr=0; rr<nr; rr++){
        for(int ss=0; ss<lsites; ss++){
          sobj s;
          for(int b=0; b<nbasis; b++)
            ((CoarseScalar *)&s)[b] = CoarseScalar(hY[(uint64_t)rr*nrows + (ss*nbasis+b)]);
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
            hX[(uint64_t)rr*N + myGsite[ss]*nbasis + b] = ComplexF(((CoarseScalar *)&s)[b]);
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
            ((CoarseScalar *)&s)[b] = CoarseScalar(hY[(uint64_t)rr*nrows + (ss*nbasis+b)]);  // Y col-major
          pokeLocalSite(s, ov, c6);
        }
      }
    }
    double t6 = usecond();
    std::cout << GridLogMessage << "DenseCoarseMatrix: apply6D " << nr << " rhs took "
              << (t6-t0)/1000.0 << " ms" << std::endl;
    std::cout << GridLogPerformance << "DenseCoarseMatrix: apply6D profile:"
              << " pack "        << (t1-t0)/1000.0
              << "  allgather "  << tprof[0]/1000.0
              << "  H2D "        << tprof[1]/1000.0
              << "  gemm+reduce "<< tprof[2]/1000.0
              << "  D2H "        << tprof[3]/1000.0
              << "  unpack "     << (t6-t5)/1000.0
              << "  ms" << std::endl;
  }
};

NAMESPACE_END(Grid);
