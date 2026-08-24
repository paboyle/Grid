/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: RecursiveSchurInverse.h

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

#include <algorithm>

NAMESPACE_BEGIN(Grid);

///////////////////////////////////////////////////////////////////////////////
// RecursiveSchurInverse: distributed dense inversion by recursive Schur
// complement over a binary rank-range tree.
//
// Contract: rank r owns global rows [rowStart[r], rowStart[r+1]) of an
// N x N matrix in rank-major ordering, and receives its rows of the
// inverse in the same layout.  Consumes GridBase collectives, GridBLAS
// and GridBLASInverse only; Eigen reference backends permit CPU unit
// testing under mpirun (Test_schur_inverse).
//
// Arithmetic is fp64 throughout; the caller rounds once into fp32 storage.
//
// Execution: SPMD full-tree walk.  Every rank makes the identical call
// sequence; data participation is ownership-gated; all collectives are
// world-wide, so no deadlock surface exists.
//
// Storage: BlockRows is column-major, ld = rows; element (i,j) at
// data[i + j*ld]; a column window is the contiguous slice at data[col0*ld].
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// My rows of a distributed dense matrix: rows x cols, column major, ld = rows.
///////////////////////////////////////////////////////////////////////////////
class BlockRows
{
public:
  deviceVector<ComplexD> data;
  int64_t                rows;
  int64_t                cols;
  int64_t                ld;

  BlockRows()
  {
    rows = 0;
    cols = 0;
    ld   = 0;
  }
  void Resize(int64_t r, int64_t c)
  {
    rows = r;
    cols = c;
    ld   = r;
    data.resize((uint64_t)r*c);
  }
  ComplexD *ColumnWindow(int64_t col0)
  {
    GRID_ASSERT( col0 >= 0 );
    GRID_ASSERT( col0 <= cols );
    return &data[(uint64_t)col0*ld];
  }
};

class RecursiveSchurInverse
{
public:
  GridBase             *grid;       // world collectives only
  int64_t               N;          // global matrix dimension
  int                   P;          // ranks
  int                   me;         // this rank
  std::vector<int64_t>  rowStart;   // P+1 entries: rank-major row ownership
  int64_t               myRow0;
  int64_t               myNrows;
  int64_t               panelBytes; // gather panel budget (DENSE_PANEL_BYTES)

  GridBLAS              BLAS;
  GridBLASInverse       INV;

  // Growth telemetry: one entry per merge node, walk order
  std::vector<double>   telNormB;      // ||B||_F = ||A11inv A12||_F
  std::vector<double>   telSratio;     // ||S||_F / ||A22||_F
  double                telLeafMaxInv; // max |(leaf inverse)_ij| over leaves

  // Phase timers/counters, reported by ReportTelemetry
  double                tMemset;       // device panel zero-fill
  double                tDeposit;      // owner rows -> panel (device kernel)
  double                tAllreduce;    // panel collective (see tBarrier)
  double                tBarrier;      // arrival skew, when DENSE_BARRIER_PROBE
  double                tRepack;       // rank-major -> panel (gather path)
  double                tGemm;         // strided gemm + synchronise
  double                tLeaf;         // leaf inversions
  double                tARmin;        // fastest single panel collective
  double                tARmax;        // slowest single panel collective
  std::vector<double>   tARall;        // every panel collective, for percentiles
  uint64_t              bytesAllreduce;
  uint64_t              nAllreduce;    // panel collectives
  uint64_t              nGatherGemm;   // GatherGemm calls
  uint64_t              nGather;       // gather-path collectives (debug gate)
  uint64_t              nGatherV;      // ... of which MPI_Allgatherv
  uint64_t              nBcast;        // ... of which Bcast-assembled

  // Persistent grow-only device panel; assembly and collectives are
  // device-resident.  Device builds require GPU-aware MPI.
  deviceVector<ComplexD> dPanelBuf;
  // Rank-major receive staging for the AllGatherV path, plus the per-panel-row
  // owner maps that drive the repack.  Grow-only, same discipline as dPanelBuf:
  // a fresh device allocation per collective is catastrophic (measured).
  deviceVector<ComplexD> dRecvBuf;
  deviceVector<int64_t>  dOwnerOff;    // panel row -> owner's first panel row
  deviceVector<int64_t>  dOwnerRows;   // panel row -> owner's row count

  // DENSE_GATHER            : transport for the panel assembly.
  //     0 = zero-fill + GlobalSumVector  (default; moves the payload twice)
  //     1 = MPI_Allgatherv               (KNOWN BROKEN here, see guard below)
  //     2 = C sequential MPI_Bcast, one per contributing rank.  Bcast has no
  //         count vector, so the zero-count shape that breaks Allgatherv
  //         cannot arise.  Costs C collectives instead of 1 and the roots do
  //         not transmit concurrently, so the byte cost is ~2*panel per rank
  //         -- the same as the allreduce.  It wins only if Bcast is faster
  //         PER BYTE than Allreduce, which is why it must be measured.
  // DENSE_GATHER_MIN_BYTES  : panels below this stay on the allreduce path.
  //                           Default 0 (always gather).  The Schur tree puts
  //                           94% of its bytes in panels with ~3.7 MB per-rank
  //                           messages, i.e. squarely bandwidth bound; only the
  //                           two deepest levels (1.2% of bytes, 16-65 KB per
  //                           rank) are latency bound.  This exists so that
  //                           tail can be excluded with an env var rather than
  //                           a rebuild, if it ever proves to matter.
  // DENSE_BARRIER_PROBE     : Barrier() before each panel collective.
  //   DEFAULTS ON, and it is an optimisation rather than instrumentation.
  //   Measured at 288 ranks, N=138240: comms 380.3 s -> 314.7 s and the whole
  //   invert 394.3 s -> 327.7 s, with the worst single collective falling
  //   2474.9 ms -> 256.2 ms and effective rate 2.40 -> 4.27 GB/s.  Convoy
  //   accounting alone predicts NO saving (barrier and collective both
  //   complete at max-over-ranks), so the mechanism is that the collective
  //   itself degrades under skewed arrival, not that skew is rebooked.
  //   Set DENSE_BARRIER_PROBE=0 to recover the old behaviour.  It also
  //   separates arrival skew (tBarrier) from transfer (tAllreduce).
  int                    useGather;
  int64_t                gatherMinBytes;
  int                    barrierProbe;
  // DENSE_GATHER_DEBUG=N : dump the descriptor of the first N AllGatherV
  // calls from rank 0.  The tiling invariant below is checked ALWAYS -- it
  // costs one O(P) host loop against a collective, and a counts array that
  // does not tile the panel is exactly the failure we are hunting.
  int                    gatherDebug;

  ///////////////////////////////////////////////////////////////////////////
  // Ownership-table validation: a proper partition of [0,N).
  // Static and communicator-free so synthetic tables unit-test directly.
  ///////////////////////////////////////////////////////////////////////////
  static void CheckRowStart(const std::vector<int64_t> &table, int64_t N)
  {
    int P = (int)table.size() - 1;
    GRID_ASSERT( P >= 1 );
    GRID_ASSERT( table[0] == 0 );
    GRID_ASSERT( table[P] == N );
    for(int r=0; r<P; r++)
    {
      GRID_ASSERT( table[r+1] >= table[r] );  // zero-row ranks permitted
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Build the ownership table from each rank's local row count: zero-fill
  // allgather (the standing comms idiom) then prefix sum.  Every rank
  // returns the identical table.
  ///////////////////////////////////////////////////////////////////////////
  static std::vector<int64_t> MakeRowStart(GridBase *g, int64_t myNrows)
  {
    int P  = g->ProcessorCount();
    int me = g->ThisRank();

    std::vector<uint64_t> counts(P, 0);
    counts[me] = (uint64_t)myNrows;
    g->GlobalSumVector(&counts[0], P);

    std::vector<int64_t> table(P+1);
    table[0] = 0;
    for(int r=0; r<P; r++)
    {
      table[r+1] = table[r] + (int64_t)counts[r];
    }
    CheckRowStart(table, table[P]);
    return table;
  }

  RecursiveSchurInverse(GridBase *g,
                        int64_t N_,
                        std::vector<int64_t> &rowStart_,
                        int64_t panelBytes_)
  {
    grid       = g;
    N          = N_;
    P          = g->ProcessorCount();
    me         = g->ThisRank();
    rowStart   = rowStart_;
    panelBytes = panelBytes_;

    GRID_ASSERT( (int)rowStart.size() == P+1 );
    CheckRowStart(rowStart, N);

    myRow0  = rowStart[me];
    myNrows = rowStart[me+1] - rowStart[me];

    telLeafMaxInv = 0.0;

    useGather      = getenv("DENSE_GATHER") ? atoi(getenv("DENSE_GATHER")) : 0;
    gatherMinBytes = getenv("DENSE_GATHER_MIN_BYTES")
                   ? atol(getenv("DENSE_GATHER_MIN_BYTES")) : 0;
    barrierProbe   = getenv("DENSE_BARRIER_PROBE") ? atoi(getenv("DENSE_BARRIER_PROBE")) : 1;
    gatherDebug    = getenv("DENSE_GATHER_DEBUG") ? atoi(getenv("DENSE_GATHER_DEBUG")) : 0;

    ///////////////////////////////////////////////////////////////////////
    // DENSE_GATHER=1 is KNOWN BROKEN on Cray MPICH and refuses to start.
    //
    // The gather asks MPI_Allgatherv to assemble a panel to which only the
    // rank sub-range [rB0,rB1) contributes; every other rank passes count 0.
    // At Schur depth d that is P/2^(d+1) contributors out of P, i.e. 9 of 288
    // at depth 4 and ~1 of 288 at depth 7.  Measured consequences at 288
    // ranks (tests/debug/Test_allgather):
    //   T5  288/288 contributing, 18.4 MB on device : sub-second.
    //   T6  144/288 contributing,  9.2 MB on device : ~53 s.
    // and in production, 9/288 contributing at 299 MB hangs and then aborts
    // inside PMPI_Allgatherv ("req != NULL", mpir_request.h:508).
    //
    // The primitive itself is sound -- Test_allgather T1-T6 pass, and T3/T4
    // verify it BITWISE against the zero-fill+GlobalSum idiom.  What fails is
    // MPI's handling of a collective in which almost every rank contributes
    // nothing.  Fixing it needs either P2P (a bisection allgather over
    // [r0,r1)) or the 2D block-cyclic layout, which removes the shape
    // altogether by giving every rank part of every sub-block.
    //
    // Fail here, in the first second, rather than 100 s and 36 nodes into a
    // job.  DENSE_GATHER_FORCE=1 proceeds anyway for debugging.
    ///////////////////////////////////////////////////////////////////////
    if ( (useGather==1) && !getenv("DENSE_GATHER_FORCE") )
    {
      std::cout << GridLogError
                << "DENSE_GATHER=1 is disabled: MPI_Allgatherv on this MPI is"
                << " pathological when most ranks contribute count 0 (see the"
                << " note at RecursiveSchurInverse.h, and Test_allgather T5 vs"
                << " T6).  Unset DENSE_GATHER, or set DENSE_GATHER_FORCE=1 to"
                << " proceed anyway." << std::endl;
      GRID_ASSERT(0 && "DENSE_GATHER=1 known broken: see comment above");
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // The communication primitive.
  //
  //   C(:, colC : colC+widthB) <-   beta  * C(:, colC : colC+widthB)
  //                               + alpha * A(:, colA : colA+widthA) * Bsub
  //
  // Bsub is the widthA x widthB sub-block of a row-distributed operand
  // owned by ranks [rB0, rB1): owner r contributes its rows of
  // B(:, colB : colB+widthB) at sub-block row offset
  // rowStart[r] - rowStart[rB0], gathered in panelBytes row-chunks by
  // device zero-fill + deposit kernel + GlobalSumVector.
  //
  // Every rank calls; non-owners of B add zeros; ranks with A.rows == 0
  // skip local compute but make every collective call.  Column offsets
  // are local buffer offsets -- non-participants pass 0.
  ///////////////////////////////////////////////////////////////////////////
  void GatherGemm(ComplexD alpha,
                  BlockRows &A, int64_t colA, int64_t widthA,
                  int rB0, int rB1,
                  BlockRows &B, int64_t colB, int64_t widthB,
                  ComplexD beta,
                  BlockRows &C, int64_t colC)
  {
    GRID_ASSERT( rB0 >= 0 );
    GRID_ASSERT( rB1 >  rB0 );
    GRID_ASSERT( rB1 <= P );

    int64_t k = rowStart[rB1] - rowStart[rB0];
    int64_t m = A.rows;
    int64_t n = widthB;
    GRID_ASSERT( widthA == k );
    GRID_ASSERT( n >= 1 );

    int owner = ( me >= rB0 ) && ( me < rB1 ) && ( B.rows > 0 );
    int64_t myOff = 0;
    if ( owner )
    {
      myOff = rowStart[me] - rowStart[rB0];
    }

    if ( m > 0 )
    {
      GRID_ASSERT( colA + widthA <= A.cols );
      GRID_ASSERT( colC + widthB <= C.cols );
      GRID_ASSERT( C.rows == m );
    }

    nGatherGemm++;
    GRID_TRACE("GatherGemm");

    if ( owner )
    {
      GRID_ASSERT( colB + widthB <= B.cols );
    }

    // COLUMN chunking, not row chunking.  Each rank's contribution to a
    // column chunk is rows_r x nchunk at ld = B.rows -- i.e. exactly the
    // contiguous window B.ColumnWindow(colB+j0) -- so the AllGatherV send
    // needs no pack.  Row chunking would slice each column and force one.
    // It also makes every chunk a full-k product, so beta applies directly
    // and the cross-chunk accumulation disappears.
    int64_t bufs = useGather ? 2 : 1;         // panel, plus receive staging
    int64_t nc = panelBytes / ( bufs * (int64_t)sizeof(ComplexD) * k );
    if ( nc < 1 ) nc = 1;
    if ( nc > n ) nc = n;
    GRID_ASSERT( k*nc < 2147483647L );        // collective counts are int

    deviceVector<ComplexD> &dPanel = dPanelBuf;
    if ( dPanel.size() < (uint64_t)k*nc ) dPanel.resize((uint64_t)k*nc);
    deviceVector<ComplexD*> ap(1);
    deviceVector<ComplexD*> bp(1);
    deviceVector<ComplexD*> cp(1);
    std::vector<ComplexD*>  ptr(1);

    // Panel-row -> owner maps for the repack.  Depend only on [rB0,rB1),
    // so they are built once per call rather than once per chunk.
    if ( useGather )
    {
      if ( dRecvBuf.size()   < (uint64_t)k*nc ) dRecvBuf.resize((uint64_t)k*nc);
      if ( dOwnerOff.size()  < (uint64_t)k )    dOwnerOff.resize((uint64_t)k);
      if ( dOwnerRows.size() < (uint64_t)k )    dOwnerRows.resize((uint64_t)k);
      std::vector<int64_t> hoff(k), hrows(k);
      for(int r=rB0; r<rB1; r++)
      {
        int64_t off_r  = rowStart[r]   - rowStart[rB0];
        int64_t rows_r = rowStart[r+1] - rowStart[r];
        for(int64_t i=0;i<rows_r;i++){ hoff[off_r+i]=off_r; hrows[off_r+i]=rows_r; }
      }
      acceleratorCopyToDevice(&hoff[0], &dOwnerOff[0], (uint64_t)k*sizeof(int64_t));
      acceleratorCopyToDevice(&hrows[0],&dOwnerRows[0],(uint64_t)k*sizeof(int64_t));
      if ( owner ) GRID_ASSERT( B.rows == rowStart[me+1]-rowStart[me] );
    }

    for(int64_t j0=0; j0<n; j0+=nc)
    {
      int64_t nchunk = std::min(nc, n-j0);
      uint64_t panelWords = (uint64_t)k*nchunk;
      uint64_t panelBytesThis = panelWords*sizeof(ComplexD);

      // The MODE (0/1/2), not a boolean: `useGather && cond` collapses 2 to 1
      // and silently dispatched DENSE_GATHER=2 onto the known-broken
      // AllGatherV path (Frontier hang, 2026-08-24).  Identical on all ranks.
      int gatherThis = ( (int64_t)panelBytesThis >= gatherMinBytes ) ? useGather : 0;

      // Arrival skew is charged to the barrier, so that tAllreduce measures
      // transfer alone.  Diagnostic only; off by default.
      if ( barrierProbe ) { double tb = -usecond(); grid->Barrier(); tBarrier += tb+usecond(); }

      double tar = -usecond();
      if ( gatherThis )
      {
        //////////////////////////////////////////////////////////////////
        // AllGatherV: move the payload once, no arithmetic, no zero fill.
        // Receive is rank major -- a concatenation of rows_r x nchunk
        // column-major blocks -- which is the correct ROW order but the
        // wrong LAYOUT, so one device repack follows.
        //////////////////////////////////////////////////////////////////
        std::vector<int> counts(P,0), displs(P,0);
        for(int r=rB0; r<rB1; r++)
        {
          counts[r] = (int)((rowStart[r+1]-rowStart[r])*nchunk);
          displs[r] = (int)((rowStart[r]  -rowStart[rB0])*nchunk);
        }
        //////////////////////////////////////////////////////////////////
        // The counts MUST tile the panel exactly: every one of the
        // k*nchunk words has exactly one owner, or MPI is handed an
        // inconsistent descriptor.  Always checked; O(P) against a
        // collective is free.
        //////////////////////////////////////////////////////////////////
        int64_t csum = 0; int nz = 0;
        for(int r=0;r<P;r++){ csum += counts[r]; if ( counts[r] ) nz++; }
        GRID_ASSERT( csum == (int64_t)panelWords );
        GRID_ASSERT( displs[rB1-1] + counts[rB1-1] == (int)panelWords );
        // Gate on the count of GATHERS, not of all collectives.  nAllreduce
        // is already in the hundreds by the first gather -- the recursion is
        // depth first and everything below the size threshold takes the
        // allreduce path -- so gating on it prints nothing at any sane value.
        nGather++;
        if ( gatherDebug && ((int)nGather <= gatherDebug) && (me==0) ) {
          std::cout << GridLogMessage << "GATHER["<<nGather-1<<"]"
                    << "  ranks ["<<rB0<<","<<rB1<<")"
                    << "  k "<<k<<"  n "<<n<<"  nchunk "<<nchunk
                    << "  panelWords "<<panelWords
                    << " ("<<panelBytesThis/1024/1024<<" MB)"
                    << "  contributors "<<nz<<"/"<<P
                    << "  maxcount "<<*std::max_element(counts.begin(),counts.end())
                    << "  rank0: count "<<counts[me]<<(owner?" owner":" non-owner")
                    << std::endl;
        }
        if ( gatherThis == 1 )
        {
          nGatherV++;
          void *send = owner ? (void *)B.ColumnWindow(colB+j0) : (void *)&dRecvBuf[0];
          grid->AllGatherV(send, counts[me],
                           (void *)&dRecvBuf[0], counts, displs, sizeof(ComplexD));
        }
        else
        {
          ////////////////////////////////////////////////////////////////
          // C sequential broadcasts into the SAME rank-major staging
          // buffer, so the repack below is shared with the AllGatherV arm.
          // The root must already hold its own block: B.ColumnWindow is
          // contiguous rows_r x nchunk, and so is its slot in dRecvBuf.
          // Every rank issues all C broadcasts in the same order, so
          // collective matching stays positional and safe.
          ////////////////////////////////////////////////////////////////
          nBcast++;
          if ( owner )
          {
            int64_t off_me = rowStart[me] - rowStart[rB0];
            acceleratorCopyDeviceToDevice((void *)B.ColumnWindow(colB+j0),
                                          (void *)&dRecvBuf[off_me*nchunk],
                                          (uint64_t)B.rows*nchunk*sizeof(ComplexD));
          }
          for(int r=rB0;r<rB1;r++)
          {
            int64_t off_r  = rowStart[r]   - rowStart[rB0];
            int64_t rows_r = rowStart[r+1] - rowStart[r];
            if ( rows_r )
              grid->Broadcast(r,(void *)&dRecvBuf[off_r*nchunk],
                              (uint64_t)rows_r*nchunk*sizeof(ComplexD));
          }
        }
        tar += usecond();

        tRepack -= usecond();
        {
          ComplexD *pan = &dPanel[0];
          ComplexD *rcv = &dRecvBuf[0];
          int64_t  *ooff  = &dOwnerOff[0];
          int64_t  *orows = &dOwnerRows[0];
          int64_t   kk = k, ncw = nchunk;
          accelerator_for(idx, panelWords, 1, {
            int64_t j = idx / kk;
            int64_t p = idx - j*kk;
            int64_t o = ooff[p];
            int64_t rr= orows[p];
            pan[idx] = rcv[(uint64_t)(o*ncw + (p-o) + j*rr)];
          });
        }
        tRepack += usecond();
      }
      else
      {
        //////////////////////////////////////////////////////////////////
        // Zero fill + deposit + GlobalSumVector.  Exact because the zero
        // fill leaves exactly one contributing rank per element, but it
        // moves the payload twice and reduces over zeros.
        //////////////////////////////////////////////////////////////////
        tar += usecond();
        tMemset -= usecond();
        acceleratorMemSet(&dPanel[0], 0, panelBytesThis);
        tMemset += usecond();
        if ( owner )
        {
          int64_t   brows = B.rows;
          int64_t   kk    = k;
          int64_t   dof   = myOff;
          ComplexD *src   = B.ColumnWindow(colB+j0);
          ComplexD *dst   = &dPanel[0];
          tDeposit -= usecond();
          accelerator_for(idx, (uint64_t)(brows*nchunk), 1, {
            int64_t j = idx / brows;
            int64_t i = idx - j*brows;
            dst[(uint64_t)(dof + i + j*kk)] = src[(uint64_t)(i + j*brows)];
          });
          tDeposit += usecond();
        }
        tar = -usecond();
        grid->GlobalSumVector(&dPanel[0], (int)panelWords);
        tar += usecond();
      }
      tAllreduce += tar;
      tARmin = std::min(tARmin, tar);
      tARmax = std::max(tARmax, tar);
      tARall.push_back(tar);
      bytesAllreduce += panelBytesThis;
      nAllreduce++;

      if ( m > 0 )
      {
        // Full-k product into this column chunk of C: beta applies directly.
        ptr[0] = A.ColumnWindow(colA);
        acceleratorCopyToDevice(&ptr[0], &ap[0], sizeof(ComplexD*));
        ptr[0] = &dPanel[0];
        acceleratorCopyToDevice(&ptr[0], &bp[0], sizeof(ComplexD*));
        ptr[0] = C.ColumnWindow(colC+j0);
        acceleratorCopyToDevice(&ptr[0], &cp[0], sizeof(ComplexD*));

        tGemm -= usecond();
        BLAS.gemmBatched(GridBLAS_OP_N, GridBLAS_OP_N,
                         (int)m, (int)nchunk, (int)k,
                         alpha,    ap, (int)A.ld,
                                   bp, (int)k,
                         beta,     cp, (int)C.ld);
        BLAS.synchronise();
        tGemm += usecond();
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Local Frobenius norm-squared of a full-height column window.
  // NO comms; callers GlobalSum the result.  Host staging, setup-scale.
  ///////////////////////////////////////////////////////////////////////////
  double FrobNorm2Local(BlockRows &X, int64_t col0, int64_t w)
  {
    if ( X.rows == 0 ) return 0.0;
    GRID_ASSERT( col0 + w <= X.cols );
    uint64_t len = (uint64_t)X.rows*w;
    std::vector<ComplexD> h(len);
    acceleratorCopyFromDevice(X.ColumnWindow(col0), &h[0], len*sizeof(ComplexD));
    // Member real()/imag(): portable across std::complex (CPU) and
    // thrust::complex (HIP), where std::norm does not resolve.
    double s = 0.0;
    for(uint64_t i=0; i<len; i++)
    {
      double re = h[i].real();
      double im = h[i].imag();
      s += re*re + im*im;
    }
    return s;
  }

  ///////////////////////////////////////////////////////////////////////////
  // dst(:, dcol0 : dcol0+w) = - src(:, 0:w).  Both operands have ld == rows
  // so full-height windows are contiguous: flat elementwise device copy.
  ///////////////////////////////////////////////////////////////////////////
  void NegateCopy(BlockRows &src, BlockRows &dst, int64_t dcol0, int64_t w)
  {
    GRID_ASSERT( src.rows == dst.rows );
    GRID_ASSERT( w <= src.cols );
    GRID_ASSERT( dcol0 + w <= dst.cols );
    if ( src.rows == 0 ) return;
    uint64_t len = (uint64_t)src.rows*w;
    ComplexD *s = &src.data[0];
    ComplexD *d = dst.ColumnWindow(dcol0);
    accelerator_for(i, len, 1, {
      d[i] = -s[i];
    });
  }

  ///////////////////////////////////////////////////////////////////////////
  // Leaf inversion: local, in place on the contiguous diagonal window.
  // No collectives.
  ///////////////////////////////////////////////////////////////////////////
  void LeafInvert(int64_t col0, int64_t width, BlockRows &Arows)
  {
    GRID_TRACE("SchurLeaf");
    GRID_ASSERT( width == Arows.rows );
    GRID_ASSERT( col0 + width <= Arows.cols );
    int64_t  w   = width;
    uint64_t len = (uint64_t)w*w;
    tLeaf -= usecond();

    deviceVector<ComplexD*> bp(1);
    std::vector<ComplexD*>  ptr(1);
    ptr[0] = Arows.ColumnWindow(col0);
    acceleratorCopyToDevice(&ptr[0], &bp[0], sizeof(ComplexD*));
    INV.inverseBatched(w, bp);

    // Telemetry: max |element| of the leaf inverse
    {
      std::vector<ComplexD> h(len);
      acceleratorCopyFromDevice(Arows.ColumnWindow(col0), &h[0], len*sizeof(ComplexD));
      double mx = 0.0;
      for(uint64_t i=0; i<len; i++)
      {
        double re = h[i].real();
        double im = h[i].imag();
        mx = std::max(mx, re*re + im*im);
      }
      telLeafMaxInv = std::max(telLeafMaxInv, std::sqrt(mx));
    }
    tLeaf += usecond();
  }

  ///////////////////////////////////////////////////////////////////////////
  // The recursion (plan 3.5 / 4B.3).  Inverts the diagonal block of the
  // rank-major matrix spanned by ranks [r0, r1), living in every member
  // rank's column window [col0, col0+width) -- IN PLACE.
  //
  // SPMD: every rank calls with IDENTICAL (r0, r1, width) and its own
  // local (col0, Arows); ranks outside [r0, r1) participate in the
  // collectives only (dummy operands, zero contributions).  The collective
  // sequence -- 5 GatherGemm calls + 3 scalar GlobalSums per merge node --
  // is identical on every rank by construction.
  //
  //   I = [r0, mid)   J = [mid, r1)   widths WI, WJ
  //   1. recurse I:      A11 -> A11inv
  //   2. B = A11inv.A12                     (I rows)
  //   3. C = A21.A11inv                     (J rows)
  //   4. S = A22 - A21.B    in place        (J rows)   [alpha=-1, beta=1]
  //   5. recurse J:      S -> Sinv
  //   6. T = Sinv.C                         (J rows)
  //   7. U = B.Sinv                         (I rows)
  //   8. X11 = A11inv + U.C  in place       (I rows)   [beta=1]
  //   9. X12 = -U, X21 = -T  local negates; X22 = Sinv already in place
  ///////////////////////////////////////////////////////////////////////////
  void SchurNode(int r0, int r1, int64_t col0, int64_t width, BlockRows &Arows)
  {
    int span = r1 - r0;
    GRID_ASSERT( span >= 1 );
    GRID_ASSERT( width == rowStart[r1] - rowStart[r0] );

    if ( span == 1 )
    {
      if ( ( me == r0 ) && ( myNrows > 0 ) )
      {
        LeafInvert(col0, width, Arows);
      }
      return;
    }

    int     mid = ( r0 + r1 ) / 2;
    int64_t WI  = rowStart[mid] - rowStart[r0];
    int64_t WJ  = rowStart[r1]  - rowStart[mid];

    // Zero-width child ranges (all ranks of a half owning no rows) are a
    // KNOWN LIMITATION: fail loudly rather than divide mysteriously.
    GRID_ASSERT( WI > 0 );
    GRID_ASSERT( WJ > 0 );

    int inI = ( me >= r0  ) && ( me < mid );
    int inJ = ( me >= mid ) && ( me < r1  );

    ComplexD one  ( 1.0,0.0);
    ComplexD mone (-1.0,0.0);
    ComplexD zero ( 0.0,0.0);

    BlockRows dummy;

    // 1. A11 -> A11inv
    SchurNode(r0, mid, col0, WI, Arows);

    // 2. B = A11inv . A12   (I rows; gather A12 from I owners)
    BlockRows Bbuf;
    if ( inI ) Bbuf.Resize(myNrows, WJ);
    {
      BlockRows &Aop = inI ? Arows : dummy;
      BlockRows &Cop = inI ? Bbuf  : dummy;
      int64_t    cA  = inI ? col0  : 0;
      GatherGemm(one, Aop, cA, WI,
                 r0, mid,
                 Arows, col0+WI, WJ,
                 zero, Cop, 0);
    }
    double nB = FrobNorm2Local(Bbuf, 0, inI ? WJ : 0);
    grid->GlobalSumVector(&nB, 1);
    telNormB.push_back(std::sqrt(nB));

    // 3. C = A21 . A11inv   (J rows; gather A11inv from I owners)
    BlockRows Cbuf;
    if ( inJ ) Cbuf.Resize(myNrows, WI);
    {
      BlockRows &Aop = inJ ? Arows : dummy;
      BlockRows &Cop = inJ ? Cbuf  : dummy;
      int64_t    cA  = inJ ? col0  : 0;
      GatherGemm(one, Aop, cA, WI,
                 r0, mid,
                 Arows, col0, WI,
                 zero, Cop, 0);
    }

    // 4. S = A22 - A21 . B  in place on my A22 window (J rows)
    double nA22 = FrobNorm2Local( inJ ? Arows : dummy, inJ ? col0+WI : 0, inJ ? WJ : 0 );
    grid->GlobalSumVector(&nA22, 1);
    {
      BlockRows &Aop = inJ ? Arows   : dummy;
      BlockRows &Cop = inJ ? Arows   : dummy;
      int64_t    cA  = inJ ? col0    : 0;
      int64_t    cC  = inJ ? col0+WI : 0;
      GatherGemm(mone, Aop, cA, WI,
                 r0, mid,
                 Bbuf, 0, WJ,
                 one, Cop, cC);
    }
    double nS = FrobNorm2Local( inJ ? Arows : dummy, inJ ? col0+WI : 0, inJ ? WJ : 0 );
    grid->GlobalSumVector(&nS, 1);
    telSratio.push_back( std::sqrt(nS) / ( std::sqrt(nA22) + 1.0e-300 ) );

    // 5. S -> Sinv
    SchurNode(mid, r1, col0+WI, WJ, Arows);

    // 6. T = Sinv . C       (J rows; gather C from J owners)
    BlockRows Tbuf;
    if ( inJ ) Tbuf.Resize(myNrows, WI);
    {
      BlockRows &Aop = inJ ? Arows   : dummy;
      BlockRows &Cop = inJ ? Tbuf    : dummy;
      int64_t    cA  = inJ ? col0+WI : 0;
      GatherGemm(one, Aop, cA, WJ,
                 mid, r1,
                 Cbuf, 0, WI,
                 zero, Cop, 0);
    }

    // 7. U = B . Sinv       (I rows; gather Sinv from J owners)
    BlockRows Ubuf;
    if ( inI ) Ubuf.Resize(myNrows, WJ);
    {
      BlockRows &Aop = inI ? Bbuf : dummy;
      BlockRows &Cop = inI ? Ubuf : dummy;
      GatherGemm(one, Aop, 0, WJ,
                 mid, r1,
                 Arows, col0+WI, WJ,
                 zero, Cop, 0);
    }

    // 8. X11 = A11inv + U . C  in place (I rows; gather C from J owners)
    {
      BlockRows &Aop = inI ? Ubuf  : dummy;
      BlockRows &Cop = inI ? Arows : dummy;
      int64_t    cC  = inI ? col0  : 0;
      GatherGemm(one, Aop, 0, WJ,
                 mid, r1,
                 Cbuf, 0, WI,
                 one, Cop, cC);
    }

    // 9. Off-diagonal signs, local
    if ( inI ) NegateCopy(Ubuf, Arows, col0+WI, WJ);
    if ( inJ ) NegateCopy(Tbuf, Arows, col0,    WI);
  }

  ///////////////////////////////////////////////////////////////////////////
  // PUBLIC ENTRY.  Arows: my rows of the rank-major N x N matrix (fp64).
  // On exit Arows holds my rows of the inverse, still fp64; the caller
  // owns the single terminal rounding into its fp32 apply storage.
  ///////////////////////////////////////////////////////////////////////////
  void Invert(BlockRows &Arows)
  {
    GRID_ASSERT( Arows.rows == myNrows );
    GRID_ASSERT( Arows.cols == N );

    telNormB.resize(0);
    telSratio.resize(0);
    telLeafMaxInv = 0.0;

    tMemset        = 0.0;
    tDeposit       = 0.0;
    tAllreduce     = 0.0;
    tBarrier       = 0.0;
    tRepack        = 0.0;
    tGemm          = 0.0;
    tLeaf          = 0.0;
    tARmin         = 1.0e30;
    tARmax         = 0.0;
    tARall.clear();
    bytesAllreduce = 0;
    nAllreduce     = 0;
    nGatherGemm    = 0;
    nGather        = 0;
    nGatherV       = 0;
    nBcast         = 0;

    SchurNode(0, P, 0, N, Arows);

    RealD mx = telLeafMaxInv;
    grid->GlobalMax(mx);
    telLeafMaxInv = mx;
  }

  // All telemetry values are globally reduced or boss-local; safe to
  // stream on every rank (Grid quiesces stdout to the boss unless
  // --debug-stdout).  NOTE: tAllreduce INCLUDES wait/imbalance -- a rank
  // arriving early books its wait here; the min/max spread across ranks
  // separates true wire time (~min) from skew (max-min).
  void ReportTelemetry(void)
  {
    for(uint64_t i=0; i<telNormB.size(); i++)
    {
      std::cout << GridLogPerformance
                << "SchurNode " << i
                << " ||B||_F "       << telNormB[i]
                << " ||S||/||A22|| " << telSratio[i]
                << std::endl;
    }
    std::cout << GridLogPerformance
              << "Schur leaves max|Ainv| " << telLeafMaxInv
              << std::endl;

    RealD armax =  tAllreduce;
    RealD armin = -tAllreduce;
    grid->GlobalMax(armax);
    grid->GlobalMax(armin);
    armin = -armin;

    std::cout << GridLogMessage << "Schur phases (boss rank, seconds):"
              << "  memset "    << tMemset/1.0e6
              << "  deposit "   << tDeposit/1.0e6
              << "  allreduce " << tAllreduce/1.0e6
              << "  barrier "   << tBarrier/1.0e6
              << "  repack "    << tRepack/1.0e6
              << "  gemm "      << tGemm/1.0e6
              << "  leaf "      << tLeaf/1.0e6
              << std::endl;
    std::cout << GridLogMessage << "Schur comms:"
              << "  [transports run: allreduce " << (nAllreduce - nGatherV - nBcast)
              << " allgatherv " << nGatherV << " bcast " << nBcast << "]"
              << ( barrierProbe ? " [barrier probe on]" : "" )
              << "  GatherGemm calls " << nGatherGemm
              << "  panel allreduces " << nAllreduce
              << "  allreduce GB "     << bytesAllreduce/1024./1024./1024.
              << "  allreduce s min/max over ranks " << armin/1.0e6
              << " / " << armax/1.0e6
              << std::endl;
    std::cout << GridLogMessage << "Schur comms per-call (boss):"
              << "  min " << (nAllreduce ? tARmin/1.0e3 : 0.0) << " ms"
              << "  avg " << (nAllreduce ? tAllreduce/nAllreduce/1.0e3 : 0.0) << " ms"
              << "  max " << tARmax/1.0e3 << " ms"
              << "   effective " << (tAllreduce>0 ? bytesAllreduce/tAllreduce*1.0e6/1.0e9 : 0.0)
              << " GB/s"
              << std::endl;
    //////////////////////////////////////////////////////////////////////
    // Distribution, not just min/avg/max.  A barrier that merely REBOOKS
    // skew shifts the median; a barrier that suppresses pathological
    // collectives shortens the tail.  p50 vs p99 separates the two from a
    // single run, which min/avg/max cannot.
    //////////////////////////////////////////////////////////////////////
    if ( tARall.size() ) {
      std::vector<double> v = tARall;
      std::sort(v.begin(),v.end());
      size_t n = v.size();
      auto pc = [&](double f){ size_t i=(size_t)(f*(n-1)); return v[i]/1.0e3; };
      double med = pc(0.50);
      // how much time is spent in calls that are gross outliers
      double tail=0.0; size_t ntail=0;
      for(size_t i=0;i<n;i++) if ( v[i] > 4.0*med*1.0e3 ) { tail+=v[i]; ntail++; }
      std::cout << GridLogMessage << "Schur comms distribution (boss):"
                << "  p50 "  << med
                << "  p90 "  << pc(0.90)
                << "  p99 "  << pc(0.99)
                << " ms   calls > 4x median: " << ntail << "/" << n
                << " carrying " << tail/1.0e6 << " s"
                << std::endl;
    }
  }
};

NAMESPACE_END(Grid);
