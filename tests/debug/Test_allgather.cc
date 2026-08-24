/*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_allgather.cc

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See the full license in the file "LICENSE" in the top level distribution
    directory
*************************************************************************************/
/*  END LEGAL */

//////////////////////////////////////////////////////////////////////////////
// Regression gate for CartesianCommunicator::AllGather / AllGatherV.
//
// CPU build under mpirun:
//   mpirun -n 1 ./Test_allgather --grid 16.16.16.32 --mpi 1.1.1.1
//   mpirun -n 2 ./Test_allgather --grid 16.16.16.32 --mpi 1.1.1.2
//   mpirun -n 3 ./Test_allgather --grid 16.16.16.48 --mpi 1.1.1.3
//   mpirun -n 4 ./Test_allgather --grid 16.16.16.32 --mpi 1.1.1.4
//
//   T1 : uniform AllGather, rank-ordered concatenation.
//   T2 : AllGatherV with non-uniform counts and displacements.
//   T3 : bitwise equivalence with the zero-fill + GlobalSumVector idiom
//        that AllGatherV is intended to replace.
//   T4 : the same equivalence for the GatherGemm panel layout -- a
//        column-major kchunk x n panel whose rows are owned in contiguous
//        rank-major ranges, requiring a pack on send and a repack on receive.
//
// T3 and T4 are the ones that matter: they certify AllGatherV as a drop-in
// for the existing idiom, bit for bit, before it is used anywhere.
//////////////////////////////////////////////////////////////////////////////

#include <Grid/Grid.h>

using namespace Grid;

static int failures = 0;

static void Report(const std::string &name, bool pass, const std::string &detail="")
{
  std::cout << GridLogMessage << "  " << name << (pass ? "  PASS" : "  ** FAIL **");
  if ( detail.size() ) std::cout << "   " << detail;
  std::cout << std::endl;
  if ( !pass ) failures++;
}

// A value that is unique to (rank, index) so that any misplacement shows up.
static ComplexD Stamp(int rank, int64_t idx)
{
  return ComplexD( 1000.0*(rank+1) + (double)idx, -(double)(rank+1) );
}

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                       GridDefaultSimd(Nd, vComplexD::Nsimd()),
                                                       GridDefaultMpi());
  const int P  = grid->ProcessorCount();
  const int me = grid->ThisRank();

  std::cout << GridLogMessage << "AllGather regression: " << P << " ranks" << std::endl;

  ////////////////////////////////////////////////////////////////////////
  // T1 : uniform AllGather
  ////////////////////////////////////////////////////////////////////////
  {
    const int64_t W = 37;                       // words per rank
    std::vector<ComplexD> mine(W), all((uint64_t)W*P);
    for(int64_t i=0;i<W;i++) mine[i] = Stamp(me,i);

    grid->AllGather((void *)&mine[0], (void *)&all[0], (uint64_t)W, sizeof(ComplexD));

    bool ok = true;
    for(int r=0;r<P;r++)
      for(int64_t i=0;i<W;i++)
        if ( all[(uint64_t)r*W + i] != Stamp(r,i) ) ok = false;
    Report("T1  uniform AllGather", ok);
  }

  ////////////////////////////////////////////////////////////////////////
  // T2 : AllGatherV, non-uniform counts
  ////////////////////////////////////////////////////////////////////////
  {
    std::vector<int> counts(P), displs(P);
    int64_t total = 0;
    for(int r=0;r<P;r++){ counts[r] = r+1; displs[r] = (int)total; total += counts[r]; }

    std::vector<ComplexD> mine(counts[me]), all(total);
    for(int64_t i=0;i<counts[me];i++) mine[i] = Stamp(me,i);

    grid->AllGatherV((void *)&mine[0], counts[me],
                     (void *)&all[0], counts, displs, sizeof(ComplexD));

    bool ok = true;
    for(int r=0;r<P;r++)
      for(int64_t i=0;i<counts[r];i++)
        if ( all[displs[r]+i] != Stamp(r,i) ) ok = false;
    Report("T2  AllGatherV non-uniform", ok);
  }

  ////////////////////////////////////////////////////////////////////////
  // T3 : bitwise equivalence with zero-fill + GlobalSumVector
  ////////////////////////////////////////////////////////////////////////
  {
    std::vector<int> counts(P), displs(P);
    int64_t total = 0;
    for(int r=0;r<P;r++){ counts[r] = 5 + (r%3); displs[r] = (int)total; total += counts[r]; }

    std::vector<ComplexD> mine(counts[me]);
    for(int64_t i=0;i<counts[me];i++) mine[i] = Stamp(me,i);

    // (a) the existing idiom
    std::vector<ComplexD> viaSum(total, ComplexD(0.0,0.0));
    for(int64_t i=0;i<counts[me];i++) viaSum[displs[me]+i] = mine[i];
    grid->GlobalSumVector(&viaSum[0], (int)total);

    // (b) the primitive
    std::vector<ComplexD> viaGather(total);
    grid->AllGatherV((void *)&mine[0], counts[me],
                     (void *)&viaGather[0], counts, displs, sizeof(ComplexD));

    bool ok = true;
    for(int64_t i=0;i<total;i++) if ( viaSum[i] != viaGather[i] ) ok = false;
    Report("T3  bitwise == zero-fill+GlobalSumVector", ok);
  }

  ////////////////////////////////////////////////////////////////////////
  // T4 : GatherGemm panel layout.
  //
  // Panel is kchunk x n, COLUMN major, ld = kchunk.  Rank r owns panel rows
  // [rowStart[r],rowStart[r+1]) and holds them locally as a rows_r x n
  // column-major block with ld = rows_r.  Its contribution to the panel is
  // therefore strided, so AllGatherV needs a pack on send (contiguous
  // rows_r x n, which the local block already is) and a repack on receive.
  ////////////////////////////////////////////////////////////////////////
  {
    const int64_t n = 6;                                   // panel columns
    std::vector<int64_t> rowStart(P+1);
    rowStart[0]=0;
    for(int r=0;r<P;r++) rowStart[r+1] = rowStart[r] + (3 + (r%4));   // uneven ownership
    const int64_t kchunk = rowStart[P];
    const int64_t myRows = rowStart[me+1]-rowStart[me];

    // local block, column major, ld = myRows
    std::vector<ComplexD> mine((uint64_t)myRows*n);
    for(int64_t j=0;j<n;j++)
      for(int64_t i=0;i<myRows;i++)
        mine[i + j*myRows] = Stamp(me, i + j*myRows);

    // (a) zero-fill + deposit + GlobalSumVector
    std::vector<ComplexD> viaSum((uint64_t)kchunk*n, ComplexD(0.0,0.0));
    for(int64_t j=0;j<n;j++)
      for(int64_t i=0;i<myRows;i++)
        viaSum[rowStart[me]+i + j*kchunk] = mine[i + j*myRows];
    grid->GlobalSumVector(&viaSum[0], (int)(kchunk*n));

    // (b) AllGatherV into rank-major order, then repack into the panel
    std::vector<int> counts(P), displs(P);
    int64_t tot=0;
    for(int r=0;r<P;r++){
      counts[r] = (int)((rowStart[r+1]-rowStart[r])*n);
      displs[r] = (int)tot;
      tot += counts[r];
    }
    std::vector<ComplexD> rankMajor(tot);
    grid->AllGatherV((void *)&mine[0], counts[me],
                     (void *)&rankMajor[0], counts, displs, sizeof(ComplexD));

    std::vector<ComplexD> viaGather((uint64_t)kchunk*n);
    for(int r=0;r<P;r++){
      int64_t rows_r = rowStart[r+1]-rowStart[r];
      for(int64_t j=0;j<n;j++)
        for(int64_t i=0;i<rows_r;i++)
          viaGather[rowStart[r]+i + j*kchunk] = rankMajor[displs[r] + i + j*rows_r];
    }

    bool ok = true;
    for(int64_t i=0;i<kchunk*n;i++) if ( viaSum[i] != viaGather[i] ) ok = false;
    Report("T4  GatherGemm panel, pack+gather+repack", ok,
           "panel "+std::to_string(kchunk)+"x"+std::to_string(n));
  }

  ////////////////////////////////////////////////////////////////////////
  // T5 : DEVICE buffers, full participation.
  // T6 : DEVICE buffers, only ranks [0,P/2) contribute -- the shape
  //      GatherGemm actually produces, where at depth d only P/2^(d+1)
  //      ranks have a non-zero count and everyone else sends zero.
  //
  // These exist because AllGatherV on device pointers aborted inside
  // PMPI_Allgatherv on Cray MPICH ("req != NULL", mpir_request.h:508) at
  // 288 ranks, while the host-buffer stages above passed.  Reproducing that
  // here costs one node and a second instead of a 36-node job.
  ////////////////////////////////////////////////////////////////////////
  // AG_WORDS : ComplexD words each CONTRIBUTING rank sends (default 4096)
  // AG_NPART : how many ranks contribute at all in T6 (default P/2).
  //   The GatherGemm shape at Schur depth d has only P/2^(d+1) contributors,
  //   so depth 7 at 288 ranks is AG_NPART=2 AG_WORDS=518400 (8.3 MB each,
  //   18.7 MB gathered).  Depth 0 is AG_NPART=144 AG_WORDS=232800.
  const int64_t AGW = getenv("AG_WORDS") ? atol(getenv("AG_WORDS")) : 4096;
  const int     AGN = getenv("AG_NPART") ? atoi(getenv("AG_NPART")) : (P+1)/2;
  // Stage skips.  T6 is the pathological one: at 288 ranks with 9 contributors
  // it takes ~26 minutes and aborts on at least one rank, which would prevent
  // T7 from ever running.  T7 is therefore ordered FIRST, and AG_T6=0 skips
  // the slow stage entirely when only the Bcast answer is wanted.
  const int     doT5 = getenv("AG_T5") ? atoi(getenv("AG_T5")) : 1;
  const int     doT6 = getenv("AG_T6") ? atoi(getenv("AG_T6")) : 1;
  const int     doT7 = getenv("AG_T7") ? atoi(getenv("AG_T7")) : 1;

  if ( doT5 ) {
    const int64_t W = AGW;                        // words per contributing rank
    std::vector<int> counts(P,0), displs(P,0);
    int64_t total=0;
    for(int r=0;r<P;r++){ counts[r]=(int)W; displs[r]=(int)total; total+=W; }

    deviceVector<ComplexD> dmine(W), dall(total);
    std::vector<ComplexD>  hmine(W), hall(total);
    for(int64_t i=0;i<W;i++) hmine[i]=Stamp(me,i);
    acceleratorCopyToDevice(&hmine[0],&dmine[0],W*sizeof(ComplexD));

    grid->AllGatherV((void *)&dmine[0], counts[me],
                     (void *)&dall[0], counts, displs, sizeof(ComplexD));

    acceleratorCopyFromDevice(&dall[0],&hall[0],total*sizeof(ComplexD));
    bool ok=true;
    for(int r=0;r<P;r++)
      for(int64_t i=0;i<W;i++)
        if ( hall[displs[r]+i] != Stamp(r,i) ) ok=false;
    Report("T5  AllGatherV on DEVICE, all "+std::to_string(P)+" contributing, "+
           std::to_string(W*16/1024)+" KB each, "+
           std::to_string(total*16/1048576)+" MB gathered", ok);
  }

  ////////////////////////////////////////////////////////////////////////
  // T7 : the SAME shape as T6, but assembled by C sequential MPI_Bcast --
  // one broadcast per contributing rank -- instead of one MPI_Allgatherv.
  //
  // This is the transport of DENSE_GATHER=2.  Bcast takes no count vector,
  // so the zero-count asymmetry that makes T6 run at ~0.18 MB/s and trip
  // mpir_request.h:508 cannot arise.  It costs C collectives rather than 1
  // and the roots do not transmit concurrently, so the byte cost is about
  // 2x that of an ideal allgather; it is worth having only if Bcast is
  // faster PER BYTE than the alternatives.  Compare the T6 and T7 times.
  ////////////////////////////////////////////////////////////////////////
  if ( doT7 ) {
    const int64_t W = AGW;
    int nroot = AGN < 1 ? 1 : (AGN > P ? P : AGN);
    int64_t total = (int64_t)nroot*W;

    deviceVector<ComplexD> dall(total>0?total:1);
    std::vector<ComplexD>  hall(total), hmine(W);
    for(int64_t i=0;i<W;i++) hmine[i]=Stamp(me,i);
    if ( me < nroot )                                  // roots seed their slot
      acceleratorCopyToDevice(&hmine[0],&dall[(int64_t)me*W],W*sizeof(ComplexD));

    for(int r=0;r<nroot;r++)
      grid->Broadcast(r,(void *)&dall[(int64_t)r*W],(uint64_t)W*sizeof(ComplexD));

    acceleratorCopyFromDevice(&dall[0],&hall[0],total*sizeof(ComplexD));
    bool ok=true;
    for(int r=0;r<nroot;r++)
      for(int64_t i=0;i<W;i++)
        if ( hall[(int64_t)r*W+i] != Stamp(r,i) ) ok=false;
    Report("T7  "+std::to_string(nroot)+" x MPI_Bcast on DEVICE, "+
           std::to_string(W*16/1024)+" KB each, "+
           std::to_string(total*16/1048576)+" MB assembled", ok);
  }

  if ( doT6 ) {
    // Only the first half of the ranks contribute; the rest send count 0.
    const int64_t W = AGW;
    int half = AGN < 1 ? 1 : (AGN > P ? P : AGN);
    std::vector<int> counts(P,0), displs(P,0);
    int64_t total=0;
    for(int r=0;r<P;r++){
      counts[r] = (r<half) ? (int)W : 0;
      displs[r] = (int)total;
      total += counts[r];
    }
    deviceVector<ComplexD> dmine(W>0?W:1), dall(total>0?total:1);
    std::vector<ComplexD>  hmine(W), hall(total);
    for(int64_t i=0;i<W;i++) hmine[i]=Stamp(me,i);
    if ( counts[me] ) acceleratorCopyToDevice(&hmine[0],&dmine[0],W*sizeof(ComplexD));

    grid->AllGatherV((void *)&dmine[0], counts[me],
                     (void *)&dall[0], counts, displs, sizeof(ComplexD));

    acceleratorCopyFromDevice(&dall[0],&hall[0],total*sizeof(ComplexD));
    bool ok=true;
    for(int r=0;r<half;r++)
      for(int64_t i=0;i<W;i++)
        if ( hall[displs[r]+i] != Stamp(r,i) ) ok=false;
    Report("T6  AllGatherV on DEVICE, "+std::to_string(half)+"/"+std::to_string(P)+
           " contributing, "+std::to_string(W*16/1024)+" KB each, "+
           std::to_string(total*16/1048576)+" MB gathered", ok);
  }

  // Reduce across ranks before reporting.  Report() prints from rank 0 only
  // and `failures` is rank-local, so a stage that passes on rank 0 and fails
  // (or aborts) elsewhere would otherwise be announced as ALL PASS -- which
  // is exactly what happened at 288 ranks with 9 contributors.
  {
    uint64_t f = failures;
    grid->GlobalSum(f);
    if ( f && !failures )
      std::cout << GridLogMessage << "  ** failures on OTHER ranks: " << f
                << " (this rank saw none) **" << std::endl;
    failures = (int)f;
  }

  std::cout << GridLogMessage << (failures ? "AllGather regression: FAILURES" 
                                           : "AllGather regression: ALL PASS") << std::endl;

  Grid_finalize();
  return failures ? 1 : 0;
}
