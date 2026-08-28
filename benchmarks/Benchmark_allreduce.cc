/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./benchmarks/Benchmark_allreduce.cc

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
// Vector all-reduce shoot-out: MPI_Allreduce vs the P2P rings, across sizes,
// so GlobalSumVector's dispatch (if any) can be decided from a table rather
// than a guess.  For each size, each type, the SAME data are reduced by:
//
//   MPI-host      GlobalSumVector on a host buffer          (today's path)
//   MPI-bare      MPI_Allreduce(MPI_IN_PLACE,..) directly on grid->communicator,
//                 host buffer: GlobalSumVector minus its FlightRecorder step log
//   MPI-dev       GlobalSumVector on a device buffer        (GPU-aware MPI;
//                                                            Cray MPICH aborts
//                                                            above ~8 MB: capped
//                                                            by BENCH_MPI_DEV_MAX_MB)
//   cart-dev      CartesianRingAllReduce, device buffer     (~sum_d 2N bytes)
//   flat-dev      RingAllReduce, device buffer              (2N(P-1)/P bytes)
//   cart-host     H2D copy + cart-dev + D2H                 (the host path when
//   flat-host     H2D copy + flat-dev + D2H                  the caller holds a
//                                                            std::vector)
//
// Every method is checked against MPI-host; the deterministic rings are also
// checked bitwise against themselves on a repeat.  Reported per method: time,
// payload rate (N bytes / t), and "wire" rate (bytes actually moved by that
// method / t), min over reps.  Rerun at 1, 4 and 36 nodes: the crossovers (if
// any) move with P as sum_d 2(P_d-1) vs 2(P-1) steps.
//
// SCALAR section (latency, min/mean of many reps): GlobalSum(RealD) -- which is
// the latency-hiding GlobalSumP2P (all shifts of a dimension posted at once,
// one completion per dimension) -- vs bare MPI_Allreduce of one double vs
// GlobalSumVector(1) vs the two rings at n=1.
//
//   mpirun -n 8 ./Benchmark_allreduce --grid 16.16.16.32 --mpi 1.1.2.4
//   env: BENCH_MIN_KB (4)  BENCH_MAX_MB (512; 0 = scalar section only)  BENCH_REPS (5)
//        BENCH_MPI_DEV_MAX_MB (4: above this MPI-dev is skipped, not attempted)
//        BENCH_SCALAR_REPS (200)
//        BENCH_SCALAR_GAP_US (0): a device kernel of about this many us between
//        scalar reps.  Back-to-back sums (gap 0) measured GlobalSumP2P at 12.8 ms
//        mean / 0.4 ms min at P=288 vs MPI 45 us (2026-08-27), yet the solver's
//        fine-smoother steps (16 ms incl. a ~10 ms matvec) cannot be paying that
//        per reduction: the solver spaces its sums by ms of GPU work.  The gap
//        reproduces the solver's condition; if P2P is fast with a gap and slow
//        without, the tight loop is the pathology, not the primitive.
//////////////////////////////////////////////////////////////////////////////

#include <Grid/Grid.h>

using namespace Grid;

template<class T> T Fill(uint64_t i,int rank){ return T(0.5*std::sin(0.01*i+0.7*rank)+1.0e-3*rank); }

template<class T> MPI_Datatype MpiType(void);
template<> MPI_Datatype MpiType<RealF>(void){ return MPI_FLOAT; }
template<> MPI_Datatype MpiType<RealD>(void){ return MPI_DOUBLE; }
template<class T> void BareAllreduce(GridCartesian *grid, T *buf, uint64_t n)
{
#ifdef GRID_COMMS_MPI3
  int ierr = MPI_Allreduce(MPI_IN_PLACE, buf, (int)n, MpiType<T>(), MPI_SUM, grid->communicator);
  GRID_ASSERT(ierr==0);
#else
  grid->GlobalSumVector(buf,(int)n);
#endif
}

template<class T> struct Stats { double tmin=1e30, tmax=0, tsum=0; int n=0; void add(double t){ tmin=std::min(tmin,t); tmax=std::max(tmax,t); tsum+=t; n++; } };

template<class T>
void Run(GridCartesian *grid, const char *tname, uint64_t nbytes_lo, uint64_t nbytes_hi, int reps, uint64_t mpiDevMax)
{
  int P  = grid->ProcessorCount();
  int me = grid->ThisRank();
  int Nd = grid->_ndimension;
  // bytes moved per rank by each method, per byte of payload
  double flatFactor = 2.0*(P-1)/P;
  double cartFactor = 0.0; for(int d=0;d<Nd;d++) if(grid->_processors[d]>1) cartFactor += 2.0*(grid->_processors[d]-1)/grid->_processors[d];

  if ( me==0 ) {
    std::cout << GridLogMessage << "==== " << tname << "  P=" << P << " grid " << grid->_processors
              << "  bytes moved per payload byte: flat " << flatFactor << "  cartesian " << cartFactor << std::endl;
    std::cout << GridLogMessage << std::setw(10) << "bytes" << std::setw(12) << "elements"
              << " | " << std::setw(9) << "MPI-host" << std::setw(9) << "MPI-bare" << std::setw(9) << "MPI-dev" << std::setw(9) << "cart-dev" << std::setw(9) << "flat-dev"
              << std::setw(10) << "cart-host" << std::setw(10) << "flat-host" << "   (ms, min of " << reps << ")"
              << " | payload GB/s: " << "cart-dev flat-dev | wire GB/s: cart-dev flat-dev | check" << std::endl;
  }
  for(uint64_t bytes=nbytes_lo; bytes<=nbytes_hi; bytes*=2){
    uint64_t n = bytes/sizeof(T);
    if ( n==0 ) continue;
    std::vector<T> h(n); for(uint64_t i=0;i<n;i++) h[i]=Fill<T>(i,me);
    std::vector<T> ref(h); grid->GlobalSumVector(&ref[0],(int)n);         // MPI-host reference (also warm-up)
    std::vector<T> out(n);
    deviceVector<T> d(n);
    Stats<T> sMh, sMb, sMd, sCd, sFd, sCh, sFh;
    double worst[7]={0,0,0,0,0,0,0}; int done[7]={1,0,1,1,1,1,1};
    auto check=[&](int k, const std::vector<T> &o){ double w=0; for(uint64_t i=0;i<n;i++) w=std::max(w,(double)std::abs(o[i]-ref[i])); worst[k]=std::max(worst[k],w); };
    int doMpiDev = (bytes <= mpiDevMax); done[1]=doMpiDev;
    for(int r=0;r<reps;r++){
      // MPI-host
      out=h; { double t0=usecond(); grid->GlobalSumVector(&out[0],(int)n); sMh.add(usecond()-t0); } check(0,out);
      // MPI-bare
      out=h; { double t0=usecond(); BareAllreduce(grid,&out[0],n); sMb.add(usecond()-t0); } check(6,out);
      // MPI-dev
      if ( doMpiDev ) { acceleratorCopyToDevice(&h[0],&d[0],n*sizeof(T)); double t0=usecond(); grid->GlobalSumVector(&d[0],(int)n); sMd.add(usecond()-t0);
        acceleratorCopyFromDevice(&d[0],&out[0],n*sizeof(T)); check(1,out); }
      // cart-dev
      { acceleratorCopyToDevice(&h[0],&d[0],n*sizeof(T)); double t0=usecond(); CartesianRingAllReduce(grid,&d[0],n); sCd.add(usecond()-t0);
        acceleratorCopyFromDevice(&d[0],&out[0],n*sizeof(T)); check(2,out); }
      // flat-dev
      { acceleratorCopyToDevice(&h[0],&d[0],n*sizeof(T)); double t0=usecond(); RingAllReduce(grid,&d[0],n); sFd.add(usecond()-t0);
        acceleratorCopyFromDevice(&d[0],&out[0],n*sizeof(T)); check(3,out); }
      // cart-host: staged
      { out=h; double t0=usecond(); acceleratorCopyToDevice(&out[0],&d[0],n*sizeof(T)); CartesianRingAllReduce(grid,&d[0],n); acceleratorCopyFromDevice(&d[0],&out[0],n*sizeof(T)); sCh.add(usecond()-t0); check(4,out); }
      // flat-host: staged
      { out=h; double t0=usecond(); acceleratorCopyToDevice(&out[0],&d[0],n*sizeof(T)); RingAllReduce(grid,&d[0],n); acceleratorCopyFromDevice(&d[0],&out[0],n*sizeof(T)); sFh.add(usecond()-t0); check(5,out); }
    }
    // rings: bitwise repeatable?
    int bitwise=1;
    { std::vector<T> a(n),b(n);
      acceleratorCopyToDevice(&h[0],&d[0],n*sizeof(T)); CartesianRingAllReduce(grid,&d[0],n); acceleratorCopyFromDevice(&d[0],&a[0],n*sizeof(T));
      acceleratorCopyToDevice(&h[0],&d[0],n*sizeof(T)); CartesianRingAllReduce(grid,&d[0],n); acceleratorCopyFromDevice(&d[0],&b[0],n*sizeof(T));
      if ( memcmp(&a[0],&b[0],n*sizeof(T)) ) bitwise=0; }
    // worst error over ranks; tolerance scales with P and magnitude
    double scale=0; for(uint64_t i=0;i<n;i++) scale=std::max(scale,(double)std::abs(ref[i]));
    RealD wmax=0; for(int k=0;k<7;k++) if(done[k]) wmax=std::max(wmax,worst[k]); grid->GlobalMax(wmax);
    RealD bw=1.0-bitwise; grid->GlobalMax(bw);
    double tol = (sizeof(T)==4 ? 1.0e-5 : 1.0e-12) * std::max(scale,1.0) * P;
    // times: use the slowest rank's min (the collective completes when the last rank does)
    auto mx=[&](double v){ RealD x=v; grid->GlobalMax(x); return (double)x; };
    double tMh=mx(sMh.tmin), tMb=mx(sMb.tmin), tMd=doMpiDev?mx(sMd.tmin):0, tCd=mx(sCd.tmin), tFd=mx(sFd.tmin), tCh=mx(sCh.tmin), tFh=mx(sFh.tmin);
    if ( me==0 ) {
      std::streamsize prec=std::cout.precision();
      std::cout << GridLogMessage << std::setw(10) << bytes << std::setw(12) << n << " | " << std::fixed << std::setprecision(3)
                << std::setw(9) << tMh/1000. << std::setw(9) << tMb/1000. << std::setw(9) << (doMpiDev ? tMd/1000. : -1.0)
                << std::setw(9) << tCd/1000. << std::setw(9) << tFd/1000. << std::setw(10) << tCh/1000. << std::setw(10) << tFh/1000.
                << "   | " << std::setprecision(2) << std::setw(8) << bytes/tCd/1.0e3 << std::setw(9) << bytes/tFd/1.0e3
                << " | " << std::setw(8) << cartFactor*bytes/tCd/1.0e3 << std::setw(9) << flatFactor*bytes/tFd/1.0e3
                << " | " << (wmax<tol ? "ok" : "**MISMATCH**") << (bw>0 ? " **ring not bitwise**" : "")
                << (doMpiDev ? "" : "  (MPI-dev skipped: > BENCH_MPI_DEV_MAX_MB)") << std::endl;
      std::cout.unsetf(std::ios::fixed); std::cout.precision(prec);
    }
  }
}

void Scalar(GridCartesian *grid, int reps, double gap_us)
{
  int me=grid->ThisRank();
  Stats<RealD> sP2P, sBare, sVec, sCart, sFlat, sP2Pdirect;
  RealD x; deviceVector<RealD> d(1); std::vector<RealD> h(1);
  // gap: a bandwidth-bound device kernel sized to take ~gap_us (calibrated once)
  uint64_t gapN = 0; deviceVector<RealD> gapbuf(1);
  if ( gap_us > 0 ) {
    gapN = 1<<20; gapbuf.resize(gapN);
    RealD *g=&gapbuf[0];
    accelerator_for(i,gapN,1,{ g[i]=1.0; });
    double t0=usecond(); accelerator_for(i,gapN,1,{ g[i]=g[i]*1.000001+1.0e-9; }); double t1=usecond();
    double per = (t1-t0)/gapN;                              // us per element at this size
    gapN = (uint64_t)std::max(1.0, gap_us/std::max(per,1.0e-6)); gapbuf.resize(gapN);
    g=&gapbuf[0]; accelerator_for(i,gapN,1,{ g[i]=1.0; });
    double t2=usecond(); accelerator_for(i,gapN,1,{ g[i]=g[i]*1.000001+1.0e-9; }); double t3=usecond();
    if ( me==0 ) std::cout << GridLogMessage << "Scalar gap kernel: " << gapN << " elements, " << (t3-t2) << " us (requested " << gap_us << ")" << std::endl;
  }
  auto gap=[&](void){ if(gapN){ RealD *g=&gapbuf[0]; accelerator_for(i,gapN,1,{ g[i]=g[i]*1.000001+1.0e-9; }); } };
  // A Barrier before each timed call so that node skew (from the gap kernel or
  // anything else) is not charged to the reduction: the timer measures the
  // collective from a synchronised start.
  for(int r=0;r<reps;r++){
    gap(); x=1.0+me;   grid->Barrier(); { double t0=usecond(); grid->GlobalSum(x);                 sP2P.add(usecond()-t0); }
    gap(); x=1.0+me;   grid->Barrier(); { double t0=usecond(); grid->GlobalSumP2P(x);              sP2Pdirect.add(usecond()-t0); }
    gap(); h[0]=1.0+me;grid->Barrier(); { double t0=usecond(); BareAllreduce(grid,&h[0],1);        sBare.add(usecond()-t0); }
    gap(); h[0]=1.0+me;grid->Barrier(); { double t0=usecond(); grid->GlobalSumVector(&h[0],1);     sVec.add(usecond()-t0); }
    gap(); acceleratorPut(d[0],h[0]); grid->Barrier(); { double t0=usecond(); CartesianRingAllReduce(grid,&d[0],1); sCart.add(usecond()-t0); }
    gap(); acceleratorPut(d[0],h[0]); grid->Barrier(); { double t0=usecond(); RingAllReduce(grid,&d[0],1);          sFlat.add(usecond()-t0); }
  }
  auto mx=[&](double v){ RealD y=v; grid->GlobalMax(y); return (double)y; };
  auto line=[&](const char *nm, Stats<RealD> &st){
    double tmin=mx(st.tmin), tmean=mx(st.tsum/st.n);
    if ( me==0 ) std::cout << GridLogMessage << "  " << std::setw(34) << std::left << nm << std::right
                           << " min " << std::setw(8) << tmin << " us   mean " << std::setw(8) << tmean << " us" << std::endl; };
  if ( me==0 ) std::cout << GridLogMessage << "==== SCALAR latency (RealD, " << reps << " reps, slowest rank, Barrier before each timed call, gap " << gap_us << " us)  P=" << grid->ProcessorCount() << std::endl;
  line("GlobalSum(RealD) = GlobalSumP2P",   sP2P);
  line("GlobalSum(RealD) = GlobalSumP2Pdirect",   sP2Pdirect);
  line("bare MPI_Allreduce(1 double)",       sBare);
  line("GlobalSumVector(double*,1) [MPI]",   sVec);
  line("CartesianRingAllReduce n=1 (device)",sCart);
  line("RingAllReduce flat n=1 (device)",    sFlat);
}

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);
  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  Scalar(grid, getenv("BENCH_SCALAR_REPS") ? atoi(getenv("BENCH_SCALAR_REPS")) : 200,
               getenv("BENCH_SCALAR_GAP_US") ? atof(getenv("BENCH_SCALAR_GAP_US")) : 0.0);
  uint64_t lo   = (getenv("BENCH_MIN_KB")         ? atol(getenv("BENCH_MIN_KB"))         : 4) * 1024ull;
  uint64_t hi   = (getenv("BENCH_MAX_MB")         ? atol(getenv("BENCH_MAX_MB"))         : 512) * 1024ull*1024ull;
  int      reps = getenv("BENCH_REPS")            ? atoi(getenv("BENCH_REPS"))            : 5;
  uint64_t mdev = (getenv("BENCH_MPI_DEV_MAX_MB") ? atol(getenv("BENCH_MPI_DEV_MAX_MB")) : 4) * 1024ull*1024ull;
  std::cout << GridLogMessage << "Benchmark_allreduce: P=" << grid->ProcessorCount() << " grid " << grid->_processors
            << "  sizes " << lo << " .. " << hi << " bytes, reps " << reps << ", MPI-dev attempted up to " << mdev << " bytes" << std::endl;
  Run<RealF>(grid,"RealF (MPI_FLOAT)", lo,hi,reps,mdev);
  Run<RealD>(grid,"RealD (MPI_DOUBLE)",lo,hi,reps,mdev);
  Grid_finalize();
}
