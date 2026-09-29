/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_sloppy_dagger.cc

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

//
// Isolate a sloppy-comms dagger defect: on a CPU/GEN multi-rank build the
// sloppy dagger halo can be wrong by ~5% of norm2, deterministically, which
// Benchmark_dwf's sloppy pass catches as a failed dagger Cshift check.
//
// The sloppy compressed-buffer pool (StencilBuffer::DeviceCommBuf) is a
// SHARED STATIC across all stencil objects, so scenarios contaminate each
// other inside one process: this program runs exactly ONE scenario per
// invocation, selected with --seq <name>:
//
//   fresh-dag       sloppy dagger is the operator's first call
//   nodag-dag       sloppy nodag x3 (one source), then sloppy dagger (same source)
//   nodag-dag-fresh nodag x3 (fresh source each), then dagger (fresh source)
//   nodag-dag-newsrc nodag x3 (one source), then dagger (different source)
//   dag-nodag-dag   dagger, nodag x3, dagger  (dagger-first history)
//
// Every call is checked against a separate never-sloppy reference
// operator (exact Dhop == the Cshift construction at 1e-31, certified by
// Benchmark_dwf's exact pass).  The last dagger's error is fingerprinted
// by slice along each MPI-decomposed direction.
//
//   OMP_NUM_THREADS=1 mpirun -n 2 ./Test_sloppy_dagger --grid 8.8.8.16 \
//        --mpi 1.1.1.2 --seq nodag-dag
//
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::string seq("nodag-dag");
  if( GridCmdOptionExists(argv,argv+argc,"--seq") )
    seq = GridCmdOptionPayload(argv,argv+argc,"--seq");

  const int Ls=8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                                   GridDefaultSimd(Nd,vComplex::Nsimd()),
                                                                   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  GridParallelRNG RNG4(UGrid); RNG4.SeedFixedIntegers({1,2,3,4});
  GridParallelRNG RNG5(FGrid); RNG5.SeedFixedIntegers({5,6,7,8});

  LatticeGaugeField Umu(UGrid);
  SU<Nc>::HotConfiguration(RNG4,Umu);

  RealD mass=0.1, M5=1.8;
  DomainWallFermionD Dref(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);  // never sloppy
  DomainWallFermionD Dw  (Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);  // under test
  Dw.SloppyComms(1);

  Coordinate mpi = GridDefaultMpi();
  LatticeFermionD src(FGrid), exact(FGrid), sloppy(FGrid), err(FGrid);

  // One call: fresh or reused source, dag or not; always checked.
  LatticeFermionD srcA(FGrid); gaussian(RNG5,srcA);
  auto call = [&](int dag, int fresh, const char *tag){
    if ( fresh ) gaussian(RNG5,src); else src = srcA;
    Dref.Dhop(src,exact,dag);
    Dw.Dhop(src,sloppy,dag);
    err = sloppy - exact;
    std::cout << GridLogMessage << "SEQ[" << seq << "] " << tag
              << (dag?" DAG   ":" NODAG ") << " norm2(err) " << norm2(err) << std::endl;
  };
  // First bad site: coordinate + exact-vs-sloppy hex words (PB: look for
  // expected/actual mismatching in the high-order half word).
  auto firstbad = [&](LatticeFermionD &ex, LatticeFermionD &sl){
    typedef LatticeFermionD::scalar_object sobj;
    Coordinate gdims = FGrid->GlobalDimensions();
    int64_t gsites = FGrid->gSites();
    for(int64_t g=0; g<gsites; g++){
      Coordinate gcoor(5);
      Lexicographic::CoorFromIndex(gcoor,g,gdims);
      sobj se, ss_;
      peekSite(se,ex,gcoor);
      peekSite(ss_,sl,gcoor);
      uint64_t *we=(uint64_t *)&se;
      uint64_t *ws=(uint64_t *)&ss_;
      int nw = sizeof(sobj)/8;
      for(int w=0;w<nw;w++){
        double de=((double *)we)[w], ds=((double *)ws)[w];
        if ( fabs(de-ds) > 1.0e-4 ) {
          std::cout << GridLogMessage << "FIRSTBAD site " << gcoor
                    << " word " << w << " (spin "<<(w/6)%4<<" col "<<(w/2)%3<<" reim "<<(w%2)<<")"
                    << "  exact " << std::hex << we[w]
                    << "  sloppy " << ws[w] << std::dec << std::endl;
          for(int w2=0;w2<8&&w2<nw;w2++)
            std::cout << GridLogMessage << "   word["<<w2<<"] exact "<<std::hex<<we[w2]
                      <<" sloppy "<<ws[w2]<<std::dec<<std::endl;
          return;
        }
      }
    }
    std::cout << GridLogMessage << "FIRSTBAD: none above 1e-4" << std::endl;
  };
  auto fingerprint = [&](void){
    for(int mu=0;mu<Nd;mu++){
      if ( mpi[mu] == 1 ) continue;
      std::vector<RealD> sn;
      sliceNorm(sn,err,mu+1);
      std::cout << GridLogMessage << "SEQ[" << seq << "] err by slice of dim " << mu << ":";
      for(int t=0;t<(int)sn.size();t++) std::cout << " " << sn[t];
      std::cout << std::endl;
    }
  };

  if        ( seq == "fresh-dag" ) {
    call(1,0,"call0");
  } else if ( seq == "nodag-dag" ) {
    call(0,0,"call0"); call(0,0,"call1"); call(0,0,"call2");
    call(1,0,"call3");
  } else if ( seq == "nodag-dag-fresh" ) {
    call(0,1,"call0"); call(0,1,"call1"); call(0,1,"call2");
    call(1,1,"call3");
  } else if ( seq == "nodag-dag-newsrc" ) {
    call(0,0,"call0"); call(0,0,"call1"); call(0,0,"call2");
    call(1,1,"call3");
  } else if ( seq == "dag-nodag-dag" ) {
    call(1,0,"call0");
    call(0,0,"call1"); call(0,0,"call2"); call(0,0,"call3");
    call(1,0,"call4");
  } else if ( seq == "noleave" ) {
    // NO exact-operator call between sloppy calls: references precomputed.
    LatticeFermionD refnodag(FGrid), refdag(FGrid);
    Dref.Dhop(srcA,refnodag,DaggerNo);
    Dref.Dhop(srcA,refdag,DaggerYes);
    for(int i=0;i<2;i++){
      Dw.Dhop(srcA,sloppy,DaggerNo);
      err = sloppy - refnodag;
      std::cout << GridLogMessage << "SEQ[noleave] call" << i << " NODAG  norm2(err) " << norm2(err) << std::endl;
      if ( i==0 ) firstbad(refnodag,sloppy);
    }
    Dw.Dhop(srcA,sloppy,DaggerYes);
    err = sloppy - refdag;
    std::cout << GridLogMessage << "SEQ[noleave] call3 DAG    norm2(err) " << norm2(err) << std::endl;
  } else if ( seq == "twoop" ) {
    // A SECOND sloppy operator shares the static pool; sequence runs on it
    // with interleaved exact references (as the honest runs had).
    DomainWallFermionD Dw2(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
    Dw2.SloppyComms(1);
    Dw.Dhop(srcA,sloppy,DaggerNo);   // prime the FIRST sloppy op's state
    for(int i=0;i<3;i++){
      Dref.Dhop(srcA,exact,DaggerNo);
      Dw2.Dhop(srcA,sloppy,DaggerNo);
      err = sloppy - exact;
      std::cout << GridLogMessage << "SEQ[twoop] call" << i << " NODAG  norm2(err) " << norm2(err) << std::endl;
    }
    Dref.Dhop(srcA,exact,DaggerYes);
    Dw2.Dhop(srcA,sloppy,DaggerYes);
    err = sloppy - exact;
    std::cout << GridLogMessage << "SEQ[twoop] call3 DAG    norm2(err) " << norm2(err) << std::endl;
  } else if ( seq == "self-exact" ) {
    // Exact call on the SAME object (Dw), opposite sense, then sloppy:
    // distinguishes per-object from shared state.
    LatticeFermionD refnodag(FGrid);
    Dref.Dhop(srcA,refnodag,DaggerNo);              // reference only, early
    Dw.SloppyComms(0);
    Dw.Dhop(srcA,sloppy,DaggerYes);                 // exact DAG on Dw itself
    Dw.SloppyComms(1);
    Dw.Dhop(srcA,sloppy,DaggerNo);                  // sloppy NODAG: sense mismatch
    err = sloppy - refnodag;
    std::cout << GridLogMessage << "SEQ[self-exact] sloppy NODAG after own exact DAG: " << norm2(err) << std::endl;
  } else if ( seq == "cold" ) {
    // NO exact Dhop anywhere before the sloppy calls.
    LatticeFermionD outn(FGrid), outd(FGrid);
    Dw.Dhop(srcA,outn,DaggerNo);
    Dw.Dhop(srcA,outd,DaggerYes);
    LatticeFermionD refnodag(FGrid), refdag(FGrid);
    Dref.Dhop(srcA,refnodag,DaggerNo);
    Dref.Dhop(srcA,refdag,DaggerYes);
    err = outn - refnodag;
    std::cout << GridLogMessage << "SEQ[cold] first-ever sloppy NODAG: " << norm2(err) << std::endl;
    err = outd - refdag;
    std::cout << GridLogMessage << "SEQ[cold] then      sloppy DAG:   " << norm2(err) << std::endl;
  } else {
    std::cout << GridLogMessage << "unknown --seq " << seq << std::endl;
  }
  fingerprint();

  Grid_finalize();
}
