/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/solver/Test_split_mobius_batched.cc

    Copyright (C) 2026

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
/////////////////////////////////////////////////////////////////////////////////////////////
// Production-size timing of MixedPrecisionConjugateGradientBatched for Mobius with the
// Hadrons default SchurDiagMooeeOperator: the same batch solved without and then with
// split inner solves (--batched-solver-split), reporting wall clock, per-rhs iterations
// and true residuals for each.
//
//   --Ls 12 --mass 0.026 --M5 1.8 --b 1.5 --c 0.5 --nbatch 4 --tol 1e-8
//   --config <NERSC file>    (omit for a hot configuration: timing only, not physics)
//   --nounsplit              (skip the reference unsplit solve)
//   --repeat N               (split solve N times; host RSS must not grow between them)
//
// MEMORY lines report host RSS (current and peak) and allocator cache sizes, maximum over
// ranks, at each phase: with --enable-unified=no every Lattice lives in host memory.
//
// Only one solution vector is kept, so the driver's own footprint is two batches of
// double red-black 5d fields; the solver adds about as much again.
/////////////////////////////////////////////////////////////////////////////////////////////
#include <Grid/Grid.h>
#include <sys/resource.h>
#ifdef __APPLE__
#include <mach/mach.h>
#endif

using namespace std;
using namespace Grid;

// Host memory of this process in GB: current resident set and its high-water mark
void HostRSS(RealD &current,RealD &peak)
{
  struct rusage ru;
  getrusage(RUSAGE_SELF,&ru);
#ifdef __APPLE__
  peak = ru.ru_maxrss/1.0e9;              // bytes on macOS
  mach_task_basic_info_data_t info;
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  task_info(mach_task_self(),MACH_TASK_BASIC_INFO,(task_info_t)&info,&count);
  current = info.resident_size/1.0e9;
#else
  peak = ru.ru_maxrss*1024.0/1.0e9;       // kilobytes on Linux
  long pages = 0;
  long resident = 0;
  FILE *f = fopen("/proc/self/statm","r");
  if ( f ) {
    if ( fscanf(f,"%ld %ld",&pages,&resident) != 2 ) {
      resident = 0;
    }
    fclose(f);
  }
  current = resident*(RealD)sysconf(_SC_PAGESIZE)/1.0e9;
#endif
}

// Largest values over ranks: host RSS now and at peak, and the allocator caches
void ReportMemory(GridBase *grid,const std::string &phase)
{
  RealD rss;
  RealD peak;
  HostRSS(rss,peak);
  RealD hostcache = MemoryManager::HostCacheBytes()/1.0e9;
  RealD devcache  = MemoryManager::DeviceCacheBytes()/1.0e9;
  grid->GlobalMax(rss);
  grid->GlobalMax(peak);
  grid->GlobalMax(hostcache);
  grid->GlobalMax(devcache);
  std::cout << GridLogMessage << "MEMORY " << phase
            << " : host RSS " << rss << " GB, peak " << peak
            << " GB; allocator cache host " << hostcache << " GB, device " << devcache
            << " GB (max over ranks)" << std::endl;
}

typedef LatticeFermionD FieldD;
typedef LatticeFermionF FieldF;

template<class T>
T CmdOption(int argc,char **argv,const std::string &name,T def)
{
  T val = def;
  if ( GridCmdOptionExists(argv,argv+argc,name) ) {
    std::stringstream ss(GridCmdOptionPayload(argv,argv+argc,name));
    ss >> val;
  }
  return val;
}

void SolveAndReport(const std::string &label,
                    MixedPrecisionConjugateGradientBatched<FieldD,FieldF> &mCG,
                    LinearOperatorBase<FieldD> &Linop_d,
                    std::vector<FieldD> &src,
                    std::vector<FieldD> &sol)
{
  int nbatch = src.size();
  for(int i=0;i<nbatch;i++){
    sol[i].Checkerboard() = src[i].Checkerboard();
    sol[i] = Zero();
  }

  std::cout << GridLogMessage << "==================================================" << std::endl;
  std::cout << GridLogMessage << label << " batched solve, nbatch " << nbatch << std::endl;
  std::cout << GridLogMessage << "==================================================" << std::endl;

  RealD t0 = usecond();
  mCG(src,sol);
  RealD t1 = usecond();

  FieldD Msol(src[0].Grid());
  Msol.Checkerboard() = src[0].Checkerboard();
  RealD worst = 0.0;
  for(int i=0;i<nbatch;i++){
    Linop_d.HermOp(sol[i],Msol);
    Msol = Msol - src[i];
    RealD resid = std::sqrt(norm2(Msol)/norm2(src[i]));
    worst = std::max(worst,resid);
    std::cout << GridLogMessage << label << ": rhs " << i << " true residual " << resid << std::endl;
  }
  std::cout << GridLogMessage << label << ": SUMMARY wall clock " << (t1-t0)/1.0e6
            << " s for " << nbatch << " rhs, " << (t1-t0)/1.0e6/nbatch
            << " s/rhs, worst true residual " << worst << std::endl;
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int         Ls        = CmdOption<int>        (argc,argv,"--Ls",12);
  RealD       mass      = CmdOption<RealD>      (argc,argv,"--mass",0.026);
  RealD       M5        = CmdOption<RealD>      (argc,argv,"--M5",1.8);
  RealD       b         = CmdOption<RealD>      (argc,argv,"--b",1.5);
  RealD       c         = CmdOption<RealD>      (argc,argv,"--c",0.5);
  int         nbatch    = CmdOption<int>        (argc,argv,"--nbatch",4);
  RealD       tol       = CmdOption<RealD>      (argc,argv,"--tol",1.0e-8);
  std::string config    = CmdOption<std::string>(argc,argv,"--config",std::string(""));
  bool        unsplit   = !GridCmdOptionExists(argv,argv+argc,"--nounsplit");
  int         repeat    = CmdOption<int>        (argc,argv,"--repeat",1);

  std::cout << GridLogMessage << "Mobius Ls " << Ls << " mass " << mass << " M5 " << M5
            << " b " << b << " c " << c << " nbatch " << nbatch << " tol " << tol << std::endl;

  GridCartesian         *UGrid_d   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid_d = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_d);
  GridCartesian         *FGrid_d   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_d);
  GridRedBlackCartesian *FrbGrid_d = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_d);

  GridCartesian         *UGrid_f   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid_f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_f);
  GridCartesian         *FGrid_f   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_f);
  GridRedBlackCartesian *FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_f);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG RNG4(UGrid_d);
  GridParallelRNG RNG5(FGrid_d);
  RNG4.SeedFixedIntegers(seeds4);
  RNG5.SeedFixedIntegers(seeds5);

  LatticeGaugeFieldD Umu_d(UGrid_d);
  LatticeGaugeFieldF Umu_f(UGrid_f);
  if ( config.size() ) {
    FieldMetaData header;
    NerscIO::readConfiguration(Umu_d,header,config);
  } else {
    std::cout << GridLogMessage << "No --config: hot configuration, timing only" << std::endl;
    SU<Nc>::HotConfiguration(RNG4,Umu_d);
  }
  precisionChange(Umu_f,Umu_d);
  ReportMemory(UGrid_d,"gauge field ready");

  // Antiperiodic in time, as in production
  WilsonImplParams params;
  params.boundary_phases[Nd-1] = -1.0;

  MobiusFermionD Dd(Umu_d,*FGrid_d,*FrbGrid_d,*UGrid_d,*UrbGrid_d,mass,M5,b,c,params);
  MobiusFermionF Df(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5,b,c,params);
  SchurDiagMooeeOperator<MobiusFermionD,FieldD> Linop_d(Dd);
  SchurDiagMooeeOperator<MobiusFermionF,FieldF> Linop_f(Df);

  std::vector<FieldD> src(nbatch,FrbGrid_d);
  std::vector<FieldD> sol(nbatch,FrbGrid_d);
  {
    FieldD tmp(FGrid_d);
    for(int i=0;i<nbatch;i++){
      random(RNG5,tmp);
      pickCheckerboard(Odd,src[i],tmp);
    }
  }

  ReportMemory(UGrid_d,"operators and sources ready");

  MixedPrecisionConjugateGradientBatched<FieldD,FieldF> mCG(tol,10000,50,10000,FrbGrid_f,Linop_f,Linop_d);

  Coordinate split     = mCG.BatchedSplit;
  bool       splitnode = mCG.BatchedSplitNode;

  if ( unsplit ) {
    mCG.BatchedSplit     = Coordinate();
    mCG.BatchedSplitNode = false;
    SolveAndReport("UNSPLIT",mCG,Linop_d,src,sol);
    ReportMemory(UGrid_d,"after unsplit solve");
  }

  mCG.BatchedSplit     = split;
  mCG.BatchedSplitNode = splitnode;
  // Repeated split solves expose allocations not released between calls
  for(int r=0;r<repeat;r++){
    SolveAndReport("SPLIT",mCG,Linop_d,src,sol);
    ReportMemory(UGrid_d,"after split solve "+std::to_string(r));
  }

  Grid_finalize();
}
