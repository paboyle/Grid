/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/debug/Test_ring_allreduce.cc

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
// Gate for RingAllReduce / CartesianRingAllReduce (P2P-only vector
// all-reduce) against GlobalSumVector (MPI_Allreduce):
//   T1 flat ring      == GlobalSumVector, RealD/ComplexD/ComplexF/RealF,
//                        n in {1, P-1, P, P+1, 1000003, 2^20} (n<P, n%P!=0)
//   T2 cartesian ring == GlobalSumVector, same sweep
//   T3 bitwise repeatable (deterministic order)
//   T4 timing at 16 MB, both rings vs GlobalSumVector
//
//   mpirun -n 4 ./Test_ring_allreduce --grid 16.16.16.32 --mpi 1.1.2.2
// (2D process grid so the cartesian variant exercises more than one ring)
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

template<class T> double Mag(const T &x){ return std::fabs((double)x); }
template<> double Mag<ComplexD>(const ComplexD &z){ return std::sqrt(z.real()*z.real()+z.imag()*z.imag()); }
template<> double Mag<ComplexF>(const ComplexF &z){ return std::sqrt((double)z.real()*z.real()+(double)z.imag()*z.imag()); }

template<class T> T Fill(uint64_t i, int rank){ return T( 0.5*std::sin(0.01*i + 0.7*rank) + 1.0e-3*rank ); }
template<> ComplexD Fill<ComplexD>(uint64_t i,int rank){ return ComplexD(0.5*std::sin(0.01*i+0.7*rank), 0.3*std::cos(0.02*i-0.1*rank)); }
template<> ComplexF Fill<ComplexF>(uint64_t i,int rank){ return ComplexF(0.5*std::sin(0.01*i+0.7*rank), 0.3*std::cos(0.02*i-0.1*rank)); }

template<class T>
void Check(const std::string &tname, GridCartesian *grid, double tol)
{
  int P  = grid->ProcessorCount();
  int me = grid->ThisRank();
  std::vector<uint64_t> sizes({1,(uint64_t)std::max(P-1,1),(uint64_t)P,(uint64_t)P+1,1000003,1<<20});
  for(uint64_t n : sizes){
    std::vector<T> h(n); for(uint64_t i=0;i<n;i++) h[i]=Fill<T>(i,me);
    std::vector<T> ref(h); grid->GlobalSumVector(&ref[0],(int)n);

    deviceVector<T> d(n);
    for(int variant=0;variant<2;variant++){
      acceleratorCopyToDevice(&h[0],&d[0],n*sizeof(T));
      if(variant==0) RingAllReduce(grid,&d[0],n);
      else           CartesianRingAllReduce(grid,&d[0],n);
      std::vector<T> out(n); acceleratorCopyFromDevice(&d[0],&out[0],n*sizeof(T));
      double worst=0.0, scale=0.0;
      for(uint64_t i=0;i<n;i++){ worst=std::max(worst,Mag<T>(out[i]-ref[i])); scale=std::max(scale,Mag<T>(ref[i])); }
      RealD w=worst; grid->GlobalMax(w);
      std::ostringstream os; os<<"n="<<n<<" worst abs "<<w<<" (scale "<<scale<<")";
      Report(std::string(variant?"T2 cartesian ":"T1 flat      ")+tname, w<tol*std::max(scale,1.0), os.str());
    }
  }
  // T3 determinism
  {
    uint64_t n=1<<18;
    std::vector<T> h(n); for(uint64_t i=0;i<n;i++) h[i]=Fill<T>(i,me);
    deviceVector<T> d(n); std::vector<T> a(n),b(n);
    acceleratorCopyToDevice(&h[0],&d[0],n*sizeof(T)); RingAllReduce(grid,&d[0],n); acceleratorCopyFromDevice(&d[0],&a[0],n*sizeof(T));
    acceleratorCopyToDevice(&h[0],&d[0],n*sizeof(T)); RingAllReduce(grid,&d[0],n); acceleratorCopyFromDevice(&d[0],&b[0],n*sizeof(T));
    RealD diff = (memcmp(&a[0],&b[0],n*sizeof(T))!=0) ? 1.0 : 0.0;
    grid->GlobalSum(diff);
    Report("T3 bitwise repeat "+tname, diff==0.0);
  }
}

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);
  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                             GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  std::cout << GridLogMessage << "Ring allreduce test: P=" << grid->ProcessorCount()
            << " processor grid " << grid->_processors << std::endl;

  Check<RealD>   ("RealD   ", grid, 1.0e-13);
  Check<ComplexD>("ComplexD", grid, 1.0e-13);
  Check<RealF>   ("RealF   ", grid, 1.0e-5);
  Check<ComplexF>("ComplexF", grid, 1.0e-5);

  // T4 timing at 16 MB of ComplexF (the dense-apply size at 12 RHS is 13.3 MB)
  {
    uint64_t n = 2*1024*1024;
    deviceVector<ComplexF> d(n); std::vector<ComplexF> h(n,ComplexF(1.0,0.0));
    for(int rep=0;rep<2;rep++){
      acceleratorCopyToDevice(&h[0],&d[0],n*sizeof(ComplexF));
      double t0=usecond(); RingAllReduce(grid,&d[0],n);          double t1=usecond();
      acceleratorCopyToDevice(&h[0],&d[0],n*sizeof(ComplexF));
      double t2=usecond(); CartesianRingAllReduce(grid,&d[0],n); double t3=usecond();
      double t4=usecond(); grid->GlobalSumVector(&h[0],(int)n);   double t5=usecond();
      if(rep) std::cout << GridLogMessage << "T4 timing 16 MB ComplexF: flat ring " << (t1-t0)/1000.
                        << " ms, cartesian ring " << (t3-t2)/1000. << " ms, GlobalSumVector(host) " << (t5-t4)/1000. << " ms" << std::endl;
    }
  }

  std::cout << GridLogMessage << (failures ? "Test_ring_allreduce: FAILURES" : "Test_ring_allreduce: ALL PASS") << std::endl;
  Grid_finalize();
  return failures ? 1 : 0;
}
