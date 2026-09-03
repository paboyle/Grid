    /*************************************************************************************

    grid` physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_cshift.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
#include <Grid/Grid.h>

using namespace Grid;

template<class LatticeObject>
void bench(GridCartesian *grid, std::string name)
{
  Coordinate latt_size   = grid->_fdimensions;
  std::cout<<"*************************************************"<<std::endl;
  std::cout<<" Benchmarking FFT of "<<name<<" on random field   "<<std::endl;
  std::cout<<"*************************************************"<<std::endl;

  // Random source: works for every Lattice type (incl. LatticeFermionD, where the
  // scalar plane-wave broadcast S=S+C does not compile).  Content is irrelevant to
  // FFT timing and plan cost; Parseval still holds -- a full forward FFT scales
  // norm2 by exactly vol for ANY input.
  GridParallelRNG RNG(grid);  RNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));
  LatticeObject S(grid);  gaussian(RNG,S);

  typedef typename LatticeObject::vector_object vobj;
  const int nrep = 10;

  // Correctness + Parseval pass (unplanned), also the WARMUP that absorbs
  // FFTW/allocator/MPI first-touch so the two timed passes below are both warm
  // and the comparison is fair (order confound removed).
  {
    LatticeObject Stilde(grid);  Stilde=S;
    FFT theFFT(grid);
    std::cout << " norm2(s) "<<norm2(Stilde)<<std::endl;
    theFFT.FFT_dim(Stilde,Stilde,0,FFT::forward); std::cout << theFFT.MFlops()<<" mflops "<<norm2(Stilde)<<std::endl;
    theFFT.FFT_dim(Stilde,Stilde,1,FFT::forward); std::cout << theFFT.MFlops()<<" mflops "<<norm2(Stilde)<<std::endl;
    theFFT.FFT_dim(Stilde,Stilde,2,FFT::forward); std::cout << theFFT.MFlops()<<" mflops "<<norm2(Stilde)<<std::endl;
    theFFT.FFT_dim(Stilde,Stilde,3,FFT::forward); std::cout << theFFT.MFlops()<<" mflops "<<norm2(Stilde)<<std::endl;
  }

  // NOTE: norm2() is a GlobalSum (MPI all-reduce over all ranks) -- it must NOT
  // sit inside the timed loops or it dominates the number (that was the 80% gap
  // in the Frontier LatticeFermionD run: 4 norm2 reductions + I/O bracketed by
  // the outer timer, NOT plan creation).  Timed loops below do FFT only.

  double t_unplanned, t_planned, t_build;

  // ---- UNPLANNED: the FFT class builds AND destroys an FFTW plan on EVERY
  //      FFT_dim call -- averaged over nrep so the recurring plan cost shows.
  {
    LatticeObject Stilde(grid);
    double tt= -usecond();
    for(int r=0;r<nrep;r++){
      Stilde=S;
      FFT theFFT(grid);
      for(int mu=0;mu<4;mu++) theFFT.FFT_dim(Stilde,Stilde,mu,FFT::forward);
    }
    tt+= usecond();
    t_unplanned = tt/1.e6/nrep;
  }

  // ---- PLANNED: PlannedFFT builds all plans once in its ctor, reused every
  //      call -- the plan-free steady state (what repeated applies pay).
  {
    LatticeObject Stilde(grid);
    double tc= -usecond();
    PlannedFFT<vobj> theFFT(grid);           // one-time plan build (all dims, fwd+bwd)
    tc+= usecond();
    t_build = tc/1.e6;
    double tt= -usecond();
    for(int r=0;r<nrep;r++){
      Stilde=S;
      for(int mu=0;mu<4;mu++) theFFT.FFT_dim(Stilde,Stilde,mu,FFT::forward);
    }
    tt+= usecond();
    t_planned = tt/1.e6/nrep;
  }

  // The decisive line: (unplanned - planned) is the per-call plan create+destroy
  // cost, since the execute path is identical.  If ~0, PlannedFFT buys nothing.
  std::cout<<" PLANCOST "<<name<<" : unplanned "<<t_unplanned<<" s/call   planned "<<t_planned
           <<" s/call   => plan create+destroy "<<(t_unplanned - t_planned)<<" s/call"
           <<"   (one-time build "<<t_build<<" s)"<<std::endl;
  std::cout<<"*************************************************"<<std::endl;

}


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  int vol = 1;
  for(int d=0;d<latt_size.size();d++){
    vol = vol * latt_size[d];
  }
  GridCartesian         GRID(latt_size,simd_layout,mpi_layout);


  bench<LatticeComplexD>(&GRID,std::string("LatticeComplexD"));
  bench<LatticeColourMatrixD>(&GRID,std::string("LatticeColourMatrixD"));
  bench<LatticeFermionD>(&GRID,std::string("LatticeFermionD"));      // Ncomp=12, the FreePropagator/Fourier-precon path
  bench<LatticePropagatorD>(&GRID,std::string("LatticePropagatorD"));
  
  Grid_finalize();
}
