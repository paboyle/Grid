    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_optimized_staple_gaugebc.cc

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@bnl.gov>
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
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>

using namespace std;
using namespace Grid;
 
int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate latt_size  = GridDefaultLatt();
  Coordinate simd_layout= GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate mpi_layout = GridDefaultMpi();
  std::cout << " mpi "<<mpi_layout<<std::endl;
  std::cout << " simd "<<simd_layout<<std::endl;
  std::cout << " latt "<<latt_size<<std::endl;
  GridCartesian GRID(latt_size,simd_layout,mpi_layout);

  GridParallelRNG   pRNG(&GRID);
  pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
  LatticeGaugeField U(&GRID);

  SU<Nc>::HotConfiguration(pRNG,U);

  //#define PRD
#ifdef PRD
  typedef PeriodicGimplD Gimpl;
#else
  typedef ConjugateGimplD Gimpl;
  std::vector<int> conj_dirs(Nd,0); conj_dirs[0]=1; conj_dirs[3]=1;
  Gimpl::setDirections(conj_dirs);
#endif

  typedef typename WilsonLoops<Gimpl>::GaugeMat GaugeMat;
  typedef typename WilsonLoops<Gimpl>::GaugeLorentz GaugeLorentz;

  int count = 0;
  double torig=0, topt=0;
     
  std::vector<GaugeMat> Umu(Nd,&GRID), U2(Nd,&GRID);
  for(int mu=0;mu<Nd;mu++){
    Umu[mu] = PeekIndex<LorentzIndex>(U,mu);
    WilsonLoops<Gimpl>::RectStapleDouble(U2[mu], Umu[mu], mu);
  }

  std::cout << GridLogMessage << "Checking optimized vs unoptimized RectStaple" << std::endl;
  for(int mu=0;mu<Nd;mu++){
    GaugeMat staple_orig(&GRID), staple_opt(&GRID), staple_U2(&GRID);
    double t0 = usecond();
    WilsonLoops<Gimpl>::RectStapleUnoptimised(staple_orig,U,mu);
    double t1 = usecond();
    WilsonLoops<Gimpl>::RectStapleOptimised(staple_opt, U2, Umu, mu);
    double t2 = usecond();
    torig += t1-t0;  topt += t2-t1;
    ++count;
    
    GaugeMat diff = staple_orig - staple_opt;
    double n = norm2(diff);
    std::cout << GridLogMessage << mu << " " << n << std::endl;
    assert(n<1e-10);
  }
  std::cout << GridLogMessage << "RectStaple timings orig: " << torig/1000/count << "ms,  optimized: " << topt/1000/count << "ms" << std::endl;
  
  Grid_finalize();
}
