    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_poisson_fft.cc

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
#include <Grid/lattice/Lattice_slice_gpu.h>

using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int N=16;
  std::vector<int> latt_size  ({N,N,N,N});
  std::vector<int> simd_layout({vComplexD::Nsimd(),1,1,1});
  std::vector<int> mpi_layout ({1,1,1,1});

  GridCartesian         GRID(latt_size,simd_layout,mpi_layout);

  LatticeComplexD  rn(&GRID);

  GridParallelRNG RNG(&GRID);
  RNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));  
  gaussian(RNG,rn);

  std::vector<TComplex> reduced_ref;
  std::vector<TComplex> reduced_gpu;
  for(int d=0;d<4;d++){
    {
      RealD t=-usecond();
      sliceSum(rn,reduced_ref,d);
      t+=usecond();
      std::cout << " sliceSum took "<<t<<" usecs"<<std::endl;
    }
    {
      RealD t=-usecond();
      sliceSumGpu(rn,reduced_gpu,d);
      t+=usecond();
      std::cout << " sliceSumGpu took "<<t<<" usecs"<<std::endl;
    }
    for(int t=0;t<reduced_ref.size();t++){
      std::cout << t<<" ref "<< reduced_ref[t] <<" opt " << reduced_gpu[t] << " diff "<<reduced_ref[t]-reduced_gpu[t]<<std::endl;
      TComplex diff = reduced_ref[t]-reduced_gpu[t];
      assert(abs(TensorRemove(diff)) < 1e-8 );
    }
  }
  Grid_finalize();
}
