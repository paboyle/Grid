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

using namespace Grid;
 ;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  int N=16;
  
  std::vector<int> latt_size  ({N,N,N,N});
  std::vector<int> simd_layout({vComplexD::Nsimd(),1,1,1});
  std::vector<int> mpi_layout ({1,1,1,1});

  int vol = 1;
  int nd  = latt_size.size();
  for(int d=0;d<nd;d++){
    vol = vol * latt_size[d];
  }

  GridCartesian         GRID(latt_size,simd_layout,mpi_layout);

  LatticeComplexD      zz(&GRID);
  LatticeInteger     coor(&GRID);
  LatticeComplexD  rn(&GRID);
  LatticeComplexD  sl(&GRID);

  zz  = ComplexD(0.0,0.0);

  GridParallelRNG RNG(&GRID);
  RNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));  
  gaussian(RNG,rn);

  RealD nn=norm2(rn);
  for(int mu=0;mu<nd;mu++){
    RealD ns=0.0;
    for(int t=0;t<latt_size[mu];t++){
      LatticeCoordinate(coor,mu);
      sl=where(coor==Integer(t),rn,zz);
      //      std::cout <<GridLogMessage<< " sl " << sl<<std::endl;
      std::cout <<GridLogMessage<<" slice "<<t<<" " << norm2(sl)<<std::endl;
      ns=ns+norm2(sl);
    }
    std::cout <<GridLogMessage <<" sliceNorm" <<mu<<" "<< nn <<" "<<ns<<" " << nn-ns<<std::endl;
  }

  Grid_finalize();
}
