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
  LatticeComplexD       C(grid);
  LatticeComplexD       coor(grid);

  ComplexD ci(0.0,1.0);
  Coordinate p({1,2,3,4});

  Coordinate latt_size   = grid->_fdimensions;
  std::cout<<"*************************************************"<<std::endl;
  std::cout<<" Benchmarking FFT of "<<name<<" on plane wave    "<<std::endl;
  std::cout<<"*************************************************"<<std::endl;
  C=Zero();
  for(int mu=0;mu<4;mu++){
    RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
    LatticeCoordinate(coor,mu);
    C = C + (TwoPiL * p[mu]) * coor;
  }
  C = exp(C*ci);

  LatticeObject S(grid);
  LatticeObject Stilde(grid);

  S=Zero();
  S = S+C;

  FFT theFFT(grid);

  Stilde=S;
  std::cout << " norm2(s) "<<norm2(Stilde)<<std::endl;
  double tt= -usecond();
  theFFT.FFT_dim(Stilde,Stilde,0,FFT::forward); std::cout << theFFT.MFlops()<<" mflops "<<norm2(Stilde)<<std::endl;
  theFFT.FFT_dim(Stilde,Stilde,1,FFT::forward); std::cout << theFFT.MFlops()<<" mflops "<<norm2(Stilde)<<std::endl;
  theFFT.FFT_dim(Stilde,Stilde,2,FFT::forward); std::cout << theFFT.MFlops()<<" mflops "<<norm2(Stilde)<<std::endl;
  theFFT.FFT_dim(Stilde,Stilde,3,FFT::forward); std::cout << theFFT.MFlops()<<" mflops "<<norm2(Stilde)<<std::endl;
  tt+= usecond();

  std::cout<<"*************************************************"<<std::endl;
  std::cout<<" FFT of "<<latt_size <<" "<<name<<" took "<<tt/1.e6<<" s"<<std::endl;
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
  bench<LatticePropagatorD>(&GRID,std::string("LatticePropagatorD"));
  
  Grid_finalize();
}
