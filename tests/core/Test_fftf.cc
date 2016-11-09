    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

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
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout( { vComplexF::Nsimd(),1,1,1});
  std::vector<int> mpi_layout  = GridDefaultMpi();

  int vol = 1;
  for(int d=0;d<latt_size.size();d++){
    vol = vol * latt_size[d];
  }
  GridCartesian        Fine(latt_size,simd_layout,mpi_layout);

  LatticeComplexF     one(&Fine);
  LatticeComplexF      zz(&Fine);
  LatticeComplexF       C(&Fine);
  LatticeComplexF  Ctilde(&Fine);
  LatticeComplexF    coor(&Fine);

  LatticeSpinMatrixF    S(&Fine);
  LatticeSpinMatrixF    Stilde(&Fine);
  
  std::vector<int> p({1,2,3,2});

  one = ComplexF(1.0,0.0);
  zz  = ComplexF(0.0,0.0);

  ComplexF ci(0.0,1.0);

  C=zero;
  for(int mu=0;mu<4;mu++){
    RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
    LatticeCoordinate(coor,mu);
    C = C - (TwoPiL * p[mu]) * coor;
  }

  C = exp(C*ci);

  S=zero;
  S = S+C;

  FFT theFFT(&Fine);

  theFFT.FFT_dim(Ctilde,C,0,FFT::forward);  C=Ctilde; std::cout << theFFT.MFlops()<<std::endl;
  theFFT.FFT_dim(Ctilde,C,1,FFT::forward);  C=Ctilde; std::cout << theFFT.MFlops()<<std::endl;
  theFFT.FFT_dim(Ctilde,C,2,FFT::forward);  C=Ctilde; std::cout << theFFT.MFlops()<<std::endl;
  theFFT.FFT_dim(Ctilde,C,3,FFT::forward);  std::cout << theFFT.MFlops()<<std::endl;

  //  C=zero;
  //  Ctilde = where(abs(Ctilde)<1.0e-10,C,Ctilde);
  TComplexF cVol;
  cVol()()() = vol;

  C=zero;
  pokeSite(cVol,C,p);
  C=C-Ctilde;
  std::cout << "diff scalar "<<norm2(C) << std::endl;

  theFFT.FFT_dim(Stilde,S,0,FFT::forward);  S=Stilde;std::cout << theFFT.MFlops()<< " "<<theFFT.USec() <<std::endl;
  theFFT.FFT_dim(Stilde,S,1,FFT::forward);  S=Stilde;std::cout << theFFT.MFlops()<< " "<<theFFT.USec() <<std::endl;
  theFFT.FFT_dim(Stilde,S,2,FFT::forward);  S=Stilde;std::cout << theFFT.MFlops()<< " "<<theFFT.USec() <<std::endl;
  theFFT.FFT_dim(Stilde,S,3,FFT::forward);std::cout << theFFT.MFlops()<<" "<<theFFT.USec() <<std::endl;

  SpinMatrixF Sp; 
  Sp = zero; Sp = Sp+cVol;

  S=zero;
  pokeSite(Sp,S,p);

  S= S-Stilde;
  std::cout << "diff FT[SpinMat] "<<norm2(S) << std::endl;

  Grid_finalize();
}
