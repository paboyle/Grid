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
  std::vector<int> simd_layout( { vComplexD::Nsimd(),1,1,1});
  std::vector<int> mpi_layout  = GridDefaultMpi();

  int vol = 1;
  for(int d=0;d<latt_size.size();d++){
    vol = vol * latt_size[d];
  }
  GridCartesian         GRID(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian RBGRID(latt_size,simd_layout,mpi_layout);

  LatticeComplexD     one(&GRID);
  LatticeComplexD      zz(&GRID);
  LatticeComplexD       C(&GRID);
  LatticeComplexD  Ctilde(&GRID);
  LatticeComplexD  Cref  (&GRID);
  LatticeComplexD  Csav  (&GRID);
  LatticeComplexD    coor(&GRID);

  LatticeSpinMatrixD    S(&GRID);
  LatticeSpinMatrixD    Stilde(&GRID);
  
  std::vector<int> p({1,2,3,2});

  one = ComplexD(1.0,0.0);
  zz  = ComplexD(0.0,0.0);

  ComplexD ci(0.0,1.0);

  C=zero;
  for(int mu=0;mu<4;mu++){
    RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
    LatticeCoordinate(coor,mu);
    C = C - (TwoPiL * p[mu]) * coor;
  }

  C = exp(C*ci);
  Csav = C;
  S=zero;
  S = S+C;

  FFT theFFT(&GRID);

  theFFT.FFT_dim(Ctilde,C,0,FFT::forward);  C=Ctilde;
  theFFT.FFT_dim(Ctilde,C,1,FFT::forward);  C=Ctilde; std::cout << theFFT.MFlops()<<" Mflops "<<std::endl;
  theFFT.FFT_dim(Ctilde,C,2,FFT::forward);  C=Ctilde;
  theFFT.FFT_dim(Ctilde,C,3,FFT::forward);  


  //  C=zero;
  //  Ctilde = where(abs(Ctilde)<1.0e-10,C,Ctilde);
  TComplexD cVol;
  cVol()()() = vol;

  Cref=zero;
  pokeSite(cVol,Cref,p);
  Cref=Cref-Ctilde;
  std::cout << "diff scalar "<<norm2(Cref) << std::endl;

  C=Csav;
  theFFT.FFT_all_dim(Ctilde,C,FFT::forward);
  theFFT.FFT_all_dim(Cref,Ctilde,FFT::backward); 

  std::cout << norm2(C) << " " << norm2(Ctilde) << " " << norm2(Cref)<< " vol " << vol<< std::endl;

  Cref= Cref - C;
  std::cout << " invertible check " << norm2(Cref)<<std::endl;

  theFFT.FFT_dim(Stilde,S,0,FFT::forward);  S=Stilde;
  theFFT.FFT_dim(Stilde,S,1,FFT::forward);  S=Stilde;std::cout << theFFT.MFlops()<<" mflops "<<std::endl;
  theFFT.FFT_dim(Stilde,S,2,FFT::forward);  S=Stilde;
  theFFT.FFT_dim(Stilde,S,3,FFT::forward);

  SpinMatrixD Sp; 
  Sp = zero; Sp = Sp+cVol;

  S=zero;
  pokeSite(Sp,S,p);

  S= S-Stilde;
  std::cout << "diff FT[SpinMat] "<<norm2(S) << std::endl;

  /*
   */
  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&GRID);
  pRNG.SeedFixedIntegers(seeds);

  LatticeGaugeFieldD Umu(&GRID);
  LatticeFermionD    src(&GRID); gaussian(pRNG,src);
  LatticeFermionD    tmp(&GRID);
  LatticeFermionD    ref(&GRID);

  SU3::ColdConfiguration(pRNG,Umu); // Unit gauge
  
  RealD mass=0.1;
  WilsonFermionD Dw(Umu,GRID,RBGRID,mass);

  Dw.M(src,tmp);

  std::cout << " src = " <<norm2(src)<<std::endl;
  std::cout << " tmp = " <<norm2(tmp)<<std::endl;

  Dw.FreePropagator(tmp,ref);

  std::cout << " ref = " <<norm2(ref)<<std::endl;
  
  ref = ref - src;
  
  std::cout << " ref-src = " <<norm2(ref)<<std::endl;

  Grid_finalize();
}
