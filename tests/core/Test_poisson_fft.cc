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
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  int N=128;
  int N2=64;
  int W=16;
  int D=8;
  std::vector<int> latt_size  ({N,N});
  std::vector<int> simd_layout({vComplexD::Nsimd(),1});
  std::vector<int> mpi_layout ({1,1});

  int vol = 1;
  int nd  = latt_size.size();
  for(int d=0;d<nd;d++){
    vol = vol * latt_size[d];
  }

  GridCartesian         GRID(latt_size,simd_layout,mpi_layout);

  LatticeComplexD     pos(&GRID);
  LatticeComplexD      zz(&GRID);
  LatticeComplexD     neg(&GRID);
  LatticeInteger     coor(&GRID);
  LatticeComplexD  Charge(&GRID);
  LatticeComplexD  ChargeTilde(&GRID);
  LatticeComplexD       V(&GRID);
  LatticeComplexD       Vtilde(&GRID);

  pos = ComplexD(1.0,0.0);
  neg = -pos;
  zz  = ComplexD(0.0,0.0);

  Charge=zero;
  
  // Parallel plate capacitor
  {
    int mu=0;
    LatticeCoordinate(coor,mu);

    Charge=where(coor==Integer(N2-D),pos,zz);
    Charge=where(coor==Integer(N2+D),neg,Charge);
  }

  {
    int mu=1;
    LatticeCoordinate(coor,mu);
    Charge=where(coor<Integer(N2-W),zz,Charge);
    Charge=where(coor>Integer(N2+W),zz,Charge);
  }

  //  std::cout << Charge <<std::endl;

  std::vector<LatticeComplexD> k(4,&GRID);
  LatticeComplexD ksq(&GRID);

  ksq=zero;
  for(int mu=0;mu<nd;mu++) {
    
    Integer L=latt_size[mu];
    
    LatticeCoordinate(coor,mu);
    LatticeCoordinate(k[mu],mu);

    k[mu] = where ( coor > (L/2), k[mu]-L, k[mu]);

    //    std::cout << k[mu]<<std::endl;

    RealD TwoPiL =  M_PI * 2.0/ L;

    k[mu] = TwoPiL * k[mu];      

    ksq = ksq + k[mu]*k[mu];

  }
  
  // D^2 V   = - rho
  // ksq Vtilde = rhoTilde
  // Vtilde = rhoTilde/Ksq
  // Fix zero of potential : Vtilde(0) = 0;
  std::vector<int> zero_mode(nd,0);
  TComplexD Tone = ComplexD(1.0,0.0);
  pokeSite(Tone,ksq,zero_mode); 

  //  std::cout << "Charge\n" << Charge <<std::endl;

  FFT theFFT(&GRID);
  theFFT.FFT_all_dim(ChargeTilde,Charge,FFT::forward); 
  //  std::cout << "Rhotilde\n" << ChargeTilde <<std::endl;

  Vtilde = ChargeTilde / ksq;
  //  std::cout << "Vtilde\n" << Vtilde <<std::endl;

  TComplexD Tzero = ComplexD(0.0,0.0);
  pokeSite(Tzero,Vtilde,zero_mode);

  theFFT.FFT_all_dim(V,Vtilde,FFT::backward); 

  std::cout << "V\n" << V <<std::endl;

  Grid_finalize();
}
