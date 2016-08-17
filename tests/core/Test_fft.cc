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

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout( { vComplexD::Nsimd(),1,1,1});
  std::vector<int> mpi_layout  = GridDefaultMpi();

  GridCartesian        Fine(latt_size,simd_layout,mpi_layout);

  LatticeComplexD     one(&Fine);
  LatticeComplexD      zz(&Fine);
  LatticeComplexD       C(&Fine);
  LatticeComplexD  Ctilde(&Fine);
  LatticeComplexD    coor(&Fine);
  
  std::vector<RealD> p({1.0,2.0,3.0,2.0});

  one = ComplexD(1.0,0.0);
  zz  = ComplexD(0.0,0.0);

  ComplexD ci(0.0,1.0);

  C=zero;
  for(int mu=0;mu<4;mu++){
    RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
    LatticeCoordinate(coor,mu);
    C = C - TwoPiL * p[mu] * coor;
  }

  std::cout << GridLogMessage<< " C " << C<<std::endl;

  C = C*ci;
  std::cout << GridLogMessage<< " C " << C<<std::endl;

  C = exp(C);
  std::cout << GridLogMessage<< " C " << C<<std::endl;

  FFT theFFT(&Fine);
  theFFT.FFT_dim(Ctilde,C,0,FFT::forward);
  std::cout << GridLogMessage<< "FT[C] " << Ctilde<<std::endl;

  C=Ctilde;
  theFFT.FFT_dim(Ctilde,C,1,FFT::forward);
  std::cout << GridLogMessage<< "FT[C] " << Ctilde<<std::endl;
  C=Ctilde;
  theFFT.FFT_dim(Ctilde,C,2,FFT::forward);
  std::cout << GridLogMessage<< "FT[C] " << Ctilde<<std::endl;
  C=Ctilde;
  theFFT.FFT_dim(Ctilde,C,3,FFT::forward);
  std::cout << GridLogMessage<< "FT[C] " << Ctilde<<std::endl;

  Grid_finalize();
}
