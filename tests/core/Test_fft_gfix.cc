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
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  std::vector<int> seeds({1,2,3,4});

  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout( { vComplex::Nsimd(),1,1,1});
  std::vector<int> mpi_layout  = GridDefaultMpi();

  int vol = 1;
  for(int d=0;d<latt_size.size();d++){
    vol = vol * latt_size[d];
  }

  GridCartesian         GRID(latt_size,simd_layout,mpi_layout);
  GridSerialRNG          sRNG;  sRNG.SeedFixedIntegers(seeds); // naughty seeding
  GridParallelRNG          pRNG(&GRID);   pRNG.SeedFixedIntegers(seeds);

  FFT theFFT(&GRID);

  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  std::cout<< "*****************************************************************" <<std::endl;
  std::cout<< "* Testing we can gauge fix steep descent a RGT of Unit gauge    *" <<std::endl;
  std::cout<< "*****************************************************************" <<std::endl;

  //  int coulomb_dir = -1;
  int coulomb_dir = Nd-1;

  LatticeGaugeField   Umu(&GRID);
  LatticeGaugeField   Urnd(&GRID);
  LatticeGaugeField   Uorg(&GRID);
  LatticeGaugeField   Utmp(&GRID);
  LatticeColourMatrix   g(&GRID); // Gauge xform

  LatticeColourMatrix   xform1(&GRID); // Gauge xform
  LatticeColourMatrix   xform2(&GRID); // Gauge xform
  LatticeColourMatrix   xform3(&GRID); // Gauge xform
  
  SU3::ColdConfiguration(pRNG,Umu); // Unit gauge
  Uorg=Umu;
  Urnd=Umu;

  SU3::RandomGaugeTransform(pRNG,Urnd,g); // Unit gauge

  Real plaq=WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
  std::cout << " Initial plaquette "<<plaq << std::endl;

  Real alpha=0.1;

  Umu = Urnd;
  FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(Umu,xform1,alpha,10000,1.0e-12, 1.0e-12,false);

  // Check the gauge xform matrices
  Utmp=Urnd;
  SU<Nc>::GaugeTransform(Utmp,xform1);
  Utmp = Utmp - Umu;
  std::cout << " Norm Difference of xformed gauge "<< norm2(Utmp) << std::endl;
  

  plaq=WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
  std::cout << " Final plaquette "<<plaq << std::endl;

  Uorg = Uorg - Umu;
  std::cout << " Norm Difference "<< norm2(Uorg) << std::endl;
  std::cout << " Norm "<< norm2(Umu) << std::endl;


  std::cout<< "*****************************************************************" <<std::endl;
  std::cout<< "* Testing Fourier accelerated fixing                            *" <<std::endl;
  std::cout<< "*****************************************************************" <<std::endl;
  Umu=Urnd;
  FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(Umu,xform2,alpha,10000,1.0e-12, 1.0e-12,true);

  Utmp=Urnd;
  SU<Nc>::GaugeTransform(Utmp,xform2);
  Utmp = Utmp - Umu;
  std::cout << " Norm Difference of xformed gauge "<< norm2(Utmp) << std::endl;


  plaq=WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
  std::cout << " Final plaquette "<<plaq << std::endl;

  std::cout<< "*****************************************************************" <<std::endl;
  std::cout<< "* Testing non-unit configuration                                *" <<std::endl;
  std::cout<< "*****************************************************************" <<std::endl;

  SU3::HotConfiguration(pRNG,Umu); // Unit gauge

  plaq=WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
  std::cout << " Initial plaquette "<<plaq << std::endl;

  FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(Umu,alpha,10000,1.0e-12, 1.0e-12,true);

  plaq=WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
  std::cout << " Final plaquette "<<plaq << std::endl;

  std::cout<< "*****************************************************************" <<std::endl;
  std::cout<< "* Testing Fourier accelerated fixing to coulomb gauge           *" <<std::endl;
  std::cout<< "*****************************************************************" <<std::endl;

  Umu=Urnd;
  SU3::HotConfiguration(pRNG,Umu); // Unit gauge

  plaq=WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
  std::cout << " Initial plaquette "<<plaq << std::endl;

  FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(Umu,xform3,alpha,10000,1.0e-12, 1.0e-12,true,coulomb_dir);

  std::cout << Umu<<std::endl;

  plaq=WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
  std::cout << " Final plaquette "<<plaq << std::endl;

  Grid_finalize();
}
