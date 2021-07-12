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

template<typename Gimpl>
void run(double alpha, bool do_fft_gfix){
  std::vector<int> seeds({1,2,3,4});
  int threads = GridThread::GetThreads();

  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout( { vComplex::Nsimd(),1,1,1});
  Coordinate mpi_layout  = GridDefaultMpi();

  int vol = 1;
  for(int d=0;d<latt_size.size();d++){
    vol = vol * latt_size[d];
  }

  GridCartesian         GRID(latt_size,simd_layout,mpi_layout);
  GridSerialRNG          sRNG;  sRNG.SeedFixedIntegers(seeds); // naughty seeding
  GridParallelRNG          pRNG(&GRID);   pRNG.SeedFixedIntegers(seeds);

  FFT theFFT(&GRID);

  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
  std::cout<<GridLogMessage << "Using alpha=" << alpha << std::endl;

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

  //#########################################################################################

  std::cout<< "*********************************************************************************************************" <<std::endl;
  std::cout<< "* Testing steepest descent fixing to Landau gauge with randomly transformed unit gauge configuration    *" <<std::endl;
  std::cout<< "*********************************************************************************************************" <<std::endl;
  
  SU<Nc>::ColdConfiguration(pRNG,Umu); // Unit gauge
  Uorg=Umu;

  Real init_plaq=WilsonLoops<Gimpl>::avgPlaquette(Umu);
  std::cout << " Initial plaquette "<< init_plaq << std::endl;

  //Apply a random gauge transformation to the unit gauge config
  Urnd=Umu;
  SU<Nc>::RandomGaugeTransform<Gimpl>(pRNG,Urnd,g);

  //Gauge fix the randomly transformed field 
  Umu = Urnd;
  FourierAcceleratedGaugeFixer<Gimpl>::SteepestDescentGaugeFix(Umu,xform1,alpha,10000,1.0e-12, 1.0e-12,false);

  // Check the gauge xform matrices
  Utmp=Urnd;
  SU<Nc>::GaugeTransform<Gimpl>(Utmp,xform1);
  Utmp = Utmp - Umu;
  std::cout << " Check the output gauge transformation matrices applied to the original field produce the xformed field "<< norm2(Utmp) << " (expect 0)" << std::endl;
  

  Real plaq=WilsonLoops<Gimpl>::avgPlaquette(Umu);
  std::cout << " Final plaquette "<<plaq << " diff " << plaq - init_plaq << " (expect 0)" << std::endl;

  Uorg = Uorg - Umu;
  std::cout << " Norm difference between a unit gauge configuration and the gauge fixed configuration "<< norm2(Uorg) << " (expect 0)" << std::endl;
  std::cout << " Norm of gauge fixed configuration "<< norm2(Umu) << std::endl;

  //#########################################################################################
  if(do_fft_gfix){
    std::cout<< "*************************************************************************************" <<std::endl;
    std::cout<< "* Testing Fourier accelerated fixing to Landau gauge with unit gauge configuration  *" <<std::endl;
    std::cout<< "*************************************************************************************" <<std::endl;
    Umu=Urnd;
    FourierAcceleratedGaugeFixer<Gimpl>::SteepestDescentGaugeFix(Umu,xform2,alpha,10000,1.0e-12, 1.0e-12,true);

    Utmp=Urnd;
    SU<Nc>::GaugeTransform<Gimpl>(Utmp,xform2);
    Utmp = Utmp - Umu;
    std::cout << " Check the output gauge transformation matrices applied to the original field produce the xformed field "<< norm2(Utmp) << " (expect 0)" << std::endl;


    plaq=WilsonLoops<Gimpl>::avgPlaquette(Umu);
    std::cout << " Final plaquette "<<plaq << " diff " << plaq - init_plaq << " (expect 0)" << std::endl;
  }
  //#########################################################################################

  std::cout<< "******************************************************************************************" <<std::endl;
  std::cout<< "* Testing steepest descent fixing to Landau gauge with random configuration             **" <<std::endl;
  std::cout<< "******************************************************************************************" <<std::endl;

  SU<Nc>::HotConfiguration(pRNG,Umu);

  init_plaq=WilsonLoops<Gimpl>::avgPlaquette(Umu);
  std::cout << " Initial plaquette "<< init_plaq << std::endl;

  FourierAcceleratedGaugeFixer<Gimpl>::SteepestDescentGaugeFix(Umu,alpha,10000,1.0e-12, 1.0e-12,false);

  plaq=WilsonLoops<Gimpl>::avgPlaquette(Umu);
  std::cout << " Final plaquette "<<plaq << " diff " << plaq - init_plaq << " (expect 0)" << std::endl;

  //#########################################################################################
  if(do_fft_gfix){
    std::cout<< "******************************************************************************************" <<std::endl;
    std::cout<< "* Testing Fourier accelerated fixing to Landau gauge with random configuration          **" <<std::endl;
    std::cout<< "******************************************************************************************" <<std::endl;

    SU<Nc>::HotConfiguration(pRNG,Umu);

    init_plaq=WilsonLoops<Gimpl>::avgPlaquette(Umu);
    std::cout << " Initial plaquette "<< init_plaq << std::endl;

    FourierAcceleratedGaugeFixer<Gimpl>::SteepestDescentGaugeFix(Umu,alpha,10000,1.0e-12, 1.0e-12,true);

    plaq=WilsonLoops<Gimpl>::avgPlaquette(Umu);
    std::cout << " Final plaquette "<<plaq << " diff " << plaq - init_plaq << " (expect 0)" << std::endl;
  }
  //#########################################################################################
  
  std::cout<< "*******************************************************************************************" <<std::endl;
  std::cout<< "* Testing steepest descent fixing to coulomb gauge with random configuration           *" <<std::endl;
  std::cout<< "*******************************************************************************************" <<std::endl;

  Umu=Urnd;
  SU<Nc>::HotConfiguration(pRNG,Umu);

  init_plaq=WilsonLoops<Gimpl>::avgPlaquette(Umu);
  std::cout << " Initial plaquette "<< init_plaq << std::endl;

  FourierAcceleratedGaugeFixer<Gimpl>::SteepestDescentGaugeFix(Umu,xform3,alpha,10000,1.0e-12, 1.0e-12,false,coulomb_dir);

  plaq=WilsonLoops<Gimpl>::avgPlaquette(Umu);
  std::cout << " Final plaquette "<<plaq << " diff " << plaq - init_plaq << " (expect 0)" << std::endl;


  //#########################################################################################
  if(do_fft_gfix){
    std::cout<< "*******************************************************************************************" <<std::endl;
    std::cout<< "* Testing Fourier accelerated fixing to coulomb gauge with random configuration           *" <<std::endl;
    std::cout<< "*******************************************************************************************" <<std::endl;

    Umu=Urnd;
    SU<Nc>::HotConfiguration(pRNG,Umu);

    init_plaq=WilsonLoops<Gimpl>::avgPlaquette(Umu);
    std::cout << " Initial plaquette "<< init_plaq << std::endl;

    FourierAcceleratedGaugeFixer<Gimpl>::SteepestDescentGaugeFix(Umu,xform3,alpha,10000,1.0e-12, 1.0e-12,true,coulomb_dir);

    plaq=WilsonLoops<Gimpl>::avgPlaquette(Umu);
    std::cout << " Final plaquette "<<plaq << " diff " << plaq - init_plaq << " (expect 0)" << std::endl;
  }
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  double alpha=0.1; //step size
  std::string gimpl = "periodic";
  bool do_fft_gfix = true; //test fourier transformed gfix as well as steepest descent
  for(int i=1;i<argc;i++){
    std::string sarg(argv[i]);
    if(sarg == "--gimpl"){
      assert(i<argc-1 && "--gimpl option requires an argument");
      gimpl = argv[i+1];
      if(gimpl != "periodic" && gimpl != "conjugate")
	assert(0 && "Invalid gimpl");
    }else if(sarg == "--no-fft-gfix"){
      std::cout << "Not doing the Fourier accelerated gauge fixing tests" << std::endl;
      do_fft_gfix = false;
    }else if(sarg == "--alpha"){
      assert(i<argc-1 && "--alpha option requires an argument");
      std::istringstream ss(argv[i+1]); ss >> alpha;
    }
  }


  if(gimpl == "periodic"){
    std::cout << GridLogMessage << "Using periodic boundary condition" << std::endl;
    run<PeriodicGimplR>(alpha, do_fft_gfix);
  }else{
    std::vector<int> conjdirs = {1,1,0,0}; //test with 2 conjugate dirs and 2 not
    std::cout << GridLogMessage << "Using complex conjugate boundary conditions in dimensions ";
    for(int i=0;i<Nd;i++)
      if(conjdirs[i])
	std::cout << i << " ";   
    std::cout << std::endl;

    ConjugateGimplR::setDirections(conjdirs);
    run<ConjugateGimplR>(alpha, do_fft_gfix);
  }
  
  Grid_finalize();
}
