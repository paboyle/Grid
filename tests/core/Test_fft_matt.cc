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

Gamma::Algebra Gmu [] = {
  Gamma::Algebra::GammaX,
  Gamma::Algebra::GammaY,
  Gamma::Algebra::GammaZ,
  Gamma::Algebra::GammaT,
  Gamma::Algebra::Gamma5
};

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
  GridRedBlackCartesian RBGRID(&GRID);

  LatticeComplexD    coor(&GRID);
  ComplexD ci(0.0,1.0);

  std::vector<int> seeds({1,2,3,4});
  GridSerialRNG          sRNG;  sRNG.SeedFixedIntegers(seeds); // naughty seeding
  GridParallelRNG          pRNG(&GRID);
  pRNG.SeedFixedIntegers(seeds);

  LatticeGaugeFieldD Umu(&GRID);
  SU<Nc>::ColdConfiguration(pRNG,Umu); // Unit gauge

  ////////////////////////////////////////////////////
  // Wilson test
  ////////////////////////////////////////////////////
  {
    LatticeFermionD    src(&GRID); gaussian(pRNG,src);
    LatticeFermionD    src_p(&GRID);
    LatticeFermionD    tmp(&GRID);
    LatticeFermionD    ref(&GRID);
    LatticeFermionD    result(&GRID);
    
    RealD mass=0.1;
    WilsonFermionD Dw(Umu,GRID,RBGRID,mass);
    
    Dw.M(src,ref);
    std::cout << "Norm src "<<norm2(src)<<std::endl;
    std::cout << "Norm Dw x src "<<norm2(ref)<<std::endl;
    {
      FFT theFFT(&GRID);

      ////////////////
      // operator in Fourier space
      ////////////////
      tmp =ref;
      theFFT.FFT_all_dim(result,tmp,FFT::forward);
      std::cout<<"FFT[ Dw x src ]  "<< norm2(result)<<std::endl;    

      tmp = src;
      theFFT.FFT_all_dim(src_p,tmp,FFT::forward);
      std::cout<<"FFT[ src      ]  "<< norm2(src_p)<<std::endl;
      
      /////////////////////////////////////////////////////////////////
      // work out the predicted FT from Fourier
      /////////////////////////////////////////////////////////////////
      auto FGrid = &GRID;
      LatticeFermionD    Kinetic(FGrid); Kinetic = Zero();
      LatticeComplexD    kmu(FGrid); 
      LatticeInteger     scoor(FGrid); 
      LatticeComplexD    sk (FGrid); sk = Zero();
      LatticeComplexD    sk2(FGrid); sk2= Zero();
      LatticeComplexD    W(FGrid); W= Zero();
      LatticeComplexD    one(FGrid); one =ComplexD(1.0,0.0);
      ComplexD ci(0.0,1.0);
    
      for(int mu=0;mu<Nd;mu++) {
	
	RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];

	LatticeCoordinate(kmu,mu);

	kmu = TwoPiL * kmu;
      
	sk2 = sk2 + 2.0*sin(kmu*0.5)*sin(kmu*0.5);
	sk  = sk  +     sin(kmu)    *sin(kmu); 
      
	// -1/2 Dw ->  1/2 gmu (eip - emip) = i sinp gmu
	Kinetic = Kinetic + sin(kmu)*ci*(Gamma(Gmu[mu])*src_p);
	
      }
    
      W = mass + sk2; 
      Kinetic = Kinetic + W * src_p;
    
      std::cout<<"Momentum space src         "<< norm2(src_p)<<std::endl;
      std::cout<<"Momentum space Dw x src    "<< norm2(Kinetic)<<std::endl;
      std::cout<<"FT[Coordinate space Dw]    "<< norm2(result)<<std::endl;
    
      result = result - Kinetic;
      std::cout<<"diff "<< norm2(result)<<std::endl;
      
    }

    std::cout << " =======================================" <<std::endl;
    std::cout << " Checking FourierFreePropagator x Dw = 1" <<std::endl;
    std::cout << " =======================================" <<std::endl;
    std::cout << "Dw src = " <<norm2(src)<<std::endl;
    std::cout << "Dw tmp = " <<norm2(tmp)<<std::endl;
    Dw.M(src,tmp);
    Dw.FreePropagator(tmp,ref,mass);

    std::cout << "Dw ref = " <<norm2(ref)<<std::endl;
    
    ref = ref - src;
    
    std::cout << "Dw ref-src = " <<norm2(ref)<<std::endl;
  }


  ////////////////////////////////////////////////////
  // Wilson prop
  ////////////////////////////////////////////////////
  {
    std::cout<<"****************************************"<<std::endl;
    std::cout << "Wilson Mom space 4d propagator \n";
    std::cout<<"****************************************"<<std::endl;

    LatticeFermionD    src(&GRID); gaussian(pRNG,src);
    LatticeFermionD    tmp(&GRID);
    LatticeFermionD    ref(&GRID);
    LatticeFermionD    diff(&GRID);

    src=Zero();
    Coordinate point(4,0); // 0,0,0,0
    SpinColourVectorD ferm;
    ferm=Zero();
    ferm()(0)(0) = ComplexD(1.0);
    pokeSite(ferm,src,point);

    RealD mass=0.1;
    WilsonFermionD Dw(Umu,GRID,RBGRID,mass);

    // Momentum space prop
    std::cout << " Solving by FFT and Feynman rules" <<std::endl;
    Dw.FreePropagator(src,ref,mass) ;

    Gamma G5(Gamma::Algebra::Gamma5);

    LatticeFermionD    result(&GRID); 
    const int sdir=0;
    
    ////////////////////////////////////////////////////////////////////////
    // Conjugate gradient on normal equations system
    ////////////////////////////////////////////////////////////////////////
    std::cout << " Solving by Conjugate Gradient (CGNE)" <<std::endl;
    Dw.Mdag(src,tmp);
    src=tmp;
    MdagMLinearOperator<WilsonFermionD,LatticeFermionD> HermOp(Dw);
    ConjugateGradient<LatticeFermionD> CG(1.0e-10,10000);
    CG(HermOp,src,result);
    
    ////////////////////////////////////////////////////////////////////////
    std::cout << " Taking difference" <<std::endl;
    std::cout << "Dw result "<<norm2(result)<<std::endl;
    std::cout << "Dw ref     "<<norm2(ref)<<std::endl;
    
    diff = ref - result;
    std::cout << "result - ref     "<<norm2(diff)<<std::endl;

    DumpSliceNorm("Slice Norm Solution ",result,Nd-1);
  }

  ////////////////////////////////////////////////////
  //Gauge invariance test
  ////////////////////////////////////////////////////
  {
    std::cout<<"****************************************"<<std::endl;
    std::cout << "Gauge invariance test \n";
    std::cout<<"****************************************"<<std::endl;
    LatticeGaugeField     U_GT(&GRID); // Gauge transformed field
    LatticeColourMatrix   g(&GRID);    // local Gauge xform matrix
    U_GT = Umu;
    // Make a random xform to teh gauge field
    SU<Nc>::RandomGaugeTransform(pRNG,U_GT,g); // Unit gauge

    LatticeFermionD    src(&GRID);
    LatticeFermionD    tmp(&GRID);
    LatticeFermionD    ref(&GRID);
    LatticeFermionD    diff(&GRID);

    // could loop over colors
    src=Zero();
    Coordinate point(4,0); // 0,0,0,0
    SpinColourVectorD ferm;
    ferm=Zero();
    ferm()(0)(0) = ComplexD(1.0);
    pokeSite(ferm,src,point);

    RealD mass=0.1;
    WilsonFermionD Dw(U_GT,GRID,RBGRID,mass);

    // Momentum space prop
    std::cout << " Solving by FFT and Feynman rules" <<std::endl;
    Dw.FreePropagator(src,ref,mass) ;

    Gamma G5(Gamma::Algebra::Gamma5);

    LatticeFermionD    result(&GRID); 
    const int sdir=0;
    
    ////////////////////////////////////////////////////////////////////////
    // Conjugate gradient on normal equations system
    ////////////////////////////////////////////////////////////////////////
    std::cout << " Solving by Conjugate Gradient (CGNE)" <<std::endl;
    Dw.Mdag(src,tmp);
    src=tmp;
    MdagMLinearOperator<WilsonFermionD,LatticeFermionD> HermOp(Dw);
    ConjugateGradient<LatticeFermionD> CG(1.0e-10,10000);
    CG(HermOp,src,result);
    
    ////////////////////////////////////////////////////////////////////////
    std::cout << " Taking difference" <<std::endl;
    std::cout << "Dw result "<<norm2(result)<<std::endl;
    std::cout << "Dw ref     "<<norm2(ref)<<std::endl;
    
    diff = ref - result;
    std::cout << "result - ref     "<<norm2(diff)<<std::endl;

    DumpSliceNorm("Slice Norm Solution ",result,Nd-1);
  }
  
  
  Grid_finalize();
}
