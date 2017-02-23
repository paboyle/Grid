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
  
  std::vector<int> p({1,3,2,3});

  one = ComplexD(1.0,0.0);
  zz  = ComplexD(0.0,0.0);

  ComplexD ci(0.0,1.0);


  std::cout<<"*************************************************"<<std::endl;
  std::cout<<"Testing Fourier from of known plane wave         "<<std::endl;
  std::cout<<"*************************************************"<<std::endl;
  C=zero;
  for(int mu=0;mu<4;mu++){
    RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
    LatticeCoordinate(coor,mu);
    C = C + (TwoPiL * p[mu]) * coor;
  }

  C = exp(C*ci);
  Csav = C;
  S=zero;
  S = S+C;

  FFT theFFT(&GRID);

  Ctilde=C;
  std::cout<<" Benchmarking FFT of LatticeComplex  "<<std::endl;
  theFFT.FFT_dim(Ctilde,Ctilde,0,FFT::forward); std::cout << theFFT.MFlops()<<" Mflops "<<std::endl;
  theFFT.FFT_dim(Ctilde,Ctilde,1,FFT::forward); std::cout << theFFT.MFlops()<<" Mflops "<<std::endl;
  theFFT.FFT_dim(Ctilde,Ctilde,2,FFT::forward); std::cout << theFFT.MFlops()<<" Mflops "<<std::endl;
  theFFT.FFT_dim(Ctilde,Ctilde,3,FFT::forward); std::cout << theFFT.MFlops()<<" Mflops "<<std::endl;

  //  C=zero;
  //  Ctilde = where(abs(Ctilde)<1.0e-10,C,Ctilde);
  TComplexD cVol;
  cVol()()() = vol;

  Cref=zero;
  pokeSite(cVol,Cref,p);
  //  std::cout <<"Ctilde "<< Ctilde <<std::endl;
  //  std::cout <<"Cref   "<< Cref <<std::endl;

  Cref=Cref-Ctilde;
  std::cout << "diff scalar "<<norm2(Cref) << std::endl;
  C=Csav;
  theFFT.FFT_all_dim(Ctilde,C,FFT::forward);
  theFFT.FFT_all_dim(Cref,Ctilde,FFT::backward); 

  std::cout << norm2(C) << " " << norm2(Ctilde) << " " << norm2(Cref)<< " vol " << vol<< std::endl;

  Cref= Cref - C;
  std::cout << " invertible check " << norm2(Cref)<<std::endl;

  Stilde=S;
  std::cout<<" Benchmarking FFT of LatticeSpinMatrix  "<<std::endl;
  theFFT.FFT_dim(Stilde,S,0,FFT::forward); std::cout << theFFT.MFlops()<<" mflops "<<std::endl;
  theFFT.FFT_dim(Stilde,S,1,FFT::forward); std::cout << theFFT.MFlops()<<" mflops "<<std::endl;
  theFFT.FFT_dim(Stilde,S,2,FFT::forward); std::cout << theFFT.MFlops()<<" mflops "<<std::endl;
  theFFT.FFT_dim(Stilde,S,3,FFT::forward); std::cout << theFFT.MFlops()<<" mflops "<<std::endl;

  SpinMatrixD Sp; 
  Sp = zero; Sp = Sp+cVol;

  S=zero;
  pokeSite(Sp,S,p);

  S= S-Stilde;
  std::cout << "diff FT[SpinMat] "<<norm2(S) << std::endl;

  /*
   */
  std::vector<int> seeds({1,2,3,4});
  GridSerialRNG          sRNG;  sRNG.SeedFixedIntegers(seeds); // naughty seeding
  GridParallelRNG          pRNG(&GRID);
  pRNG.SeedFixedIntegers(seeds);

  LatticeGaugeFieldD Umu(&GRID);

  SU3::ColdConfiguration(pRNG,Umu); // Unit gauge
  //  Umu=zero;
  ////////////////////////////////////////////////////
  // Wilson test
  ////////////////////////////////////////////////////
  {
    LatticeFermionD    src(&GRID); gaussian(pRNG,src);
    LatticeFermionD    tmp(&GRID);
    LatticeFermionD    ref(&GRID);
    
    RealD mass=0.01;
    WilsonFermionD Dw(Umu,GRID,RBGRID,mass);
    
    Dw.M(src,tmp);

    std::cout << "Dw src = " <<norm2(src)<<std::endl;
    std::cout << "Dw tmp = " <<norm2(tmp)<<std::endl;
    
    Dw.FreePropagator(tmp,ref,mass);

    std::cout << "Dw ref = " <<norm2(ref)<<std::endl;
    
    ref = ref - src;
    
    std::cout << "Dw ref-src = " <<norm2(ref)<<std::endl;
  }

  ////////////////////////////////////////////////////
  // Dwf matrix
  ////////////////////////////////////////////////////
  {
    std::cout<<"****************************************"<<std::endl;
    std::cout<<"Testing Fourier representation of Ddwf"<<std::endl;
    std::cout<<"****************************************"<<std::endl;
    
    const int Ls=16;
    const int sdir=0;
    RealD mass=0.01;
    RealD M5  =1.0;
    Gamma G5(Gamma::Algebra::Gamma5);

    GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,&GRID);
    GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,&GRID);

    std::cout<<"Making Ddwf"<<std::endl;
    DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,GRID,RBGRID,mass,M5);

    GridParallelRNG          RNG5(FGrid);      RNG5.SeedFixedIntegers(seeds);
    LatticeFermionD    src5(FGrid); gaussian(RNG5,src5);
    LatticeFermionD    src5_p(FGrid); 
    LatticeFermionD    result5(FGrid); 
    LatticeFermionD    ref5(FGrid); 
    LatticeFermionD    tmp5(FGrid); 

    /////////////////////////////////////////////////////////////////
    // result5 is the non pert operator in 4d mom space
    /////////////////////////////////////////////////////////////////
    Ddwf.M(src5,tmp5); 
    ref5 = tmp5;
    
    FFT theFFT5(FGrid);
      
    theFFT5.FFT_dim(result5,tmp5,1,FFT::forward); tmp5 = result5;
    theFFT5.FFT_dim(result5,tmp5,2,FFT::forward); tmp5 = result5;
    theFFT5.FFT_dim(result5,tmp5,3,FFT::forward); tmp5 = result5;
    theFFT5.FFT_dim(result5,tmp5,4,FFT::forward); result5 = result5*ComplexD(::sqrt(1.0/vol),0.0);
    
    std::cout<<"Fourier xformed Ddwf"<<std::endl;
    
    tmp5 = src5;
    theFFT5.FFT_dim(src5_p,tmp5,1,FFT::forward); tmp5 = src5_p;
    theFFT5.FFT_dim(src5_p,tmp5,2,FFT::forward); tmp5 = src5_p;
    theFFT5.FFT_dim(src5_p,tmp5,3,FFT::forward); tmp5 = src5_p;
    theFFT5.FFT_dim(src5_p,tmp5,4,FFT::forward); src5_p = src5_p*ComplexD(::sqrt(1.0/vol),0.0);

    std::cout<<"Fourier xformed src5"<<std::endl;
      
    /////////////////////////////////////////////////////////////////
    // work out the predicted from Fourier
    /////////////////////////////////////////////////////////////////
    Gamma::Algebra Gmu [] = {
      Gamma::Algebra::GammaX,
      Gamma::Algebra::GammaY,
      Gamma::Algebra::GammaZ,
      Gamma::Algebra::GammaT,
      Gamma::Algebra::Gamma5
    };
    LatticeFermionD    Kinetic(FGrid); Kinetic = zero;
    LatticeComplexD    kmu(FGrid); 
    LatticeInteger     scoor(FGrid); 
    LatticeComplexD    sk (FGrid); sk = zero;
    LatticeComplexD    sk2(FGrid); sk2= zero;
    LatticeComplexD    W(FGrid); W= zero;
    //      LatticeComplexD    a(FGrid); a= zero;
    LatticeComplexD    one(FGrid); one =ComplexD(1.0,0.0);
    ComplexD ci(0.0,1.0);
    
    for(int mu=0;mu<Nd;mu++) {
      
      LatticeCoordinate(kmu,mu+1);
      
      RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
      
      kmu = TwoPiL * kmu;
      
      sk2 = sk2 + 2.0*sin(kmu*0.5)*sin(kmu*0.5);
      sk  = sk  +     sin(kmu)    *sin(kmu); 
      
      // -1/2 Dw ->  1/2 gmu (eip - emip) = i sinp gmu
      Kinetic = Kinetic + sin(kmu)*ci*(Gamma(Gmu[mu])*src5_p);
      
    }
    
    // NB implicit sum over mu
    //
    // 1-1/2 Dw = 1 - 1/2 ( eip+emip)
    //          = - 1/2 (ei - 2 +  emi) 
    //          = - 1/4 2 (eih - eimh)(eih - eimh) 
    //          = 2 sink/2 ink/2 = sk2
    
    W = one - M5 + sk2; 
    Kinetic = Kinetic + W * src5_p;
    
    LatticeCoordinate(scoor,sdir);
    
    tmp5 = Cshift(src5_p,sdir,+1);
    tmp5 = (tmp5 - G5*tmp5)*0.5;
    tmp5 = where(scoor==Integer(Ls-1),mass*tmp5,-tmp5);
    Kinetic = Kinetic + tmp5;
    
    tmp5 = Cshift(src5_p,sdir,-1);
    tmp5 = (tmp5 + G5*tmp5)*0.5;
    tmp5 = where(scoor==Integer(0),mass*tmp5,-tmp5);
    Kinetic = Kinetic + tmp5;
    
    std::cout<<"Momentum space Ddwf  "<< norm2(Kinetic)<<std::endl;
    std::cout<<"Stencil Ddwf         "<< norm2(result5)<<std::endl;
    
    result5 = result5 - Kinetic;
    std::cout<<"diff "<< norm2(result5)<<std::endl;
    
  }

  ////////////////////////////////////////////////////
  // Dwf prop
  ////////////////////////////////////////////////////
  {
    std::cout<<"****************************************"<<std::endl;
    std::cout << "Testing Ddwf Ht Mom space 4d propagator \n";
    std::cout<<"****************************************"<<std::endl;

    LatticeFermionD    src(&GRID); gaussian(pRNG,src);
    LatticeFermionD    tmp(&GRID);
    LatticeFermionD    ref(&GRID);
    LatticeFermionD    diff(&GRID);

    std::vector<int> point(4,0);
    src=zero;
    SpinColourVectorD ferm; gaussian(sRNG,ferm);
    pokeSite(ferm,src,point);

    const int Ls=32;
    GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,&GRID);
    GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,&GRID);

    RealD mass=0.01;
    RealD M5  =0.8;
    DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,GRID,RBGRID,mass,M5);

    // Momentum space prop
    std::cout << " Solving by FFT and Feynman rules" <<std::endl;
    Ddwf.FreePropagator(src,ref,mass) ;

    Gamma G5(Gamma::Algebra::Gamma5);

    LatticeFermionD    src5(FGrid); src5=zero;
    LatticeFermionD    tmp5(FGrid); 
    LatticeFermionD    result5(FGrid); result5=zero;
    LatticeFermionD    result4(&GRID); 
    const int sdir=0;

    ////////////////////////////////////////////////////////////////////////
    // Domain wall physical field source
    ////////////////////////////////////////////////////////////////////////
    /*
	chi_5[0]   = chiralProjectPlus(chi);
	chi_5[Ls-1]= chiralProjectMinus(chi);
    */      
    tmp =   (src + G5*src)*0.5;      InsertSlice(tmp,src5,   0,sdir);
    tmp =   (src - G5*src)*0.5;      InsertSlice(tmp,src5,Ls-1,sdir);
    
    ////////////////////////////////////////////////////////////////////////
    // Conjugate gradient on normal equations system
    ////////////////////////////////////////////////////////////////////////
    std::cout << " Solving by Conjugate Gradient (CGNE)" <<std::endl;
    Ddwf.Mdag(src5,tmp5);
    src5=tmp5;
    MdagMLinearOperator<DomainWallFermionD,LatticeFermionD> HermOp(Ddwf);
    ConjugateGradient<LatticeFermionD> CG(1.0e-16,10000);
    CG(HermOp,src5,result5);
    
    ////////////////////////////////////////////////////////////////////////
    // Domain wall physical field propagator
    ////////////////////////////////////////////////////////////////////////
    /*
      psi  = chiralProjectMinus(psi_5[0]);
      psi += chiralProjectPlus(psi_5[Ls-1]);
    */
    ExtractSlice(tmp,result5,0   ,sdir);   result4 =         (tmp-G5*tmp)*0.5;
    ExtractSlice(tmp,result5,Ls-1,sdir);   result4 = result4+(tmp+G5*tmp)*0.5;
    
    std::cout << " Taking difference" <<std::endl;
    std::cout << "Ddwf result4 "<<norm2(result4)<<std::endl;
    std::cout << "Ddwf ref     "<<norm2(ref)<<std::endl;
    
    diff = ref - result4;
    std::cout << "result - ref     "<<norm2(diff)<<std::endl;

  }



  ////////////////////////////////////////////////////
  // Dwf prop
  ////////////////////////////////////////////////////
  {
    std::cout<<"****************************************"<<std::endl;
    std::cout << "Testing Dov Ht Mom space 4d propagator \n";
    std::cout<<"****************************************"<<std::endl;

    LatticeFermionD    src(&GRID); gaussian(pRNG,src);
    LatticeFermionD    tmp(&GRID);
    LatticeFermionD    ref(&GRID);
    LatticeFermionD    diff(&GRID);

    std::vector<int> point(4,0);
    src=zero;
    SpinColourVectorD ferm; gaussian(sRNG,ferm);
    pokeSite(ferm,src,point);

    const int Ls=48;
    GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,&GRID);
    GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,&GRID);

    RealD mass=0.01;
    RealD M5  =0.8;

    OverlapWilsonCayleyTanhFermionD Dov(Umu,*FGrid,*FrbGrid,GRID,RBGRID,mass,M5,1.0);

    // Momentum space prop
    std::cout << " Solving by FFT and Feynman rules" <<std::endl;
    Dov.FreePropagator(src,ref,mass) ;

    Gamma G5(Gamma::Algebra::Gamma5);

    LatticeFermionD    src5(FGrid); src5=zero;
    LatticeFermionD    tmp5(FGrid); 
    LatticeFermionD    result5(FGrid); result5=zero;
    LatticeFermionD    result4(&GRID); 
    const int sdir=0;

    ////////////////////////////////////////////////////////////////////////
    // Domain wall physical field source; need D_minus
    ////////////////////////////////////////////////////////////////////////
    /*
	chi_5[0]   = chiralProjectPlus(chi);
	chi_5[Ls-1]= chiralProjectMinus(chi);
    */      
    tmp =   (src + G5*src)*0.5;      InsertSlice(tmp,src5,   0,sdir);
    tmp =   (src - G5*src)*0.5;      InsertSlice(tmp,src5,Ls-1,sdir);

    
    ////////////////////////////////////////////////////////////////////////
    // Conjugate gradient on normal equations system
    ////////////////////////////////////////////////////////////////////////
    std::cout << " Solving by Conjugate Gradient (CGNE)" <<std::endl;
    Dov.Dminus(src5,tmp5);
    src5=tmp5;
    Dov.Mdag(src5,tmp5);
    src5=tmp5;
    MdagMLinearOperator<OverlapWilsonCayleyTanhFermionD,LatticeFermionD> HermOp(Dov);
    ConjugateGradient<LatticeFermionD> CG(1.0e-16,10000);
    CG(HermOp,src5,result5);
    
    ////////////////////////////////////////////////////////////////////////
    // Domain wall physical field propagator
    ////////////////////////////////////////////////////////////////////////
    /*
      psi  = chiralProjectMinus(psi_5[0]);
      psi += chiralProjectPlus(psi_5[Ls-1]);
    */
    ExtractSlice(tmp,result5,0   ,sdir);   result4 =         (tmp-G5*tmp)*0.5;
    ExtractSlice(tmp,result5,Ls-1,sdir);   result4 = result4+(tmp+G5*tmp)*0.5;
    
    std::cout << " Taking difference" <<std::endl;
    std::cout << "Dov result4 "<<norm2(result4)<<std::endl;
    std::cout << "Dov ref     "<<norm2(ref)<<std::endl;
    
    diff = ref - result4;
    std::cout << "result - ref     "<<norm2(diff)<<std::endl;

  }

  {
    /*
     * 
    typedef GaugeImplTypes<vComplexD, 1> QEDGimplTypesD;
    typedef Photon<QEDGimplTypesD>       QEDGaction;

    QEDGaction Maxwell(QEDGaction::FEYNMAN_L);
    QEDGaction::GaugeField Prop(&GRID);Prop=zero;
    QEDGaction::GaugeField Source(&GRID);Source=zero;

    Maxwell.FreePropagator (Source,Prop);
    std::cout << " MaxwellFree propagator\n";
    */
  }
  Grid_finalize();
}
