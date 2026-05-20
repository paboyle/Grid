/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/core/Test_planned_fft.cc

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

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  int vol = 1;
  for(int d=0;d<latt_size.size();d++) vol *= latt_size[d];

  GridCartesian         GRID(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian RBGRID(&GRID);

  LatticeComplexD     one(&GRID);
  LatticeComplexD      zz(&GRID);
  LatticeComplexD       C(&GRID);
  LatticeComplexD  Ctilde(&GRID);
  LatticeComplexD  Cref  (&GRID);
  LatticeComplexD  Csav  (&GRID);
  LatticeComplexD    coor(&GRID);

  LatticeSpinMatrixD    S(&GRID);
  LatticeSpinMatrixD    Stilde(&GRID);

  Coordinate p({1,3,2,3});

  one = ComplexD(1.0,0.0);
  zz  = ComplexD(0.0,0.0);
  ComplexD ci(0.0,1.0);

  std::cout<<"*************************************************"<<std::endl;
  std::cout<<"Testing Fourier form of known plane wave         "<<std::endl;
  std::cout<<"*************************************************"<<std::endl;
  C=Zero();
  for(int mu=0;mu<4;mu++){
    RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
    LatticeCoordinate(coor,mu);
    C = C + (TwoPiL * p[mu]) * coor;
  }
  C = exp(C*ci);
  Csav = C;
  S=Zero();
  S = S+C;

  // PlannedFFT is templated on the lattice element type (vector_object), not the Lattice<> itself.
  PlannedFFT<LatticeComplexD::vector_object>    theFFT(&GRID);
  PlannedFFT<LatticeSpinMatrixD::vector_object> theFFT_spin(&GRID);

  Ctilde=C;
  std::cout<<" Benchmarking PlannedFFT of LatticeComplex  "<<std::endl;
  theFFT.FFT_dim(Ctilde,Ctilde,0,FFTbase::forward); std::cout << theFFT.MFlops()<<" Mflops "<<std::endl;
  theFFT.FFT_dim(Ctilde,Ctilde,1,FFTbase::forward); std::cout << theFFT.MFlops()<<" Mflops "<<std::endl;
  theFFT.FFT_dim(Ctilde,Ctilde,2,FFTbase::forward); std::cout << theFFT.MFlops()<<" Mflops "<<std::endl;
  theFFT.FFT_dim(Ctilde,Ctilde,3,FFTbase::forward); std::cout << theFFT.MFlops()<<" Mflops "<<std::endl;

  TComplexD cVol;
  cVol()()() = vol;

  Cref=Zero();
  pokeSite(cVol,Cref,p);

  Cref=Cref-Ctilde;
  std::cout << "diff scalar "<<norm2(Cref) << std::endl;

  C=Csav;
  theFFT.FFT_all_dim(Ctilde,C,FFTbase::forward);
  theFFT.FFT_all_dim(Cref,Ctilde,FFTbase::backward);

  std::cout << norm2(C) << " " << norm2(Ctilde) << " " << norm2(Cref)<< " vol " << vol<< std::endl;

  Cref= Cref - C;
  std::cout << " invertible check " << norm2(Cref)<<std::endl;

  Stilde=S;
  std::cout<<" Benchmarking PlannedFFT of LatticeSpinMatrix  "<<std::endl;
  theFFT_spin.FFT_dim(Stilde,Stilde,0,FFTbase::forward); std::cout << theFFT_spin.MFlops()<<" mflops "<<std::endl;
  theFFT_spin.FFT_dim(Stilde,Stilde,1,FFTbase::forward); std::cout << theFFT_spin.MFlops()<<" mflops "<<std::endl;
  theFFT_spin.FFT_dim(Stilde,Stilde,2,FFTbase::forward); std::cout << theFFT_spin.MFlops()<<" mflops "<<std::endl;
  theFFT_spin.FFT_dim(Stilde,Stilde,3,FFTbase::forward); std::cout << theFFT_spin.MFlops()<<" mflops "<<std::endl;

  SpinMatrixD Sp;
  Sp = Zero(); Sp = Sp+cVol;

  S=Zero();
  pokeSite(Sp,S,p);

  S= S-Stilde;
  std::cout << "diff FT[SpinMat] "<<norm2(S) << std::endl;

  std::vector<int> seeds({1,2,3,4});
  GridSerialRNG          sRNG;  sRNG.SeedFixedIntegers(seeds);
  GridParallelRNG          pRNG(&GRID);
  pRNG.SeedFixedIntegers(seeds);

  LatticeGaugeFieldD Umu(&GRID);
  SU<Nc>::ColdConfiguration(pRNG,Umu);

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
  // Dwf matrix — verify Fourier representation using PlannedFFT<LatticeFermionD>
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

    DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,GRID,RBGRID,mass,M5);

    GridParallelRNG          RNG5(FGrid);      RNG5.SeedFixedIntegers(seeds);
    LatticeFermionD    src5(FGrid); gaussian(RNG5,src5);
    LatticeFermionD    src5_p(FGrid);
    LatticeFermionD    result5(FGrid);
    LatticeFermionD    ref5(FGrid);
    LatticeFermionD    tmp5(FGrid);

    Ddwf.M(src5,tmp5);
    ref5 = tmp5;

    PlannedFFT<LatticeFermionD::vector_object> theFFT5(FGrid);

    theFFT5.FFT_dim(result5,tmp5,1,FFTbase::forward); tmp5 = result5;
    std::cout<<"Fourier xformed Ddwf 1 "<<norm2(result5)<<std::endl;
    theFFT5.FFT_dim(result5,tmp5,2,FFTbase::forward); tmp5 = result5;
    std::cout<<"Fourier xformed Ddwf 2 "<<norm2(result5)<<std::endl;
    theFFT5.FFT_dim(result5,tmp5,3,FFTbase::forward); tmp5 = result5;
    std::cout<<"Fourier xformed Ddwf 3 "<<norm2(result5)<<std::endl;
    theFFT5.FFT_dim(result5,tmp5,4,FFTbase::forward);
    std::cout<<"Fourier xformed Ddwf 4 "<<norm2(result5)<<std::endl;
    result5 = result5*ComplexD(::sqrt(1.0/vol),0.0);

    std::cout<<"Fourier xformed Ddwf "<<norm2(result5)<<std::endl;

    tmp5 = src5;
    theFFT5.FFT_dim(src5_p,tmp5,1,FFTbase::forward); tmp5 = src5_p;
    theFFT5.FFT_dim(src5_p,tmp5,2,FFTbase::forward); tmp5 = src5_p;
    theFFT5.FFT_dim(src5_p,tmp5,3,FFTbase::forward); tmp5 = src5_p;
    theFFT5.FFT_dim(src5_p,tmp5,4,FFTbase::forward); src5_p = src5_p*ComplexD(::sqrt(1.0/vol),0.0);

    std::cout<<"Fourier xformed src5"<< norm2(src5)<<" -> "<<norm2(src5_p)<<std::endl;

    Gamma::Algebra Gmu [] = {
      Gamma::Algebra::GammaX,
      Gamma::Algebra::GammaY,
      Gamma::Algebra::GammaZ,
      Gamma::Algebra::GammaT,
      Gamma::Algebra::Gamma5
    };
    LatticeFermionD    Kinetic(FGrid); Kinetic = Zero();
    LatticeComplexD    kmu(FGrid);
    LatticeInteger     scoor(FGrid);
    LatticeComplexD    sk (FGrid); sk = Zero();
    LatticeComplexD    sk2(FGrid); sk2= Zero();
    LatticeComplexD    W(FGrid); W= Zero();
    LatticeComplexD    one5(FGrid); one5 =ComplexD(1.0,0.0);

    for(int mu=0;mu<Nd;mu++) {
      LatticeCoordinate(kmu,mu+1);
      RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
      kmu = TwoPiL * kmu;
      sk2 = sk2 + 2.0*sin(kmu*0.5)*sin(kmu*0.5);
      sk  = sk  +     sin(kmu)    *sin(kmu);
      Kinetic = Kinetic + sin(kmu)*ci*(Gamma(Gmu[mu])*src5_p);
    }
    std::cout << " src5    "<<norm2(src5_p)<<std::endl;
    std::cout << " Kinetic "<<norm2(Kinetic)<<std::endl;

    W = one5 - M5 + sk2;
    std::cout << " W "<<norm2(W)<<std::endl;
    Kinetic = Kinetic + W * src5_p;
    std::cout << " Kinetic "<<norm2(Kinetic)<<std::endl;

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
    GRID_ASSERT(norm2(result5)<1.0e-4);
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

    Coordinate point(4,0);
    src=Zero();
    SpinColourVectorD ferm; gaussian(sRNG,ferm);
    pokeSite(ferm,src,point);

    const int Ls=32;
    GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,&GRID);
    GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,&GRID);

    RealD mass=0.01;
    RealD M5  =0.8;
    DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,GRID,RBGRID,mass,M5);

    std::cout << " Solving by FFT and Feynman rules" <<std::endl;
    bool fiveD = false;
    Ddwf.FreePropagator(src,ref,mass,fiveD);

    Gamma G5(Gamma::Algebra::Gamma5);

    LatticeFermionD    src5(FGrid); src5=Zero();
    LatticeFermionD    tmp5(FGrid);
    LatticeFermionD    result5(FGrid); result5=Zero();
    LatticeFermionD    result4(&GRID);
    const int sdir=0;

    tmp =   (src + G5*src)*0.5;      InsertSlice(tmp,src5,   0,sdir);
    tmp =   (src - G5*src)*0.5;      InsertSlice(tmp,src5,Ls-1,sdir);

    std::cout << " Solving by Conjugate Gradient (CGNE)" <<std::endl;
    Ddwf.Mdag(src5,tmp5);
    src5=tmp5;
    MdagMLinearOperator<DomainWallFermionD,LatticeFermionD> HermOp(Ddwf);
    ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
    CG(HermOp,src5,result5);

    ExtractSlice(tmp,result5,0   ,sdir);   result4 =         (tmp-G5*tmp)*0.5;
    ExtractSlice(tmp,result5,Ls-1,sdir);   result4 = result4+(tmp+G5*tmp)*0.5;

    std::cout << " Taking difference" <<std::endl;
    std::cout << "Ddwf result4 "<<norm2(result4)<<std::endl;
    std::cout << "Ddwf ref     "<<norm2(ref)<<std::endl;

    diff = ref - result4;
    std::cout << "result - ref     "<<norm2(diff)<<std::endl;
    GRID_ASSERT(norm2(diff)<1.0e-4);
  }

  Grid_finalize();
}
