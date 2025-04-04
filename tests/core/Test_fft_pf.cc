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

  ComplexD ci(0.0,1.0);

  std::vector<int> seeds({1,2,3,4});
  GridSerialRNG          sRNG;  sRNG.SeedFixedIntegers(seeds); // naughty seeding
  GridParallelRNG          pRNG(&GRID);
  pRNG.SeedFixedIntegers(seeds);

  LatticeGaugeFieldD Umu(&GRID);

  SU<Nc>::ColdConfiguration(pRNG,Umu); // Unit gauge

  ////////////////////////////////////////////////////
  // PF prop
  ////////////////////////////////////////////////////
  LatticeFermionD    src(&GRID);

  gaussian(pRNG,src);
#if 1
    Coordinate point(4,0);
    src=Zero();
    SpinColourVectorD ferm; gaussian(sRNG,ferm);
    pokeSite(ferm,src,point);
#endif
  
  {
    std::cout<<"****************************************"<<std::endl;
    std::cout << "Testing OverlapWilsonPartialFractionTanhFermionD Hw kernel Mom space 4d propagator \n";
    std::cout<<"****************************************"<<std::endl;

    //    LatticeFermionD    src(&GRID); gaussian(pRNG,src);
    LatticeFermionD    tmp(&GRID);
    LatticeFermionD    ref(&GRID);
    LatticeFermionD    diff(&GRID);

    const int Ls=48+1;
    GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,&GRID);
    GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,&GRID);

    RealD mass=0.1;
    RealD M5  =0.8;
    OverlapWilsonPartialFractionTanhFermionD Dov(Umu,*FGrid,*FrbGrid,GRID,RBGRID,mass,M5,1.0);

    // Momentum space prop
    std::cout << " Solving by FFT and Feynman rules" <<std::endl;
    bool fiveD = false; //calculate 4d free propagator

    std::cout << " Free propagator " <<std::endl;
    Dov.FreePropagator(src,ref,mass) ;
    std::cout << " Free propagator norm "<< norm2(ref) <<std::endl;

    Gamma G5(Gamma::Algebra::Gamma5);

    LatticeFermionD    src5(FGrid); src5=Zero();
    LatticeFermionD    tmp5(FGrid); 
    LatticeFermionD    result5(FGrid); result5=Zero();
    LatticeFermionD    result4(&GRID); 
    const int sdir=0;

    ////////////////////////////////////////////////////////////////////////
    // Import
    ////////////////////////////////////////////////////////////////////////
    std::cout << " Free propagator Import "<< norm2(src) <<std::endl;
    Dov.ImportPhysicalFermionSource  (src,src5);
    std::cout << " Free propagator Imported "<< norm2(src5) <<std::endl;
    
    ////////////////////////////////////////////////////////////////////////
    // Conjugate gradient on normal equations system
    ////////////////////////////////////////////////////////////////////////
    std::cout << " Solving by Conjugate Gradient (CGNE)" <<std::endl;
    Dov.Mdag(src5,tmp5);
    src5=tmp5;
    MdagMLinearOperator<OverlapWilsonPartialFractionTanhFermionD,LatticeFermionD> HermOp(Dov);
    ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
    CG(HermOp,src5,result5);
    std::cout << " Solved by Conjugate Gradient (CGNE)" <<std::endl;
    ////////////////////////////////////////////////////////////////////////
    // Domain wall physical field propagator
    ////////////////////////////////////////////////////////////////////////
    Dov.ExportPhysicalFermionSolution(result5,result4);

    // From DWF4d.pdf :
    //
    // Dov_pf = 2/(1-m) D_cayley_ovlap  [ Page 43 ]
    // Dinv_cayley_ovlap = 2/(1-m) Dinv_pf 
    // Dinv_cayley_surface =1/(1-m) ( Dinv_cayley_ovlap - 1 ) =>  2/(1-m)^2 Dinv_pf - 1/(1-m) * src   [ Eq.2.67 ]

    RealD scale = 2.0/(1.0-mass)/(1.0-mass);
    result4 = result4 * scale;
    result4 = result4 - src*(1.0/(1.0-mass)); // Subtract contact term
    DumpSliceNorm("Src",src);
    DumpSliceNorm("Grid",result4);
    DumpSliceNorm("Fourier",ref);

    std::cout << "Dov result4 "<<norm2(result4)<<std::endl;
    std::cout << "Dov ref     "<<norm2(ref)<<std::endl;

    diff = result4- ref;
    DumpSliceNorm("diff ",diff);
    
  }
  
  ////////////////////////////////////////////////////
  // Dwf prop
  ////////////////////////////////////////////////////
  {
    std::cout<<"****************************************"<<std::endl;
    std::cout << "Testing OverlapWilsonCayleyTanhFermionD space 4d propagator \n";
    std::cout<<"****************************************"<<std::endl;

    LatticeFermionD    tmp(&GRID);
    LatticeFermionD    ref(&GRID);
    LatticeFermionD    diff(&GRID);

    const int Ls=48;
    GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,&GRID);
    GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,&GRID);

    RealD mass=0.1;
    RealD M5  =0.8;

    OverlapWilsonCayleyTanhFermionD Dov(Umu,*FGrid,*FrbGrid,GRID,RBGRID,mass,M5,1.0);

    // Momentum space prop
    std::cout << " Solving by FFT and Feynman rules" <<std::endl;
    Dov.FreePropagator(src,ref,mass) ;

    Gamma G5(Gamma::Algebra::Gamma5);

    LatticeFermionD    src5(FGrid); src5=Zero();
    LatticeFermionD    tmp5(FGrid); 
    LatticeFermionD    result5(FGrid); result5=Zero();
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
    DumpSliceNorm("Grid",result4);
    DumpSliceNorm("Fourier",ref);
    diff = ref - result4;
    std::cout << "result - ref     "<<norm2(diff)<<std::endl;
    
    DumpSliceNorm("diff",diff);

  }

  
  {
    std::cout<<"****************************************"<<std::endl;
    std::cout<<"Testing OverlapWilsonPartialFractionTanhFermionD Hw kernel Mom space 4d propagator with q\n";
    std::cout<<"****************************************"<<std::endl;

    //    LatticeFermionD    src(&GRID); gaussian(pRNG,src);
    LatticeFermionD    tmp(&GRID);
    LatticeFermionD    ref(&GRID);
    LatticeFermionD    diff(&GRID);

    const int Ls=48+1;
    GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,&GRID);
    GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,&GRID);

    RealD mass=0.1;
    RealD M5  =0.8;
    OverlapWilsonPartialFractionTanhFermionD Dov(Umu,*FGrid,*FrbGrid,GRID,RBGRID,mass,M5,1.0);
    std::vector<RealD> qmu({1.0,0.0,0.0,0.0});
    Dov.set_qmu(qmu);

    // Momentum space prop
    std::cout << " Solving by FFT and Feynman rules" <<std::endl;
    bool fiveD = false; //calculate 4d free propagator

    std::cout << " Free propagator " <<std::endl;
    Dov.FreePropagator(src,ref,mass) ;
    std::cout << " Free propagator norm "<< norm2(ref) <<std::endl;

    Gamma G5(Gamma::Algebra::Gamma5);

    LatticeFermionD    src5(FGrid); src5=Zero();
    LatticeFermionD    tmp5(FGrid); 
    LatticeFermionD    result5(FGrid); result5=Zero();
    LatticeFermionD    result4(&GRID); 
    const int sdir=0;

    ////////////////////////////////////////////////////////////////////////
    // Import
    ////////////////////////////////////////////////////////////////////////
    std::cout << " Free propagator Import "<< norm2(src) <<std::endl;
    Dov.ImportPhysicalFermionSource  (src,src5);
    std::cout << " Free propagator Imported "<< norm2(src5) <<std::endl;
    
    ////////////////////////////////////////////////////////////////////////
    // Conjugate gradient on normal equations system
    ////////////////////////////////////////////////////////////////////////
    std::cout << " Solving by Conjugate Gradient (CGNE)" <<std::endl;
    Dov.Mdag(src5,tmp5);
    src5=tmp5;
    MdagMLinearOperator<OverlapWilsonPartialFractionTanhFermionD,LatticeFermionD> HermOp(Dov);
    ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
    CG(HermOp,src5,result5);
    ////////////////////////////////////////////////////////////////////////
    // Domain wall physical field propagator
    ////////////////////////////////////////////////////////////////////////
    Dov.ExportPhysicalFermionSolution(result5,result4);

    // From DWF4d.pdf :
    //
    // Dov_pf = 2/(1-m) D_cayley_ovlap  [ Page 43 ]
    // Dinv_cayley_ovlap = 2/(1-m) Dinv_pf 
    // Dinv_cayley_surface =1/(1-m) ( Dinv_cayley_ovlap - 1 ) =>  2/(1-m)^2 Dinv_pf - 1/(1-m) * src   [ Eq.2.67 ]

    RealD scale = 2.0/(1.0-mass)/(1.0-mass);
    result4 = result4 * scale;
    result4 = result4 - src*(1.0/(1.0-mass)); // Subtract contact term
    DumpSliceNorm("Src",src);
    DumpSliceNorm("Grid",result4);
    DumpSliceNorm("Fourier",ref);

    std::cout << "Dov result4 "<<norm2(result4)<<std::endl;
    std::cout << "Dov ref     "<<norm2(ref)<<std::endl;

    diff = result4- ref;
    DumpSliceNorm("diff ",diff);
    
  }

  
  Grid_finalize();
}
