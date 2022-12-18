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
 ;

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

    RealD mass=0.01;
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

  
  Grid_finalize();
}
