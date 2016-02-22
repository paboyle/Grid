    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_zmm.cc

    Copyright (C) 2015

Author: paboyle <paboyle@ph.ed.ac.uk>

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
#include <Grid.h>
#include <PerfCount.h>
#include <simd/Avx512Asm.h>


using namespace Grid;
using namespace Grid::QCD;
void WilsonDslashAvx512(void *ptr1,void *ptr2,void *ptr3);
void WilsonDslashAvx512F(void *ptr1,void *ptr2,void *ptr3);
void TimesIAvx512F(void *ptr1,void *ptr3);
void TimesIAvx512(void *ptr1,void *ptr3);



int main(int argc,char **argv)
{
  Grid_init(&argc,&argv);

  
  std::vector<int> latt4 = GridDefaultLatt();
  const int Ls=16;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt4, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  int threads = GridThread::GetThreads();

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});

  GridSerialRNG sRNG; sRNG.SeedFixedIntegers(seeds4);

  vColourMatrixD mat;
  vHalfSpinColourVectorD vec;
  vHalfSpinColourVectorD matvec;
  vHalfSpinColourVectorD ref;
  vComplexD err;

  random(sRNG,mat);
  random(sRNG,vec);

  ref = mat*vec;
  
  WilsonDslashAvx512((void *)&vec, (void *)&mat,(void *)&matvec);

  ref = ref - matvec;
  err = TensorRemove(innerProduct(ref,ref));
  std::cout <<"Double SU3 x 2spin diff   "<< Reduce(err)<<std::endl;

  vColourMatrixF matF;
  vHalfSpinColourVectorF vecF;
  vHalfSpinColourVectorF matvecF;
  vHalfSpinColourVectorF refF;
  vComplexF errF;

  random(sRNG,matF);
  random(sRNG,vecF);

  refF = matF*vecF;

  WilsonDslashAvx512F((void *)&vecF, (void *)&matF,(void *)&matvecF);
  
  refF = refF-matvecF;
  errF = TensorRemove(innerProduct(refF,refF));
  std::cout <<"Single SU3 x 2spin diff   "<< Reduce(errF)<<std::endl;

  TimesIAvx512F((void *)&vecF,(void *)&matvecF);
  refF = timesI(vecF)-matvecF;
  errF = TensorRemove(innerProduct(refF,refF));
  std::cout <<" timesI single diff  "<< Reduce(errF)<<std::endl;

  TimesIAvx512((void *)&vec,(void *)&matvec);
  
  ref = timesI(vec)-matvec;
  err = TensorRemove(innerProduct(ref,ref));
  std::cout <<" timesI double diff  "<< Reduce(err)<<std::endl;

  LatticeFermion src (FGrid);
  LatticeFermion srce(FrbGrid);

  LatticeFermion resulto(FrbGrid); resulto=zero;
  LatticeFermion resulta(FrbGrid); resulta=zero;
  LatticeFermion diff(FrbGrid); 
  LatticeGaugeField Umu(UGrid);

#if 1
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  random(RNG5,src);
  random(RNG4,Umu);
#else
  int mmu=3;
  std::vector<LatticeColourMatrix> U(4,UGrid);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
    if ( mu!=mmu ) U[mu] = zero;
    if ( mu==mmu ) U[mu] = 1.0;
    PokeIndex<LorentzIndex>(Umu,U[mu],mu);
  }
#endif
 pickCheckerboard(Even,srce,src);

  RealD mass=0.1;
  RealD M5  =1.8;
  DomainWallFermionR Dw(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

  std::cout<<GridLogMessage << "Calling Dw"<<std::endl;
  int ncall=50;
  double t0=usecond();
  for(int i=0;i<ncall;i++){
    Dw.DhopOE(srce,resulto,0);
  }
  double t1=usecond();

  double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
  double flops=1344*volume/2;
  
  std::cout<<GridLogMessage << "Called Dw"<<std::endl;
  std::cout<<GridLogMessage << "norm result "<< norm2(resulto)<<std::endl;
  std::cout<<GridLogMessage << "mflop/s =   "<< flops*ncall/(t1-t0)<<std::endl;

  QCD::WilsonFermion5DStatic::AsmOptDslash=1;
  t0=usecond();
  for(int i=0;i<ncall;i++){
    Dw.DhopOE(srce,resulta,0);
  }
  t1=usecond();


  for(int i=0;i< PerformanceCounter::NumTypes(); i++ ){
    Dw.DhopOE(srce,resulta,0);
    PerformanceCounter Counter(i);
    Counter.Start();
    Dw.DhopOE(srce,resulta,0);
    Counter.Stop();
    Counter.Report();
  }
  resulta = (-0.5) * resulta;

  std::cout<<GridLogMessage << "Called Asm Dw"<<std::endl;
  std::cout<<GridLogMessage << "norm result "<< norm2(resulta)<<std::endl;
  std::cout<<GridLogMessage << "mflop/s =   "<< flops*ncall/(t1-t0)<<std::endl;
  diff = resulto-resulta;
  std::cout<<GridLogMessage << "diff "<< norm2(diff)<<std::endl;

}

#undef VLOAD
#undef VSTORE
#undef VMUL
#undef VMADD
#undef ZEND1
#undef ZEND2
#undef ZLOAD
#undef ZMUL
#undef ZMADD

#define VZERO(A) VZEROd(A)
#define VTIMESI(A,B,C) VTIMESId(A,B,C)
#define VTIMESMINUSI(A,B,C) VTIMESMINUSId(A,B,C)

#define VLOAD(OFF,PTR,DEST)       VLOADd(OFF,PTR,DEST)
#define VSTORE(OFF,PTR,SRC)       VSTOREd(OFF,PTR,SRC)
#define VMUL(Uri,Uir,Chi,UChi,Z)  VMULd(Uri,Uir,Chi,UChi,Z)
#define VMADD(Uri,Uir,Chi,UChi,Z) VMADDd(Uri,Uir,Chi,UChi,Z)
#define ZEND1(A,B,C)               ZEND1d(A,B,C)
#define ZEND2(A,B,C)               ZEND2d(A,B,C)
#define ZLOAD(A,B,C,D)            ZLOADd(A,B,C,D)
#define ZMUL(A,B,C,D,E)           ZMULd(A,B,C,D,E)
#define ZMADD(A,B,C,D,E)          ZMADDd(A,B,C,D,E)
#define ZMULMEM2SP(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr) ZMULMEM2SPd(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr) 
#define ZMADDMEM2SP(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr) ZMADDMEM2SPd(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr) 

#define zz Z0

void TimesIAvx512(void *ptr1,void *ptr3)
{
  __asm__ ("mov     $0xAAAA, %%eax "  : : :"%eax");
  __asm__ ("kmov    %%eax, %%k6 " : : :);
  __asm__ ("knot     %%k6, %%k7 " : : :);


  MASK_REGS;
  
  LOAD_CHI(ptr1);
  
  __asm__ (
  VZERO(zz)
  VTIMESI(Chi_00,UChi_00,zz)
  VTIMESI(Chi_01,UChi_01,zz)
  VTIMESI(Chi_02,UChi_02,zz)
  VTIMESI(Chi_10,UChi_10,zz)
  VTIMESI(Chi_11,UChi_11,zz)
  VTIMESI(Chi_12,UChi_12,zz)
  );

  SAVE_UCHI(ptr3);
}

void WilsonDslashAvx512(void *ptr1,void *ptr2,void *ptr3)
{
  int return_address;
  // prototype computed goto to eliminate ABI save restore on call/return in
  // generated assembly.
  static void * table[] = { &&save, &&mult };

  MASK_REGS;

  LOAD_CHI(ptr1);

  return_address = 0;
  goto mult;

 save:
  SAVE_UCHI(ptr3);
  return;

 mult:
  MULT_2SPIN(ptr2);
  goto *table[return_address];

}

#undef VLOAD
#undef VSTORE
#undef VMUL
#undef VMADD
#undef ZEND1
#undef ZEND2
#undef ZLOAD
#undef ZMUL
#undef ZMADD
#undef VZERO
#undef VTIMESI
#undef VTIMESI0
#undef VTIMESI1
#undef VTIMESI2
#undef VTIMESMINUSI
#undef ZMULMEM2SP
#undef ZMADDMEM2SP
#define VZERO(A) VZEROf(A)
#define VMOV(A,B) VMOVf(A,B)
#define VADD(A,B,C) VADDf(A,B,C)
#define VSUB(A,B,C) VSUBf(A,B,C)
#define VTIMESI(A,B,C) VTIMESIf(A,B,C)
#define VTIMESMINUSI(A,B,C) VTIMESMINUSIf(A,B,C)

#define VLOAD(OFF,PTR,DEST)       VLOADf(OFF,PTR,DEST)
#define VSTORE(OFF,PTR,SRC)       VSTOREf(OFF,PTR,SRC)
#define VMUL(Uri,Uir,Chi,UChi,Z)  VMULf(Uri,Uir,Chi,UChi,Z)
#define VMADD(Uri,Uir,Chi,UChi,Z) VMADDf(Uri,Uir,Chi,UChi,Z)
#define ZEND1(A,B,C)               ZEND1f(A,B,C)
#define ZEND2(A,B,C)               ZEND2f(A,B,C)
#define ZLOAD(A,B,C,D)            ZLOADf(A,B,C,D)
#define ZMUL(A,B,C,D,E)           ZMULf(A,B,C,D,E)
#define ZMADD(A,B,C,D,E)          ZMADDf(A,B,C,D,E)
#define ZMULMEM2SP(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr)  ZMULMEM2SPf(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr) 
#define ZMADDMEM2SP(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr) ZMADDMEM2SPf(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr) 

void TimesIAvx512F(void *ptr1,void *ptr3)
{
  MASK_REGS;
  
  LOAD_CHI(ptr1);
  __asm__ (
  VZERO(zz)
  VTIMESI(Chi_00,UChi_00,zz)
  VTIMESI(Chi_01,UChi_01,zz)
  VTIMESI(Chi_02,UChi_02,zz)
  VTIMESI(Chi_10,UChi_10,zz)
  VTIMESI(Chi_11,UChi_11,zz)
  VTIMESI(Chi_12,UChi_12,zz)
	   );
  SAVE_UCHI(ptr3);
}

void WilsonDslashAvx512F(void *ptr1,void *ptr2,void *ptr3)
{
  MASK_REGS;

  LOAD_CHI(ptr1);

  MULT_2SPIN(ptr2);

  SAVE_UCHI(ptr3);

  return;
}

