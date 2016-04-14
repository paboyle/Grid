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
#include <simd/Intel512wilson.h>


using namespace Grid;
using namespace Grid::QCD;

void ZmulF(void *ptr1,void *ptr2,void *ptr3);
void Zmul(void *ptr1,void *ptr2,void *ptr3);
void WilsonDslashAvx512(void *ptr1,void *ptr2,void *ptr3);
void WilsonDslashAvx512F(void *ptr1,void *ptr2,void *ptr3);
void TimesIAvx512F(void *ptr1,void *ptr3);
void TimesIAvx512(void *ptr1,void *ptr3);
void TimesMinusIAvx512F(void *ptr1,void *ptr3);
void TimesMinusIAvx512(void *ptr1,void *ptr3);



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
  vHalfSpinColourVectorD vec1;
  vHalfSpinColourVectorD vec2;
  vHalfSpinColourVectorD vec3;
 
  vHalfSpinColourVectorD matvec;
  vHalfSpinColourVectorD ref;
  vComplexD err;

  random(sRNG,vec1); 
  vec1 = std::complex<double>(0.1,3.0);
  random(sRNG,vec2);
  vec2=2.0;
  random(sRNG,vec3);

  //std::cout << "Zmul  vec1"<<vec1<<" &vec1 "<<& vec1<<std::endl;
  //std::cout << "Zmul  vec2"<<vec2<<" &vec2 "<<& vec2<<std::endl;
  //std::cout << "Zmul  vec3"<<vec3<<" &vec3 "<<& vec3<<std::endl;
  for(int sp=0;sp<2;sp++){
  for(int co=0;co<3;co++){
  ref()(sp)(co) = vec1()(sp)(co)*vec2()(sp)(co);
  }}

  Zmul((void *)&vec1,(void *)&vec2,(void *)&vec3);
  //std::cout << "Zmul  vec3"<<vec3<<" &vec3 "<<& vec3<<std::endl;
  //std::cout << "Zmul \n\t ref "<<ref<<"\n\t vec3"<<vec3 <<std::endl;
  ref = ref - vec3;
  err = TensorRemove(innerProduct(ref,ref));
  std::cout <<"Zmul diff   "<< Reduce(err)<<std::endl;

  random(sRNG,mat);
  mat = zero;
  mat()()(0,0) = 1.0;
  random(sRNG,vec);

  ref = mat*vec;
  
  WilsonDslashAvx512((void *)&vec, (void *)&mat,(void *)&matvec);

  //std::cout << ref   <<std::endl;
  //std::cout << matvec<<std::endl;
  ref = ref - matvec;
  err = TensorRemove(innerProduct(ref,ref));
  std::cout <<"Double SU3 x 2spin diff   "<< Reduce(err)<<std::endl;
  vColourMatrixF matF;
  vHalfSpinColourVectorF vec1F;
  vHalfSpinColourVectorF vec2F;
  vHalfSpinColourVectorF vec3F;
  vHalfSpinColourVectorF vecF;
  vHalfSpinColourVectorF matvecF;
  vHalfSpinColourVectorF refF;
  vComplexF errF;

  random(sRNG,matF);
  matF = zero;
  matF()()(0,0)=1.0;
  random(sRNG,vecF);

  refF = matF*vecF;

  WilsonDslashAvx512F((void *)&vecF, (void *)&matF,(void *)&matvecF);
  //std::cout << refF   <<std::endl;
  //std::cout << matvecF<<std::endl;
 
  refF = refF-matvecF;
  errF = TensorRemove(innerProduct(refF,refF));
  std::cout <<"Single SU3 x 2spin diff   "<< Reduce(errF)<<std::endl;

  TimesIAvx512F((void *)&vecF,(void *)&matvecF);
  //std::cout << timesI(vecF)<<std::endl;
  //std::cout << matvecF<<std::endl;
  refF = timesI(vecF)-matvecF;
  errF = TensorRemove(innerProduct(refF,refF));
  std::cout <<" timesI single diff  "<< Reduce(errF)<<std::endl;

  TimesIAvx512((void *)&vec,(void *)&matvec);
  //std::cout << timesI(vec)<<std::endl;
  //std::cout << matvec<<std::endl;
 
  ref = timesI(vec)-matvec;
  err = TensorRemove(innerProduct(ref,ref));
  std::cout <<" timesI double diff  "<< Reduce(err)<<std::endl;

  TimesMinusIAvx512F((void *)&vecF,(void *)&matvecF);
  //std::cout << timesMinusI(vecF)<<std::endl;
  //std::cout << matvecF<<std::endl;
  refF = timesMinusI(vecF)-matvecF;
  errF = TensorRemove(innerProduct(refF,refF));
  std::cout <<" timesMinusI single diff  "<< Reduce(errF)<<std::endl;

  TimesMinusIAvx512((void *)&vec,(void *)&matvec);
  //std::cout << timesMinusI(vec)<<std::endl;
  //std::cout << matvec<<std::endl;

  ref = timesMinusI(vec)-matvec;
  err = TensorRemove(innerProduct(ref,ref));
  std::cout <<" timesMinusI double diff  "<< Reduce(err)<<std::endl;


  LatticeFermion src (FGrid);
  LatticeFermion tmp (FGrid);
  LatticeFermion srce(FrbGrid);

  LatticeFermion resulto(FrbGrid); resulto=zero;
  LatticeFermion resulta(FrbGrid); resulta=zero;
  LatticeFermion diff(FrbGrid); 
  LatticeGaugeField Umu(UGrid);


  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  random(RNG5,src);
#if 1
  random(RNG4,Umu);
#else
  int mmu=2;
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

#if 1
  for(int i=0;i< PerformanceCounter::NumTypes(); i++ ){
    Dw.DhopOE(srce,resulta,0);
    PerformanceCounter Counter(i);
    Counter.Start();
    Dw.DhopOE(srce,resulta,0);
    Counter.Stop();
    Counter.Report();
  }
#endif
  //resulta = (-0.5) * resulta;

  std::cout<<GridLogMessage << "Called Asm Dw"<<std::endl;
  std::cout<<GridLogMessage << "norm result "<< norm2(resulta)<<std::endl;
  std::cout<<GridLogMessage << "mflop/s =   "<< flops*ncall/(t1-t0)<<std::endl;
  diff = resulto-resulta;
  std::cout<<GridLogMessage << "diff "<< norm2(diff)<<std::endl;
  std::cout<<std::endl;
#if 0
  std::cout<<"=========== result Grid ============="<<std::endl;
  std::cout<<std::endl;
  tmp = zero;
  setCheckerboard(tmp,resulto);
  std::cout<<tmp<<std::endl;
  std::cout<<std::endl;
  std::cout<<"=========== result ASM ============="<<std::endl;
  std::cout<<std::endl;
  tmp = zero;
  setCheckerboard(tmp,resulta);
  std::cout<<tmp<<std::endl;
#endif
}

#include <simd/Intel512double.h>

#define zz Z0


void Zmul(void *ptr1,void *ptr2,void *ptr3)
{
  __asm__ ("mov     $0xAAAA, %%eax "  : : :"%eax");
  __asm__ ("kmovw    %%eax, %%k6 " : : :);
  __asm__ ("mov     $0x5555, %%eax "  : : :"%eax");
  __asm__ ("kmovw    %%eax, %%k7 " : : :);

#define CC result_00
  LOAD64(%r9,ptr1);
  LOAD64(%r8,ptr2);
  LOAD64(%r10,ptr3)
  __asm__ (
  VLOAD(0,%r8,CC)
  ZLOAD(0,%r9,Chi_00,Z0) 
  ZMUL(Chi_00,Z0,CC,UChi_00,Z1)
  //VSTORE(0,%r10,UChi_00)
  //VSTORE(1,%r10,Z1)
  ZEND1(UChi_00,Z1,Z0)
  //VSTORE(2,%r10,UChi_00)
  ZEND2(UChi_00,Z1,Z0)
  //VSTORE(3,%r10,UChi_00)
  VSTORE(0,%r10,UChi_00)
  VLOAD(1,%r8,CC)
  ZLOAD(1,%r9,Chi_01,Z0) 
  ZMUL(Chi_01,Z0,CC,UChi_01,Z1)
  ZEND1(UChi_01,Z1,Z0)
  ZEND2(UChi_01,Z1,Z0)
  VSTORE(1,%r10,UChi_01)
  VLOAD(2,%r8,CC)
  ZLOAD(2,%r9,Chi_02,Z0) 
  ZMUL(Chi_02,Z0,CC,UChi_02,Z1)
  ZEND1(UChi_02,Z1,Z0)
  ZEND2(UChi_02,Z1,Z0)
  VSTORE(2,%r10,UChi_02)
  VLOAD(3,%r8,CC)
  ZLOAD(3,%r9,Chi_10,Z0) 
  ZMUL(Chi_10,Z0,CC,UChi_10,Z1)
  ZEND1(UChi_10,Z1,Z0)
  ZEND2(UChi_10,Z1,Z0)
  VSTORE(3,%r10,UChi_10)
  VLOAD(4,%r8,CC)
  ZLOAD(4,%r9,Chi_11,Z0) 
  ZMUL(Chi_11,Z0,CC,UChi_11,Z1)
  ZEND1(UChi_11,Z1,Z0)
  ZEND2(UChi_11,Z1,Z0)
  VSTORE(4,%r10,UChi_11)
  VLOAD(5,%r8,CC)
  ZLOAD(5,%r9,Chi_12,Z0) 
  ZMUL(Chi_12,Z0,CC,UChi_12,Z1)
  ZEND1(UChi_12,Z1,Z0)
  ZEND2(UChi_12,Z1,Z0)
  VSTORE(5,%r10,UChi_12)
  );
}
void TimesMinusIAvx512(void *ptr1,void *ptr3)
{
  __asm__ ("mov     $0xAAAA, %%eax "  : : :"%eax");
  __asm__ ("kmovw    %%eax, %%k6 " : : :);
  __asm__ ("mov     $0x5555, %%eax "  : : :"%eax");
  __asm__ ("kmovw    %%eax, %%k7 " : : :);

  MASK_REGS;

  LOAD_CHI(ptr1);

  __asm__ (
  VZERO(zz)
  VTIMESMINUSI(Chi_00,UChi_00,zz)
  VTIMESMINUSI(Chi_01,UChi_01,zz)
  VTIMESMINUSI(Chi_02,UChi_02,zz)
  VTIMESMINUSI(Chi_10,UChi_10,zz)
  VTIMESMINUSI(Chi_11,UChi_11,zz)
  VTIMESMINUSI(Chi_12,UChi_12,zz)
  );

  SAVE_UCHI(ptr3);
}

void TimesIAvx512(void *ptr1,void *ptr3)
{
  __asm__ ("mov     $0xAAAA, %%eax "  : : :"%eax");
  __asm__ ("kmovw    %%eax, %%k6 " : : :);
  __asm__ ("mov     $0x5555, %%eax "  : : :"%eax");
  __asm__ ("kmovw    %%eax, %%k7 " : : :);

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

#include <simd/Intel512single.h>

void ZmulF(void *ptr1,void *ptr2,void *ptr3)
{
  __asm__ ("mov     $0xAAAA, %%eax "  : : :"%eax");
  __asm__ ("kmovw    %%eax, %%k6 " : : :);
  __asm__ ("mov     $0x5555, %%eax "  : : :"%eax");
  __asm__ ("kmovw    %%eax, %%k7 " : : :);
  MASK_REGS;
  ZLOAD(0,ptr1,Chi_00,Z0);
  ZLOAD(1,ptr1,Chi_01,Z1);
  ZLOAD(2,ptr1,Chi_02,Z2);
  ZLOAD(3,ptr1,Chi_10,Z3);
  ZLOAD(4,ptr1,Chi_11,Z4);
  ZLOAD(5,ptr1,Chi_12,Z5);

  VLOAD(0,ptr2,Chi_20);
  VLOAD(1,ptr2,Chi_21);
  VLOAD(2,ptr2,Chi_22);
  VLOAD(3,ptr2,Chi_30);  
  VLOAD(4,ptr2,Chi_31);  
  VLOAD(5,ptr2,Chi_32);  

  ZMUL(Chi_00,Z0,Chi_20,UChi_00,UChi_20);
  ZMUL(Chi_01,Z1,Chi_21,UChi_01,UChi_21);
  ZMUL(Chi_02,Z2,Chi_22,UChi_02,UChi_22);
  ZMUL(Chi_10,Z3,Chi_23,UChi_10,UChi_30);
  ZMUL(Chi_11,Z4,Chi_24,UChi_11,UChi_31);
  ZMUL(Chi_12,Z5,Chi_25,UChi_12,UChi_32);
  
  ZEND1(UChi_00,UChi_20,Z0);
  ZEND1(UChi_01,UChi_21,Z1);
  ZEND1(UChi_02,UChi_22,Z2);
  ZEND1(UChi_10,UChi_30,Z3);
  ZEND1(UChi_11,UChi_31,Z4);
  ZEND1(UChi_12,UChi_32,Z5);

  ZEND2(UChi_00,UChi_20,Z0);
  ZEND2(UChi_01,UChi_21,Z1);
  ZEND2(UChi_02,UChi_22,Z2);
  ZEND2(UChi_10,UChi_30,Z3);
  ZEND2(UChi_11,UChi_31,Z4);
  ZEND2(UChi_12,UChi_32,Z5);

  SAVE_UCHI(ptr3); 
}

void TimesMinusIAvx512F(void *ptr1,void *ptr3)
{
  MASK_REGS;

  LOAD_CHI(ptr1);
  __asm__ (
  VZERO(zz)
  VTIMESMINUSI(Chi_00,UChi_00,zz)
  VTIMESMINUSI(Chi_01,UChi_01,zz)
  VTIMESMINUSI(Chi_02,UChi_02,zz)
  VTIMESMINUSI(Chi_10,UChi_10,zz)
  VTIMESMINUSI(Chi_11,UChi_11,zz)
  VTIMESMINUSI(Chi_12,UChi_12,zz)
           );
  SAVE_UCHI(ptr3);
}

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

  MULT_ADDSUB_2SPIN(ptr2);
  //MULT_2SPIN(ptr2);

  SAVE_UCHI(ptr3);

  return;
}

