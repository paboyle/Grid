    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/CayleyFermion5D.cc

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
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


#include <Grid/Grid.h>


namespace Grid {
namespace QCD {  /*
   * Dense matrix versions of routines
   */
template<class Impl>
void CayleyFermion5D<Impl>::MooeeInvDag (const FermionField &psi, FermionField &chi)
{
  this->MooeeInternal(psi,chi,DaggerYes,InverseYes);
}
  
template<class Impl>
void CayleyFermion5D<Impl>::MooeeInv(const FermionField &psi, FermionField &chi)
{
  this->MooeeInternal(psi,chi,DaggerNo,InverseYes);
}
template<class Impl>  
void CayleyFermion5D<Impl>::M5D(const FermionField &psi,
				const FermionField &phi, 
				FermionField &chi,
				std::vector<Coeff_t> &lower,
				std::vector<Coeff_t> &diag,
				std::vector<Coeff_t> &upper)
{
  GridBase *grid=psi._grid;
  int Ls   = this->Ls;
  int LLs  = grid->_rdimensions[0];
  const int nsimd= Simd::Nsimd();

  Vector<iSinglet<Simd> > u(LLs);
  Vector<iSinglet<Simd> > l(LLs);
  Vector<iSinglet<Simd> > d(LLs);

  assert(Ls/LLs==nsimd);
  assert(phi.checkerboard == psi.checkerboard);

  chi.checkerboard=psi.checkerboard;

  // just directly address via type pun
  typedef typename Simd::scalar_type scalar_type;
  scalar_type * u_p = (scalar_type *)&u[0];
  scalar_type * l_p = (scalar_type *)&l[0];
  scalar_type * d_p = (scalar_type *)&d[0];

  for(int o=0;o<LLs;o++){ // outer
  for(int i=0;i<nsimd;i++){ //inner
    int s  = o+i*LLs;
    int ss = o*nsimd+i;
    u_p[ss] = upper[s];
    l_p[ss] = lower[s];
    d_p[ss] = diag[s];
  }}


  M5Dcalls++;
  M5Dtime-=usecond();

  assert(Nc==3);

PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss+=LLs){ // adds LLs
#if 0
      alignas(64) SiteHalfSpinor hp;
      alignas(64) SiteHalfSpinor hm;
      alignas(64) SiteSpinor fp;
      alignas(64) SiteSpinor fm;

      for(int v=0;v<LLs;v++){

	int vp=(v+1)%LLs;
	int vm=(v+LLs-1)%LLs;

	spProj5m(hp,psi[ss+vp]);
	spProj5p(hm,psi[ss+vm]);

	if ( vp<=v ) rotate(hp,hp,1);
	if ( vm>=v ) rotate(hm,hm,nsimd-1);
	
	hp=0.5*hp;
        hm=0.5*hm;

	spRecon5m(fp,hp);
	spRecon5p(fm,hm);

	chi[ss+v] = d[v]*phi[ss+v];
	chi[ss+v] = chi[ss+v]     +u[v]*fp;
	chi[ss+v] = chi[ss+v]     +l[v]*fm;

      }
#else
      for(int v=0;v<LLs;v++){

	vprefetch(psi[ss+v+LLs]);

	int vp= (v==LLs-1) ? 0     : v+1;
	int vm= (v==0    ) ? LLs-1 : v-1;
	
	Simd hp_00 = psi[ss+vp]()(2)(0); 
	Simd hp_01 = psi[ss+vp]()(2)(1); 
	Simd hp_02 = psi[ss+vp]()(2)(2); 
	Simd hp_10 = psi[ss+vp]()(3)(0); 
	Simd hp_11 = psi[ss+vp]()(3)(1); 
	Simd hp_12 = psi[ss+vp]()(3)(2); 
	
	Simd hm_00 = psi[ss+vm]()(0)(0); 
	Simd hm_01 = psi[ss+vm]()(0)(1); 
	Simd hm_02 = psi[ss+vm]()(0)(2); 
	Simd hm_10 = psi[ss+vm]()(1)(0); 
	Simd hm_11 = psi[ss+vm]()(1)(1); 
	Simd hm_12 = psi[ss+vm]()(1)(2); 

	if ( vp<=v ) {
	  hp_00.v = Optimization::Rotate::tRotate<2>(hp_00.v);
	  hp_01.v = Optimization::Rotate::tRotate<2>(hp_01.v);
	  hp_02.v = Optimization::Rotate::tRotate<2>(hp_02.v);
	  hp_10.v = Optimization::Rotate::tRotate<2>(hp_10.v);
	  hp_11.v = Optimization::Rotate::tRotate<2>(hp_11.v);
	  hp_12.v = Optimization::Rotate::tRotate<2>(hp_12.v);
	}
	if ( vm>=v ) {
	  hm_00.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_00.v);
	  hm_01.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_01.v);
	  hm_02.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_02.v);
	  hm_10.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_10.v);
	  hm_11.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_11.v);
	  hm_12.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_12.v);
	}

	// Can force these to real arithmetic and save 2x.
	Simd p_00  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(0)(0))  + switcheroo<Coeff_t>::mult(l[v]()()(),hm_00); 
	Simd p_01  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(0)(1))  + switcheroo<Coeff_t>::mult(l[v]()()(),hm_01); 
	Simd p_02  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(0)(2))  + switcheroo<Coeff_t>::mult(l[v]()()(),hm_02); 
	Simd p_10  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(1)(0))  + switcheroo<Coeff_t>::mult(l[v]()()(),hm_10); 
	Simd p_11  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(1)(1))  + switcheroo<Coeff_t>::mult(l[v]()()(),hm_11); 
	Simd p_12  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(1)(2))  + switcheroo<Coeff_t>::mult(l[v]()()(),hm_12); 
	Simd p_20  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(2)(0))  + switcheroo<Coeff_t>::mult(u[v]()()(),hp_00); 
	Simd p_21  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(2)(1))  + switcheroo<Coeff_t>::mult(u[v]()()(),hp_01); 
	Simd p_22  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(2)(2))  + switcheroo<Coeff_t>::mult(u[v]()()(),hp_02);  
	Simd p_30  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(3)(0))  + switcheroo<Coeff_t>::mult(u[v]()()(),hp_10); 
	Simd p_31  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(3)(1))  + switcheroo<Coeff_t>::mult(u[v]()()(),hp_11); 
	Simd p_32  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(3)(2))  + switcheroo<Coeff_t>::mult(u[v]()()(),hp_12); 

	vstream(chi[ss+v]()(0)(0),p_00);
	vstream(chi[ss+v]()(0)(1),p_01);
	vstream(chi[ss+v]()(0)(2),p_02);
	vstream(chi[ss+v]()(1)(0),p_10);
	vstream(chi[ss+v]()(1)(1),p_11);
	vstream(chi[ss+v]()(1)(2),p_12);
	vstream(chi[ss+v]()(2)(0),p_20);
	vstream(chi[ss+v]()(2)(1),p_21);
	vstream(chi[ss+v]()(2)(2),p_22);
	vstream(chi[ss+v]()(3)(0),p_30);
	vstream(chi[ss+v]()(3)(1),p_31);
	vstream(chi[ss+v]()(3)(2),p_32);

      }
#endif
  }
  M5Dtime+=usecond();
}

template<class Impl>  
void CayleyFermion5D<Impl>::M5Ddag(const FermionField &psi,
				   const FermionField &phi, 
				   FermionField &chi,
				   std::vector<Coeff_t> &lower,
				   std::vector<Coeff_t> &diag,
				   std::vector<Coeff_t> &upper)
{
  GridBase *grid=psi._grid;
  int Ls   = this->Ls;
  int LLs  = grid->_rdimensions[0];
  int nsimd= Simd::Nsimd();

  Vector<iSinglet<Simd> > u(LLs);
  Vector<iSinglet<Simd> > l(LLs);
  Vector<iSinglet<Simd> > d(LLs);

  assert(Ls/LLs==nsimd);
  assert(phi.checkerboard == psi.checkerboard);

  chi.checkerboard=psi.checkerboard;

  // just directly address via type pun
  typedef typename Simd::scalar_type scalar_type;
  scalar_type * u_p = (scalar_type *)&u[0];
  scalar_type * l_p = (scalar_type *)&l[0];
  scalar_type * d_p = (scalar_type *)&d[0];

  for(int o=0;o<LLs;o++){ // outer
  for(int i=0;i<nsimd;i++){ //inner
    int s  = o+i*LLs;
    int ss = o*nsimd+i;
    u_p[ss] = upper[s];
    l_p[ss] = lower[s];
    d_p[ss] = diag[s];
  }}

  M5Dcalls++;
  M5Dtime-=usecond();
PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss+=LLs){ // adds LLs
#if 0
    alignas(64) SiteHalfSpinor hp;
    alignas(64) SiteHalfSpinor hm;
    alignas(64) SiteSpinor fp;
    alignas(64) SiteSpinor fm;

    for(int v=0;v<LLs;v++){

      int vp=(v+1)%LLs;
      int vm=(v+LLs-1)%LLs;

      spProj5p(hp,psi[ss+vp]);
      spProj5m(hm,psi[ss+vm]);

      if ( vp<=v ) rotate(hp,hp,1);
      if ( vm>=v ) rotate(hm,hm,nsimd-1);
      
      hp=hp*0.5;
      hm=hm*0.5;
      spRecon5p(fp,hp);
      spRecon5m(fm,hm);

      chi[ss+v] = d[v]*phi[ss+v]+u[v]*fp;
      chi[ss+v] = chi[ss+v]     +l[v]*fm;

    }
#else
      for(int v=0;v<LLs;v++){

	vprefetch(psi[ss+v+LLs]);

	int vp= (v==LLs-1) ? 0     : v+1;
	int vm= (v==0    ) ? LLs-1 : v-1;
	
	Simd hp_00 = psi[ss+vp]()(0)(0); 
	Simd hp_01 = psi[ss+vp]()(0)(1); 
	Simd hp_02 = psi[ss+vp]()(0)(2); 
	Simd hp_10 = psi[ss+vp]()(1)(0); 
	Simd hp_11 = psi[ss+vp]()(1)(1); 
	Simd hp_12 = psi[ss+vp]()(1)(2); 
	
	Simd hm_00 = psi[ss+vm]()(2)(0); 
	Simd hm_01 = psi[ss+vm]()(2)(1); 
	Simd hm_02 = psi[ss+vm]()(2)(2); 
	Simd hm_10 = psi[ss+vm]()(3)(0); 
	Simd hm_11 = psi[ss+vm]()(3)(1); 
	Simd hm_12 = psi[ss+vm]()(3)(2); 

	if ( vp<=v ) {
	  hp_00.v = Optimization::Rotate::tRotate<2>(hp_00.v);
	  hp_01.v = Optimization::Rotate::tRotate<2>(hp_01.v);
	  hp_02.v = Optimization::Rotate::tRotate<2>(hp_02.v);
	  hp_10.v = Optimization::Rotate::tRotate<2>(hp_10.v);
	  hp_11.v = Optimization::Rotate::tRotate<2>(hp_11.v);
	  hp_12.v = Optimization::Rotate::tRotate<2>(hp_12.v);
	}
	if ( vm>=v ) {
	  hm_00.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_00.v);
	  hm_01.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_01.v);
	  hm_02.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_02.v);
	  hm_10.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_10.v);
	  hm_11.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_11.v);
	  hm_12.v = Optimization::Rotate::tRotate<2*Simd::Nsimd()-2>(hm_12.v);
	}

	Simd p_00  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(0)(0))  + switcheroo<Coeff_t>::mult(u[v]()()(),hp_00); 
	Simd p_01  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(0)(1))  + switcheroo<Coeff_t>::mult(u[v]()()(),hp_01); 
	Simd p_02  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(0)(2))  + switcheroo<Coeff_t>::mult(u[v]()()(),hp_02); 
	Simd p_10  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(1)(0))  + switcheroo<Coeff_t>::mult(u[v]()()(),hp_10); 
	Simd p_11  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(1)(1))  + switcheroo<Coeff_t>::mult(u[v]()()(),hp_11); 
	Simd p_12  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(1)(2))  + switcheroo<Coeff_t>::mult(u[v]()()(),hp_12); 

	Simd p_20  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(2)(0))  + switcheroo<Coeff_t>::mult(l[v]()()(),hm_00); 
	Simd p_21  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(2)(1))  + switcheroo<Coeff_t>::mult(l[v]()()(),hm_01); 
	Simd p_22  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(2)(2))  + switcheroo<Coeff_t>::mult(l[v]()()(),hm_02);  
	Simd p_30  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(3)(0))  + switcheroo<Coeff_t>::mult(l[v]()()(),hm_10); 
	Simd p_31  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(3)(1))  + switcheroo<Coeff_t>::mult(l[v]()()(),hm_11); 
	Simd p_32  = switcheroo<Coeff_t>::mult(d[v]()()(), phi[ss+v]()(3)(2))  + switcheroo<Coeff_t>::mult(l[v]()()(),hm_12); 

	vstream(chi[ss+v]()(0)(0),p_00);
	vstream(chi[ss+v]()(0)(1),p_01);
	vstream(chi[ss+v]()(0)(2),p_02);
	vstream(chi[ss+v]()(1)(0),p_10);
	vstream(chi[ss+v]()(1)(1),p_11);
	vstream(chi[ss+v]()(1)(2),p_12);
	vstream(chi[ss+v]()(2)(0),p_20);
	vstream(chi[ss+v]()(2)(1),p_21);
	vstream(chi[ss+v]()(2)(2),p_22);
	vstream(chi[ss+v]()(3)(0),p_30);
	vstream(chi[ss+v]()(3)(1),p_31);
	vstream(chi[ss+v]()(3)(2),p_32);
      }
#endif
  }
  M5Dtime+=usecond();
}


#ifdef AVX512 
#include <simd/Intel512common.h>
#include <simd/Intel512avx.h>
#include <simd/Intel512single.h>
#endif 

template<class Impl>
void CayleyFermion5D<Impl>::MooeeInternalAsm(const FermionField &psi, FermionField &chi,
					     int LLs, int site,
					     Vector<iSinglet<Simd> > &Matp,
					     Vector<iSinglet<Simd> > &Matm)
{
#ifndef AVX512
  {
  SiteHalfSpinor BcastP;
  SiteHalfSpinor BcastM;
  SiteHalfSpinor SiteChiP;
  SiteHalfSpinor SiteChiM;

  // Ls*Ls * 2 * 12 * vol flops
  for(int s1=0;s1<LLs;s1++){ 
    for(int s2=0;s2<LLs;s2++){ 
      for(int  l=0; l<Simd::Nsimd();l++){ // simd lane

        int s=s2+l*LLs;
	int lex=s2+LLs*site;
	
	if ( s2==0 && l==0) {
	  SiteChiP=zero;
	  SiteChiM=zero;
	}
	
	for(int sp=0;sp<2;sp++){
        for(int co=0;co<Nc;co++){
	  vbroadcast(BcastP()(sp  )(co),psi[lex]()(sp)(co),l);
	}}
	for(int sp=0;sp<2;sp++){
        for(int co=0;co<Nc;co++){
	  vbroadcast(BcastM()(sp  )(co),psi[lex]()(sp+2)(co),l);
	}}

	for(int sp=0;sp<2;sp++){
        for(int co=0;co<Nc;co++){
	  SiteChiP()(sp)(co)=real_madd(Matp[LLs*s+s1]()()(),BcastP()(sp)(co),SiteChiP()(sp)(co)); // 1100 us.
	  SiteChiM()(sp)(co)=real_madd(Matm[LLs*s+s1]()()(),BcastM()(sp)(co),SiteChiM()(sp)(co)); // each found by commenting out
	}}

    }}
    {
      int lex = s1+LLs*site;
      for(int sp=0;sp<2;sp++){
      for(int co=0;co<Nc;co++){
	vstream(chi[lex]()(sp)(co), SiteChiP()(sp)(co));
	vstream(chi[lex]()(sp+2)(co), SiteChiM()(sp)(co));
      }}
    }
  }

  }
#else
  {
  // pointers
    //  MASK_REGS;
#define Chi_00 %%zmm1
#define Chi_01 %%zmm2
#define Chi_02 %%zmm3
#define Chi_10 %%zmm4
#define Chi_11 %%zmm5
#define Chi_12 %%zmm6
#define Chi_20 %%zmm7
#define Chi_21 %%zmm8
#define Chi_22 %%zmm9
#define Chi_30 %%zmm10
#define Chi_31 %%zmm11
#define Chi_32 %%zmm12

#define BCAST0   %%zmm13
#define BCAST1   %%zmm14
#define BCAST2   %%zmm15
#define BCAST3   %%zmm16
#define BCAST4   %%zmm17
#define BCAST5   %%zmm18
#define BCAST6   %%zmm19
#define BCAST7   %%zmm20
#define BCAST8   %%zmm21
#define BCAST9   %%zmm22
#define BCAST10  %%zmm23
#define BCAST11  %%zmm24

  int incr=LLs*LLs*sizeof(iSinglet<Simd>);
  for(int s1=0;s1<LLs;s1++){ 
    for(int s2=0;s2<LLs;s2++){ 
      int lex=s2+LLs*site;
      uint64_t a0 = (uint64_t)&Matp[LLs*s2+s1]; // should be cacheable
      uint64_t a1 = (uint64_t)&Matm[LLs*s2+s1];
      uint64_t a2 = (uint64_t)&psi[lex];
      for(int  l=0; l<Simd::Nsimd();l++){ // simd lane
	if ( (s2+l)==0 ) {
	  asm (
  	           VPREFETCH1(0,%2)  	     VPREFETCH1(0,%1)
  	           VPREFETCH1(12,%2)  	     VPREFETCH1(13,%2)
  	           VPREFETCH1(14,%2)  	     VPREFETCH1(15,%2)         
		   VBCASTCDUP(0,%2,BCAST0)   
		   VBCASTCDUP(1,%2,BCAST1)   
		   VBCASTCDUP(2,%2,BCAST2)   
		   VBCASTCDUP(3,%2,BCAST3)   
		   VBCASTCDUP(4,%2,BCAST4)     VMULMEM (0,%0,BCAST0,Chi_00)
		   VBCASTCDUP(5,%2,BCAST5)     VMULMEM (0,%0,BCAST1,Chi_01)
		   VBCASTCDUP(6,%2,BCAST6)     VMULMEM (0,%0,BCAST2,Chi_02)
		   VBCASTCDUP(7,%2,BCAST7)     VMULMEM (0,%0,BCAST3,Chi_10)
		   VBCASTCDUP(8,%2,BCAST8)     VMULMEM (0,%0,BCAST4,Chi_11)
		   VBCASTCDUP(9,%2,BCAST9)     VMULMEM (0,%0,BCAST5,Chi_12)
		   VBCASTCDUP(10,%2,BCAST10)   VMULMEM (0,%1,BCAST6,Chi_20)
		   VBCASTCDUP(11,%2,BCAST11)   VMULMEM (0,%1,BCAST7,Chi_21)
		   VMULMEM (0,%1,BCAST8,Chi_22)         
		   VMULMEM (0,%1,BCAST9,Chi_30)
		   VMULMEM (0,%1,BCAST10,Chi_31)       
		   VMULMEM (0,%1,BCAST11,Chi_32)
		   : : "r" (a0), "r" (a1), "r" (a2)  );
	} else { 
	  asm (
		   VBCASTCDUP(0,%2,BCAST0)   VMADDMEM (0,%0,BCAST0,Chi_00)
		   VBCASTCDUP(1,%2,BCAST1)   VMADDMEM (0,%0,BCAST1,Chi_01)
		   VBCASTCDUP(2,%2,BCAST2)   VMADDMEM (0,%0,BCAST2,Chi_02)
		   VBCASTCDUP(3,%2,BCAST3)   VMADDMEM (0,%0,BCAST3,Chi_10)
		   VBCASTCDUP(4,%2,BCAST4)   VMADDMEM (0,%0,BCAST4,Chi_11)
		   VBCASTCDUP(5,%2,BCAST5)   VMADDMEM (0,%0,BCAST5,Chi_12)
		   VBCASTCDUP(6,%2,BCAST6)   VMADDMEM (0,%1,BCAST6,Chi_20)
		   VBCASTCDUP(7,%2,BCAST7)   VMADDMEM (0,%1,BCAST7,Chi_21)
		   VBCASTCDUP(8,%2,BCAST8)   VMADDMEM (0,%1,BCAST8,Chi_22)
		   VBCASTCDUP(9,%2,BCAST9)   VMADDMEM (0,%1,BCAST9,Chi_30)
		   VBCASTCDUP(10,%2,BCAST10)  VMADDMEM (0,%1,BCAST10,Chi_31)
		   VBCASTCDUP(11,%2,BCAST11)  VMADDMEM (0,%1,BCAST11,Chi_32) 
		   : : "r" (a0), "r" (a1), "r" (a2)  );
	}
	a0 = a0+incr;
	a1 = a1+incr;
	a2 = a2+sizeof(Simd::scalar_type);
      }}
    {
      int lexa = s1+LLs*site;
      asm (
	       VSTORE(0,%0,Chi_00) VSTORE(1 ,%0,Chi_01)  VSTORE(2 ,%0,Chi_02)		
	       VSTORE(3,%0,Chi_10) VSTORE(4 ,%0,Chi_11)  VSTORE(5 ,%0,Chi_12)		
	       VSTORE(6,%0,Chi_20) VSTORE(7 ,%0,Chi_21)  VSTORE(8 ,%0,Chi_22)		
	       VSTORE(9,%0,Chi_30) VSTORE(10,%0,Chi_31)  VSTORE(11,%0,Chi_32)		
	       : : "r" ((uint64_t)&chi[lexa]) : "memory" );

    }
  }
  }
#undef Chi_00
#undef Chi_01
#undef Chi_02
#undef Chi_10
#undef Chi_11
#undef Chi_12
#undef Chi_20
#undef Chi_21
#undef Chi_22
#undef Chi_30
#undef Chi_31
#undef Chi_32

#undef BCAST0
#undef BCAST1
#undef BCAST2
#undef BCAST3
#undef BCAST4
#undef BCAST5
#undef BCAST6
#undef BCAST7
#undef BCAST8
#undef BCAST9
#undef BCAST10
#undef BCAST11
#endif
};

  // Z-mobius version
template<class Impl>
void CayleyFermion5D<Impl>::MooeeInternalZAsm(const FermionField &psi, FermionField &chi,
					     int LLs, int site, Vector<iSinglet<Simd> > &Matp, Vector<iSinglet<Simd> > &Matm)
{
#ifndef AVX512
  {
  SiteHalfSpinor BcastP;
  SiteHalfSpinor BcastM;
  SiteHalfSpinor SiteChiP;
  SiteHalfSpinor SiteChiM;

  // Ls*Ls * 2 * 12 * vol flops
  for(int s1=0;s1<LLs;s1++){ 
    for(int s2=0;s2<LLs;s2++){ 
      for(int  l=0; l<Simd::Nsimd();l++){ // simd lane

        int s=s2+l*LLs;
	int lex=s2+LLs*site;
	
	if ( s2==0 && l==0) {
	  SiteChiP=zero;
	  SiteChiM=zero;
	}
	
	for(int sp=0;sp<2;sp++){
        for(int co=0;co<Nc;co++){
	  vbroadcast(BcastP()(sp  )(co),psi[lex]()(sp)(co),l);
	}}
	for(int sp=0;sp<2;sp++){
        for(int co=0;co<Nc;co++){
	  vbroadcast(BcastM()(sp  )(co),psi[lex]()(sp+2)(co),l);
	}}

	for(int sp=0;sp<2;sp++){
        for(int co=0;co<Nc;co++){
	  SiteChiP()(sp)(co)=SiteChiP()(sp)(co)+ Matp[LLs*s+s1]()()()*BcastP()(sp)(co); 
	  SiteChiM()(sp)(co)=SiteChiM()(sp)(co)+ Matm[LLs*s+s1]()()()*BcastM()(sp)(co); 
	}}


    }}
    {
      int lex = s1+LLs*site;
      for(int sp=0;sp<2;sp++){
      for(int co=0;co<Nc;co++){
	vstream(chi[lex]()(sp)(co), SiteChiP()(sp)(co));
	vstream(chi[lex]()(sp+2)(co), SiteChiM()(sp)(co));
      }}
    }
  }

  }
#else
  {
  // pointers
  //  MASK_REGS;
#define Chi_00 %zmm0
#define Chi_01 %zmm1
#define Chi_02 %zmm2
#define Chi_10 %zmm3
#define Chi_11 %zmm4
#define Chi_12 %zmm5
#define Chi_20 %zmm6
#define Chi_21 %zmm7
#define Chi_22 %zmm8
#define Chi_30 %zmm9
#define Chi_31 %zmm10
#define Chi_32 %zmm11
#define pChi_00 %%zmm0
#define pChi_01 %%zmm1
#define pChi_02 %%zmm2
#define pChi_10 %%zmm3
#define pChi_11 %%zmm4
#define pChi_12 %%zmm5
#define pChi_20 %%zmm6
#define pChi_21 %%zmm7
#define pChi_22 %%zmm8
#define pChi_30 %%zmm9
#define pChi_31 %%zmm10
#define pChi_32 %%zmm11

#define BCAST_00   %zmm12
#define  SHUF_00   %zmm13
#define BCAST_01   %zmm14
#define  SHUF_01   %zmm15
#define BCAST_02   %zmm16
#define  SHUF_02   %zmm17
#define BCAST_10   %zmm18
#define  SHUF_10   %zmm19
#define BCAST_11   %zmm20
#define  SHUF_11   %zmm21
#define BCAST_12   %zmm22
#define  SHUF_12   %zmm23

#define Mp  %zmm24
#define Mps %zmm25
#define Mm  %zmm26
#define Mms %zmm27
#define N 8
  int incr=LLs*LLs*sizeof(iSinglet<Simd>);
  for(int s1=0;s1<LLs;s1++){ 
    for(int s2=0;s2<LLs;s2++){ 
      int lex=s2+LLs*site;
      uint64_t a0 = (uint64_t)&Matp[LLs*s2+s1]; // should be cacheable
      uint64_t a1 = (uint64_t)&Matm[LLs*s2+s1];
      uint64_t a2 = (uint64_t)&psi[lex];
      for(int  l=0; l<Simd::Nsimd();l++){ // simd lane
	if ( (s2+l)==0 ) {
	  LOAD64(%r8,a0);
	  LOAD64(%r9,a1);
	  LOAD64(%r10,a2);
	  asm (
	       VLOAD(0,%r8,Mp)// i r
	       VLOAD(0,%r9,Mm)
	       VSHUF(Mp,Mps)  // r i 
	       VSHUF(Mm,Mms)
	       VPREFETCH1(12,%r10)  	     VPREFETCH1(13,%r10)
	       VPREFETCH1(14,%r10)  	     VPREFETCH1(15,%r10)         

	       VMULIDUP(0*N,%r10,Mps,Chi_00)
	       VMULIDUP(1*N,%r10,Mps,Chi_01)
	       VMULIDUP(2*N,%r10,Mps,Chi_02)
	       VMULIDUP(3*N,%r10,Mps,Chi_10)
	       VMULIDUP(4*N,%r10,Mps,Chi_11)
	       VMULIDUP(5*N,%r10,Mps,Chi_12)

	       VMULIDUP(6*N ,%r10,Mms,Chi_20)
	       VMULIDUP(7*N ,%r10,Mms,Chi_21)
	       VMULIDUP(8*N ,%r10,Mms,Chi_22)
	       VMULIDUP(9*N ,%r10,Mms,Chi_30)
	       VMULIDUP(10*N,%r10,Mms,Chi_31)
	       VMULIDUP(11*N,%r10,Mms,Chi_32)

	       VMADDSUBRDUP(0*N,%r10,Mp,Chi_00)
	       VMADDSUBRDUP(1*N,%r10,Mp,Chi_01)
	       VMADDSUBRDUP(2*N,%r10,Mp,Chi_02)
	       VMADDSUBRDUP(3*N,%r10,Mp,Chi_10)
	       VMADDSUBRDUP(4*N,%r10,Mp,Chi_11)
	       VMADDSUBRDUP(5*N,%r10,Mp,Chi_12)

	       VMADDSUBRDUP(6*N ,%r10,Mm,Chi_20)
	       VMADDSUBRDUP(7*N ,%r10,Mm,Chi_21)
	       VMADDSUBRDUP(8*N ,%r10,Mm,Chi_22)
	       VMADDSUBRDUP(9*N ,%r10,Mm,Chi_30)
	       VMADDSUBRDUP(10*N,%r10,Mm,Chi_31)
	       VMADDSUBRDUP(11*N,%r10,Mm,Chi_32)
	       );
	} else { 
	  LOAD64(%r8,a0);
	  LOAD64(%r9,a1);
	  LOAD64(%r10,a2);
	  asm (
	       VLOAD(0,%r8,Mp)
	       VSHUF(Mp,Mps)

	       VLOAD(0,%r9,Mm)
	       VSHUF(Mm,Mms)

	       VMADDSUBIDUP(0*N,%r10,Mps,Chi_00) //  Mri * Pii +- Cir
	       VMADDSUBIDUP(1*N,%r10,Mps,Chi_01)
	       VMADDSUBIDUP(2*N,%r10,Mps,Chi_02)
	       VMADDSUBIDUP(3*N,%r10,Mps,Chi_10)
	       VMADDSUBIDUP(4*N,%r10,Mps,Chi_11)
	       VMADDSUBIDUP(5*N,%r10,Mps,Chi_12)

	       VMADDSUBIDUP(6 *N,%r10,Mms,Chi_20)
	       VMADDSUBIDUP(7 *N,%r10,Mms,Chi_21)
	       VMADDSUBIDUP(8 *N,%r10,Mms,Chi_22)
	       VMADDSUBIDUP(9 *N,%r10,Mms,Chi_30)
	       VMADDSUBIDUP(10*N,%r10,Mms,Chi_31)
	       VMADDSUBIDUP(11*N,%r10,Mms,Chi_32)

	       VMADDSUBRDUP(0*N,%r10,Mp,Chi_00) //  Cir = Mir * Prr +- ( Mri * Pii +- Cir) 
	       VMADDSUBRDUP(1*N,%r10,Mp,Chi_01) //  Ci = MiPr + Ci + MrPi ;    Cr = MrPr - ( MiPi - Cr)
	       VMADDSUBRDUP(2*N,%r10,Mp,Chi_02)
	       VMADDSUBRDUP(3*N,%r10,Mp,Chi_10)
	       VMADDSUBRDUP(4*N,%r10,Mp,Chi_11)
	       VMADDSUBRDUP(5*N,%r10,Mp,Chi_12)

	       VMADDSUBRDUP(6 *N,%r10,Mm,Chi_20)
	       VMADDSUBRDUP(7 *N,%r10,Mm,Chi_21)
	       VMADDSUBRDUP(8 *N,%r10,Mm,Chi_22)
	       VMADDSUBRDUP(9 *N,%r10,Mm,Chi_30)
	       VMADDSUBRDUP(10*N,%r10,Mm,Chi_31)
	       VMADDSUBRDUP(11*N,%r10,Mm,Chi_32)
	       );
	}
	a0 = a0+incr;
	a1 = a1+incr;
	a2 = a2+sizeof(Simd::scalar_type);
      }}
    {
      int lexa = s1+LLs*site;
      /*
      SiteSpinor tmp;
      asm (
	       VSTORE(0,%0,pChi_00) VSTORE(1 ,%0,pChi_01)  VSTORE(2 ,%0,pChi_02)		
	       VSTORE(3,%0,pChi_10) VSTORE(4 ,%0,pChi_11)  VSTORE(5 ,%0,pChi_12)		
	       VSTORE(6,%0,pChi_20) VSTORE(7 ,%0,pChi_21)  VSTORE(8 ,%0,pChi_22)		
	       VSTORE(9,%0,pChi_30) VSTORE(10,%0,pChi_31)  VSTORE(11,%0,pChi_32)		
	       : : "r" ((uint64_t)&tmp) : "memory" );
      */

      asm (
	       VSTORE(0,%0,pChi_00) VSTORE(1 ,%0,pChi_01)  VSTORE(2 ,%0,pChi_02)		
	       VSTORE(3,%0,pChi_10) VSTORE(4 ,%0,pChi_11)  VSTORE(5 ,%0,pChi_12)		
	       VSTORE(6,%0,pChi_20) VSTORE(7 ,%0,pChi_21)  VSTORE(8 ,%0,pChi_22)		
	       VSTORE(9,%0,pChi_30) VSTORE(10,%0,pChi_31)  VSTORE(11,%0,pChi_32)		
	       : : "r" ((uint64_t)&chi[lexa]) : "memory" );

      //      if ( 1 || (site==0) ) { 
      //	std::cout<<site << " s1 "<<s1<<"\n\t"<<tmp << "\n't" << chi[lexa] <<"\n\t"<<tmp-chi[lexa]<<std::endl;
      //      }
    }
  }
  }
#undef Chi_00
#undef Chi_01
#undef Chi_02
#undef Chi_10
#undef Chi_11
#undef Chi_12
#undef Chi_20
#undef Chi_21
#undef Chi_22
#undef Chi_30
#undef Chi_31
#undef Chi_32

#undef BCAST0
#undef BCAST1
#undef BCAST2
#undef BCAST3
#undef BCAST4
#undef BCAST5
#undef BCAST6
#undef BCAST7
#undef BCAST8
#undef BCAST9
#undef BCAST10
#undef BCAST11

#endif
};


template<class Impl>
void CayleyFermion5D<Impl>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv)
{
  int Ls=this->Ls;
  int LLs = psi._grid->_rdimensions[0];
  int vol = psi._grid->oSites()/LLs;

  chi.checkerboard=psi.checkerboard;
  
  Vector<iSinglet<Simd> >  Matp;
  Vector<iSinglet<Simd> >  Matm;
  Vector<iSinglet<Simd> >  *_Matp;
  Vector<iSinglet<Simd> >  *_Matm;
  
  //  MooeeInternalCompute(dag,inv,Matp,Matm);
  if ( inv && dag ) { 
    _Matp = &MatpInvDag;
    _Matm = &MatmInvDag;
  }
  if ( inv && (!dag) ) { 
    _Matp = &MatpInv;
    _Matm = &MatmInv;
  } 
  if ( !inv ) {
    MooeeInternalCompute(dag,inv,Matp,Matm);
    _Matp = &Matp;
    _Matm = &Matm;
  }
  assert(_Matp->size()==Ls*LLs);

  MooeeInvCalls++;
  MooeeInvTime-=usecond();

  if ( switcheroo<Coeff_t>::iscomplex() ) {
  PARALLEL_FOR_LOOP
    for(auto site=0;site<vol;site++){
      MooeeInternalZAsm(psi,chi,LLs,site,*_Matp,*_Matm);
    }
  } else { 
  PARALLEL_FOR_LOOP
    for(auto site=0;site<vol;site++){
      MooeeInternalAsm(psi,chi,LLs,site,*_Matp,*_Matm);
    }
  }
  MooeeInvTime+=usecond();
}

INSTANTIATE_DPERP(DomainWallVec5dImplD);
INSTANTIATE_DPERP(DomainWallVec5dImplF);
INSTANTIATE_DPERP(ZDomainWallVec5dImplD);
INSTANTIATE_DPERP(ZDomainWallVec5dImplF);

template void CayleyFermion5D<DomainWallVec5dImplF>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<DomainWallVec5dImplD>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<ZDomainWallVec5dImplF>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<ZDomainWallVec5dImplD>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);


}}
