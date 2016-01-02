    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/WilsonKernelsAsm.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#if defined(AVX512) || defined (IMCI)

#include <simd/Avx512Asm.h>

#undef VLOAD
#undef VSTORE
#undef VMUL
#undef VMADD
#undef ZEND
#undef ZLOAD
#undef ZMUL
#undef ZMADD
#undef VZERO
#undef VTIMESI
#undef VTIMESMINUSI

#define VZERO(A)                  VZEROf(A)
#define VMOV(A,B)                 VMOVf(A,B)
#define VLOAD(OFF,PTR,DEST)       VLOADf(OFF,PTR,DEST)
#define VSTORE(OFF,PTR,SRC)       VSTOREf(OFF,PTR,SRC)

#define VADD(A,B,C)               VADDf(A,B,C)
#define VSUB(A,B,C)               VSUBf(A,B,C)
#define VMUL(Uri,Uir,Chi,UChi,Z)  VMULf(Uri,Uir,Chi,UChi,Z)
#define VMADD(Uri,Uir,Chi,UChi,Z) VMADDf(Uri,Uir,Chi,UChi,Z)

#define VTIMESI(A,B,C)            VTIMESIf(A,B,C)
#define VTIMESMINUSI(A,B,C)       VTIMESMINUSIf(A,B,C)
#define VACCTIMESI(A,B,C)         VACCTIMESIf(A,B,C)
#define VACCTIMESMINUSI(A,B,C)    VACCTIMESMINUSIf(A,B,C)

#define VTIMESI0(A,B,C)            VTIMESI0f(A,B,C)
#define VTIMESMINUSI0(A,B,C)       VTIMESMINUSI0f(A,B,C)
#define VACCTIMESI0(A,B,C)         VACCTIMESI0f(A,B,C)
#define VACCTIMESMINUSI0(A,B,C)    VACCTIMESMINUSI0f(A,B,C)

#define VTIMESI1(A,B,C)            VTIMESI1f(A,B,C)
#define VTIMESMINUSI1(A,B,C)       VTIMESMINUSI1f(A,B,C)
#define VACCTIMESI1(A,B,C)         VACCTIMESI1f(A,B,C)
#define VACCTIMESMINUSI1(A,B,C)    VACCTIMESMINUSI1f(A,B,C)

#define VTIMESI2(A,B,C)            VTIMESI2f(A,B,C)
#define VTIMESMINUSI2(A,B,C)       VTIMESMINUSI2f(A,B,C)
#define VACCTIMESI2(A,B,C)         VACCTIMESI2f(A,B,C)
#define VACCTIMESMINUSI2(A,B,C)    VACCTIMESMINUSI2f(A,B,C)

#define VACCTIMESI1MEM(A,ACC,O,P) VACCTIMESI1MEMf(A,ACC,O,P)
#define VACCTIMESI2MEM(A,ACC,O,P) VACCTIMESI2MEMf(A,ACC,O,P)
#define VACCTIMESMINUSI1MEM(A,ACC,O,P) VACCTIMESMINUSI1MEMf(A,ACC,O,P)
#define VACCTIMESMINUSI2MEM(A,ACC,O,P) VACCTIMESMINUSI2MEMf(A,ACC,O,P)

#define VPERM0(A,B)               VPERM0f(A,B)
#define VPERM1(A,B)               VPERM1f(A,B)
#define VPERM2(A,B)               VPERM2f(A,B)
#define VPERM3(A,B)               VPERM3f(A,B)
#define VSHUFMEM(OFF,A,DEST)      VSHUFMEMf(OFF,A,DEST)

#define ZEND1(A,B,C)               ZEND1f(A,B,C)
#define ZEND2(A,B,C)               ZEND2f(A,B,C)
#define ZLOAD(A,B,C,D)            ZLOADf(A,B,C,D)
#define ZMUL(A,B,C,D,E)           ZMULf(A,B,C,D,E)
#define ZMADD(A,B,C,D,E)          ZMADDf(A,B,C,D,E)

#define ZMUL(A,B,C,D,E)           ZMULf(A,B,C,D,E)
#define ZMADD(A,B,C,D,E)          ZMADDf(A,B,C,D,E)

#define VADDMEM(O,A,B,C)            VADDMEMf(O,A,B,C)
#define VSUBMEM(O,A,B,C)            VSUBMEMf(O,A,B,C)

#define ZMULMEM2SP(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr)  ZMULMEM2SPf(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr) 
#define ZMADDMEM2SP(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr) ZMADDMEM2SPf(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr) 

namespace Grid {
namespace QCD {

template<class Impl>
void WilsonKernels<Impl >::DiracOptAsmDhopSite(StencilImpl &st,DoubledGaugeField &U,
						   std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
					       int ss,int sU,const FermionField &in, FermionField &out,uint64_t *timers)
{
  uint64_t  now;
  uint64_t first ;
  int offset,local,perm, ptype;
  const SiteHalfSpinor *pbuf = & buf[0];
  const SiteSpinor   *plocal = & in._odata[0];
  void *pf;
  int osites = in._grid->oSites();

  
  StencilEntry *SE;

  //#define STAMP(i) timers[i] = __rdtsc() ; 
#define STAMP(i) //timers[i] = __rdtsc() ; 

  MASK_REGS;

  first = __rdtsc();

  SE=st.GetEntry(ptype,Xm,ss);

#if 0
  if (SE->_is_local) pf=(void *)&plocal[SE->_offset];
  else               pf=(void *)&pbuf[SE->_offset];

  LOAD64(%r9,pf);
  __asm__( 
	  VPREFETCH(0,%r9)
	  VPREFETCH(1,%r9)
	  VPREFETCH(2,%r9)
	  VPREFETCH(3,%r9)
	  VPREFETCH(4,%r9)
	  VPREFETCH(5,%r9)
	  VPREFETCH(6,%r9)
	  VPREFETCH(7,%r9)
	  VPREFETCH(8,%r9)
	  VPREFETCH(9,%r9)
	  VPREFETCH(10,%r9)
	  VPREFETCH(11,%r9) );
#endif

  // Xm
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;

  // Prefetch
  SE=st.GetEntry(ptype,Ym,ss);
  if (SE->_is_local) pf=(void *)&plocal[SE->_offset];
  else               pf=(void *)&pbuf[SE->_offset];
  
  if ( local ) {
    XM_PROJMEM(&plocal[offset]);
    if ( perm) {
      PERMUTE_DIR3; // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI(&pbuf[offset]);
  }
  {
    MULT_2SPIN_DIR_PFXM(Xm,pf);
  }
  XM_RECON;

  // Ym
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;

  // Prefetch
  SE=st.GetEntry(ptype,Zm,ss);
  if (SE->_is_local) pf=(void *)&plocal[SE->_offset];
  else               pf=(void *)&pbuf[SE->_offset];
  
  if ( local ) {
    YM_PROJMEM(&plocal[offset]);
    if ( perm) {
      PERMUTE_DIR2; // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI(&pbuf[offset]);
  }
  {
    MULT_2SPIN_DIR_PFYM(Ym,pf);
  }
  YM_RECON_ACCUM;

  // Zm
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;

  // Prefetch
  SE=st.GetEntry(ptype,Tm,ss);
  if (SE->_is_local) pf=(void *)&plocal[SE->_offset];
  else               pf=(void *)&pbuf[SE->_offset];

  if ( local ) {
    ZM_PROJMEM(&plocal[offset]);
    if ( perm) {
      PERMUTE_DIR1; // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI(&pbuf[offset]);
  }
  {
    MULT_2SPIN_DIR_PFZM(Zm,pf);
  }
  ZM_RECON_ACCUM;

  // Tm
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  SE=st.GetEntry(ptype,Tp,ss);
  if (SE->_is_local) pf=(void *)&plocal[SE->_offset];
  else               pf=(void *)&pbuf[SE->_offset];


  if ( local ) {
    TM_PROJMEM(&plocal[offset]);
    if ( perm) {
      PERMUTE_DIR0; // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI(&pbuf[offset]);
  }
  {
    MULT_2SPIN_DIR_PFTM(Tm,pf);
  }
  TM_RECON_ACCUM;

  // Tp
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;

  // Prefetch
  SE=st.GetEntry(ptype,Zp,ss);
  if (SE->_is_local) pf=(void *)&plocal[SE->_offset];
  else               pf=(void *)&pbuf[SE->_offset];
  
  if ( local ) {
    TP_PROJMEM(&plocal[offset]);
    if ( perm) {
      PERMUTE_DIR0; // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI(&pbuf[offset]);
  }
  {
    MULT_2SPIN_DIR_PFTP(Tp,pf);
  }
  TP_RECON_ACCUM;

  // Zp
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;

  // Prefetch
  SE=st.GetEntry(ptype,Yp,ss);
  if (SE->_is_local) pf=(void *)&plocal[SE->_offset];
  else               pf=(void *)&pbuf[SE->_offset];

  if ( local ) {
    ZP_PROJMEM(&plocal[offset]);
    if ( perm) {
      PERMUTE_DIR1; // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI(&pbuf[offset]);
  }
  {
    MULT_2SPIN_DIR_PFZP(Zp,pf);
  }
  ZP_RECON_ACCUM;


  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;

  // Prefetch
  SE=st.GetEntry(ptype,Xp,ss);
  if (SE->_is_local) pf=(void *)&plocal[SE->_offset];
  else               pf=(void *)&pbuf[SE->_offset];
  
  if ( local ) {
    YP_PROJMEM(&plocal[offset]);
    if ( perm) {
      PERMUTE_DIR2; // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI(&pbuf[offset]);
  }
  {
    MULT_2SPIN_DIR_PFYP(Yp,pf);
  }
  YP_RECON_ACCUM;

  // Xp
  perm   = SE->_permute;
  offset = SE->_offset;
  local  = SE->_is_local;
    
  //  PREFETCH_R(A);

  // Prefetch
  SE=st.GetEntry(ptype,Xm,(ss+1)%osites);
  if (SE->_is_local) pf=(void *)&plocal[SE->_offset];
  else               pf=(void *)&pbuf[SE->_offset];

  if ( local ) {
    XP_PROJMEM(&plocal[offset]);
    if ( perm) {
      PERMUTE_DIR3; // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI(&pbuf[offset]);
  }
  {
    MULT_2SPIN_DIR_PFXP(Xp,pf);
  }
  XP_RECON_ACCUM;

 debug:
  SAVE_RESULT(&out._odata[ss]);

}

  template class WilsonKernels<WilsonImplF>;		
  template class WilsonKernels<WilsonImplD>; 

}}
#endif
