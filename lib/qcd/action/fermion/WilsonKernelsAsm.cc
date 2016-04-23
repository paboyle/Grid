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
#if defined(AVX512) 
//#if defined (IMCI)

#include <simd/Intel512wilson.h>

#include <simd/Intel512single.h>


namespace Grid {
namespace QCD {

template<class Impl>
void WilsonKernels<Impl >::DiracOptAsmDhopSite(StencilImpl &st,DoubledGaugeField &U,
						   std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
					       int ss,int sU,const FermionField &in, FermionField &out)
{
  uint64_t  now;
  uint64_t first ;
  int offset,local,perm, ptype;
  const SiteHalfSpinor *pbuf = & buf[0];
  const SiteSpinor   *plocal = & in._odata[0];
  void *pf;
  int osites = in._grid->oSites();

  
  StencilEntry *SE;

  //#define STAMP(i) timers[i] = cyclecount() ; 
#define STAMP(i) //timers[i] = cyclecount() ; 

  MASK_REGS;

  first = cyclecount();

  SE=st.GetEntry(ptype,Xm,ss);

  // Xm
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;

  // Prefetch
  SE=st.GetEntry(ptype,Ym,ss);
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
    MULT_2SPIN_DIR_PFXM(Xm,pf);
  }
  XP_RECON;

  // Ym
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;

  // Prefetch
  SE=st.GetEntry(ptype,Zm,ss);
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
    MULT_2SPIN_DIR_PFYM(Ym,pf);
  }
  YP_RECON_ACCUM;

  // Zm
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;

  // Prefetch
  SE=st.GetEntry(ptype,Tm,ss);
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
    MULT_2SPIN_DIR_PFZM(Zm,pf);
  }
  ZP_RECON_ACCUM;

  // Tm
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  SE=st.GetEntry(ptype,Tp,ss);
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
    MULT_2SPIN_DIR_PFTM(Tm,pf);
  }
  TP_RECON_ACCUM;

  // Tp
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;

  // Prefetch
  SE=st.GetEntry(ptype,Zp,ss);
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
    MULT_2SPIN_DIR_PFTP(Tp,pf);
  }
  TM_RECON_ACCUM;

  // Zp
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;

  // Prefetch
  SE=st.GetEntry(ptype,Yp,ss);
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
    MULT_2SPIN_DIR_PFZP(Zp,pf);
  }
  ZM_RECON_ACCUM;


  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;

  // Prefetch
  SE=st.GetEntry(ptype,Xp,ss);
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
    MULT_2SPIN_DIR_PFYP(Yp,pf);
  }
  YM_RECON_ACCUM;

  // Xp
  perm   = SE->_permute;
  offset = SE->_offset;
  local  = SE->_is_local;
    
  // Prefetch
  SE=st.GetEntry(ptype,Xm,(ss+1)%osites);
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
    MULT_2SPIN_DIR_PFXP(Xp,pf);
  }
  XM_RECON_ACCUM;

 debug:
  SAVE_RESULT(&out._odata[ss]);

}

  template class WilsonKernels<WilsonImplF>;		
  template class WilsonKernels<WilsonImplD>; 
  template class WilsonKernels<GparityWilsonImplF>;
  template class WilsonKernels<GparityWilsonImplD>;
  template class WilsonKernels<DomainWallRedBlack5dImplF>;
  template class WilsonKernels<DomainWallRedBlack5dImplD>;
}}
#endif
