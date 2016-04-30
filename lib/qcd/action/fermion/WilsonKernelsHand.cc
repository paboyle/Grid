    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/WilsonKernelsHand.cc

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

#define REGISTER

#define LOAD_CHIMU \
  const SiteSpinor & ref (in._odata[offset]);	\
    Chimu_00=ref()(0)(0);\
    Chimu_01=ref()(0)(1);\
    Chimu_02=ref()(0)(2);\
    Chimu_10=ref()(1)(0);\
    Chimu_11=ref()(1)(1);\
    Chimu_12=ref()(1)(2);\
    Chimu_20=ref()(2)(0);\
    Chimu_21=ref()(2)(1);\
    Chimu_22=ref()(2)(2);\
    Chimu_30=ref()(3)(0);\
    Chimu_31=ref()(3)(1);\
    Chimu_32=ref()(3)(2);

#define LOAD_CHI\
  const SiteHalfSpinor &ref(buf[offset]);	\
    Chi_00 = ref()(0)(0);\
    Chi_01 = ref()(0)(1);\
    Chi_02 = ref()(0)(2);\
    Chi_10 = ref()(1)(0);\
    Chi_11 = ref()(1)(1);\
    Chi_12 = ref()(1)(2);

// To splat or not to splat depends on the implementation
#define MULT_2SPIN(A)\
   auto & ref(U._odata[sU](A));	\
   Impl::loadLinkElement(U_00,ref()(0,0));	\
   Impl::loadLinkElement(U_10,ref()(1,0));	\
   Impl::loadLinkElement(U_20,ref()(2,0));	\
   Impl::loadLinkElement(U_01,ref()(0,1));	\
   Impl::loadLinkElement(U_11,ref()(1,1));	\
   Impl::loadLinkElement(U_21,ref()(2,1));	\
    UChi_00 = U_00*Chi_00;\
    UChi_10 = U_00*Chi_10;\
    UChi_01 = U_10*Chi_00;\
    UChi_11 = U_10*Chi_10;\
    UChi_02 = U_20*Chi_00;\
    UChi_12 = U_20*Chi_10;\
    UChi_00+= U_01*Chi_01;\
    UChi_10+= U_01*Chi_11;\
    UChi_01+= U_11*Chi_01;\
    UChi_11+= U_11*Chi_11;\
    UChi_02+= U_21*Chi_01;\
    UChi_12+= U_21*Chi_11;\
    Impl::loadLinkElement(U_00,ref()(0,2));	\
    Impl::loadLinkElement(U_10,ref()(1,2));	\
    Impl::loadLinkElement(U_20,ref()(2,2));	\
    UChi_00+= U_00*Chi_02;\
    UChi_10+= U_00*Chi_12;\
    UChi_01+= U_10*Chi_02;\
    UChi_11+= U_10*Chi_12;\
    UChi_02+= U_20*Chi_02;\
    UChi_12+= U_20*Chi_12;


#define PERMUTE_DIR(dir)			\
      permute##dir(Chi_00,Chi_00);\
      permute##dir(Chi_01,Chi_01);\
      permute##dir(Chi_02,Chi_02);\
      permute##dir(Chi_10,Chi_10);\
      permute##dir(Chi_11,Chi_11);\
      permute##dir(Chi_12,Chi_12);

//      hspin(0)=fspin(0)+timesI(fspin(3));
//      hspin(1)=fspin(1)+timesI(fspin(2));
#define XP_PROJ \
    Chi_00 = Chimu_00+timesI(Chimu_30);\
    Chi_01 = Chimu_01+timesI(Chimu_31);\
    Chi_02 = Chimu_02+timesI(Chimu_32);\
    Chi_10 = Chimu_10+timesI(Chimu_20);\
    Chi_11 = Chimu_11+timesI(Chimu_21);\
    Chi_12 = Chimu_12+timesI(Chimu_22);

#define YP_PROJ \
    Chi_00 = Chimu_00-Chimu_30;\
    Chi_01 = Chimu_01-Chimu_31;\
    Chi_02 = Chimu_02-Chimu_32;\
    Chi_10 = Chimu_10+Chimu_20;\
    Chi_11 = Chimu_11+Chimu_21;\
    Chi_12 = Chimu_12+Chimu_22;

#define ZP_PROJ \
  Chi_00 = Chimu_00+timesI(Chimu_20);		\
  Chi_01 = Chimu_01+timesI(Chimu_21);		\
  Chi_02 = Chimu_02+timesI(Chimu_22);		\
  Chi_10 = Chimu_10-timesI(Chimu_30);		\
  Chi_11 = Chimu_11-timesI(Chimu_31);		\
  Chi_12 = Chimu_12-timesI(Chimu_32);

#define TP_PROJ \
  Chi_00 = Chimu_00+Chimu_20;		\
  Chi_01 = Chimu_01+Chimu_21;		\
  Chi_02 = Chimu_02+Chimu_22;		\
  Chi_10 = Chimu_10+Chimu_30;		\
  Chi_11 = Chimu_11+Chimu_31;		\
  Chi_12 = Chimu_12+Chimu_32;


//      hspin(0)=fspin(0)-timesI(fspin(3));
//      hspin(1)=fspin(1)-timesI(fspin(2));
#define XM_PROJ \
    Chi_00 = Chimu_00-timesI(Chimu_30);\
    Chi_01 = Chimu_01-timesI(Chimu_31);\
    Chi_02 = Chimu_02-timesI(Chimu_32);\
    Chi_10 = Chimu_10-timesI(Chimu_20);\
    Chi_11 = Chimu_11-timesI(Chimu_21);\
    Chi_12 = Chimu_12-timesI(Chimu_22);

#define YM_PROJ \
    Chi_00 = Chimu_00+Chimu_30;\
    Chi_01 = Chimu_01+Chimu_31;\
    Chi_02 = Chimu_02+Chimu_32;\
    Chi_10 = Chimu_10-Chimu_20;\
    Chi_11 = Chimu_11-Chimu_21;\
    Chi_12 = Chimu_12-Chimu_22;

#define ZM_PROJ \
  Chi_00 = Chimu_00-timesI(Chimu_20);		\
  Chi_01 = Chimu_01-timesI(Chimu_21);		\
  Chi_02 = Chimu_02-timesI(Chimu_22);		\
  Chi_10 = Chimu_10+timesI(Chimu_30);		\
  Chi_11 = Chimu_11+timesI(Chimu_31);		\
  Chi_12 = Chimu_12+timesI(Chimu_32);

#define TM_PROJ \
  Chi_00 = Chimu_00-Chimu_20;		\
  Chi_01 = Chimu_01-Chimu_21;		\
  Chi_02 = Chimu_02-Chimu_22;		\
  Chi_10 = Chimu_10-Chimu_30;		\
  Chi_11 = Chimu_11-Chimu_31;		\
  Chi_12 = Chimu_12-Chimu_32;

//      fspin(0)=hspin(0);
//      fspin(1)=hspin(1);
//      fspin(2)=timesMinusI(hspin(1));
//      fspin(3)=timesMinusI(hspin(0));
#define XP_RECON\
  result_00 = UChi_00;\
  result_01 = UChi_01;\
  result_02 = UChi_02;\
  result_10 = UChi_10;\
  result_11 = UChi_11;\
  result_12 = UChi_12;\
  result_20 = timesMinusI(UChi_10);\
  result_21 = timesMinusI(UChi_11);\
  result_22 = timesMinusI(UChi_12);\
  result_30 = timesMinusI(UChi_00);\
  result_31 = timesMinusI(UChi_01);\
  result_32 = timesMinusI(UChi_02);

#define XP_RECON_ACCUM\
  result_00+=UChi_00;\
  result_01+=UChi_01;\
  result_02+=UChi_02;\
  result_10+=UChi_10;\
  result_11+=UChi_11;\
  result_12+=UChi_12;\
  result_20-=timesI(UChi_10);\
  result_21-=timesI(UChi_11);\
  result_22-=timesI(UChi_12);\
  result_30-=timesI(UChi_00);\
  result_31-=timesI(UChi_01);\
  result_32-=timesI(UChi_02);

#define XM_RECON\
  result_00 = UChi_00;\
  result_01 = UChi_01;\
  result_02 = UChi_02;\
  result_10 = UChi_10;\
  result_11 = UChi_11;\
  result_12 = UChi_12;\
  result_20 = timesI(UChi_10);\
  result_21 = timesI(UChi_11);\
  result_22 = timesI(UChi_12);\
  result_30 = timesI(UChi_00);\
  result_31 = timesI(UChi_01);\
  result_32 = timesI(UChi_02);

#define XM_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20+= timesI(UChi_10);\
  result_21+= timesI(UChi_11);\
  result_22+= timesI(UChi_12);\
  result_30+= timesI(UChi_00);\
  result_31+= timesI(UChi_01);\
  result_32+= timesI(UChi_02);

#define YP_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20+= UChi_10;\
  result_21+= UChi_11;\
  result_22+= UChi_12;\
  result_30-= UChi_00;\
  result_31-= UChi_01;\
  result_32-= UChi_02;

#define YM_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20-= UChi_10;\
  result_21-= UChi_11;\
  result_22-= UChi_12;\
  result_30+= UChi_00;\
  result_31+= UChi_01;\
  result_32+= UChi_02;

#define ZP_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20-= timesI(UChi_00);			\
  result_21-= timesI(UChi_01);			\
  result_22-= timesI(UChi_02);			\
  result_30+= timesI(UChi_10);			\
  result_31+= timesI(UChi_11);			\
  result_32+= timesI(UChi_12);

#define ZM_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20+= timesI(UChi_00);			\
  result_21+= timesI(UChi_01);			\
  result_22+= timesI(UChi_02);			\
  result_30-= timesI(UChi_10);			\
  result_31-= timesI(UChi_11);			\
  result_32-= timesI(UChi_12);

#define TP_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20+= UChi_00;			\
  result_21+= UChi_01;			\
  result_22+= UChi_02;			\
  result_30+= UChi_10;			\
  result_31+= UChi_11;			\
  result_32+= UChi_12;

#define TM_RECON_ACCUM\
  result_00+= UChi_00;\
  result_01+= UChi_01;\
  result_02+= UChi_02;\
  result_10+= UChi_10;\
  result_11+= UChi_11;\
  result_12+= UChi_12;\
  result_20-= UChi_00;	\
  result_21-= UChi_01;	\
  result_22-= UChi_02;	\
  result_30-= UChi_10;	\
  result_31-= UChi_11;	\
  result_32-= UChi_12;

namespace Grid {
namespace QCD {


template<class Impl>
int WilsonKernels<Impl >::DiracOptHandDhopSite(StencilImpl &st,DoubledGaugeField &U,
					       std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
					       int ss,int sU,const FermionField &in, FermionField &out)
{
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  REGISTER Simd result_00; // 12 regs on knc
  REGISTER Simd result_01;
  REGISTER Simd result_02;

  REGISTER Simd result_10;
  REGISTER Simd result_11;
  REGISTER Simd result_12;

  REGISTER Simd result_20;
  REGISTER Simd result_21;
  REGISTER Simd result_22;

  REGISTER Simd result_30;
  REGISTER Simd result_31;
  REGISTER Simd result_32; // 20 left

  REGISTER Simd Chi_00;    // two spinor; 6 regs
  REGISTER Simd Chi_01;
  REGISTER Simd Chi_02;

  REGISTER Simd Chi_10;
  REGISTER Simd Chi_11;
  REGISTER Simd Chi_12;   // 14 left

  REGISTER Simd UChi_00;  // two spinor; 6 regs
  REGISTER Simd UChi_01;
  REGISTER Simd UChi_02;

  REGISTER Simd UChi_10;
  REGISTER Simd UChi_11;
  REGISTER Simd UChi_12;  // 8 left

  REGISTER Simd U_00;  // two rows of U matrix
  REGISTER Simd U_10;
  REGISTER Simd U_20;  
  REGISTER Simd U_01;
  REGISTER Simd U_11;
  REGISTER Simd U_21;  // 2 reg left.

#define Chimu_00 Chi_00
#define Chimu_01 Chi_01
#define Chimu_02 Chi_02
#define Chimu_10 Chi_10
#define Chimu_11 Chi_11
#define Chimu_12 Chi_12
#define Chimu_20 UChi_00
#define Chimu_21 UChi_01
#define Chimu_22 UChi_02
#define Chimu_30 UChi_10
#define Chimu_31 UChi_11
#define Chimu_32 UChi_12


  int offset,local,perm, ptype;
  StencilEntry *SE;

  // Xp
  SE=st.GetEntry(ptype,Xp,ss);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHIMU;
    XM_PROJ;
    if ( perm) {
      PERMUTE_DIR(3); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Xp);
  }
  XM_RECON;
  
  // Yp
  SE=st.GetEntry(ptype,Yp,ss);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHIMU;
    YM_PROJ;
    if ( perm) {
      PERMUTE_DIR(2); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Yp);
  }
  YM_RECON_ACCUM;


  // Zp
  SE=st.GetEntry(ptype,Zp,ss);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHIMU;
    ZM_PROJ;
    if ( perm) {
      PERMUTE_DIR(1); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Zp);
  }
  ZM_RECON_ACCUM;

  // Tp
  SE=st.GetEntry(ptype,Tp,ss);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHIMU;
    TM_PROJ;
    if ( perm) {
      PERMUTE_DIR(0); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Tp);
  }
  TM_RECON_ACCUM;
  
  // Xm
  SE=st.GetEntry(ptype,Xm,ss);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHIMU;
    XP_PROJ;
    if ( perm) {
      PERMUTE_DIR(3); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Xm);
  }
  XP_RECON_ACCUM;
  
  
  // Ym
  SE=st.GetEntry(ptype,Ym,ss);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHIMU;
    YP_PROJ;
    if ( perm) {
      PERMUTE_DIR(2); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Ym);
  }
  YP_RECON_ACCUM;

  // Zm
  SE=st.GetEntry(ptype,Zm,ss);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHIMU;
    ZP_PROJ;
    if ( perm) {
      PERMUTE_DIR(1); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Zm);
  }
  ZP_RECON_ACCUM;

  // Tm
  SE=st.GetEntry(ptype,Tm,ss);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHIMU;
    TP_PROJ;
    if ( perm) {
      PERMUTE_DIR(0); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Tm);
  }
  TP_RECON_ACCUM;

  {
    SiteSpinor & ref (out._odata[ss]);
    vstream(ref()(0)(0),result_00);
    vstream(ref()(0)(1),result_01);
    vstream(ref()(0)(2),result_02);
    vstream(ref()(1)(0),result_10);
    vstream(ref()(1)(1),result_11);
    vstream(ref()(1)(2),result_12);
    vstream(ref()(2)(0),result_20);
    vstream(ref()(2)(1),result_21);
    vstream(ref()(2)(2),result_22);
    vstream(ref()(3)(0),result_30);
    vstream(ref()(3)(1),result_31);
    vstream(ref()(3)(2),result_32);
  }
  return 0;
}


template<class Impl>
int WilsonKernels<Impl >::DiracOptHandDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
					       std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
					       int ss,int sU,const FermionField &in, FermionField &out)
{
  //  std::cout << "Hand op Dhop "<<std::endl;
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  REGISTER Simd result_00; // 12 regs on knc
  REGISTER Simd result_01;
  REGISTER Simd result_02;
  
  REGISTER Simd result_10;
  REGISTER Simd result_11;
  REGISTER Simd result_12;

  REGISTER Simd result_20;
  REGISTER Simd result_21;
  REGISTER Simd result_22;

  REGISTER Simd result_30;
  REGISTER Simd result_31;
  REGISTER Simd result_32; // 20 left

  REGISTER Simd Chi_00;    // two spinor; 6 regs
  REGISTER Simd Chi_01;
  REGISTER Simd Chi_02;

  REGISTER Simd Chi_10;
  REGISTER Simd Chi_11;
  REGISTER Simd Chi_12;   // 14 left

  REGISTER Simd UChi_00;  // two spinor; 6 regs
  REGISTER Simd UChi_01;
  REGISTER Simd UChi_02;

  REGISTER Simd UChi_10;
  REGISTER Simd UChi_11;
  REGISTER Simd UChi_12;  // 8 left

  REGISTER Simd U_00;  // two rows of U matrix
  REGISTER Simd U_10;
  REGISTER Simd U_20;  
  REGISTER Simd U_01;
  REGISTER Simd U_11;
  REGISTER Simd U_21;  // 2 reg left.

#define Chimu_00 Chi_00
#define Chimu_01 Chi_01
#define Chimu_02 Chi_02
#define Chimu_10 Chi_10
#define Chimu_11 Chi_11
#define Chimu_12 Chi_12
#define Chimu_20 UChi_00
#define Chimu_21 UChi_01
#define Chimu_22 UChi_02
#define Chimu_30 UChi_10
#define Chimu_31 UChi_11
#define Chimu_32 UChi_12


  StencilEntry *SE;
  int offset,local,perm, ptype;
  
  // Xp
  SE=st.GetEntry(ptype,Xp,ss);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHIMU;
    XP_PROJ;
    if ( perm) {
      PERMUTE_DIR(3); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI;
  }

  {
    MULT_2SPIN(Xp);
  }
  XP_RECON;

  // Yp
  SE=st.GetEntry(ptype,Yp,ss);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHIMU;
    YP_PROJ;
    if ( perm) {
      PERMUTE_DIR(2); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Yp);
  }
  YP_RECON_ACCUM;


  // Zp
  SE=st.GetEntry(ptype,Zp,ss);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHIMU;
    ZP_PROJ;
    if ( perm) {
      PERMUTE_DIR(1); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Zp);
  }
  ZP_RECON_ACCUM;

  // Tp
  SE=st.GetEntry(ptype,Tp,ss);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHIMU;
    TP_PROJ;
    if ( perm) {
      PERMUTE_DIR(0); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Tp);
  }
  TP_RECON_ACCUM;
  
  // Xm
  SE=st.GetEntry(ptype,Xm,ss);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHIMU;
    XM_PROJ;
    if ( perm) {
      PERMUTE_DIR(3); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Xm);
  }
  XM_RECON_ACCUM;
  
  // Ym
  SE=st.GetEntry(ptype,Ym,ss);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;
  
  if ( local ) {
    LOAD_CHIMU;
    YM_PROJ;
    if ( perm) {
      PERMUTE_DIR(2); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Ym);
  }
  YM_RECON_ACCUM;

  // Zm
  SE=st.GetEntry(ptype,Zm,ss);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;

  if ( local ) {
    LOAD_CHIMU;
    ZM_PROJ;
    if ( perm) {
      PERMUTE_DIR(1); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Zm);
  }
  ZM_RECON_ACCUM;

  // Tm
  SE=st.GetEntry(ptype,Tm,ss);
  offset = SE->_offset;
  local  = SE->_is_local;
  perm   = SE->_permute;

  if ( local ) {
    LOAD_CHIMU;
    TM_PROJ;
    if ( perm) {
      PERMUTE_DIR(0); // T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
    }
  } else { 
    LOAD_CHI;
  }
  {
    MULT_2SPIN(Tm);
  }
  TM_RECON_ACCUM;

  {
    SiteSpinor & ref (out._odata[ss]);
    vstream(ref()(0)(0),result_00);
    vstream(ref()(0)(1),result_01);
    vstream(ref()(0)(2),result_02);
    vstream(ref()(1)(0),result_10);
    vstream(ref()(1)(1),result_11);
    vstream(ref()(1)(2),result_12);
    vstream(ref()(2)(0),result_20);
    vstream(ref()(2)(1),result_21);
    vstream(ref()(2)(2),result_22);
    vstream(ref()(3)(0),result_30);
    vstream(ref()(3)(1),result_31);
    vstream(ref()(3)(2),result_32);
  }
  return 0;
}


  ////////////////////////////////////////////////
  // Specialise Gparity to simple implementation
  ////////////////////////////////////////////////
template<>
int WilsonKernels<GparityWilsonImplF>::DiracOptHandDhopSite(StencilImpl &st,DoubledGaugeField &U,
							     std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
							     int sF,int sU,const FermionField &in, FermionField &out)
{
  DiracOptDhopSite(st,U,buf,sF,sU,in,out); // returns void, will template override for Wilson Nc=3
  return 0;
  
}

template<>
int WilsonKernels<GparityWilsonImplF>::DiracOptHandDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
								std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
								int sF,int sU,const FermionField &in, FermionField &out)
{
  DiracOptDhopSiteDag(st,U,buf,sF,sU,in,out); // will template override for Wilson Nc=3
  return 0;
}

template<>
int WilsonKernels<GparityWilsonImplD>::DiracOptHandDhopSite(StencilImpl &st,DoubledGaugeField &U,
							     std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
							     int sF,int sU,const FermionField &in, FermionField &out)
{
  DiracOptDhopSite(st,U,buf,sF,sU,in,out); // will template override for Wilson Nc=3
  return 0;
}

template<>
int WilsonKernels<GparityWilsonImplD>::DiracOptHandDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
								std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
								int sF,int sU,const FermionField &in, FermionField &out)
{
  DiracOptDhopSiteDag(st,U,buf,sF,sU,in,out); // will template override for Wilson Nc=3
  return 0;
}



  //////////////
/*
template<>
int WilsonKernels<DomainWallRedBlack5dImplF>::DiracOptHandDhopSite(StencilImpl &st,DoubledGaugeField &U,
							     std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
							     int sF,int sU,const FermionField &in, FermionField &out)
{
  DiracOptDhopSite(st,U,buf,sF,sU,in,out); // returns void, will template override for Wilson Nc=3
  return 0;
  
}

template<>
int WilsonKernels<DomainWallRedBlack5dImplF>::DiracOptHandDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
								std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
								int sF,int sU,const FermionField &in, FermionField &out)
{
  DiracOptDhopSiteDag(st,U,buf,sF,sU,in,out); // will template override for Wilson Nc=3
  return 0;
}

template<>
int WilsonKernels<DomainWallRedBlack5dImplD>::DiracOptHandDhopSite(StencilImpl &st,DoubledGaugeField &U,
							     std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
							     int sF,int sU,const FermionField &in, FermionField &out)
{
  DiracOptDhopSite(st,U,buf,sF,sU,in,out); // will template override for Wilson Nc=3
  return 0;
}

template<>
int WilsonKernels<DomainWallRedBlack5dImplD>::DiracOptHandDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
								std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
								int sF,int sU,const FermionField &in, FermionField &out)
{
  DiracOptDhopSiteDag(st,U,buf,sF,sU,in,out); // will template override for Wilson Nc=3
  return 0;
}

*/

template int WilsonKernels<WilsonImplF>::DiracOptHandDhopSite(StencilImpl &st,DoubledGaugeField &U,
							       std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
							       int ss,int sU,const FermionField &in, FermionField &out);
template int WilsonKernels<WilsonImplD>::DiracOptHandDhopSite(StencilImpl &st,DoubledGaugeField &U,
							       std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
							       int ss,int sU,const FermionField &in, FermionField &out);
template int WilsonKernels<WilsonImplF>::DiracOptHandDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
								  std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
								  int ss,int sU,const FermionField &in, FermionField &out);
template int WilsonKernels<WilsonImplD>::DiracOptHandDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
								  std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
								  int ss,int sU,const FermionField &in, FermionField &out);


template int WilsonKernels<GparityWilsonImplF>::DiracOptHandDhopSite(StencilImpl &st,DoubledGaugeField &U,
								      std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
								      int ss,int sU,const FermionField &in, FermionField &out);
template int WilsonKernels<GparityWilsonImplD>::DiracOptHandDhopSite(StencilImpl &st,DoubledGaugeField &U,
								      std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
								      int ss,int sU,const FermionField &in, FermionField &out);
template int WilsonKernels<GparityWilsonImplF>::DiracOptHandDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
									 std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
									 int ss,int sU,const FermionField &in, FermionField &out);
template int WilsonKernels<GparityWilsonImplD>::DiracOptHandDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
									 std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
									 int ss,int sU,const FermionField &in, FermionField &out);




template int WilsonKernels<DomainWallRedBlack5dImplF>::DiracOptHandDhopSite(StencilImpl &st,DoubledGaugeField &U,
								      std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
								      int ss,int sU,const FermionField &in, FermionField &out);
template int WilsonKernels<DomainWallRedBlack5dImplD>::DiracOptHandDhopSite(StencilImpl &st,DoubledGaugeField &U,
								      std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
								      int ss,int sU,const FermionField &in, FermionField &out);
template int WilsonKernels<DomainWallRedBlack5dImplF>::DiracOptHandDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
									 std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
									 int ss,int sU,const FermionField &in, FermionField &out);
template int WilsonKernels<DomainWallRedBlack5dImplD>::DiracOptHandDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
									 std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
									 int ss,int sU,const FermionField &in, FermionField &out);


}}
