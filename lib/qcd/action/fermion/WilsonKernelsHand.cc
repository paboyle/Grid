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
#include <Grid/qcd/action/fermion/FermionCore.h>

#define REGISTER

#define LOAD_CHIMU_BODY(F)			\
  Chimu_00=ref(F)(0)(0);			\
  Chimu_01=ref(F)(0)(1);			\
  Chimu_02=ref(F)(0)(2);			\
  Chimu_10=ref(F)(1)(0);			\
  Chimu_11=ref(F)(1)(1);			\
  Chimu_12=ref(F)(1)(2);			\
  Chimu_20=ref(F)(2)(0);			\
  Chimu_21=ref(F)(2)(1);			\
  Chimu_22=ref(F)(2)(2);			\
  Chimu_30=ref(F)(3)(0);			\
  Chimu_31=ref(F)(3)(1);			\
  Chimu_32=ref(F)(3)(2)

#define LOAD_CHIMU(DIR,F,PERM)						\
  { const SiteSpinor & ref (in._odata[offset]); LOAD_CHIMU_BODY(F); }

#define LOAD_CHI_BODY(F)				\
    Chi_00 = ref(F)(0)(0);\
    Chi_01 = ref(F)(0)(1);\
    Chi_02 = ref(F)(0)(2);\
    Chi_10 = ref(F)(1)(0);\
    Chi_11 = ref(F)(1)(1);\
    Chi_12 = ref(F)(1)(2)

#define LOAD_CHI(DIR,F,PERM)					\
  {const SiteHalfSpinor &ref(buf[offset]); LOAD_CHI_BODY(F); }


//G-parity implementations using in-place intrinsic ops

//1l 1h -> 1h 1l
//0l 0h , 1h 1l -> 0l 1h 0h,1l
//0h,1l -> 1l,0h
//if( (distance == 1 && !perm_will_occur) || (distance == -1 && perm_will_occur) )
//Pulled fermion through forwards face, GPBC on upper component
//Need 0= 0l 1h   1= 1l 0h
//else if( (distance == -1 && !perm) || (distance == 1 && perm) )
//Pulled fermion through backwards face, GPBC on lower component
//Need 0= 1l 0h   1= 0l 1h

//1l 1h -> 1h 1l
//0l 0h , 1h 1l -> 0l 1h 0h,1l
#define DO_TWIST_0L_1H(INTO,S,C,F, PERM, tmp1, tmp2, tmp3)			\
  permute##PERM(tmp1, ref(1)(S)(C));				\
  exchange##PERM(tmp2,tmp3, ref(0)(S)(C), tmp1);		\
  INTO = tmp2;

//0l 0h -> 0h 0l
//1l 1h, 0h 0l -> 1l 0h, 1h 0l
#define DO_TWIST_1L_0H(INTO,S,C,F, PERM, tmp1, tmp2, tmp3)			\
  permute##PERM(tmp1, ref(0)(S)(C));				\
  exchange##PERM(tmp2,tmp3, ref(1)(S)(C), tmp1);		\
  INTO = tmp2;




#define LOAD_CHI_SETUP(DIR,F)						\
  g = F;								\
  direction = st._directions[DIR];				\
  distance = st._distances[DIR];				\
  sl = st._grid->_simd_layout[direction];			\
  inplace_twist = 0;						\
  if(SE->_around_the_world && this->Params.twists[DIR % 4]){		\
    if(sl == 1){							\
      g = (F+1) % 2;							\
    }else{								\
      inplace_twist = 1;						\
    }									\
  }  

#define LOAD_CHIMU_GPARITY_INPLACE_TWIST(DIR,F,PERM)			\
  { const SiteSpinor &ref(in._odata[offset]);				\
    LOAD_CHI_SETUP(DIR,F);						\
    if(!inplace_twist){							\
      LOAD_CHIMU_BODY(g);						\
    }else{								\
      if(  ( F==0 && ((distance == 1 && !perm) || (distance == -1 && perm)) ) || \
	   ( F==1 && ((distance == -1 && !perm) || (distance == 1 && perm)) ) ){ \
	DO_TWIST_0L_1H(Chimu_00,0,0,F,PERM,  U_00,U_01,U_10);		\
	DO_TWIST_0L_1H(Chimu_01,0,1,F,PERM,  U_11,U_20,U_21);		\
	DO_TWIST_0L_1H(Chimu_02,0,2,F,PERM,  U_00,U_01,U_10);		\
	DO_TWIST_0L_1H(Chimu_10,1,0,F,PERM,  U_11,U_20,U_21);		\
	DO_TWIST_0L_1H(Chimu_11,1,1,F,PERM,  U_00,U_01,U_10);		\
	DO_TWIST_0L_1H(Chimu_12,1,2,F,PERM,  U_11,U_20,U_21);		\
	DO_TWIST_0L_1H(Chimu_20,2,0,F,PERM,  U_00,U_01,U_10);		\
	DO_TWIST_0L_1H(Chimu_21,2,1,F,PERM,  U_11,U_20,U_21);		\
	DO_TWIST_0L_1H(Chimu_22,2,2,F,PERM,  U_00,U_01,U_10);		\
	DO_TWIST_0L_1H(Chimu_30,3,0,F,PERM,  U_11,U_20,U_21);		\
	DO_TWIST_0L_1H(Chimu_31,3,1,F,PERM,  U_00,U_01,U_10);		\
	DO_TWIST_0L_1H(Chimu_32,3,2,F,PERM,  U_11,U_20,U_21);		\
      }else{								\
	DO_TWIST_1L_0H(Chimu_00,0,0,F,PERM,  U_00,U_01,U_10);		\
	DO_TWIST_1L_0H(Chimu_01,0,1,F,PERM,  U_11,U_20,U_21);		\
	DO_TWIST_1L_0H(Chimu_02,0,2,F,PERM,  U_00,U_01,U_10);		\
	DO_TWIST_1L_0H(Chimu_10,1,0,F,PERM,  U_11,U_20,U_21);		\
	DO_TWIST_1L_0H(Chimu_11,1,1,F,PERM,  U_00,U_01,U_10);		\
	DO_TWIST_1L_0H(Chimu_12,1,2,F,PERM,  U_11,U_20,U_21);		\
	DO_TWIST_1L_0H(Chimu_20,2,0,F,PERM,  U_00,U_01,U_10);		\
	DO_TWIST_1L_0H(Chimu_21,2,1,F,PERM,  U_11,U_20,U_21);		\
	DO_TWIST_1L_0H(Chimu_22,2,2,F,PERM,  U_00,U_01,U_10);		\
	DO_TWIST_1L_0H(Chimu_30,3,0,F,PERM,  U_11,U_20,U_21);		\
	DO_TWIST_1L_0H(Chimu_31,3,1,F,PERM,  U_00,U_01,U_10);		\
	DO_TWIST_1L_0H(Chimu_32,3,2,F,PERM,  U_11,U_20,U_21);		\
      } \
    } \
  }


#define LOAD_CHI_GPARITY_INPLACE_TWIST(DIR,F,PERM)				\
  { const SiteHalfSpinor &ref(buf[offset]);				\
    LOAD_CHI_SETUP(DIR,F);						\
    if(!inplace_twist){							\
      LOAD_CHI_BODY(g);							\
    }else{								\
      if(  ( F==0 && ((distance == 1 && !perm) || (distance == -1 && perm)) ) || \
	   ( F==1 && ((distance == -1 && !perm) || (distance == 1 && perm)) ) ){ \
	DO_TWIST_0L_1H(Chi_00,0,0,F,PERM,  U_00,U_01,U_10);			\
	DO_TWIST_0L_1H(Chi_01,0,1,F,PERM,  U_11,U_20,U_21);			\
	DO_TWIST_0L_1H(Chi_02,0,2,F,PERM,  UChi_00,UChi_01,UChi_02);		\
	DO_TWIST_0L_1H(Chi_10,1,0,F,PERM,  UChi_10,UChi_11,UChi_12);		\
	DO_TWIST_0L_1H(Chi_11,1,1,F,PERM,  U_00,U_01,U_10);			\
	DO_TWIST_0L_1H(Chi_12,1,2,F,PERM,  U_11,U_20,U_21);			\
      }else{								\
	DO_TWIST_1L_0H(Chi_00,0,0,F,PERM,  U_00,U_01,U_10);			\
	DO_TWIST_1L_0H(Chi_01,0,1,F,PERM,  U_11,U_20,U_21);			\
	DO_TWIST_1L_0H(Chi_02,0,2,F,PERM,  UChi_00,UChi_01,UChi_02);		\
	DO_TWIST_1L_0H(Chi_10,1,0,F,PERM,  UChi_10,UChi_11,UChi_12);		\
	DO_TWIST_1L_0H(Chi_11,1,1,F,PERM,  U_00,U_01,U_10);			\
	DO_TWIST_1L_0H(Chi_12,1,2,F,PERM,  U_11,U_20,U_21);			\
      }									\
    }									\
  }


#define LOAD_CHI_GPARITY(DIR,F,PERM) LOAD_CHI_GPARITY_INPLACE_TWIST(DIR,F,PERM)
#define LOAD_CHIMU_GPARITY(DIR,F,PERM) LOAD_CHIMU_GPARITY_INPLACE_TWIST(DIR,F,PERM)

// To splat or not to splat depends on the implementation
#define MULT_2SPIN_BODY \
  Impl::loadLinkElement(U_00,ref()(0,0));	\
  Impl::loadLinkElement(U_10,ref()(1,0));	\
  Impl::loadLinkElement(U_20,ref()(2,0));	\
  Impl::loadLinkElement(U_01,ref()(0,1));	\
  Impl::loadLinkElement(U_11,ref()(1,1));	\
  Impl::loadLinkElement(U_21,ref()(2,1));	\
  UChi_00 = U_00*Chi_00;			\
  UChi_10 = U_00*Chi_10;			\
  UChi_01 = U_10*Chi_00;			\
  UChi_11 = U_10*Chi_10;			\
  UChi_02 = U_20*Chi_00;			\
  UChi_12 = U_20*Chi_10;			\
  UChi_00+= U_01*Chi_01;			\
  UChi_10+= U_01*Chi_11;			\
  UChi_01+= U_11*Chi_01;			\
  UChi_11+= U_11*Chi_11;			\
  UChi_02+= U_21*Chi_01;			\
  UChi_12+= U_21*Chi_11;			\
  Impl::loadLinkElement(U_00,ref()(0,2));	\
  Impl::loadLinkElement(U_10,ref()(1,2));	\
  Impl::loadLinkElement(U_20,ref()(2,2));	\
  UChi_00+= U_00*Chi_02;			\
  UChi_10+= U_00*Chi_12;			\
  UChi_01+= U_10*Chi_02;			\
  UChi_11+= U_10*Chi_12;			\
  UChi_02+= U_20*Chi_02;			\
  UChi_12+= U_20*Chi_12


#define MULT_2SPIN(A,F)					\
  {auto & ref(U._odata[sU](A)); MULT_2SPIN_BODY; }

#define MULT_2SPIN_GPARITY(A,F)				\
  {auto & ref(U._odata[sU](F)(A)); MULT_2SPIN_BODY; }


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

#define HAND_STENCIL_LEG(PROJ,PERM,DIR,RECON,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL) \
  SE=st.GetEntry(ptype,DIR,ss);			\
  offset = SE->_offset;				\
  local  = SE->_is_local;			\
  perm   = SE->_permute;			\
  if ( local ) {				\
    LOAD_CHIMU_IMPL(DIR,F,PERM);			\
    PROJ;					\
    if ( perm) {				\
      PERMUTE_DIR(PERM);			\
    }						\
  } else {					\
    LOAD_CHI_IMPL(DIR,F,PERM);			\
  }						\
  MULT_2SPIN_IMPL(DIR,F);			\
  RECON;					


#define HAND_STENCIL_LEG_INT(PROJ,PERM,DIR,RECON,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL)	\
  SE=st.GetEntry(ptype,DIR,ss);			\
  offset = SE->_offset;				\
  local  = SE->_is_local;			\
  perm   = SE->_permute;			\
  if ( local ) {				\
    LOAD_CHIMU_IMPL(DIR,F,PERM);			\
    PROJ;					\
    if ( perm) {				\
      PERMUTE_DIR(PERM);			\
    }						\
  } else if ( st.same_node[DIR] ) {		\
    LOAD_CHI_IMPL(DIR,F,PERM);			\
  }						\
  if (local || st.same_node[DIR] ) {		\
    MULT_2SPIN_IMPL(DIR,F);			\
    RECON;					\
  }

#define HAND_STENCIL_LEG_EXT(PROJ,PERM,DIR,RECON,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL)	\
  SE=st.GetEntry(ptype,DIR,ss);			\
  offset = SE->_offset;				\
  if((!SE->_is_local)&&(!st.same_node[DIR]) ) {	\
    LOAD_CHI_IMPL(DIR,F,PERM);			\
    MULT_2SPIN_IMPL(DIR,F);			\
    RECON;					\
    nmu++;					\
  }

#define HAND_RESULT(ss,F)			\
  {						\
    SiteSpinor & ref (out._odata[ss]);		\
    vstream(ref(F)(0)(0),result_00);		\
    vstream(ref(F)(0)(1),result_01);		\
    vstream(ref(F)(0)(2),result_02);		\
    vstream(ref(F)(1)(0),result_10);		\
    vstream(ref(F)(1)(1),result_11);		\
    vstream(ref(F)(1)(2),result_12);		\
    vstream(ref(F)(2)(0),result_20);		\
    vstream(ref(F)(2)(1),result_21);		\
    vstream(ref(F)(2)(2),result_22);		\
    vstream(ref(F)(3)(0),result_30);		\
    vstream(ref(F)(3)(1),result_31);		\
    vstream(ref(F)(3)(2),result_32);		\
  }

#define HAND_RESULT_EXT(ss,F)			\
  if (nmu){					\
    SiteSpinor & ref (out._odata[ss]);		\
    ref(F)(0)(0)+=result_00;		\
    ref(F)(0)(1)+=result_01;		\
    ref(F)(0)(2)+=result_02;		\
    ref(F)(1)(0)+=result_10;		\
    ref(F)(1)(1)+=result_11;		\
    ref(F)(1)(2)+=result_12;		\
    ref(F)(2)(0)+=result_20;		\
    ref(F)(2)(1)+=result_21;		\
    ref(F)(2)(2)+=result_22;		\
    ref(F)(3)(0)+=result_30;		\
    ref(F)(3)(1)+=result_31;		\
    ref(F)(3)(2)+=result_32;		\
  }


#define HAND_DECLARATIONS(a)			\
  Simd result_00;				\
  Simd result_01;				\
  Simd result_02;				\
  Simd result_10;				\
  Simd result_11;				\
  Simd result_12;				\
  Simd result_20;				\
  Simd result_21;				\
  Simd result_22;				\
  Simd result_30;				\
  Simd result_31;				\
  Simd result_32;				\
  Simd Chi_00;					\
  Simd Chi_01;					\
  Simd Chi_02;					\
  Simd Chi_10;					\
  Simd Chi_11;					\
  Simd Chi_12;					\
  Simd UChi_00;					\
  Simd UChi_01;					\
  Simd UChi_02;					\
  Simd UChi_10;					\
  Simd UChi_11;					\
  Simd UChi_12;					\
  Simd U_00;					\
  Simd U_10;					\
  Simd U_20;					\
  Simd U_01;					\
  Simd U_11;					\
  Simd U_21;

#define ZERO_RESULT				\
  result_00=zero;				\
  result_01=zero;				\
  result_02=zero;				\
  result_10=zero;				\
  result_11=zero;				\
  result_12=zero;				\
  result_20=zero;				\
  result_21=zero;				\
  result_22=zero;				\
  result_30=zero;				\
  result_31=zero;				\
  result_32=zero;			

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

namespace Grid {
namespace QCD {

template<class Impl> void 
WilsonKernels<Impl>::HandDhopSite(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor  *buf,
					  int ss,int sU,const FermionField &in, FermionField &out)
{
// T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  HAND_DECLARATIONS(ignore);

  int offset,local,perm, ptype;
  StencilEntry *SE;

#define HAND_DOP_SITE(F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL) \
  HAND_STENCIL_LEG(XM_PROJ,3,Xp,XM_RECON,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG(YM_PROJ,2,Yp,YM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL);	\
  HAND_STENCIL_LEG(ZM_PROJ,1,Zp,ZM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG(TM_PROJ,0,Tp,TM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG(XP_PROJ,3,Xm,XP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG(YP_PROJ,2,Ym,YP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG(ZP_PROJ,1,Zm,ZP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG(TP_PROJ,0,Tm,TP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_RESULT(ss,F)

  HAND_DOP_SITE(, LOAD_CHI,LOAD_CHIMU,MULT_2SPIN);
}

template<class Impl>
void WilsonKernels<Impl>::HandDhopSiteDag(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
						  int ss,int sU,const FermionField &in, FermionField &out)
{
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  HAND_DECLARATIONS(ignore);

  StencilEntry *SE;
  int offset,local,perm, ptype;

#define HAND_DOP_SITE_DAG(F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL) \
  HAND_STENCIL_LEG(XP_PROJ,3,Xp,XP_RECON,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG(YP_PROJ,2,Yp,YP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG(ZP_PROJ,1,Zp,ZP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG(TP_PROJ,0,Tp,TP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG(XM_PROJ,3,Xm,XM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG(YM_PROJ,2,Ym,YM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG(ZM_PROJ,1,Zm,ZM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG(TM_PROJ,0,Tm,TM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_RESULT(ss,F)

  HAND_DOP_SITE_DAG(, LOAD_CHI,LOAD_CHIMU,MULT_2SPIN);
}

template<class Impl> void 
WilsonKernels<Impl>::HandDhopSiteInt(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor  *buf,
					  int ss,int sU,const FermionField &in, FermionField &out)
{
// T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  HAND_DECLARATIONS(ignore);

  int offset,local,perm, ptype;
  StencilEntry *SE;

#define HAND_DOP_SITE_INT(F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL) \
  ZERO_RESULT; \
  HAND_STENCIL_LEG_INT(XM_PROJ,3,Xp,XM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_INT(YM_PROJ,2,Yp,YM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_INT(ZM_PROJ,1,Zp,ZM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_INT(TM_PROJ,0,Tp,TM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_INT(XP_PROJ,3,Xm,XP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_INT(YP_PROJ,2,Ym,YP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_INT(ZP_PROJ,1,Zm,ZP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_INT(TP_PROJ,0,Tm,TP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_RESULT(ss,F)

  HAND_DOP_SITE_INT(, LOAD_CHI,LOAD_CHIMU,MULT_2SPIN);
}

template<class Impl>
void WilsonKernels<Impl>::HandDhopSiteDagInt(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
						  int ss,int sU,const FermionField &in, FermionField &out)
{
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  HAND_DECLARATIONS(ignore);

  StencilEntry *SE;
  int offset,local,perm, ptype;

#define HAND_DOP_SITE_DAG_INT(F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL)				\
  ZERO_RESULT;							\
  HAND_STENCIL_LEG_INT(XP_PROJ,3,Xp,XP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL);		\
  HAND_STENCIL_LEG_INT(YP_PROJ,2,Yp,YP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL);		\
  HAND_STENCIL_LEG_INT(ZP_PROJ,1,Zp,ZP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL);		\
  HAND_STENCIL_LEG_INT(TP_PROJ,0,Tp,TP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL);		\
  HAND_STENCIL_LEG_INT(XM_PROJ,3,Xm,XM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL);		\
  HAND_STENCIL_LEG_INT(YM_PROJ,2,Ym,YM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL);		\
  HAND_STENCIL_LEG_INT(ZM_PROJ,1,Zm,ZM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL);		\
  HAND_STENCIL_LEG_INT(TM_PROJ,0,Tm,TM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL);		\
  HAND_RESULT(ss,F)
  
  HAND_DOP_SITE_DAG_INT(, LOAD_CHI,LOAD_CHIMU,MULT_2SPIN);
}

template<class Impl> void 
WilsonKernels<Impl>::HandDhopSiteExt(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor  *buf,
					  int ss,int sU,const FermionField &in, FermionField &out)
{
// T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  HAND_DECLARATIONS(ignore);

  int offset,local,perm, ptype;
  StencilEntry *SE;
  int nmu=0;

#define HAND_DOP_SITE_EXT(F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL) \
  ZERO_RESULT; \
  HAND_STENCIL_LEG_EXT(XM_PROJ,3,Xp,XM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_EXT(YM_PROJ,2,Yp,YM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_EXT(ZM_PROJ,1,Zp,ZM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_EXT(TM_PROJ,0,Tp,TM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_EXT(XP_PROJ,3,Xm,XP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_EXT(YP_PROJ,2,Ym,YP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_EXT(ZP_PROJ,1,Zm,ZP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_EXT(TP_PROJ,0,Tm,TP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_RESULT_EXT(ss,F)

  HAND_DOP_SITE_EXT(, LOAD_CHI,LOAD_CHIMU,MULT_2SPIN);
}

template<class Impl>
void WilsonKernels<Impl>::HandDhopSiteDagExt(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor *buf,
						  int ss,int sU,const FermionField &in, FermionField &out)
{
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  HAND_DECLARATIONS(ignore);

  StencilEntry *SE;
  int offset,local,perm, ptype;
  int nmu=0;

#define HAND_DOP_SITE_DAG_EXT(F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL) \
  ZERO_RESULT; \
  HAND_STENCIL_LEG_EXT(XP_PROJ,3,Xp,XP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_EXT(YP_PROJ,2,Yp,YP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_EXT(ZP_PROJ,1,Zp,ZP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_EXT(TP_PROJ,0,Tp,TP_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_EXT(XM_PROJ,3,Xm,XM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_EXT(YM_PROJ,2,Ym,YM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_EXT(ZM_PROJ,1,Zm,ZM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_STENCIL_LEG_EXT(TM_PROJ,0,Tm,TM_RECON_ACCUM,F,LOAD_CHI_IMPL,LOAD_CHIMU_IMPL,MULT_2SPIN_IMPL); \
  HAND_RESULT_EXT(ss,F)

  HAND_DOP_SITE_DAG_EXT(, LOAD_CHI,LOAD_CHIMU,MULT_2SPIN);
}

  ////////////////////////////////////////////////
  // Specialise Gparity to simple implementation
  ////////////////////////////////////////////////
#define HAND_SPECIALISE_EMPTY(IMPL)					\
  template<> void							\
  WilsonKernels<IMPL>::HandDhopSite(StencilImpl &st,			\
				    LebesgueOrder &lo,			\
				    DoubledGaugeField &U,		\
				    SiteHalfSpinor *buf,		\
				    int sF,int sU,			\
				    const FermionField &in,		\
				    FermionField &out){ assert(0); }	\
  template<> void							\
  WilsonKernels<IMPL>::HandDhopSiteDag(StencilImpl &st,			\
				    LebesgueOrder &lo,			\
				    DoubledGaugeField &U,		\
				    SiteHalfSpinor *buf,		\
				    int sF,int sU,			\
				    const FermionField &in,		\
				    FermionField &out){ assert(0); }	\
  template<> void							\
  WilsonKernels<IMPL>::HandDhopSiteInt(StencilImpl &st,			\
				    LebesgueOrder &lo,			\
				    DoubledGaugeField &U,		\
				    SiteHalfSpinor *buf,		\
				    int sF,int sU,			\
				    const FermionField &in,		\
				    FermionField &out){ assert(0); }	\
  template<> void							\
  WilsonKernels<IMPL>::HandDhopSiteExt(StencilImpl &st,			\
				    LebesgueOrder &lo,			\
				    DoubledGaugeField &U,		\
				    SiteHalfSpinor *buf,		\
				    int sF,int sU,			\
				    const FermionField &in,		\
				    FermionField &out){ assert(0); }	\
  template<> void							\
  WilsonKernels<IMPL>::HandDhopSiteDagInt(StencilImpl &st,	       	\
				    LebesgueOrder &lo,			\
				    DoubledGaugeField &U,		\
				    SiteHalfSpinor *buf,		\
				    int sF,int sU,			\
				    const FermionField &in,		\
				    FermionField &out){ assert(0); }	\
  template<> void							\
  WilsonKernels<IMPL>::HandDhopSiteDagExt(StencilImpl &st,	       	\
				    LebesgueOrder &lo,			\
				    DoubledGaugeField &U,		\
				    SiteHalfSpinor *buf,		\
				    int sF,int sU,			\
				    const FermionField &in,		\
				    FermionField &out){ assert(0); }	\



#define HAND_SPECIALISE_GPARITY(IMPL)					\
  template<> void							\
  WilsonKernels<IMPL>::HandDhopSite(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor  *buf, \
				    int ss,int sU,const FermionField &in, FermionField &out) \
  {									\
    typedef IMPL Impl;							\
    typedef typename Simd::scalar_type S;				\
    typedef typename Simd::vector_type V;				\
									\
    HAND_DECLARATIONS(ignore);						\
									\
    int offset,local,perm, ptype, g, direction, distance, sl, inplace_twist; \
    StencilEntry *SE;							\
    HAND_DOP_SITE(0, LOAD_CHI_GPARITY,LOAD_CHIMU_GPARITY,MULT_2SPIN_GPARITY); \
    HAND_DOP_SITE(1, LOAD_CHI_GPARITY,LOAD_CHIMU_GPARITY,MULT_2SPIN_GPARITY); \
  }									\
									\
  template<>								\
  void WilsonKernels<IMPL>::HandDhopSiteDag(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor *buf, \
					    int ss,int sU,const FermionField &in, FermionField &out) \
  {									\
    typedef IMPL Impl;							\
    typedef typename Simd::scalar_type S;				\
    typedef typename Simd::vector_type V;				\
									\
    HAND_DECLARATIONS(ignore);						\
									\
    StencilEntry *SE;							\
    int offset,local,perm, ptype, g, direction, distance, sl, inplace_twist;					\
    HAND_DOP_SITE_DAG(0, LOAD_CHI_GPARITY,LOAD_CHIMU_GPARITY,MULT_2SPIN_GPARITY); \
    HAND_DOP_SITE_DAG(1, LOAD_CHI_GPARITY,LOAD_CHIMU_GPARITY,MULT_2SPIN_GPARITY); \
  }									\
									\
  template<> void							\
  WilsonKernels<IMPL>::HandDhopSiteInt(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor  *buf, \
						     int ss,int sU,const FermionField &in, FermionField &out) \
  {									\
    typedef IMPL Impl;							\
    typedef typename Simd::scalar_type S;				\
    typedef typename Simd::vector_type V;				\
									\
    HAND_DECLARATIONS(ignore);						\
									\
    int offset,local,perm, ptype, g, direction, distance, sl, inplace_twist;					\
    StencilEntry *SE;							\
    HAND_DOP_SITE_INT(0, LOAD_CHI_GPARITY,LOAD_CHIMU_GPARITY,MULT_2SPIN_GPARITY); \
    HAND_DOP_SITE_INT(1, LOAD_CHI_GPARITY,LOAD_CHIMU_GPARITY,MULT_2SPIN_GPARITY); \
  }									\
									\
  template<>								\
  void WilsonKernels<IMPL>::HandDhopSiteDagInt(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor *buf, \
							     int ss,int sU,const FermionField &in, FermionField &out) \
  {									\
    typedef IMPL Impl;							\
    typedef typename Simd::scalar_type S;				\
    typedef typename Simd::vector_type V;				\
									\
    HAND_DECLARATIONS(ignore);						\
									\
    StencilEntry *SE;							\
    int offset,local,perm, ptype, g, direction, distance, sl, inplace_twist; \
    HAND_DOP_SITE_DAG_INT(0, LOAD_CHI_GPARITY,LOAD_CHIMU_GPARITY,MULT_2SPIN_GPARITY); \
    HAND_DOP_SITE_DAG_INT(1, LOAD_CHI_GPARITY,LOAD_CHIMU_GPARITY,MULT_2SPIN_GPARITY); \
  }									\
									\
  template<> void							\
  WilsonKernels<IMPL>::HandDhopSiteExt(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor  *buf, \
						     int ss,int sU,const FermionField &in, FermionField &out) \
  {									\
    typedef IMPL Impl;							\
    typedef typename Simd::scalar_type S;				\
    typedef typename Simd::vector_type V;				\
									\
    HAND_DECLARATIONS(ignore);						\
									\
    int offset,local,perm, ptype, g, direction, distance, sl, inplace_twist; \
    StencilEntry *SE;							\
    int nmu=0;								\
    HAND_DOP_SITE_EXT(0, LOAD_CHI_GPARITY,LOAD_CHIMU_GPARITY,MULT_2SPIN_GPARITY); \
    nmu = 0;								\
    HAND_DOP_SITE_EXT(1, LOAD_CHI_GPARITY,LOAD_CHIMU_GPARITY,MULT_2SPIN_GPARITY); \
  }									\
  template<>								\
  void WilsonKernels<IMPL>::HandDhopSiteDagExt(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor *buf, \
							     int ss,int sU,const FermionField &in, FermionField &out) \
  {									\
    typedef IMPL Impl;							\
    typedef typename Simd::scalar_type S;				\
    typedef typename Simd::vector_type V;				\
									\
    HAND_DECLARATIONS(ignore);						\
									\
    StencilEntry *SE;							\
    int offset,local,perm, ptype, g, direction, distance, sl, inplace_twist; \
    int nmu=0;								\
    HAND_DOP_SITE_DAG_EXT(0, LOAD_CHI_GPARITY,LOAD_CHIMU_GPARITY,MULT_2SPIN_GPARITY); \
    nmu = 0;								\
    HAND_DOP_SITE_DAG_EXT(1, LOAD_CHI_GPARITY,LOAD_CHIMU_GPARITY,MULT_2SPIN_GPARITY); \
  }


HAND_SPECIALISE_GPARITY(GparityWilsonImplF);
HAND_SPECIALISE_GPARITY(GparityWilsonImplD);
HAND_SPECIALISE_GPARITY(GparityWilsonImplFH);
HAND_SPECIALISE_GPARITY(GparityWilsonImplDF);










  
////////////// Wilson ; uses this implementation /////////////////////

#define INSTANTIATE_THEM(A) \
template void WilsonKernels<A>::HandDhopSite(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor *buf,\
					     int ss,int sU,const FermionField &in, FermionField &out); \
template void WilsonKernels<A>::HandDhopSiteDag(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor *buf, \
						int ss,int sU,const FermionField &in, FermionField &out);\
template void WilsonKernels<A>::HandDhopSiteInt(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor *buf,\
						int ss,int sU,const FermionField &in, FermionField &out); \
template void WilsonKernels<A>::HandDhopSiteDagInt(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor *buf, \
						   int ss,int sU,const FermionField &in, FermionField &out); \
template void WilsonKernels<A>::HandDhopSiteExt(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor *buf,\
						int ss,int sU,const FermionField &in, FermionField &out); \
template void WilsonKernels<A>::HandDhopSiteDagExt(StencilImpl &st,LebesgueOrder &lo,DoubledGaugeField &U,SiteHalfSpinor *buf, \
						   int ss,int sU,const FermionField &in, FermionField &out); 

INSTANTIATE_THEM(WilsonImplF);
INSTANTIATE_THEM(WilsonImplD);
INSTANTIATE_THEM(ZWilsonImplF);
INSTANTIATE_THEM(ZWilsonImplD);
INSTANTIATE_THEM(GparityWilsonImplF);
INSTANTIATE_THEM(GparityWilsonImplD);
INSTANTIATE_THEM(DomainWallVec5dImplF);
INSTANTIATE_THEM(DomainWallVec5dImplD);
INSTANTIATE_THEM(ZDomainWallVec5dImplF);
INSTANTIATE_THEM(ZDomainWallVec5dImplD);
INSTANTIATE_THEM(WilsonImplFH);
INSTANTIATE_THEM(WilsonImplDF);
INSTANTIATE_THEM(ZWilsonImplFH);
INSTANTIATE_THEM(ZWilsonImplDF);
INSTANTIATE_THEM(GparityWilsonImplFH);
INSTANTIATE_THEM(GparityWilsonImplDF);
INSTANTIATE_THEM(DomainWallVec5dImplFH);
INSTANTIATE_THEM(DomainWallVec5dImplDF);
INSTANTIATE_THEM(ZDomainWallVec5dImplFH);
INSTANTIATE_THEM(ZDomainWallVec5dImplDF);

}}
