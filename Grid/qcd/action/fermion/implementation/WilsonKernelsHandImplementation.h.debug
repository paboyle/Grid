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

#pragma once

#include <Grid/qcd/action/fermion/FermionCore.h>


#undef LOAD_CHIMU
#undef LOAD_CHI
#undef MULT_2SPIN
#undef PERMUTE_DIR
#undef XP_PROJ
#undef YP_PROJ
#undef ZP_PROJ
#undef TP_PROJ
#undef XM_PROJ
#undef YM_PROJ
#undef ZM_PROJ
#undef TM_PROJ
#undef XP_RECON
#undef XP_RECON_ACCUM
#undef XM_RECON
#undef XM_RECON_ACCUM
#undef YP_RECON_ACCUM
#undef YM_RECON_ACCUM
#undef ZP_RECON_ACCUM
#undef ZM_RECON_ACCUM
#undef TP_RECON_ACCUM
#undef TM_RECON_ACCUM
#undef ZERO_RESULT
#undef Chimu_00
#undef Chimu_01
#undef Chimu_02
#undef Chimu_10
#undef Chimu_11
#undef Chimu_12
#undef Chimu_20
#undef Chimu_21
#undef Chimu_22
#undef Chimu_30
#undef Chimu_31
#undef Chimu_32
#undef HAND_STENCIL_LEG
#undef HAND_STENCIL_LEG_INT
#undef HAND_STENCIL_LEG_EXT
#undef HAND_RESULT
#undef HAND_RESULT_INT
#undef HAND_RESULT_EXT

#define REGISTER

#define LOAD_CHIMU \
  {const SiteSpinor & ref (in[offset]);	\
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
    Chimu_32=ref()(3)(2);\
    std::cout << std::endl << "DEBUG -- LOAD_CHIMU" << std::endl; \
    std::cout << "Chimu_00 -- " <<  Chimu_00 << std::endl; \
    std::cout << "Chimu_01 -- " <<  Chimu_01 << std::endl; \
    std::cout << "Chimu_02 -- " <<  Chimu_02 << std::endl; \
    std::cout << "Chimu_10 -- " <<  Chimu_10 << std::endl; \
    std::cout << "Chimu_11 -- " <<  Chimu_11 << std::endl; \
    std::cout << "Chimu_12 -- " <<  Chimu_12 << std::endl; \
    std::cout << "Chimu_20 -- " <<  Chimu_20 << std::endl; \
    std::cout << "Chimu_21 -- " <<  Chimu_21 << std::endl; \
    std::cout << "Chimu_22 -- " <<  Chimu_22 << std::endl; \
    std::cout << "Chimu_30 -- " <<  Chimu_30 << std::endl; \
    std::cout << "Chimu_31 -- " <<  Chimu_31 << std::endl; \
    std::cout << "Chimu_32 -- " <<  Chimu_32 << std::endl; \
}

#define LOAD_CHI\
  {const SiteHalfSpinor &ref(buf[offset]);	\
    Chi_00 = ref()(0)(0);\
    Chi_01 = ref()(0)(1);\
    Chi_02 = ref()(0)(2);\
    Chi_10 = ref()(1)(0);\
    Chi_11 = ref()(1)(1);\
    Chi_12 = ref()(1)(2);\
    std::cout << std::endl << "DEBUG -- LOAD_CHI" << std::endl; \
    std::cout << "Chi_00 -- " <<  Chi_00 << std::endl; \
    std::cout << "Chi_01 -- " <<  Chi_01 << std::endl; \
    std::cout << "Chi_02 -- " <<  Chi_02 << std::endl; \
    std::cout << "Chi_10 -- " <<  Chi_10 << std::endl; \
    std::cout << "Chi_11 -- " <<  Chi_11 << std::endl; \
    std::cout << "Chi_12 -- " <<  Chi_12 << std::endl; \
  }

// To splat or not to splat depends on the implementation
#define MULT_2SPIN(A)\
  {auto & ref(U[sU](A));			\
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
    UChi_12+= U_20*Chi_12;\
    std::cout << std::endl << "DEBUG -- MULT_2SPIN" << std::endl; \
    std::cout << "UChi_00 -- " <<  UChi_00 << std::endl; \
    std::cout << "UChi_01 -- " <<  UChi_01 << std::endl; \
    std::cout << "UChi_02 -- " <<  UChi_02 << std::endl; \
    std::cout << "UChi_10 -- " <<  UChi_10 << std::endl; \
    std::cout << "UChi_11 -- " <<  UChi_11 << std::endl; \
    std::cout << "UChi_12 -- " <<  UChi_12 << std::endl; \
    }


#define PERMUTE_DIR(dir)			\
std::cout << std::endl << "DEBUG -- PERM PRE" << std::endl; \
std::cout << "Chi_00 -- " <<  Chi_00 << std::endl; \
std::cout << "Chi_01 -- " <<  Chi_01 << std::endl; \
std::cout << "Chi_02 -- " <<  Chi_02 << std::endl; \
std::cout << "Chi_10 -- " <<  Chi_10 << std::endl; \
std::cout << "Chi_11 -- " <<  Chi_11 << std::endl; \
std::cout << "Chi_12 -- " <<  Chi_12 << std::endl; \
      permute##dir(Chi_00,Chi_00);\
      permute##dir(Chi_01,Chi_01);\
      permute##dir(Chi_02,Chi_02);\
      permute##dir(Chi_10,Chi_10);\
      permute##dir(Chi_11,Chi_11);\
      permute##dir(Chi_12,Chi_12);\
  std::cout << std::endl << "DEBUG -- PERM POST" << std::endl; \
  std::cout << "Chi_00 -- " <<  Chi_00 << std::endl; \
  std::cout << "Chi_01 -- " <<  Chi_01 << std::endl; \
  std::cout << "Chi_02 -- " <<  Chi_02 << std::endl; \
  std::cout << "Chi_10 -- " <<  Chi_10 << std::endl; \
  std::cout << "Chi_11 -- " <<  Chi_11 << std::endl; \
  std::cout << "Chi_12 -- " <<  Chi_12 << std::endl;

//      hspin(0)=fspin(0)+timesI(fspin(3));
//      hspin(1)=fspin(1)+timesI(fspin(2));
#define XP_PROJ \
    Chi_00 = Chimu_00+timesI(Chimu_30);\
    Chi_01 = Chimu_01+timesI(Chimu_31);\
    Chi_02 = Chimu_02+timesI(Chimu_32);\
    Chi_10 = Chimu_10+timesI(Chimu_20);\
    Chi_11 = Chimu_11+timesI(Chimu_21);\
    Chi_12 = Chimu_12+timesI(Chimu_22);\
    std::cout << std::endl << "DEBUG -- XP_PROJ" << std::endl; \
    std::cout << "Chi_00 -- " <<  Chi_00 << std::endl; \
    std::cout << "Chi_01 -- " <<  Chi_01 << std::endl; \
    std::cout << "Chi_02 -- " <<  Chi_02 << std::endl; \
    std::cout << "Chi_10 -- " <<  Chi_10 << std::endl; \
    std::cout << "Chi_11 -- " <<  Chi_11 << std::endl; \
    std::cout << "Chi_12 -- " <<  Chi_12 << std::endl;

#define YP_PROJ \
    Chi_00 = Chimu_00-Chimu_30;\
    Chi_01 = Chimu_01-Chimu_31;\
    Chi_02 = Chimu_02-Chimu_32;\
    Chi_10 = Chimu_10+Chimu_20;\
    Chi_11 = Chimu_11+Chimu_21;\
    Chi_12 = Chimu_12+Chimu_22;\
    std::cout << std::endl << "DEBUG -- YP_PROJ" << std::endl; \
    std::cout << "Chi_00 -- " <<  Chi_00 << std::endl; \
    std::cout << "Chi_01 -- " <<  Chi_01 << std::endl; \
    std::cout << "Chi_02 -- " <<  Chi_02 << std::endl; \
    std::cout << "Chi_10 -- " <<  Chi_10 << std::endl; \
    std::cout << "Chi_11 -- " <<  Chi_11 << std::endl; \
    std::cout << "Chi_12 -- " <<  Chi_12 << std::endl;

#define ZP_PROJ \
  Chi_00 = Chimu_00+timesI(Chimu_20);		\
  Chi_01 = Chimu_01+timesI(Chimu_21);		\
  Chi_02 = Chimu_02+timesI(Chimu_22);		\
  Chi_10 = Chimu_10-timesI(Chimu_30);		\
  Chi_11 = Chimu_11-timesI(Chimu_31);		\
  Chi_12 = Chimu_12-timesI(Chimu_32);\
  std::cout << std::endl << "DEBUG -- ZP_PROJ" << std::endl; \
  std::cout << "Chi_00 -- " <<  Chi_00 << std::endl; \
  std::cout << "Chi_01 -- " <<  Chi_01 << std::endl; \
  std::cout << "Chi_02 -- " <<  Chi_02 << std::endl; \
  std::cout << "Chi_10 -- " <<  Chi_10 << std::endl; \
  std::cout << "Chi_11 -- " <<  Chi_11 << std::endl; \
  std::cout << "Chi_12 -- " <<  Chi_12 << std::endl;

#define TP_PROJ \
  Chi_00 = Chimu_00+Chimu_20;		\
  Chi_01 = Chimu_01+Chimu_21;		\
  Chi_02 = Chimu_02+Chimu_22;		\
  Chi_10 = Chimu_10+Chimu_30;		\
  Chi_11 = Chimu_11+Chimu_31;		\
  Chi_12 = Chimu_12+Chimu_32;\
  std::cout << std::endl << "DEBUG -- TP_PROJ" << std::endl; \
  std::cout << "Chi_00 -- " <<  Chi_00 << std::endl; \
  std::cout << "Chi_01 -- " <<  Chi_01 << std::endl; \
  std::cout << "Chi_02 -- " <<  Chi_02 << std::endl; \
  std::cout << "Chi_10 -- " <<  Chi_10 << std::endl; \
  std::cout << "Chi_11 -- " <<  Chi_11 << std::endl; \
  std::cout << "Chi_12 -- " <<  Chi_12 << std::endl;


//      hspin(0)=fspin(0)-timesI(fspin(3));
//      hspin(1)=fspin(1)-timesI(fspin(2));
#define XM_PROJ \
    Chi_00 = Chimu_00-timesI(Chimu_30);\
    Chi_01 = Chimu_01-timesI(Chimu_31);\
    Chi_02 = Chimu_02-timesI(Chimu_32);\
    Chi_10 = Chimu_10-timesI(Chimu_20);\
    Chi_11 = Chimu_11-timesI(Chimu_21);\
    Chi_12 = Chimu_12-timesI(Chimu_22);\
    std::cout << std::endl << "DEBUG -- XM_PROJ" << std::endl; \
    std::cout << "Chi_00 -- " <<  Chi_00 << std::endl; \
    std::cout << "Chi_01 -- " <<  Chi_01 << std::endl; \
    std::cout << "Chi_02 -- " <<  Chi_02 << std::endl; \
    std::cout << "Chi_10 -- " <<  Chi_10 << std::endl; \
    std::cout << "Chi_11 -- " <<  Chi_11 << std::endl; \
    std::cout << "Chi_12 -- " <<  Chi_12 << std::endl;

#define YM_PROJ \
    Chi_00 = Chimu_00+Chimu_30;\
    Chi_01 = Chimu_01+Chimu_31;\
    Chi_02 = Chimu_02+Chimu_32;\
    Chi_10 = Chimu_10-Chimu_20;\
    Chi_11 = Chimu_11-Chimu_21;\
    Chi_12 = Chimu_12-Chimu_22;\
    std::cout << std::endl << "DEBUG -- YM_PROJ" << std::endl; \
    std::cout << "Chi_00 -- " <<  Chi_00 << std::endl; \
    std::cout << "Chi_01 -- " <<  Chi_01 << std::endl; \
    std::cout << "Chi_02 -- " <<  Chi_02 << std::endl; \
    std::cout << "Chi_10 -- " <<  Chi_10 << std::endl; \
    std::cout << "Chi_11 -- " <<  Chi_11 << std::endl; \
    std::cout << "Chi_12 -- " <<  Chi_12 << std::endl;

#define ZM_PROJ \
  Chi_00 = Chimu_00-timesI(Chimu_20);		\
  Chi_01 = Chimu_01-timesI(Chimu_21);		\
  Chi_02 = Chimu_02-timesI(Chimu_22);		\
  Chi_10 = Chimu_10+timesI(Chimu_30);		\
  Chi_11 = Chimu_11+timesI(Chimu_31);		\
  Chi_12 = Chimu_12+timesI(Chimu_32);\
  std::cout << std::endl << "DEBUG -- ZM_PROJ" << std::endl; \
  std::cout << "Chi_00 -- " <<  Chi_00 << std::endl; \
  std::cout << "Chi_01 -- " <<  Chi_01 << std::endl; \
  std::cout << "Chi_02 -- " <<  Chi_02 << std::endl; \
  std::cout << "Chi_10 -- " <<  Chi_10 << std::endl; \
  std::cout << "Chi_11 -- " <<  Chi_11 << std::endl; \
  std::cout << "Chi_12 -- " <<  Chi_12 << std::endl;

#define TM_PROJ \
  Chi_00 = Chimu_00-Chimu_20;		\
  Chi_01 = Chimu_01-Chimu_21;		\
  Chi_02 = Chimu_02-Chimu_22;		\
  Chi_10 = Chimu_10-Chimu_30;		\
  Chi_11 = Chimu_11-Chimu_31;		\
  Chi_12 = Chimu_12-Chimu_32;\
  std::cout << std::endl << "DEBUG -- TM_PROJ" << std::endl; \
  std::cout << "Chi_00 -- " <<  Chi_00 << std::endl; \
  std::cout << "Chi_01 -- " <<  Chi_01 << std::endl; \
  std::cout << "Chi_02 -- " <<  Chi_02 << std::endl; \
  std::cout << "Chi_10 -- " <<  Chi_10 << std::endl; \
  std::cout << "Chi_11 -- " <<  Chi_11 << std::endl; \
  std::cout << "Chi_12 -- " <<  Chi_12 << std::endl;

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
  result_32 = timesMinusI(UChi_02);\
  std::cout << std::endl << "DEBUG -- XP_RECON" << std::endl; \
  std::cout << "result_00 -- " <<  result_00 << std::endl; \
  std::cout << "result_01 -- " <<  result_01 << std::endl; \
  std::cout << "result_02 -- " <<  result_02 << std::endl; \
  std::cout << "result_10 -- " <<  result_10 << std::endl; \
  std::cout << "result_11 -- " <<  result_11 << std::endl; \
  std::cout << "result_12 -- " <<  result_12 << std::endl; \
  std::cout << "result_20 -- " <<  result_20 << std::endl; \
  std::cout << "result_21 -- " <<  result_21 << std::endl; \
  std::cout << "result_22 -- " <<  result_22 << std::endl; \
  std::cout << "result_30 -- " <<  result_30 << std::endl; \
  std::cout << "result_31 -- " <<  result_31 << std::endl; \
  std::cout << "result_32 -- " <<  result_32 << std::endl;

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
  result_32-=timesI(UChi_02);\
  std::cout << std::endl << "DEBUG -- XP_RECON_ACCUM" << std::endl; \
  std::cout << "result_00 -- " <<  result_00 << std::endl; \
  std::cout << "result_01 -- " <<  result_01 << std::endl; \
  std::cout << "result_02 -- " <<  result_02 << std::endl; \
  std::cout << "result_10 -- " <<  result_10 << std::endl; \
  std::cout << "result_11 -- " <<  result_11 << std::endl; \
  std::cout << "result_12 -- " <<  result_12 << std::endl; \
  std::cout << "result_20 -- " <<  result_20 << std::endl; \
  std::cout << "result_21 -- " <<  result_21 << std::endl; \
  std::cout << "result_22 -- " <<  result_22 << std::endl; \
  std::cout << "result_30 -- " <<  result_30 << std::endl; \
  std::cout << "result_31 -- " <<  result_31 << std::endl; \
  std::cout << "result_32 -- " <<  result_32 << std::endl;

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
  result_32 = timesI(UChi_02);\
  std::cout << std::endl << "DEBUG -- XM_RECON" << std::endl; \
  std::cout << "result_00 -- " <<  result_00 << std::endl; \
  std::cout << "result_01 -- " <<  result_01 << std::endl; \
  std::cout << "result_02 -- " <<  result_02 << std::endl; \
  std::cout << "result_10 -- " <<  result_10 << std::endl; \
  std::cout << "result_11 -- " <<  result_11 << std::endl; \
  std::cout << "result_12 -- " <<  result_12 << std::endl; \
  std::cout << "result_20 -- " <<  result_20 << std::endl; \
  std::cout << "result_21 -- " <<  result_21 << std::endl; \
  std::cout << "result_22 -- " <<  result_22 << std::endl; \
  std::cout << "result_30 -- " <<  result_30 << std::endl; \
  std::cout << "result_31 -- " <<  result_31 << std::endl; \
  std::cout << "result_32 -- " <<  result_32 << std::endl;

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
  result_32+= timesI(UChi_02);\
  std::cout << std::endl << "DEBUG -- XM_RECON_ACCUM" << std::endl; \
  std::cout << "result_00 -- " <<  result_00 << std::endl; \
  std::cout << "result_01 -- " <<  result_01 << std::endl; \
  std::cout << "result_02 -- " <<  result_02 << std::endl; \
  std::cout << "result_10 -- " <<  result_10 << std::endl; \
  std::cout << "result_11 -- " <<  result_11 << std::endl; \
  std::cout << "result_12 -- " <<  result_12 << std::endl; \
  std::cout << "result_20 -- " <<  result_20 << std::endl; \
  std::cout << "result_21 -- " <<  result_21 << std::endl; \
  std::cout << "result_22 -- " <<  result_22 << std::endl; \
  std::cout << "result_30 -- " <<  result_30 << std::endl; \
  std::cout << "result_31 -- " <<  result_31 << std::endl; \
  std::cout << "result_32 -- " <<  result_32 << std::endl;

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
  result_32-= UChi_02;\
  std::cout << std::endl << "DEBUG -- YP_RECON_ACCUM" << std::endl; \
  std::cout << "result_00 -- " <<  result_00 << std::endl; \
  std::cout << "result_01 -- " <<  result_01 << std::endl; \
  std::cout << "result_02 -- " <<  result_02 << std::endl; \
  std::cout << "result_10 -- " <<  result_10 << std::endl; \
  std::cout << "result_11 -- " <<  result_11 << std::endl; \
  std::cout << "result_12 -- " <<  result_12 << std::endl; \
  std::cout << "result_20 -- " <<  result_20 << std::endl; \
  std::cout << "result_21 -- " <<  result_21 << std::endl; \
  std::cout << "result_22 -- " <<  result_22 << std::endl; \
  std::cout << "result_30 -- " <<  result_30 << std::endl; \
  std::cout << "result_31 -- " <<  result_31 << std::endl; \
  std::cout << "result_32 -- " <<  result_32 << std::endl;

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
  result_32+= UChi_02;\
  std::cout << std::endl << "DEBUG -- YM_RECON_ACCUM" << std::endl; \
  std::cout << "result_00 -- " <<  result_00 << std::endl; \
  std::cout << "result_01 -- " <<  result_01 << std::endl; \
  std::cout << "result_02 -- " <<  result_02 << std::endl; \
  std::cout << "result_10 -- " <<  result_10 << std::endl; \
  std::cout << "result_11 -- " <<  result_11 << std::endl; \
  std::cout << "result_12 -- " <<  result_12 << std::endl; \
  std::cout << "result_20 -- " <<  result_20 << std::endl; \
  std::cout << "result_21 -- " <<  result_21 << std::endl; \
  std::cout << "result_22 -- " <<  result_22 << std::endl; \
  std::cout << "result_30 -- " <<  result_30 << std::endl; \
  std::cout << "result_31 -- " <<  result_31 << std::endl; \
  std::cout << "result_32 -- " <<  result_32 << std::endl;

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
  result_32+= timesI(UChi_12);\
  std::cout << std::endl << "DEBUG -- ZP_RECON_ACCUM" << std::endl; \
  std::cout << "result_00 -- " <<  result_00 << std::endl; \
  std::cout << "result_01 -- " <<  result_01 << std::endl; \
  std::cout << "result_02 -- " <<  result_02 << std::endl; \
  std::cout << "result_10 -- " <<  result_10 << std::endl; \
  std::cout << "result_11 -- " <<  result_11 << std::endl; \
  std::cout << "result_12 -- " <<  result_12 << std::endl; \
  std::cout << "result_20 -- " <<  result_20 << std::endl; \
  std::cout << "result_21 -- " <<  result_21 << std::endl; \
  std::cout << "result_22 -- " <<  result_22 << std::endl; \
  std::cout << "result_30 -- " <<  result_30 << std::endl; \
  std::cout << "result_31 -- " <<  result_31 << std::endl; \
  std::cout << "result_32 -- " <<  result_32 << std::endl;

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
  result_32-= timesI(UChi_12);\
  std::cout << std::endl << "DEBUG -- ZM_RECON_ACCUM" << std::endl; \
  std::cout << "result_00 -- " <<  result_00 << std::endl; \
  std::cout << "result_01 -- " <<  result_01 << std::endl; \
  std::cout << "result_02 -- " <<  result_02 << std::endl; \
  std::cout << "result_10 -- " <<  result_10 << std::endl; \
  std::cout << "result_11 -- " <<  result_11 << std::endl; \
  std::cout << "result_12 -- " <<  result_12 << std::endl; \
  std::cout << "result_20 -- " <<  result_20 << std::endl; \
  std::cout << "result_21 -- " <<  result_21 << std::endl; \
  std::cout << "result_22 -- " <<  result_22 << std::endl; \
  std::cout << "result_30 -- " <<  result_30 << std::endl; \
  std::cout << "result_31 -- " <<  result_31 << std::endl; \
  std::cout << "result_32 -- " <<  result_32 << std::endl;

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
  result_32+= UChi_12;\
  std::cout << std::endl << "DEBUG -- TP_RECON_ACCUM" << std::endl; \
  std::cout << "result_00 -- " <<  result_00 << std::endl; \
  std::cout << "result_01 -- " <<  result_01 << std::endl; \
  std::cout << "result_02 -- " <<  result_02 << std::endl; \
  std::cout << "result_10 -- " <<  result_10 << std::endl; \
  std::cout << "result_11 -- " <<  result_11 << std::endl; \
  std::cout << "result_12 -- " <<  result_12 << std::endl; \
  std::cout << "result_20 -- " <<  result_20 << std::endl; \
  std::cout << "result_21 -- " <<  result_21 << std::endl; \
  std::cout << "result_22 -- " <<  result_22 << std::endl; \
  std::cout << "result_30 -- " <<  result_30 << std::endl; \
  std::cout << "result_31 -- " <<  result_31 << std::endl; \
  std::cout << "result_32 -- " <<  result_32 << std::endl;

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
  result_32-= UChi_12;\
  std::cout << std::endl << "DEBUG -- TM_RECON_ACCUM" << std::endl; \
  std::cout << "result_00 -- " <<  result_00 << std::endl; \
  std::cout << "result_01 -- " <<  result_01 << std::endl; \
  std::cout << "result_02 -- " <<  result_02 << std::endl; \
  std::cout << "result_10 -- " <<  result_10 << std::endl; \
  std::cout << "result_11 -- " <<  result_11 << std::endl; \
  std::cout << "result_12 -- " <<  result_12 << std::endl; \
  std::cout << "result_20 -- " <<  result_20 << std::endl; \
  std::cout << "result_21 -- " <<  result_21 << std::endl; \
  std::cout << "result_22 -- " <<  result_22 << std::endl; \
  std::cout << "result_30 -- " <<  result_30 << std::endl; \
  std::cout << "result_31 -- " <<  result_31 << std::endl; \
  std::cout << "result_32 -- " <<  result_32 << std::endl;

#define HAND_STENCIL_LEG(PROJ,PERM,DIR,RECON)	\
  SE=st.GetEntry(ptype,DIR,ss);			\
  offset = SE->_offset;				\
  local  = SE->_is_local;			\
  perm   = SE->_permute;			\
  if ( local ) {				\
    LOAD_CHIMU;					\
    PROJ;					\
    if ( perm) {				\
      PERMUTE_DIR(PERM);			\
    }						\
  } else {					\
    LOAD_CHI;					\
  }						\
  MULT_2SPIN(DIR);				\
  RECON;

#define HAND_STENCIL_LEG_INT(PROJ,PERM,DIR,RECON)	\
  SE=st.GetEntry(ptype,DIR,ss);			\
  offset = SE->_offset;				\
  local  = SE->_is_local;			\
  perm   = SE->_permute;			\
  if ( local ) {				\
    LOAD_CHIMU;					\
    PROJ;					\
    if ( perm) {				\
      PERMUTE_DIR(PERM);			\
    }						\
  } else if ( st.same_node[DIR] ) {		\
    LOAD_CHI;					\
  }						\
  if (local || st.same_node[DIR] ) {		\
    MULT_2SPIN(DIR);				\
    RECON;					\
  }

#define HAND_STENCIL_LEG_EXT(PROJ,PERM,DIR,RECON)	\
  SE=st.GetEntry(ptype,DIR,ss);			\
  offset = SE->_offset;				\
  if((!SE->_is_local)&&(!st.same_node[DIR]) ) {	\
    LOAD_CHI;					\
    MULT_2SPIN(DIR);				\
    RECON;					\
    nmu++;					\
  }

#define HAND_RESULT(ss)				\
  {						\
    SiteSpinor & ref (out[ss]);		\
    vstream(ref()(0)(0),result_00);		\
    vstream(ref()(0)(1),result_01);		\
    vstream(ref()(0)(2),result_02);		\
    vstream(ref()(1)(0),result_10);		\
    vstream(ref()(1)(1),result_11);		\
    vstream(ref()(1)(2),result_12);		\
    vstream(ref()(2)(0),result_20);		\
    vstream(ref()(2)(1),result_21);		\
    vstream(ref()(2)(2),result_22);		\
    vstream(ref()(3)(0),result_30);		\
    vstream(ref()(3)(1),result_31);		\
    vstream(ref()(3)(2),result_32);		\
    std::cout << std::endl << "DEBUG -- RESULT" << std::endl; \
    std::cout << "result_00 -- " <<  result_00 << std::endl; \
    std::cout << "result_01 -- " <<  result_01 << std::endl; \
    std::cout << "result_02 -- " <<  result_02 << std::endl; \
    std::cout << "result_10 -- " <<  result_10 << std::endl; \
    std::cout << "result_11 -- " <<  result_11 << std::endl; \
    std::cout << "result_12 -- " <<  result_12 << std::endl; \
    std::cout << "result_20 -- " <<  result_20 << std::endl; \
    std::cout << "result_21 -- " <<  result_21 << std::endl; \
    std::cout << "result_22 -- " <<  result_22 << std::endl; \
    std::cout << "result_30 -- " <<  result_30 << std::endl; \
    std::cout << "result_31 -- " <<  result_31 << std::endl; \
    std::cout << "result_32 -- " <<  result_32 << std::endl;\
  }

#define HAND_RESULT_EXT(ss)			\
  if (nmu){					\
    SiteSpinor & ref (out[ss]);		\
    ref()(0)(0)+=result_00;		\
    ref()(0)(1)+=result_01;		\
    ref()(0)(2)+=result_02;		\
    ref()(1)(0)+=result_10;		\
    ref()(1)(1)+=result_11;		\
    ref()(1)(2)+=result_12;		\
    ref()(2)(0)+=result_20;		\
    ref()(2)(1)+=result_21;		\
    ref()(2)(2)+=result_22;		\
    ref()(3)(0)+=result_30;		\
    ref()(3)(1)+=result_31;		\
    ref()(3)(2)+=result_32;		\
    std::cout << std::endl << "DEBUG -- RESULT EXT" << std::endl; \
    std::cout << "result_00 -- " <<  result_00 << std::endl; \
    std::cout << "result_01 -- " <<  result_01 << std::endl; \
    std::cout << "result_02 -- " <<  result_02 << std::endl; \
    std::cout << "result_10 -- " <<  result_10 << std::endl; \
    std::cout << "result_11 -- " <<  result_11 << std::endl; \
    std::cout << "result_12 -- " <<  result_12 << std::endl; \
    std::cout << "result_20 -- " <<  result_20 << std::endl; \
    std::cout << "result_21 -- " <<  result_21 << std::endl; \
    std::cout << "result_22 -- " <<  result_22 << std::endl; \
    std::cout << "result_30 -- " <<  result_30 << std::endl; \
    std::cout << "result_31 -- " <<  result_31 << std::endl; \
    std::cout << "result_32 -- " <<  result_32 << std::endl;\
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
  Simd U_21;\
  Simd debugreg;\
  svbool_t pg1;        \
  pg1 = svptrue_b64();        \

#define ZERO_RESULT				\
  result_00=Zero();				\
  result_01=Zero();				\
  result_02=Zero();				\
  result_10=Zero();				\
  result_11=Zero();				\
  result_12=Zero();				\
  result_20=Zero();				\
  result_21=Zero();				\
  result_22=Zero();				\
  result_30=Zero();				\
  result_31=Zero();				\
  result_32=Zero();

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

NAMESPACE_BEGIN(Grid);

template<class Impl> void
WilsonKernels<Impl>::HandDhopSite(StencilView &st, DoubledGaugeFieldView &U,SiteHalfSpinor  *buf,
				  int ss,int sU,const FermionFieldView &in, FermionFieldView &out)
{
// T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  HAND_DECLARATIONS(ignore);

  int offset,local,perm, ptype;
  StencilEntry *SE;

  HAND_STENCIL_LEG(XM_PROJ,3,Xp,XM_RECON);
  HAND_STENCIL_LEG(YM_PROJ,2,Yp,YM_RECON_ACCUM);
  HAND_STENCIL_LEG(ZM_PROJ,1,Zp,ZM_RECON_ACCUM);
  HAND_STENCIL_LEG(TM_PROJ,0,Tp,TM_RECON_ACCUM);
  HAND_STENCIL_LEG(XP_PROJ,3,Xm,XP_RECON_ACCUM);
  HAND_STENCIL_LEG(YP_PROJ,2,Ym,YP_RECON_ACCUM);
  HAND_STENCIL_LEG(ZP_PROJ,1,Zm,ZP_RECON_ACCUM);
  HAND_STENCIL_LEG(TP_PROJ,0,Tm,TP_RECON_ACCUM);
  HAND_RESULT(ss);
}

template<class Impl>
void WilsonKernels<Impl>::HandDhopSiteDag(StencilView &st,DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
					  int ss,int sU,const FermionFieldView &in, FermionFieldView &out)
{
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  HAND_DECLARATIONS(ignore);

  StencilEntry *SE;
  int offset,local,perm, ptype;

  HAND_STENCIL_LEG(XP_PROJ,3,Xp,XP_RECON);
  HAND_STENCIL_LEG(YP_PROJ,2,Yp,YP_RECON_ACCUM);
  HAND_STENCIL_LEG(ZP_PROJ,1,Zp,ZP_RECON_ACCUM);
  HAND_STENCIL_LEG(TP_PROJ,0,Tp,TP_RECON_ACCUM);
  HAND_STENCIL_LEG(XM_PROJ,3,Xm,XM_RECON_ACCUM);
  HAND_STENCIL_LEG(YM_PROJ,2,Ym,YM_RECON_ACCUM);
  HAND_STENCIL_LEG(ZM_PROJ,1,Zm,ZM_RECON_ACCUM);
  HAND_STENCIL_LEG(TM_PROJ,0,Tm,TM_RECON_ACCUM);
  HAND_RESULT(ss);
}

template<class Impl> void
WilsonKernels<Impl>::HandDhopSiteInt(StencilView &st,DoubledGaugeFieldView &U,SiteHalfSpinor  *buf,
					  int ss,int sU,const FermionFieldView &in, FermionFieldView &out)
{
// T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  HAND_DECLARATIONS(ignore);

  int offset,local,perm, ptype;
  StencilEntry *SE;
  ZERO_RESULT;
  HAND_STENCIL_LEG_INT(XM_PROJ,3,Xp,XM_RECON_ACCUM);
  HAND_STENCIL_LEG_INT(YM_PROJ,2,Yp,YM_RECON_ACCUM);
  HAND_STENCIL_LEG_INT(ZM_PROJ,1,Zp,ZM_RECON_ACCUM);
  HAND_STENCIL_LEG_INT(TM_PROJ,0,Tp,TM_RECON_ACCUM);
  HAND_STENCIL_LEG_INT(XP_PROJ,3,Xm,XP_RECON_ACCUM);
  HAND_STENCIL_LEG_INT(YP_PROJ,2,Ym,YP_RECON_ACCUM);
  HAND_STENCIL_LEG_INT(ZP_PROJ,1,Zm,ZP_RECON_ACCUM);
  HAND_STENCIL_LEG_INT(TP_PROJ,0,Tm,TP_RECON_ACCUM);
  HAND_RESULT(ss);
}

template<class Impl>
void WilsonKernels<Impl>::HandDhopSiteDagInt(StencilView &st,DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
						  int ss,int sU,const FermionFieldView &in, FermionFieldView &out)
{
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  HAND_DECLARATIONS(ignore);

  StencilEntry *SE;
  int offset,local,perm, ptype;
  ZERO_RESULT;
  HAND_STENCIL_LEG_INT(XP_PROJ,3,Xp,XP_RECON_ACCUM);
  HAND_STENCIL_LEG_INT(YP_PROJ,2,Yp,YP_RECON_ACCUM);
  HAND_STENCIL_LEG_INT(ZP_PROJ,1,Zp,ZP_RECON_ACCUM);
  HAND_STENCIL_LEG_INT(TP_PROJ,0,Tp,TP_RECON_ACCUM);
  HAND_STENCIL_LEG_INT(XM_PROJ,3,Xm,XM_RECON_ACCUM);
  HAND_STENCIL_LEG_INT(YM_PROJ,2,Ym,YM_RECON_ACCUM);
  HAND_STENCIL_LEG_INT(ZM_PROJ,1,Zm,ZM_RECON_ACCUM);
  HAND_STENCIL_LEG_INT(TM_PROJ,0,Tm,TM_RECON_ACCUM);
  HAND_RESULT(ss);
}

template<class Impl> void
WilsonKernels<Impl>::HandDhopSiteExt(StencilView &st,DoubledGaugeFieldView &U,SiteHalfSpinor  *buf,
					  int ss,int sU,const FermionFieldView &in, FermionFieldView &out)
{
// T==0, Z==1, Y==2, Z==3 expect 1,2,2,2 simd layout etc...
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  HAND_DECLARATIONS(ignore);

  int offset, ptype;
  StencilEntry *SE;
  int nmu=0;
  ZERO_RESULT;
  HAND_STENCIL_LEG_EXT(XM_PROJ,3,Xp,XM_RECON_ACCUM);
  HAND_STENCIL_LEG_EXT(YM_PROJ,2,Yp,YM_RECON_ACCUM);
  HAND_STENCIL_LEG_EXT(ZM_PROJ,1,Zp,ZM_RECON_ACCUM);
  HAND_STENCIL_LEG_EXT(TM_PROJ,0,Tp,TM_RECON_ACCUM);
  HAND_STENCIL_LEG_EXT(XP_PROJ,3,Xm,XP_RECON_ACCUM);
  HAND_STENCIL_LEG_EXT(YP_PROJ,2,Ym,YP_RECON_ACCUM);
  HAND_STENCIL_LEG_EXT(ZP_PROJ,1,Zm,ZP_RECON_ACCUM);
  HAND_STENCIL_LEG_EXT(TP_PROJ,0,Tm,TP_RECON_ACCUM);
  HAND_RESULT_EXT(ss);
}

template<class Impl>
void WilsonKernels<Impl>::HandDhopSiteDagExt(StencilView &st,DoubledGaugeFieldView &U,SiteHalfSpinor *buf,
						  int ss,int sU,const FermionFieldView &in, FermionFieldView &out)
{
  typedef typename Simd::scalar_type S;
  typedef typename Simd::vector_type V;

  HAND_DECLARATIONS(ignore);

  StencilEntry *SE;
  int offset, ptype;
  int nmu=0;
  ZERO_RESULT;
  HAND_STENCIL_LEG_EXT(XP_PROJ,3,Xp,XP_RECON_ACCUM);
  HAND_STENCIL_LEG_EXT(YP_PROJ,2,Yp,YP_RECON_ACCUM);
  HAND_STENCIL_LEG_EXT(ZP_PROJ,1,Zp,ZP_RECON_ACCUM);
  HAND_STENCIL_LEG_EXT(TP_PROJ,0,Tp,TP_RECON_ACCUM);
  HAND_STENCIL_LEG_EXT(XM_PROJ,3,Xm,XM_RECON_ACCUM);
  HAND_STENCIL_LEG_EXT(YM_PROJ,2,Ym,YM_RECON_ACCUM);
  HAND_STENCIL_LEG_EXT(ZM_PROJ,1,Zm,ZM_RECON_ACCUM);
  HAND_STENCIL_LEG_EXT(TM_PROJ,0,Tm,TM_RECON_ACCUM);
  HAND_RESULT_EXT(ss);
}

////////////// Wilson ; uses this implementation /////////////////////

NAMESPACE_END(Grid);
#undef LOAD_CHIMU
#undef LOAD_CHI
#undef MULT_2SPIN
#undef PERMUTE_DIR
#undef XP_PROJ
#undef YP_PROJ
#undef ZP_PROJ
#undef TP_PROJ
#undef XM_PROJ
#undef YM_PROJ
#undef ZM_PROJ
#undef TM_PROJ
#undef XP_RECON
#undef XP_RECON_ACCUM
#undef XM_RECON
#undef XM_RECON_ACCUM
#undef YP_RECON_ACCUM
#undef YM_RECON_ACCUM
#undef ZP_RECON_ACCUM
#undef ZM_RECON_ACCUM
#undef TP_RECON_ACCUM
#undef TM_RECON_ACCUM
#undef ZERO_RESULT
#undef Chimu_00
#undef Chimu_01
#undef Chimu_02
#undef Chimu_10
#undef Chimu_11
#undef Chimu_12
#undef Chimu_20
#undef Chimu_21
#undef Chimu_22
#undef Chimu_30
#undef Chimu_31
#undef Chimu_32
#undef HAND_STENCIL_LEG
#undef HAND_STENCIL_LEG_INT
#undef HAND_STENCIL_LEG_EXT
#undef HAND_RESULT
#undef HAND_RESULT_INT
#undef HAND_RESULT_EXT
