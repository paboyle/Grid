    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/StaggerdKernelsHand.cc

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
#include <simd/Intel512common.h>
#include <simd/Intel512avx.h>
#include <simd/Intel512single.h>

// Interleave operations from two directions
// This looks just like a 2 spin multiply and reuse same sequence from the Wilson
// Kernel. But the spin index becomes a mu index instead.
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

#define UChi_00 %zmm12
#define UChi_01 %zmm13
#define UChi_02 %zmm14
#define UChi_10 %zmm15
#define UChi_11 %zmm16
#define UChi_12 %zmm17
#define UChi_20 %zmm18
#define UChi_21 %zmm19
#define UChi_22 %zmm20
#define UChi_30 %zmm21
#define UChi_31 %zmm22
#define UChi_32 %zmm23

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

#define pUChi_00 %%zmm12
#define pUChi_01 %%zmm13
#define pUChi_02 %%zmm14
#define pUChi_10 %%zmm15
#define pUChi_11 %%zmm16
#define pUChi_12 %%zmm17
#define pUChi_20 %%zmm18
#define pUChi_21 %%zmm19
#define pUChi_22 %%zmm20
#define pUChi_30 %%zmm21
#define pUChi_31 %%zmm22
#define pUChi_32 %%zmm23

#define T0 %zmm24
#define T1 %zmm25
#define T2 %zmm26
#define T3 %zmm27

#define MULT_ADD_LS(g0,g1,g2,g3)					\
  asm ( "movq %0, %%r8 \n\t"					\
	"movq %1, %%r9 \n\t"						\
        "movq %2, %%r10 \n\t"						\
        "movq %3, %%r11 \n\t" :  : "r"(g0), "r"(g1), "r"(g2), "r"(g3) : "%r8","%r9","%r10","%r11" );\
  asm (									\
  VSHUF(Chi_00,T0)      VSHUF(Chi_10,T1)				\
  VSHUF(Chi_20,T2)      VSHUF(Chi_30,T3)				\
  VMADDSUBIDUP(0,%r8,T0,UChi_00) VMADDSUBIDUP(0,%r9,T1,UChi_10)		\
  VMADDSUBIDUP(3,%r8,T0,UChi_01) VMADDSUBIDUP(3,%r9,T1,UChi_11)		\
  VMADDSUBIDUP(6,%r8,T0,UChi_02) VMADDSUBIDUP(6,%r9,T1,UChi_12)		\
  VMADDSUBIDUP(0,%r10,T2,UChi_20) VMADDSUBIDUP(0,%r11,T3,UChi_30)		\
  VMADDSUBIDUP(3,%r10,T2,UChi_21) VMADDSUBIDUP(3,%r11,T3,UChi_31)		\
  VMADDSUBIDUP(6,%r10,T2,UChi_22) VMADDSUBIDUP(6,%r11,T3,UChi_32)		\
  VMADDSUBRDUP(0,%r8,Chi_00,UChi_00) VMADDSUBRDUP(0,%r9,Chi_10,UChi_10) \
  VMADDSUBRDUP(3,%r8,Chi_00,UChi_01) VMADDSUBRDUP(3,%r9,Chi_10,UChi_11) \
  VMADDSUBRDUP(6,%r8,Chi_00,UChi_02) VMADDSUBRDUP(6,%r9,Chi_10,UChi_12) \
  VMADDSUBRDUP(0,%r10,Chi_20,UChi_20) VMADDSUBRDUP(0,%r11,Chi_30,UChi_30) \
  VMADDSUBRDUP(3,%r10,Chi_20,UChi_21) VMADDSUBRDUP(3,%r11,Chi_30,UChi_31) \
  VMADDSUBRDUP(6,%r10,Chi_20,UChi_22) VMADDSUBRDUP(6,%r11,Chi_30,UChi_32) \
  VSHUF(Chi_01,T0)	  VSHUF(Chi_11,T1)				\
  VSHUF(Chi_21,T2)	  VSHUF(Chi_31,T3)				\
  VMADDSUBIDUP(1,%r8,T0,UChi_00)     VMADDSUBIDUP(1,%r9,T1,UChi_10)	\
  VMADDSUBIDUP(4,%r8,T0,UChi_01)     VMADDSUBIDUP(4,%r9,T1,UChi_11)	\
  VMADDSUBIDUP(7,%r8,T0,UChi_02)     VMADDSUBIDUP(7,%r9,T1,UChi_12)	\
  VMADDSUBIDUP(1,%r10,T2,UChi_20)     VMADDSUBIDUP(1,%r11,T3,UChi_30)	\
  VMADDSUBIDUP(4,%r10,T2,UChi_21)     VMADDSUBIDUP(4,%r11,T3,UChi_31)	\
  VMADDSUBIDUP(7,%r10,T2,UChi_22)     VMADDSUBIDUP(7,%r11,T3,UChi_32)	\
  VMADDSUBRDUP(1,%r8,Chi_01,UChi_00) VMADDSUBRDUP(1,%r9,Chi_11,UChi_10) \
  VMADDSUBRDUP(4,%r8,Chi_01,UChi_01) VMADDSUBRDUP(4,%r9,Chi_11,UChi_11) \
  VMADDSUBRDUP(7,%r8,Chi_01,UChi_02) VMADDSUBRDUP(7,%r9,Chi_11,UChi_12) \
  VMADDSUBRDUP(1,%r10,Chi_21,UChi_20) VMADDSUBRDUP(1,%r11,Chi_31,UChi_30) \
  VMADDSUBRDUP(4,%r10,Chi_21,UChi_21) VMADDSUBRDUP(4,%r11,Chi_31,UChi_31) \
  VMADDSUBRDUP(7,%r10,Chi_21,UChi_22) VMADDSUBRDUP(7,%r11,Chi_31,UChi_32) \
  VSHUF(Chi_02,T0)    VSHUF(Chi_12,T1)					\
  VSHUF(Chi_22,T2)    VSHUF(Chi_32,T3)					\
  VMADDSUBIDUP(2,%r8,T0,UChi_00)     VMADDSUBIDUP(2,%r9,T1,UChi_10)     \
  VMADDSUBIDUP(5,%r8,T0,UChi_01)     VMADDSUBIDUP(5,%r9,T1,UChi_11)     \
  VMADDSUBIDUP(8,%r8,T0,UChi_02)     VMADDSUBIDUP(8,%r9,T1,UChi_12)     \
  VMADDSUBIDUP(2,%r10,T2,UChi_20)     VMADDSUBIDUP(2,%r11,T3,UChi_30)     \
  VMADDSUBIDUP(5,%r10,T2,UChi_21)     VMADDSUBIDUP(5,%r11,T3,UChi_31)     \
  VMADDSUBIDUP(8,%r10,T2,UChi_22)     VMADDSUBIDUP(8,%r11,T3,UChi_32)     \
  VMADDSUBRDUP(2,%r8,Chi_02,UChi_00) VMADDSUBRDUP(2,%r9,Chi_12,UChi_10) \
  VMADDSUBRDUP(5,%r8,Chi_02,UChi_01) VMADDSUBRDUP(5,%r9,Chi_12,UChi_11) \
  VMADDSUBRDUP(8,%r8,Chi_02,UChi_02) VMADDSUBRDUP(8,%r9,Chi_12,UChi_12) \
  VMADDSUBRDUP(2,%r10,Chi_22,UChi_20) VMADDSUBRDUP(2,%r11,Chi_32,UChi_30) \
  VMADDSUBRDUP(5,%r10,Chi_22,UChi_21) VMADDSUBRDUP(5,%r11,Chi_32,UChi_31) \
  VMADDSUBRDUP(8,%r10,Chi_22,UChi_22) VMADDSUBRDUP(8,%r11,Chi_32,UChi_32) );

#define MULT_LS(g0,g1,g2,g3)					\
  asm ( "movq %0, %%r8 \n\t"					\
	"movq %1, %%r9 \n\t"						\
        "movq %2, %%r10 \n\t"						\
        "movq %3, %%r11 \n\t" :  : "r"(g0), "r"(g1), "r"(g2), "r"(g3) : "%r8","%r9","%r10","%r11" );\
  asm (									\
  VSHUF(Chi_00,T0)      VSHUF(Chi_10,T1)				\
  VSHUF(Chi_20,T2)      VSHUF(Chi_30,T3)				\
  VMULIDUP(0,%r8,T0,UChi_00) VMULIDUP(0,%r9,T1,UChi_10)		\
  VMULIDUP(3,%r8,T0,UChi_01) VMULIDUP(3,%r9,T1,UChi_11)		\
  VMULIDUP(6,%r8,T0,UChi_02) VMULIDUP(6,%r9,T1,UChi_12)		\
  VMULIDUP(0,%r10,T2,UChi_20) VMULIDUP(0,%r11,T3,UChi_30)		\
  VMULIDUP(3,%r10,T2,UChi_21) VMULIDUP(3,%r11,T3,UChi_31)		\
  VMULIDUP(6,%r10,T2,UChi_22) VMULIDUP(6,%r11,T3,UChi_32)		\
  VMADDSUBRDUP(0,%r8,Chi_00,UChi_00) VMADDSUBRDUP(0,%r9,Chi_10,UChi_10) \
  VMADDSUBRDUP(3,%r8,Chi_00,UChi_01) VMADDSUBRDUP(3,%r9,Chi_10,UChi_11) \
  VMADDSUBRDUP(6,%r8,Chi_00,UChi_02) VMADDSUBRDUP(6,%r9,Chi_10,UChi_12) \
  VMADDSUBRDUP(0,%r10,Chi_20,UChi_20) VMADDSUBRDUP(0,%r11,Chi_30,UChi_30) \
  VMADDSUBRDUP(3,%r10,Chi_20,UChi_21) VMADDSUBRDUP(3,%r11,Chi_30,UChi_31) \
  VMADDSUBRDUP(6,%r10,Chi_20,UChi_22) VMADDSUBRDUP(6,%r11,Chi_30,UChi_32) \
  VSHUF(Chi_01,T0)	  VSHUF(Chi_11,T1)				\
  VSHUF(Chi_21,T2)	  VSHUF(Chi_31,T3)				\
  VMADDSUBIDUP(1,%r8,T0,UChi_00)     VMADDSUBIDUP(1,%r9,T1,UChi_10)	\
  VMADDSUBIDUP(4,%r8,T0,UChi_01)     VMADDSUBIDUP(4,%r9,T1,UChi_11)	\
  VMADDSUBIDUP(7,%r8,T0,UChi_02)     VMADDSUBIDUP(7,%r9,T1,UChi_12)	\
  VMADDSUBIDUP(1,%r10,T2,UChi_20)     VMADDSUBIDUP(1,%r11,T3,UChi_30)	\
  VMADDSUBIDUP(4,%r10,T2,UChi_21)     VMADDSUBIDUP(4,%r11,T3,UChi_31)	\
  VMADDSUBIDUP(7,%r10,T2,UChi_22)     VMADDSUBIDUP(7,%r11,T3,UChi_32)	\
  VMADDSUBRDUP(1,%r8,Chi_01,UChi_00) VMADDSUBRDUP(1,%r9,Chi_11,UChi_10) \
  VMADDSUBRDUP(4,%r8,Chi_01,UChi_01) VMADDSUBRDUP(4,%r9,Chi_11,UChi_11) \
  VMADDSUBRDUP(7,%r8,Chi_01,UChi_02) VMADDSUBRDUP(7,%r9,Chi_11,UChi_12) \
  VMADDSUBRDUP(1,%r10,Chi_21,UChi_20) VMADDSUBRDUP(1,%r11,Chi_31,UChi_30) \
  VMADDSUBRDUP(4,%r10,Chi_21,UChi_21) VMADDSUBRDUP(4,%r11,Chi_31,UChi_31) \
  VMADDSUBRDUP(7,%r10,Chi_21,UChi_22) VMADDSUBRDUP(7,%r11,Chi_31,UChi_32) \
  VSHUF(Chi_02,T0)    VSHUF(Chi_12,T1)					\
  VSHUF(Chi_22,T2)    VSHUF(Chi_32,T3)					\
  VMADDSUBIDUP(2,%r8,T0,UChi_00)     VMADDSUBIDUP(2,%r9,T1,UChi_10)     \
  VMADDSUBIDUP(5,%r8,T0,UChi_01)     VMADDSUBIDUP(5,%r9,T1,UChi_11)     \
  VMADDSUBIDUP(8,%r8,T0,UChi_02)     VMADDSUBIDUP(8,%r9,T1,UChi_12)     \
  VMADDSUBIDUP(2,%r10,T2,UChi_20)     VMADDSUBIDUP(2,%r11,T3,UChi_30)     \
  VMADDSUBIDUP(5,%r10,T2,UChi_21)     VMADDSUBIDUP(5,%r11,T3,UChi_31)     \
  VMADDSUBIDUP(8,%r10,T2,UChi_22)     VMADDSUBIDUP(8,%r11,T3,UChi_32)     \
  VMADDSUBRDUP(2,%r8,Chi_02,UChi_00) VMADDSUBRDUP(2,%r9,Chi_12,UChi_10) \
  VMADDSUBRDUP(5,%r8,Chi_02,UChi_01) VMADDSUBRDUP(5,%r9,Chi_12,UChi_11) \
  VMADDSUBRDUP(8,%r8,Chi_02,UChi_02) VMADDSUBRDUP(8,%r9,Chi_12,UChi_12) \
  VMADDSUBRDUP(2,%r10,Chi_22,UChi_20) VMADDSUBRDUP(2,%r11,Chi_32,UChi_30) \
  VMADDSUBRDUP(5,%r10,Chi_22,UChi_21) VMADDSUBRDUP(5,%r11,Chi_32,UChi_31) \
  VMADDSUBRDUP(8,%r10,Chi_22,UChi_22) VMADDSUBRDUP(8,%r11,Chi_32,UChi_32) );

#define LOAD_CHI(a0,a1,a2,a3)				\
  asm (									\
       "movq %0, %%r8 \n\t"						\
       VLOAD(0,%%r8,pChi_00)					\
       VLOAD(1,%%r8,pChi_01)					\
       VLOAD(2,%%r8,pChi_02)					\
       : : "r" (a0) : "%r8" );						\
  asm (									\
       "movq %0, %%r8 \n\t"						\
       VLOAD(0,%%r8,pChi_10)					\
       VLOAD(1,%%r8,pChi_11)					\
       VLOAD(2,%%r8,pChi_12)					\
       : : "r" (a1) : "%r8" );						\
  asm (									\
       "movq %0, %%r8 \n\t"						\
       VLOAD(0,%%r8,pChi_20)					\
       VLOAD(1,%%r8,pChi_21)					\
       VLOAD(2,%%r8,pChi_22)					\
       : : "r" (a2) : "%r8" );						\
  asm (									\
       "movq %0, %%r8 \n\t"						\
       VLOAD(0,%%r8,pChi_30)					\
       VLOAD(1,%%r8,pChi_31)					\
       VLOAD(2,%%r8,pChi_32)					\
       : : "r" (a3) : "%r8" );						

#define PF_CHI(a0)				\
  asm (									\
       "movq %0, %%r8 \n\t"						\
       VPREFETCH1(0,%%r8)					\
       VPREFETCH1(1,%%r8)					\
       VPREFETCH1(2,%%r8)						\
       : : "r" (a0) : "%r8" );						\
  

#define REDUCE(out)					\
  asm (							\
  VADD(UChi_00,UChi_10,UChi_00)				\
  VADD(UChi_01,UChi_11,UChi_01)				\
  VADD(UChi_02,UChi_12,UChi_02)				\
  VADD(UChi_30,UChi_20,UChi_30)				\
  VADD(UChi_31,UChi_21,UChi_31)				\
  VADD(UChi_32,UChi_22,UChi_32)				\
  VADD(UChi_00,UChi_30,UChi_00)				\
  VADD(UChi_01,UChi_31,UChi_01)				\
  VADD(UChi_02,UChi_32,UChi_02)				);	\
  asm (							\
  VSTORE(0,%0,pUChi_00)					\
  VSTORE(1,%0,pUChi_01)					\
  VSTORE(2,%0,pUChi_02)					\
  : : "r" (out) : "memory" );

#define PERMUTE_DIR(dir)			\
      permute##dir(Chi_0,Chi_0);\
      permute##dir(Chi_1,Chi_1);\
      permute##dir(Chi_2,Chi_2);

namespace Grid {
namespace QCD {

  // This is the single precision 5th direction vectorised kernel
template <class Impl>
void StaggeredKernels<Impl>::DhopSiteAsm(StencilImpl &st, LebesgueOrder &lo, 
					      DoubledGaugeField &U,
					      DoubledGaugeField &UUU,
					      SiteSpinor *buf, int sF,
					      int sU, const FermionField &in, SiteSpinor &out) 
{
  uint64_t gauge0,gauge1,gauge2,gauge3;
  uint64_t addr0,addr1,addr2,addr3;

  int o0,o1,o2,o3; // offsets
  int l0,l1,l2,l3; // local 
  int p0,p1,p2,p3; // perm
  int ptype;
  StencilEntry *SE0;
  StencilEntry *SE1;
  StencilEntry *SE2;
  StencilEntry *SE3;

  // Xp, Yp, Zp, Tp
#define PREPARE(X,Y,Z,T,skew,UU)			\
  SE0=st.GetEntry(ptype,X+skew,sF);			\
  o0 = SE0->_offset;					\
  l0 = SE0->_is_local;					\
  addr0 = l0 ?  (uint64_t) &in._odata[o0] : (uint64_t) &buf[o0];	\
  PF_CHI(addr0);							\
									\
  SE1=st.GetEntry(ptype,Y+skew,sF);			\
  o1 = SE1->_offset;					\
  l1 = SE1->_is_local;					\
  addr1 = l1 ?  (uint64_t) &in._odata[o1] : (uint64_t) &buf[o1];	\
  PF_CHI(addr1);							\
									\
  SE2=st.GetEntry(ptype,Z+skew,sF);			\
  o2 = SE2->_offset;					\
  l2 = SE2->_is_local;					\
  addr2 = l2 ?  (uint64_t) &in._odata[o2] : (uint64_t) &buf[o2];	\
  PF_CHI(addr2);							\
									\
  SE3=st.GetEntry(ptype,T+skew,sF);			\
  o3 = SE3->_offset;					\
  l3 = SE3->_is_local;					\
  addr3 = l3 ?  (uint64_t) &in._odata[o3] : (uint64_t) &buf[o3];	\
  PF_CHI(addr3);							\
  									\
  gauge0 =(uint64_t)&UU._odata[sU]( X ); \
  gauge1 =(uint64_t)&UU._odata[sU]( Y ); \
  gauge2 =(uint64_t)&UU._odata[sU]( Z ); \
  gauge3 =(uint64_t)&UU._odata[sU]( T ); 

  PREPARE(Xp,Yp,Zp,Tp,0,U);
  LOAD_CHI(addr0,addr1,addr2,addr3);
  MULT_LS(gauge0,gauge1,gauge2,gauge3);  

  PREPARE(Xm,Ym,Zm,Tm,0,U);
  LOAD_CHI(addr0,addr1,addr2,addr3);
  MULT_ADD_LS(gauge0,gauge1,gauge2,gauge3);  

  PREPARE(Xp,Yp,Zp,Tp,8,UUU);
  LOAD_CHI(addr0,addr1,addr2,addr3);
  MULT_ADD_LS(gauge0,gauge1,gauge2,gauge3);

  PREPARE(Xm,Ym,Zm,Tm,8,UUU);
  LOAD_CHI(addr0,addr1,addr2,addr3);
  MULT_ADD_LS(gauge0,gauge1,gauge2,gauge3);

  addr0 = (uint64_t) &out;
  REDUCE(addr0);
}

FermOpStaggeredTemplateInstantiate(StaggeredKernels);
FermOpStaggeredVec5dTemplateInstantiate(StaggeredKernels);

}}

