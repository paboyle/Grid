/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/simd/Avx512Asm.h

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
#ifndef GRID_ASM_INTEL_512_QCD_H
#define GRID_ASM_INTEL_512_QCD_H

//////////////////////////////////////////////////////////////////////////////////////////
// Register allocations for Wilson Kernel are precision indept
//////////////////////////////////////////////////////////////////////////////////////////
#define result_00 %zmm0 
#define result_01 %zmm1
#define result_02 %zmm2
  
#define result_10 %zmm3
#define result_11 %zmm4
#define result_12 %zmm5

#define result_20 %zmm6
#define result_21 %zmm7
#define result_22 %zmm8

#define result_30 %zmm9
#define result_31 %zmm10
#define result_32 %zmm11

#define Chi_00 %zmm12  
#define Chi_01 %zmm13
#define Chi_02 %zmm14

#define Chi_10 %zmm15
#define Chi_11 %zmm16
#define Chi_12 %zmm17  

#define UChi_00 %zmm18 
#define UChi_01 %zmm19
#define UChi_02 %zmm20

#define UChi_10 %zmm21
#define UChi_11 %zmm22
#define UChi_12 %zmm23 

#define Uir %zmm24 
#define Uri %zmm25  
#define T1 %zmm24
#define T2 %zmm25

#define Z0 %zmm26
#define Z1 %zmm27
#define Z2 %zmm28
#define Z3 %zmm29
#define Z4 %zmm30
#define Z5 %zmm31

#define TMP Chi_00

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

#include <simd/Intel512common.h>
#include <simd/Intel512avx.h>

//////////////////////////////////////////////////////////////////
// Macros used to build wilson kernel -- can rationalise and simplify
// a little as some duplication developed during trying different
// variants during optimisation. Could cut back to only those used.
//////////////////////////////////////////////////////////////////

//  const SiteSpinor * ptr = & in._odata[offset];	
#define LOAD_CHIMU(PTR)	 LOAD_CHIMUi(PTR)
#define LOAD_CHI(PTR)	 LOAD64(%r8,PTR) __asm__ ( LOAD_CHIi );
#define SAVE_UCHI(PTR)	 SAVE_UCHIi(PTR)
#define SAVE_CHI(PTR)	 SAVE_CHIi(PTR)
#define SAVE_RESULT(PTR) SAVE_RESULTi(PTR)

#define LOAD_CHIMUi \
	   LOAD_CHIMU01i	\
	   LOAD_CHIMU23i	);


#define LOAD_CHIMU01i\
	   VLOAD(0,%r8,Chimu_00)		\
	   VLOAD(1,%r8,Chimu_01)		\
	   VLOAD(2,%r8,Chimu_02)		\
	   VLOAD(3,%r8,Chimu_10)		\
	   VLOAD(4,%r8,Chimu_11)		\
	   VLOAD(5,%r8,Chimu_12)		

#define LOAD_CHIMU23i\
	   VLOAD(6,%r8,Chimu_20)		\
	   VLOAD(7,%r8,Chimu_21)		\
	   VLOAD(8,%r8,Chimu_22)		\
	   VLOAD(9,%r8,Chimu_30)		\
	   VLOAD(10,%r8,Chimu_31)		\
	   VLOAD(11,%r8,Chimu_32)		

#define SHUF_CHIMU23i\
	   VSHUFMEM(6,%r8,Chimu_20)		\
	   VSHUFMEM(7,%r8,Chimu_21)		\
	   VSHUFMEM(8,%r8,Chimu_22)		\
	   VSHUFMEM(9,%r8,Chimu_30)		\
	   VSHUFMEM(10,%r8,Chimu_31)		\
	   VSHUFMEM(11,%r8,Chimu_32)		


//  const SiteHalfSpinor *ptr = &buf[offset];	

#define LOAD_CHIi				\
  VLOAD(0,%r8,Chi_00)					\
  VLOAD(1,%r8,Chi_01)					\
  VLOAD(2,%r8,Chi_02)					\
  VLOAD(3,%r8,Chi_10)					\
  VLOAD(4,%r8,Chi_11)					\
  VLOAD(5,%r8,Chi_12)	
	

#define SAVE_UCHIi(PTR)				\
  LOAD64(%r8,PTR)				\
  __asm__ (					\
  VSTORE(0,%r8,UChi_00)				\
  VSTORE(1,%r8,UChi_01)				\
  VSTORE(2,%r8,UChi_02)				\
  VSTORE(3,%r8,UChi_10)				\
  VSTORE(4,%r8,UChi_11)				\
  VSTORE(5,%r8,UChi_12)				\
						);

#define SAVE_CHIi(PTR)				\
  LOAD64(%r8,PTR)				\
  __asm__ (					\
  VSTORE(0,%r8,Chi_00)				\
  VSTORE(1,%r8,Chi_01)				\
  VSTORE(2,%r8,Chi_02)				\
  VSTORE(3,%r8,Chi_10)				\
  VSTORE(4,%r8,Chi_11)				\
  VSTORE(5,%r8,Chi_12)				\
						);

#define SAVE_RESULTi(PTR)\
	   LOAD64(%r8,PTR)			\
  __asm__ (					\
	   VSTORE(0,%r8,result_00)		\
	   VSTORE(1,%r8,result_01)		\
	   VSTORE(2,%r8,result_02)		\
	   VSTORE(3,%r8,result_10)		\
	   VSTORE(4,%r8,result_11)		\
	   VSTORE(5,%r8,result_12)		\
	   VSTORE(6,%r8,result_20)		\
	   VSTORE(7,%r8,result_21)		\
	   VSTORE(8,%r8,result_22)		\
	   VSTORE(9,%r8,result_30)		\
	   VSTORE(10,%r8,result_31)		\
	   VSTORE(11,%r8,result_32) 		\
						);

#define MULT_2SPIN_DIR_PFXP(A,p) MULT_2SPIN_PFXP(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFYP(A,p) MULT_2SPIN_PFYP(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFZP(A,p) MULT_2SPIN_PFZP(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFTP(A,p) MULT_2SPIN_PFTP(&U._odata[sU](A),p)

#define MULT_2SPIN_DIR_PFXM(A,p) MULT_2SPIN_PFXM(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFYM(A,p) MULT_2SPIN_PFYM(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFZM(A,p) MULT_2SPIN_PFZM(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFTM(A,p) MULT_2SPIN_PFTM(&U._odata[sU](A),p)

#define MULT_2SPIN_PFXM(ptr,pf) MULT_2SPIN(ptr,pf)
#define MULT_2SPIN_PFYM(ptr,pf) MULT_2SPIN(ptr,pf)
#define MULT_2SPIN_PFZM(ptr,pf) MULT_2SPIN(ptr,pf)
#define MULT_2SPIN_PFTM(ptr,pf) MULT_2SPIN(ptr,pf)
#define MULT_2SPIN_PFTP(ptr,pf) MULT_2SPIN(ptr,pf)
#define MULT_2SPIN_PFZP(ptr,pf) MULT_2SPIN(ptr,pf)
#define MULT_2SPIN_PFYP(ptr,pf) MULT_2SPIN(ptr,pf)
#define MULT_2SPIN_PFXP(ptr,pf) MULT_2SPIN(ptr,pf)

//////////////////////////////////////////////////////////////////
// Dirac algebra
//////////////////////////////////////////////////////////////////

//      hspin(0)=fspin(0)+timesI(fspin(3));
//      hspin(1)=fspin(1)+timesI(fspin(2));
#define XP_PROJMEM(PTR) \
  LOAD64(%r8,PTR)							\
  __asm__ (								\
	   LOAD_CHIi						\
	   SHUF_CHIMU23i						\
	   VACCTIMESI1(Chi_00,Chi_00,Chimu_30)		\
	   VACCTIMESI1(Chi_01,Chi_01,Chimu_31)		\
	   VACCTIMESI1(Chi_02,Chi_02,Chimu_32)		\
	   VACCTIMESI1(Chi_10,Chi_10,Chimu_20)		\
	   VACCTIMESI1(Chi_11,Chi_11,Chimu_21)		\
	   VACCTIMESI1(Chi_12,Chi_12,Chimu_22)		\
	   VACCTIMESI2(Chi_00,Chi_00,Chimu_30)		\
	   VACCTIMESI2(Chi_01,Chi_01,Chimu_31)		\
	   VACCTIMESI2(Chi_02,Chi_02,Chimu_32)		\
	   VACCTIMESI2(Chi_10,Chi_10,Chimu_20)		\
	   VACCTIMESI2(Chi_11,Chi_11,Chimu_21)		\
	   VACCTIMESI2(Chi_12,Chi_12,Chimu_22)		);


#define YP_PROJMEM(ptr) \
  LOAD64(%r8,ptr)		\
  __asm__ (					\
  LOAD_CHIMU01i					\
  VSUBMEM(9,%r8 ,Chimu_00,Chi_00)		\
  VSUBMEM(10,%r8,Chimu_01,Chi_01)		\
  VSUBMEM(11,%r8,Chimu_02,Chi_02)		\
  VADDMEM(6,%r8,Chimu_10,Chi_10)		\
  VADDMEM(7,%r8,Chimu_11,Chi_11)		\
  VADDMEM(8,%r8,Chimu_12,Chi_12)		);

#define ZP_PROJMEM(PTR) \
  LOAD64(%r8,PTR)							\
  __asm__ (								\
	   LOAD_CHIi						\
	   SHUF_CHIMU23i						\
	   VACCTIMESI1(Chi_00,Chi_00,Chimu_20)				\
	   VACCTIMESI1(Chi_01,Chi_01,Chimu_21)		   	        \
	   VACCTIMESI1(Chi_02,Chi_02,Chimu_22)				\
	   VACCTIMESMINUSI1(Chi_10,Chi_10,Chimu_30)			\
	   VACCTIMESMINUSI1(Chi_11,Chi_11,Chimu_31)			\
	   VACCTIMESMINUSI1(Chi_12,Chi_12,Chimu_32)			\
	   VACCTIMESI2(Chi_00,Chi_00,Chimu_20)				\
	   VACCTIMESI2(Chi_01,Chi_01,Chimu_21)				\
	   VACCTIMESI2(Chi_02,Chi_02,Chimu_22)				\
	   VACCTIMESMINUSI2(Chi_10,Chi_10,Chimu_30)		\
	   VACCTIMESMINUSI2(Chi_11,Chi_11,Chimu_31)		\
	   VACCTIMESMINUSI2(Chi_12,Chi_12,Chimu_32)	);


#define TP_PROJMEM(ptr)				\
  LOAD64(%r8,ptr)				\
  __asm__ (					\
	   LOAD_CHIMU01i			\
	   VADDMEM(6,%r8 ,Chimu_00,Chi_00)	\
	   VADDMEM(7,%r8,Chimu_01,Chi_01)	\
	   VADDMEM(8,%r8,Chimu_02,Chi_02)	\
	   VADDMEM(9,%r8,Chimu_10,Chi_10)	\
	   VADDMEM(10,%r8,Chimu_11,Chi_11)	\
	   VADDMEM(11,%r8,Chimu_12,Chi_12)	);

//      hspin(0)=fspin(0)-timesI(fspin(3))
//      hspin(1)=fspin(1)-timesI(fspin(2))

#define XM_PROJMEM(PTR) \
  LOAD64(%r8,PTR)\
  __asm__ (								\
	   SHUF_CHIMU23i						\
	   LOAD_CHIi \
	   VACCTIMESMINUSI1(Chi_00,Chi_00,Chimu_30)\
	   VACCTIMESMINUSI1(Chi_01,Chi_01,Chimu_31)\
	   VACCTIMESMINUSI1(Chi_02,Chi_02,Chimu_32)\
	   VACCTIMESMINUSI1(Chi_10,Chi_10,Chimu_20)\
	   VACCTIMESMINUSI1(Chi_11,Chi_11,Chimu_21)\
	   VACCTIMESMINUSI1(Chi_12,Chi_12,Chimu_22)\
	   VACCTIMESMINUSI2(Chi_00,Chi_00,Chimu_30)\
	   VACCTIMESMINUSI2(Chi_01,Chi_01,Chimu_31)\
	   VACCTIMESMINUSI2(Chi_02,Chi_02,Chimu_32)\
	   VACCTIMESMINUSI2(Chi_10,Chi_10,Chimu_20)\
	   VACCTIMESMINUSI2(Chi_11,Chi_11,Chimu_21)\
	   VACCTIMESMINUSI2(Chi_12,Chi_12,Chimu_22) );

#define YM_PROJMEM(ptr)				\
  LOAD64(%r8,ptr)				\
  __asm__ (					\
  LOAD_CHIMU01i					\
  VADDMEM(9,%r8 ,Chimu_00,Chi_00)		\
  VADDMEM(10,%r8,Chimu_01,Chi_01)		\
  VADDMEM(11,%r8,Chimu_02,Chi_02)		\
  VSUBMEM(6,%r8,Chimu_10,Chi_10)		\
  VSUBMEM(7,%r8,Chimu_11,Chi_11)		\
  VSUBMEM(8,%r8,Chimu_12,Chi_12)			);

#define ZM_PROJMEM(PTR) \
  LOAD64(%r8,PTR)							\
  __asm__ (								\
	   SHUF_CHIMU23i						\
           LOAD_CHIi \
	   VACCTIMESMINUSI1(Chi_00,Chi_00,Chimu_20)\
	   VACCTIMESMINUSI1(Chi_01,Chi_01,Chimu_21)\
	   VACCTIMESMINUSI1(Chi_02,Chi_02,Chimu_22)\
	   VACCTIMESI1(Chi_10,Chi_10,Chimu_30)\
	   VACCTIMESI1(Chi_11,Chi_11,Chimu_31)\
	   VACCTIMESI1(Chi_12,Chi_12,Chimu_32)\
	   VACCTIMESMINUSI2(Chi_00,Chi_00,Chimu_20)\
	   VACCTIMESMINUSI2(Chi_01,Chi_01,Chimu_21)\
	   VACCTIMESMINUSI2(Chi_02,Chi_02,Chimu_22)\
	   VACCTIMESI2(Chi_10,Chi_10,Chimu_30)\
	   VACCTIMESI2(Chi_11,Chi_11,Chimu_31)\
	   VACCTIMESI2(Chi_12,Chi_12,Chimu_32) );

#define TM_PROJMEM(ptr)				\
  LOAD64(%r8,ptr)				\
  __asm__ (					\
  LOAD_CHIMU01i					\
  VSUBMEM(6,%r8,Chimu_00,Chi_00)		\
  VSUBMEM(7,%r8,Chimu_01,Chi_01)		\
  VSUBMEM(8,%r8,Chimu_02,Chi_02)		\
  VSUBMEM(9,%r8,Chimu_10,Chi_10)		\
  VSUBMEM(10,%r8,Chimu_11,Chi_11)		\
  VSUBMEM(11,%r8,Chimu_12,Chi_12)		);

//      fspin(0)=hspin(0)
//      fspin(1)=hspin(1)
//      fspin(2)=timesMinusI(hspin(1))
//      fspin(3)=timesMinusI(hspin(0))
#define XP_RECON __asm__ (			\
			  VZERO(TMP)		\
			  VTIMESMINUSI0(UChi_00,result_30,TMP)	\
			  VTIMESMINUSI0(UChi_10,result_20,TMP)	\
			  VTIMESMINUSI0(UChi_01,result_31,TMP)	\
			  VTIMESMINUSI0(UChi_11,result_21,TMP)	\
			  VTIMESMINUSI0(UChi_02,result_32,TMP)   \
			  VTIMESMINUSI0(UChi_12,result_22,TMP)	\
			  VMOV(UChi_00,result_00)	\
			  VMOV(UChi_10,result_10)	\
			  VMOV(UChi_01,result_01)	\
			  VMOV(UChi_11,result_11)	\
			  VMOV(UChi_02,result_02)	\
			  VMOV(UChi_12,result_12)	\
			  VTIMESMINUSI1(UChi_10,result_20,TMP)	\
			  VTIMESMINUSI1(UChi_11,result_21,TMP)	\
			  VTIMESMINUSI1(UChi_12,result_22,TMP)	\
			  VTIMESMINUSI1(UChi_00,result_30,TMP)	\
			  VTIMESMINUSI1(UChi_01,result_31,TMP)	\
			  VTIMESMINUSI1(UChi_02,result_32,TMP)   \
			  VTIMESMINUSI2(UChi_10,result_20,TMP)	\
			  VTIMESMINUSI2(UChi_11,result_21,TMP)	\
			  VTIMESMINUSI2(UChi_12,result_22,TMP)	\
			  VTIMESMINUSI2(UChi_00,result_30,TMP)	\
			  VTIMESMINUSI2(UChi_01,result_31,TMP)	\
			  VTIMESMINUSI2(UChi_02,result_32,TMP)   \
						);
  // NB could save 6 ops using addsub => 12 cycles
#define XP_RECON_ACCUM __asm__ ( \
  VZERO(TMP)\
  VACCTIMESMINUSI0(UChi_00,result_30,Z3)\
  VACCTIMESMINUSI0(UChi_10,result_20,Z0)\
  VACCTIMESMINUSI0(UChi_01,result_31,Z4)\
  VACCTIMESMINUSI0(UChi_11,result_21,Z1)\
  VACCTIMESMINUSI0(UChi_02,result_32,Z5)\
  VACCTIMESMINUSI0(UChi_12,result_22,Z2)\
  VADD(UChi_00,result_00,result_00)\
  VADD(UChi_10,result_10,result_10)\
  VADD(UChi_01,result_01,result_01)\
  VADD(UChi_11,result_11,result_11)\
  VADD(UChi_02,result_02,result_02)\
  VADD(UChi_12,result_12,result_12)\
  VACCTIMESMINUSI1(UChi_00,result_30,Z3)\
  VACCTIMESMINUSI1(UChi_10,result_20,Z0)\
  VACCTIMESMINUSI1(UChi_01,result_31,Z4)\
  VACCTIMESMINUSI1(UChi_11,result_21,Z1)\
  VACCTIMESMINUSI1(UChi_02,result_32,Z5)\
  VACCTIMESMINUSI1(UChi_12,result_22,Z2)\
  VACCTIMESMINUSI2(UChi_10,result_20,Z0)\
  VACCTIMESMINUSI2(UChi_11,result_21,Z1)\
  VACCTIMESMINUSI2(UChi_12,result_22,Z2)\
  VACCTIMESMINUSI2(UChi_00,result_30,Z3)\
  VACCTIMESMINUSI2(UChi_01,result_31,Z4)\
  VACCTIMESMINUSI2(UChi_02,result_32,Z5)\
				 );

#define XM_RECON __asm__ ( \
  VZERO(TMP)\
  VTIMESI0(UChi_00,result_30,TMP)\
  VTIMESI0(UChi_10,result_20,TMP)\
  VTIMESI0(UChi_01,result_31,TMP)\
  VTIMESI0(UChi_11,result_21,TMP)\
  VTIMESI0(UChi_02,result_32,TMP)\
  VTIMESI0(UChi_12,result_22,TMP)\
  VMOV(UChi_00,result_00)\
  VMOV(UChi_10,result_10)\
  VMOV(UChi_01,result_01)\
  VMOV(UChi_11,result_11)\
  VMOV(UChi_02,result_02)\
  VMOV(UChi_12,result_12)\
  VTIMESI1(UChi_00,result_30,TMP)\
  VTIMESI1(UChi_10,result_20,TMP)\
  VTIMESI1(UChi_01,result_31,TMP)\
  VTIMESI1(UChi_11,result_21,TMP)\
  VTIMESI1(UChi_02,result_32,TMP)\
  VTIMESI1(UChi_12,result_22,TMP)\
  VTIMESI2(UChi_10,result_20,TMP)\
  VTIMESI2(UChi_11,result_21,TMP)\
  VTIMESI2(UChi_12,result_22,TMP)\
  VTIMESI2(UChi_00,result_30,TMP)\
  VTIMESI2(UChi_01,result_31,TMP)\
  VTIMESI2(UChi_02,result_32,TMP)\
			   );

#define XM_RECON_ACCUM __asm__ ( \
  VACCTIMESI0(UChi_10,result_20,Z0)\
  VACCTIMESI0(UChi_00,result_30,Z3)\
  VACCTIMESI0(UChi_11,result_21,Z1)\
  VACCTIMESI0(UChi_01,result_31,Z4)\
  VACCTIMESI0(UChi_12,result_22,Z2)\
  VACCTIMESI0(UChi_02,result_32,Z5)\
  \
  VADD(UChi_10,result_10,result_10)\
  VADD(UChi_00,result_00,result_00)\
  VADD(UChi_11,result_11,result_11)\
  VADD(UChi_01,result_01,result_01)\
  VADD(UChi_12,result_12,result_12)\
  VADD(UChi_02,result_02,result_02)\
  \
  VACCTIMESI1(UChi_10,result_20,Z0)\
  VACCTIMESI1(UChi_00,result_30,Z3)\
  VACCTIMESI1(UChi_11,result_21,Z1)\
  VACCTIMESI1(UChi_01,result_31,Z4)\
  VACCTIMESI1(UChi_12,result_22,Z2)\
  VACCTIMESI1(UChi_02,result_32,Z5)\
  VACCTIMESI2(UChi_10,result_20,Z0)\
  VACCTIMESI2(UChi_11,result_21,Z1)\
  VACCTIMESI2(UChi_12,result_22,Z2)\
  VACCTIMESI2(UChi_00,result_30,Z3)\
  VACCTIMESI2(UChi_01,result_31,Z4)\
  VACCTIMESI2(UChi_02,result_32,Z5)\
				 );

#define YP_RECON_ACCUM __asm__ ( \
  VADD(UChi_00,result_00,result_00)\
  VADD(UChi_10,result_10,result_10)\
  VADD(UChi_01,result_01,result_01)\
  VADD(UChi_11,result_11,result_11)\
  VADD(UChi_02,result_02,result_02)\
  VADD(UChi_12,result_12,result_12)\
  VADD(UChi_10,result_20,result_20)\
  VADD(UChi_11,result_21,result_21)\
  VADD(UChi_12,result_22,result_22)\
  VSUB(UChi_00,result_30,result_30)\
  VSUB(UChi_01,result_31,result_31)\
  VSUB(UChi_02,result_32,result_32) );

#define YM_RECON_ACCUM __asm__ ( \
  VADD(UChi_00,result_00,result_00)\
  VADD(UChi_10,result_10,result_10)\
  VADD(UChi_01,result_01,result_01)\
  VADD(UChi_11,result_11,result_11)\
  VADD(UChi_02,result_02,result_02)\
  VADD(UChi_12,result_12,result_12)\
  VSUB(UChi_10,result_20,result_20)\
  VSUB(UChi_11,result_21,result_21)\
  VSUB(UChi_12,result_22,result_22)\
  VADD(UChi_00,result_30,result_30)\
  VADD(UChi_01,result_31,result_31)\
  VADD(UChi_02,result_32,result_32) );

#define ZP_RECON_ACCUM __asm__ ( \
  VACCTIMESMINUSI0(UChi_00,result_20,Z0)\
  VACCTIMESI0(UChi_10,result_30,Z3)\
  VACCTIMESMINUSI0(UChi_01,result_21,Z1)\
  VACCTIMESI0(UChi_11,result_31,Z4)\
  VACCTIMESMINUSI0(UChi_02,result_22,Z2)\
  VACCTIMESI0(UChi_12,result_32,Z5)\
  VADD(UChi_00,result_00,result_00)\
  VADD(UChi_10,result_10,result_10)\
  VADD(UChi_01,result_01,result_01)\
  VADD(UChi_11,result_11,result_11)\
  VADD(UChi_02,result_02,result_02)\
  VADD(UChi_12,result_12,result_12)\
  VACCTIMESMINUSI1(UChi_00,result_20,Z0)\
  VACCTIMESI1(UChi_10,result_30,Z3)\
  VACCTIMESMINUSI1(UChi_01,result_21,Z1)\
  VACCTIMESI1(UChi_11,result_31,Z4)\
  VACCTIMESMINUSI1(UChi_02,result_22,Z2)\
  VACCTIMESI1(UChi_12,result_32,Z5)\
  VACCTIMESMINUSI2(UChi_00,result_20,Z0)\
  VACCTIMESMINUSI2(UChi_01,result_21,Z1)\
  VACCTIMESMINUSI2(UChi_02,result_22,Z2)\
  VACCTIMESI2(UChi_10,result_30,Z3)\
  VACCTIMESI2(UChi_11,result_31,Z4)\
  VACCTIMESI2(UChi_12,result_32,Z5)\
				 );

#define ZM_RECON_ACCUM __asm__ ( \
  VACCTIMESI0(UChi_00,result_20,Z0)\
  VACCTIMESMINUSI0(UChi_10,result_30,Z3)\
  VACCTIMESI0(UChi_01,result_21,Z1)\
  VACCTIMESMINUSI0(UChi_11,result_31,Z4)\
  VACCTIMESI0(UChi_02,result_22,Z2)\
  VACCTIMESMINUSI0(UChi_12,result_32,Z5)\
  VADD(UChi_00,result_00,result_00)\
  VADD(UChi_10,result_10,result_10)\
  VADD(UChi_01,result_01,result_01)\
  VADD(UChi_11,result_11,result_11)\
  VADD(UChi_02,result_02,result_02)\
  VADD(UChi_12,result_12,result_12)\
  VACCTIMESI1(UChi_00,result_20,Z0)\
  VACCTIMESMINUSI1(UChi_10,result_30,Z3)\
  VACCTIMESI1(UChi_01,result_21,Z1)\
  VACCTIMESMINUSI1(UChi_11,result_31,Z4)\
  VACCTIMESI1(UChi_02,result_22,Z2)\
  VACCTIMESMINUSI1(UChi_12,result_32,Z5)\
  VACCTIMESI2(UChi_00,result_20,Z0)\
  VACCTIMESI2(UChi_01,result_21,Z1)\
  VACCTIMESI2(UChi_02,result_22,Z2)\
  VACCTIMESMINUSI2(UChi_10,result_30,Z3)\
  VACCTIMESMINUSI2(UChi_11,result_31,Z4)\
  VACCTIMESMINUSI2(UChi_12,result_32,Z5)\
				 );

#define TP_RECON_ACCUM __asm__ ( \
  VADD(UChi_00,result_00,result_00)\
  VADD(UChi_10,result_10,result_10)\
  VADD(UChi_01,result_01,result_01)\
  VADD(UChi_11,result_11,result_11)\
  VADD(UChi_02,result_02,result_02)\
  VADD(UChi_12,result_12,result_12)\
  VADD(UChi_00,result_20,result_20)\
  VADD(UChi_10,result_30,result_30)\
  VADD(UChi_01,result_21,result_21)\
  VADD(UChi_11,result_31,result_31)\
  VADD(UChi_02,result_22,result_22)\
  VADD(UChi_12,result_32,result_32) );

#define TM_RECON_ACCUM __asm__ ( \
  VADD(UChi_00,result_00,result_00)\
  VADD(UChi_10,result_10,result_10)\
  VADD(UChi_01,result_01,result_01)\
  VADD(UChi_11,result_11,result_11)\
  VADD(UChi_02,result_02,result_02)\
  VADD(UChi_12,result_12,result_12)\
  VSUB(UChi_00,result_20,result_20)\
  VSUB(UChi_10,result_30,result_30)\
  VSUB(UChi_01,result_21,result_21)\
  VSUB(UChi_11,result_31,result_31)\
  VSUB(UChi_02,result_22,result_22)\
  VSUB(UChi_12,result_32,result_32) );

#define PREFETCH_CHIMU(A) 
/*
  LOAD64(%r9,A)						\
	   __asm__ (						\
  VPREFETCHG(0,%r9)\
  VPREFETCHG(1,%r9)\
  VPREFETCHG(2,%r9)\
  VPREFETCHG(3,%r9)\
  VPREFETCHG(4,%r9)\
  VPREFETCHG(5,%r9)\
  VPREFETCHG(6,%r9)\
  VPREFETCHG(7,%r9)\
  VPREFETCHG(8,%r9)\
  VPREFETCHG(9,%r9)\
  VPREFETCHG(10,%r9)\
  VPREFETCHG(11,%r9));
*/
#define PERMUTE_DIR0 __asm__ ( 	\
  VPERM0(Chi_00,Chi_00)	\
  VPERM0(Chi_01,Chi_01)	\
  VPERM0(Chi_02,Chi_02)	\
  VPERM0(Chi_10,Chi_10)	\
  VPERM0(Chi_11,Chi_11)	\
  VPERM0(Chi_12,Chi_12) );

#define PERMUTE_DIR1 __asm__ (	\
  VPERM1(Chi_00,Chi_00)	\
  VPERM1(Chi_01,Chi_01)	\
  VPERM1(Chi_02,Chi_02)	\
  VPERM1(Chi_10,Chi_10)	\
  VPERM1(Chi_11,Chi_11)	\
  VPERM1(Chi_12,Chi_12));

#define PERMUTE_DIR2 __asm__ (	\
  VPERM2(Chi_00,Chi_00)	\
  VPERM2(Chi_01,Chi_01)	\
  VPERM2(Chi_02,Chi_02)	\
  VPERM2(Chi_10,Chi_10)	\
  VPERM2(Chi_11,Chi_11)	\
  VPERM2(Chi_12,Chi_12) );

#define PERMUTE_DIR3 __asm__ (	\
  VPERM3(Chi_00,Chi_00)	\
  VPERM3(Chi_01,Chi_01)	\
  VPERM3(Chi_02,Chi_02)	\
  VPERM3(Chi_10,Chi_10)	\
  VPERM3(Chi_11,Chi_11)	\
  VPERM3(Chi_12,Chi_12) );


#define MULT_ADDSUB_2SPIN(ptr,pf)					\
  LOAD64(%r8,ptr)						\
  LOAD64(%r9,pf)						\
	   __asm__ (						\
	   VPREFETCH2(9,%r8)				   VPREFETCH2(10,%r8)					   \
	   VPREFETCH2(11,%r8)					   \
	   VPREFETCH2(12,%r8)					   \
	   VPREFETCH2(13,%r8)					   \
	   VPREFETCH2(14,%r8)					   \
	   VPREFETCH2(15,%r8)					   \
	   VPREFETCH2(16,%r8)					   \
	   VPREFETCH2(17,%r8)					   \
	   VSHUF(Chi_00,T1)				\
	   VMOVIDUP(0,%r8,Z0 )					\
           VMOVIDUP(3,%r8,Z1 )					\
           VMOVIDUP(6,%r8,Z2 )	          VSHUF(Chi_10,T2)		\
	   /*6*/							\
           VMUL(Z0,T1,UChi_00)            VMOVRDUP(0,%r8,Z3 )	\
           VMUL(Z0,T2,UChi_10)            VMOVRDUP(3,%r8,Z4 )	\
           VMUL(Z1,T1,UChi_01)            VMOVRDUP(6,%r8,Z5 )	\
           VMUL(Z1,T2,UChi_11)            VMOVIDUP(1,%r8,Z0 )	\
           VMUL(Z2,T1,UChi_02)            VMOVIDUP(4,%r8,Z1 )	\
           VMUL(Z2,T2,UChi_12)            VMOVIDUP(7,%r8,Z2 )	\
	   VPREFETCHG(0,%r9)					   \
	   VPREFETCHG(1,%r9)					   \
	   VPREFETCHG(2,%r9)					   \
	   VPREFETCHG(3,%r9)					   \
	   /*18*/						\
           VMADDSUB(Z3,Chi_00,UChi_00)    VSHUF(Chi_01,T1)	\
           VMADDSUB(Z3,Chi_10,UChi_10)				\
           VMADDSUB(Z4,Chi_00,UChi_01)    VMOVRDUP(1,%r8,Z3 )	\
           VMADDSUB(Z4,Chi_10,UChi_11)    VSHUF(Chi_11,T2)	\
           VMADDSUB(Z5,Chi_00,UChi_02)    VMOVRDUP(4,%r8,Z4 )	\
           VMADDSUB(Z5,Chi_10,UChi_12)				\
	   VPREFETCHG(4,%r9)					   \
	   VPREFETCHG(5,%r9)					   \
	   VPREFETCHG(6,%r9)					   \
	   VPREFETCHG(7,%r9)					   \
	   /*28*/						\
           VMADDSUB(Z0,T1,UChi_00)        VMOVRDUP(7,%r8,Z5 )	\
           VMADDSUB(Z0,T2,UChi_10)				\
           VMADDSUB(Z1,T1,UChi_01)        VMOVIDUP(2,%r8,Z0 )	\
           VMADDSUB(Z1,T2,UChi_11)				\
           VMADDSUB(Z2,T1,UChi_02)        VMOVIDUP(5,%r8,Z1 )	\
           VMADDSUB(Z2,T2,UChi_12)        VMOVIDUP(8,%r8,Z2 )	\
	   VPREFETCH2(12,%r9)					   \
	   VPREFETCH2(13,%r9)					   \
	   VPREFETCH2(14,%r9)					   \
	   VPREFETCH2(15,%r9)					   \
	   VPREFETCH2(16,%r9)					   \
	   VPREFETCH2(17,%r9)					   \
	   VPREFETCH2(18,%r9)					   \
	   VPREFETCH2(19,%r9)					   \
	   VPREFETCH2(20,%r9)					   \
	   VPREFETCH2(21,%r9)					   \
	   VPREFETCH2(22,%r9)					   \
	   VPREFETCH2(23,%r9)					   \
           /*38*/						\
           VMADDSUB(Z3,Chi_01,UChi_00)    VSHUF(Chi_02,T1)	\
           VMADDSUB(Z3,Chi_11,UChi_10)				\
           VMADDSUB(Z4,Chi_01,UChi_01)    VMOVRDUP(2,%r8,Z3 )	\
           VMADDSUB(Z4,Chi_11,UChi_11)    VSHUF(Chi_12,T2)	\
           VMADDSUB(Z5,Chi_01,UChi_02)    VMOVRDUP(5,%r8,Z4 )	\
           VMADDSUB(Z5,Chi_11,UChi_12)				\
	   VPREFETCHG(9,%r8)				   \
	   VPREFETCHG(10,%r8)					   \
	   VPREFETCHG(11,%r8)					   \
	   VPREFETCHG(12,%r8)					   \
	   VPREFETCHG(13,%r8)					   \
	   VPREFETCHG(14,%r8)					   \
	   VPREFETCHG(15,%r8)					   \
	   VPREFETCHG(16,%r8)					   \
	   VPREFETCHG(17,%r8)					   \
	   /*48*/						\
           VMADDSUB(Z0,T1,UChi_00)        VMOVRDUP(8,%r8,Z5 ) \
           VMADDSUB(Z0,T2,UChi_10)			      \
           VMADDSUB(Z1,T1,UChi_01)			      \
           VMADDSUB(Z1,T2,UChi_11)			      \
           VMADDSUB(Z2,T1,UChi_02)			      \
           VMADDSUB(Z2,T2,UChi_12)			      \
	   VPREFETCHG(8,%r9)					   \
	   VPREFETCHG(9,%r9)					   \
	   VPREFETCHG(10,%r9)					   \
	   VPREFETCHG(11,%r9)					   \
	   /*55*/					      \
           VMADDSUB(Z3,Chi_02,UChi_00)			      \
           VMADDSUB(Z3,Chi_12,UChi_10)			      \
           VMADDSUB(Z4,Chi_02,UChi_01)			      \
           VMADDSUB(Z4,Chi_12,UChi_11)			      \
           VMADDSUB(Z5,Chi_02,UChi_02)			      \
           VMADDSUB(Z5,Chi_12,UChi_12)			      \
	   /*61 insns*/							);


#define MULT_ADDSUB_2SPIN_LS(ptr,pf)				   \
  LOAD64(%r8,ptr)						   \
  LOAD64(%r9,pf)						   \
  __asm__ (							   \
           VSHUF(Chi_00,T1)      VSHUF(Chi_10,T2)		   \
           VMULIDUP(0,%r8,T1,UChi_00) VMULIDUP(0,%r8,T2,UChi_10)   \
           VMULIDUP(3,%r8,T1,UChi_01) VMULIDUP(3,%r8,T2,UChi_11)   \
           VMULIDUP(6,%r8,T1,UChi_02) VMULIDUP(6,%r8,T2,UChi_12)   \
	   VPREFETCHG(0,%r9)					   \
	   VPREFETCHG(1,%r9)					   \
	   VPREFETCHG(2,%r9)					   \
	   VPREFETCHG(3,%r9)					   \
	   /*8*/						   \
           VSHUF(Chi_01,T1)	  VSHUF(Chi_11,T2)	       	   \
           VMADDSUBRDUP(0,%r8,Chi_00,UChi_00) VMADDSUBRDUP(0,%r8,Chi_10,UChi_10) \
           VMADDSUBRDUP(3,%r8,Chi_00,UChi_01) VMADDSUBRDUP(3,%r8,Chi_10,UChi_11) \
           VMADDSUBRDUP(6,%r8,Chi_00,UChi_02) VMADDSUBRDUP(6,%r8,Chi_10,UChi_12) \
	   VPREFETCHG(4,%r9)					   \
	   VPREFETCHG(5,%r9)					   \
	   VPREFETCHG(6,%r9)					   \
	   VPREFETCHG(7,%r9)					   \
	   /*16*/					  	   \
           VMADDSUBIDUP(1,%r8,T1,UChi_00)     VMADDSUBIDUP(1,%r8,T2,UChi_10)	   \
           VMADDSUBIDUP(4,%r8,T1,UChi_01)     VMADDSUBIDUP(4,%r8,T2,UChi_11) \
           VMADDSUBIDUP(7,%r8,T1,UChi_02)     VMADDSUBIDUP(7,%r8,T2,UChi_12) \
	   VPREFETCHG(8,%r9)					   \
	   VPREFETCHG(9,%r9)					   \
	   VPREFETCHG(10,%r9)					   \
	   VPREFETCHG(11,%r9)					   \
           /*22*/						   \
           VSHUF(Chi_02,T1)    VSHUF(Chi_12,T2)	                   \
           VMADDSUBRDUP(1,%r8,Chi_01,UChi_00) VMADDSUBRDUP(1,%r8,Chi_11,UChi_10) \
           VMADDSUBRDUP(4,%r8,Chi_01,UChi_01) VMADDSUBRDUP(4,%r8,Chi_11,UChi_11) \
           VMADDSUBRDUP(7,%r8,Chi_01,UChi_02) VMADDSUBRDUP(7,%r8,Chi_11,UChi_12) \
	   VPREFETCH2(12,%r9)					   \
	   VPREFETCH2(13,%r9)					   \
	   VPREFETCH2(14,%r9)					   \
	   VPREFETCH2(15,%r9)					   \
	   /*30*/						   \
           VMADDSUBIDUP(2,%r8,T1,UChi_00)     VMADDSUBIDUP(2,%r8,T2,UChi_10)	   \
           VMADDSUBIDUP(5,%r8,T1,UChi_01)     VMADDSUBIDUP(5,%r8,T2,UChi_11)     \
	   VPREFETCH2(16,%r9)					   \
	   VPREFETCH2(17,%r9)					   \
	   VPREFETCH2(18,%r9)					   \
	   VPREFETCH2(19,%r9)					   \
           VMADDSUBIDUP(8,%r8,T1,UChi_02)     VMADDSUBIDUP(8,%r8,T2,UChi_12)     \
	   /*36*/					           \
           VMADDSUBRDUP(2,%r8,Chi_02,UChi_00) VMADDSUBRDUP(2,%r8,Chi_12,UChi_10) \
           VMADDSUBRDUP(5,%r8,Chi_02,UChi_01) VMADDSUBRDUP(5,%r8,Chi_12,UChi_11) \
           VMADDSUBRDUP(8,%r8,Chi_02,UChi_02) VMADDSUBRDUP(8,%r8,Chi_12,UChi_12) \
	   VPREFETCH2(20,%r9)					   \
	   VPREFETCH2(21,%r9)					   \
	   VPREFETCH2(22,%r9)					   \
	   VPREFETCH2(23,%r9)					   \
	   VPREFETCHG(2,%r8)					   \
	   VPREFETCHG(3,%r8)					   \
	   VPREFETCH2(4,%r8)					   \
	   VPREFETCH2(5,%r8)					   \
	   /*42 insns*/						);

#define MULT_ADDSUB_2SPIN_LSNOPF(ptr,pf)				   \
  LOAD64(%r8,ptr)						   \
  LOAD64(%r9,pf)						   \
  __asm__ (							   \
           VSHUF(Chi_00,T1)      VSHUF(Chi_10,T2)		   \
           VMULIDUP(0,%r8,T1,UChi_00) VMULIDUP(0,%r8,T2,UChi_10)   \
           VMULIDUP(3,%r8,T1,UChi_01) VMULIDUP(3,%r8,T2,UChi_11)   \
           VMULIDUP(6,%r8,T1,UChi_02) VMULIDUP(6,%r8,T2,UChi_12)   \
	   /*8*/						   \
           VSHUF(Chi_01,T1)	  VSHUF(Chi_11,T2)	       	   \
           VMADDSUBRDUP(0,%r8,Chi_00,UChi_00) VMADDSUBRDUP(0,%r8,Chi_10,UChi_10) \
           VMADDSUBRDUP(3,%r8,Chi_00,UChi_01) VMADDSUBRDUP(3,%r8,Chi_10,UChi_11) \
           VMADDSUBRDUP(6,%r8,Chi_00,UChi_02) VMADDSUBRDUP(6,%r8,Chi_10,UChi_12) \
	   /*16*/					  	   \
           VMADDSUBIDUP(1,%r8,T1,UChi_00)     VMADDSUBIDUP(1,%r8,T2,UChi_10)	   \
           VMADDSUBIDUP(4,%r8,T1,UChi_01)     VMADDSUBIDUP(4,%r8,T2,UChi_11) \
           VMADDSUBIDUP(7,%r8,T1,UChi_02)     VMADDSUBIDUP(7,%r8,T2,UChi_12) \
           /*22*/						   \
           VSHUF(Chi_02,T1)    VSHUF(Chi_12,T2)	                   \
           VMADDSUBRDUP(1,%r8,Chi_01,UChi_00) VMADDSUBRDUP(1,%r8,Chi_11,UChi_10) \
           VMADDSUBRDUP(4,%r8,Chi_01,UChi_01) VMADDSUBRDUP(4,%r8,Chi_11,UChi_11) \
           VMADDSUBRDUP(7,%r8,Chi_01,UChi_02) VMADDSUBRDUP(7,%r8,Chi_11,UChi_12) \
	   /*30*/						   \
           VMADDSUBIDUP(2,%r8,T1,UChi_00)     VMADDSUBIDUP(2,%r8,T2,UChi_10)	   \
           VMADDSUBIDUP(5,%r8,T1,UChi_01)     VMADDSUBIDUP(5,%r8,T2,UChi_11)     \
           VMADDSUBIDUP(8,%r8,T1,UChi_02)     VMADDSUBIDUP(8,%r8,T2,UChi_12)     \
	   /*36*/					           \
           VMADDSUBRDUP(2,%r8,Chi_02,UChi_00) VMADDSUBRDUP(2,%r8,Chi_12,UChi_10) \
           VMADDSUBRDUP(5,%r8,Chi_02,UChi_01) VMADDSUBRDUP(5,%r8,Chi_12,UChi_11) \
           VMADDSUBRDUP(8,%r8,Chi_02,UChi_02) VMADDSUBRDUP(8,%r8,Chi_12,UChi_12) \
	   /*	   VPREFETCHG(2,%r8)*/				   \
	   /*	   VPREFETCHG(3,%r8)*/				   \
	   /*42 insns*/						);


#define Z6 Chi_00
#define MULT_ADDSUB_2SPIN_NEW(ptr,pf)			       \
  LOAD64(%r8,ptr)					       \
  __asm__ (							  \
   VSHUFMEM(0,%r8,Z0)					          \
   VRDUP(Chi_00,T1)           VIDUP(Chi_00,Chi_00)	          \
   VRDUP(Chi_10,T2)           VIDUP(Chi_10,Chi_10)		  \
   VMUL(Z0,Chi_00,Z1)         VMUL(Z0,Chi_10,Z2)		  \
   VSHUFMEM(3,%r8,Z0)						  \
   VMUL(Z0,Chi_00,Z3)         VMUL(Z0,Chi_10,Z4)		  \
   VSHUFMEM(6,%r8,Z0)						  \
   VMUL(Z0,Chi_00,Z5)         VMUL(Z0,Chi_10,Z6)		  \
   VMULMEM(0,%r8,T1,UChi_00)  VMULMEM(0,%r8,T2,UChi_10)		  \
   VMULMEM(3,%r8,T1,UChi_01)  VMULMEM(3,%r8,T2,UChi_11)		  \
   VMULMEM(6,%r8,T1,UChi_02)  VMULMEM(6,%r8,T2,UChi_12)		  \
   /*11 cycles*/						  \
   VSHUFMEM(1,%r8,Z0)						  \
   VRDUP(Chi_01,T1)           VIDUP(Chi_01,Chi_01)		  \
   VRDUP(Chi_11,T2)           VIDUP(Chi_11,Chi_11)		  \
   VMADD(Z0,Chi_01,Z1)        VMADD(Z0,Chi_11,Z2)		  \
   VSHUFMEM(4,%r8,Z0)						  \
   VMADD(Z0,Chi_01,Z3)        VMADD(Z0,Chi_11,Z4)		  \
   VSHUFMEM(7,%r8,Z0)						  \
   VMADD(Z0,Chi_01,Z5)        VMADD(Z0,Chi_11,Z6)		  \
   VMADDMEM(1,%r8,T1,UChi_00) VMADDMEM(1,%r8,T2,UChi_10)	  \
   VMADDMEM(4,%r8,T1,UChi_01) VMADDMEM(4,%r8,T2,UChi_11)	  \
   VMADDMEM(7,%r8,T1,UChi_02) VMADDMEM(7,%r8,T2,UChi_12)	  \
   /*22 cycles*/						  \
   VSHUFMEM(2,%r8,Z0)						  \
   VRDUP(Chi_02,T1)        VIDUP(Chi_02,Chi_02)			  \
   VRDUP(Chi_12,T2)        VIDUP(Chi_12,Chi_12)			  \
   VMADD(Z0,Chi_02,Z1)        VMADD(Z0,Chi_12,Z2)		  \
   VSHUFMEM(5,%r8,Z0)						  \
   VMADD(Z0,Chi_02,Z3)        VMADD(Z0,Chi_12,Z4)		  \
   VSHUFMEM(8,%r8,Z0)						  \
   VMADD(Z0,Chi_02,Z5)        VMADD(Z0,Chi_12,Z6)		  \
   /*33 cycles*/						  \
   VMADDSUBMEM(2,%r8,T1,Z1)   VMADDSUBMEM(2,%r8,T2,Z2)		  \
   VMADDSUBMEM(5,%r8,T1,Z3)   VMADDSUBMEM(5,%r8,T2,Z4)	          \
   VMADDSUBMEM(8,%r8,T1,Z5)   VMADDSUBMEM(8,%r8,T2,Z6)	       \
  /*stall*/						       \
  /*stall*/						       \
  /*stall*/						       \
  VADD(Z1,UChi_00,UChi_00)   VADD(Z2,UChi_10,UChi_10)	       \
  VADD(Z3,UChi_01,UChi_01)   VADD(Z4,UChi_11,UChi_11)	       \
  VADD(Z5,UChi_02,UChi_02)   VADD(Z6,UChi_12,UChi_12) )


#endif
