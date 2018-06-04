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

#ifdef AVX512
#include <simd/Intel512common.h>
#include <simd/Intel512avx.h>
#endif

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

#define Z00 %zmm26
#define Z10 %zmm27
#define Z0 Z00
#define Z1 %zmm28
#define Z2 %zmm29

#define Z3 %zmm30
#define Z4 %zmm31
#define Z5 Chi_31
#define Z6 Chi_32

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

#define MULT_ADD_XYZTa(g0,g1)					\
  asm ( "movq %0, %%r8 \n\t"					\
	"movq %1, %%r9 \n\t"	 :  : "r"(g0), "r"(g1) : "%r8","%r9");\
	   __asm__ (						\
	   VSHUF(Chi_00,T0)				\
	   VSHUF(Chi_10,T1)						\
	   VMOVIDUP(0,%r8,Z0 )						\
           VMOVIDUP(3,%r8,Z1 )						\
           VMOVIDUP(6,%r8,Z2 )						\
           VMADDSUB(Z0,T0,UChi_00)					\
	   VMADDSUB(Z1,T0,UChi_01)					\
	   VMADDSUB(Z2,T0,UChi_02)					\
									\
	   VMOVIDUP(0,%r9,Z0 )						\
           VMOVIDUP(3,%r9,Z1 )						\
           VMOVIDUP(6,%r9,Z2 )						\
           VMADDSUB(Z0,T1,UChi_10)					\
           VMADDSUB(Z1,T1,UChi_11)            \
           VMADDSUB(Z2,T1,UChi_12)            \
	   							\
								\
	   VMOVRDUP(0,%r8,Z3 )					\
	   VMOVRDUP(3,%r8,Z4 )					\
	   VMOVRDUP(6,%r8,Z5 )					\
           VMADDSUB(Z3,Chi_00,UChi_00)/*rr * ir = ri rr*/	\
           VMADDSUB(Z4,Chi_00,UChi_01)				\
           VMADDSUB(Z5,Chi_00,UChi_02)				\
								\
	   VMOVRDUP(0,%r9,Z3 )					\
	   VMOVRDUP(3,%r9,Z4 )					\
	   VMOVRDUP(6,%r9,Z5 )					\
           VMADDSUB(Z3,Chi_10,UChi_10)				\
           VMADDSUB(Z4,Chi_10,UChi_11)\
           VMADDSUB(Z5,Chi_10,UChi_12)				\
	   							\
								\
	   VMOVIDUP(1,%r8,Z0 )					\
	   VMOVIDUP(4,%r8,Z1 )					\
	   VMOVIDUP(7,%r8,Z2 )					\
	   VSHUF(Chi_01,T0)					\
           VMADDSUB(Z0,T0,UChi_00)				\
           VMADDSUB(Z1,T0,UChi_01)				\
           VMADDSUB(Z2,T0,UChi_02)				\
								\
	   VMOVIDUP(1,%r9,Z0 )					\
	   VMOVIDUP(4,%r9,Z1 )					\
	   VMOVIDUP(7,%r9,Z2 )					\
	   VSHUF(Chi_11,T1)					\
           VMADDSUB(Z0,T1,UChi_10)				\
           VMADDSUB(Z1,T1,UChi_11)				\
           VMADDSUB(Z2,T1,UChi_12)				\
								\
	   VMOVRDUP(1,%r8,Z3 )					\
	   VMOVRDUP(4,%r8,Z4 )					\
	   VMOVRDUP(7,%r8,Z5 )					\
           VMADDSUB(Z3,Chi_01,UChi_00)				\
           VMADDSUB(Z4,Chi_01,UChi_01)				\
           VMADDSUB(Z5,Chi_01,UChi_02)				\
								\
	   VMOVRDUP(1,%r9,Z3 )					\
	   VMOVRDUP(4,%r9,Z4 )					\
	   VMOVRDUP(7,%r9,Z5 )					\
           VMADDSUB(Z3,Chi_11,UChi_10)				\
           VMADDSUB(Z4,Chi_11,UChi_11)				\
           VMADDSUB(Z5,Chi_11,UChi_12)				\
	   							\
	   VSHUF(Chi_02,T0)					\
	   VSHUF(Chi_12,T1)					\
	   VMOVIDUP(2,%r8,Z0 )					\
	   VMOVIDUP(5,%r8,Z1 )					\
	   VMOVIDUP(8,%r8,Z2 )					\
           VMADDSUB(Z0,T0,UChi_00)				\
           VMADDSUB(Z1,T0,UChi_01)			      \
           VMADDSUB(Z2,T0,UChi_02)			      \
	   VMOVIDUP(2,%r9,Z0 )					\
	   VMOVIDUP(5,%r9,Z1 )					\
	   VMOVIDUP(8,%r9,Z2 )					\
           VMADDSUB(Z0,T1,UChi_10)			      \
           VMADDSUB(Z1,T1,UChi_11)			      \
           VMADDSUB(Z2,T1,UChi_12)			      \
	   /*55*/					      \
	   VMOVRDUP(2,%r8,Z3 )		  \
	   VMOVRDUP(5,%r8,Z4 )					\
	   VMOVRDUP(8,%r8,Z5 )				      \
           VMADDSUB(Z3,Chi_02,UChi_00)			      \
           VMADDSUB(Z4,Chi_02,UChi_01)			      \
           VMADDSUB(Z5,Chi_02,UChi_02)			      \
	   VMOVRDUP(2,%r9,Z3 )		  \
	   VMOVRDUP(5,%r9,Z4 )					\
	   VMOVRDUP(8,%r9,Z5 )				      \
           VMADDSUB(Z3,Chi_12,UChi_10)			      \
           VMADDSUB(Z4,Chi_12,UChi_11)			      \
           VMADDSUB(Z5,Chi_12,UChi_12)			      \
	   /*61 insns*/							);

#define MULT_ADD_XYZT(g0,g1)					\
  asm ( "movq %0, %%r8 \n\t"					\
	"movq %1, %%r9 \n\t"	 :  : "r"(g0), "r"(g1) : "%r8","%r9");\
  __asm__ (							  \
  VSHUFMEM(0,%r8,Z00)		   VSHUFMEM(0,%r9,Z10)			\
  VRDUP(Chi_00,T0)           VIDUP(Chi_00,Chi_00)	          \
   VRDUP(Chi_10,T1)           VIDUP(Chi_10,Chi_10)		  \
   VMUL(Z00,Chi_00,Z1)        VMUL(Z10,Chi_10,Z2)		  \
   VSHUFMEM(3,%r8,Z00)	      VSHUFMEM(3,%r9,Z10)		  \
   VMUL(Z00,Chi_00,Z3)        VMUL(Z10,Chi_10,Z4)		  \
   VSHUFMEM(6,%r8,Z00)	      VSHUFMEM(6,%r9,Z10)		  \
   VMUL(Z00,Chi_00,Z5)        VMUL(Z10,Chi_10,Z6)		  \
   VMADDMEM(0,%r8,T0,UChi_00)  VMADDMEM(0,%r9,T1,UChi_10)		  \
   VMADDMEM(3,%r8,T0,UChi_01)  VMADDMEM(3,%r9,T1,UChi_11)		  \
   VMADDMEM(6,%r8,T0,UChi_02)  VMADDMEM(6,%r9,T1,UChi_12)		  \
   VSHUFMEM(1,%r8,Z00)	      VSHUFMEM(1,%r9,Z10)		  \
   VRDUP(Chi_01,T0)           VIDUP(Chi_01,Chi_01)		  \
   VRDUP(Chi_11,T1)           VIDUP(Chi_11,Chi_11)		  \
   VMADD(Z00,Chi_01,Z1)       VMADD(Z10,Chi_11,Z2)		  \
   VSHUFMEM(4,%r8,Z00)	      VSHUFMEM(4,%r9,Z10)		  \
   VMADD(Z00,Chi_01,Z3)       VMADD(Z10,Chi_11,Z4)		  \
   VSHUFMEM(7,%r8,Z00)	      VSHUFMEM(7,%r9,Z10)		  \
   VMADD(Z00,Chi_01,Z5)       VMADD(Z10,Chi_11,Z6)		  \
   VMADDMEM(1,%r8,T0,UChi_00) VMADDMEM(1,%r9,T1,UChi_10)	  \
   VMADDMEM(4,%r8,T0,UChi_01) VMADDMEM(4,%r9,T1,UChi_11)	  \
   VMADDMEM(7,%r8,T0,UChi_02) VMADDMEM(7,%r9,T1,UChi_12)	  \
   VSHUFMEM(2,%r8,Z00)	      VSHUFMEM(2,%r9,Z10)			\
   VRDUP(Chi_02,T0)           VIDUP(Chi_02,Chi_02)			\
   VRDUP(Chi_12,T1)           VIDUP(Chi_12,Chi_12)			\
   VMADD(Z00,Chi_02,Z1)       VMADD(Z10,Chi_12,Z2)		  \
   VSHUFMEM(5,%r8,Z00)	      VSHUFMEM(5,%r9,Z10)		  \
   VMADD(Z00,Chi_02,Z3)       VMADD(Z10,Chi_12,Z4)		  \
   VSHUFMEM(8,%r8,Z00)	      VSHUFMEM(8,%r9,Z10)		  \
   VMADD(Z00,Chi_02,Z5)       VMADD(Z10,Chi_12,Z6)		  \
   VMADDSUBMEM(2,%r8,T0,Z1)   VMADDSUBMEM(2,%r9,T1,Z2)		  \
   VMADDSUBMEM(5,%r8,T0,Z3)   VMADDSUBMEM(5,%r9,T1,Z4)	          \
   VMADDSUBMEM(8,%r8,T0,Z5)   VMADDSUBMEM(8,%r9,T1,Z6)	       \
   VADD(Z1,UChi_00,UChi_00)   VADD(Z2,UChi_10,UChi_10)	       \
   VADD(Z3,UChi_01,UChi_01)   VADD(Z4,UChi_11,UChi_11)	       \
   VADD(Z5,UChi_02,UChi_02)   VADD(Z6,UChi_12,UChi_12) );

#define MULT_XYZT(g0,g1)					\
    asm ( "movq %0, %%r8 \n\t"						\
	"movq %1, %%r9 \n\t" :  : "r"(g0), "r"(g1) : "%r8","%r9" ); \
	   __asm__ (						\
	   VSHUF(Chi_00,T0)				\
	   VSHUF(Chi_10,T1)						\
	   VMOVIDUP(0,%r8,Z0 )						\
           VMOVIDUP(3,%r8,Z1 )						\
           VMOVIDUP(6,%r8,Z2 )						\
	   /*6*/							\
           VMUL(Z0,T0,UChi_00)            \
           VMUL(Z1,T0,UChi_01)            \
           VMUL(Z2,T0,UChi_02)            \
	   VMOVIDUP(0,%r9,Z0 )						\
           VMOVIDUP(3,%r9,Z1 )						\
           VMOVIDUP(6,%r9,Z2 )						\
           VMUL(Z0,T1,UChi_10)            \
           VMUL(Z1,T1,UChi_11)            \
           VMUL(Z2,T1,UChi_12)            \
	   VMOVRDUP(0,%r8,Z3 )					\
	   VMOVRDUP(3,%r8,Z4 )					\
	   VMOVRDUP(6,%r8,Z5 )					\
	   /*18*/						\
           VMADDSUB(Z3,Chi_00,UChi_00)				\
           VMADDSUB(Z4,Chi_00,UChi_01)\
           VMADDSUB(Z5,Chi_00,UChi_02) \
	   VMOVRDUP(0,%r9,Z3 )					\
	   VMOVRDUP(3,%r9,Z4 )					\
	   VMOVRDUP(6,%r9,Z5 )					\
           VMADDSUB(Z3,Chi_10,UChi_10)				\
           VMADDSUB(Z4,Chi_10,UChi_11)\
           VMADDSUB(Z5,Chi_10,UChi_12)				\
	   VMOVIDUP(1,%r8,Z0 )					\
	   VMOVIDUP(4,%r8,Z1 )					\
	   VMOVIDUP(7,%r8,Z2 )					\
	   /*28*/						\
	   VSHUF(Chi_01,T0)					\
           VMADDSUB(Z0,T0,UChi_00)      \
           VMADDSUB(Z1,T0,UChi_01)       \
           VMADDSUB(Z2,T0,UChi_02)        \
	   VMOVIDUP(1,%r9,Z0 )					\
	   VMOVIDUP(4,%r9,Z1 )					\
	   VMOVIDUP(7,%r9,Z2 )					\
	   VSHUF(Chi_11,T1)					\
           VMADDSUB(Z0,T1,UChi_10)				\
           VMADDSUB(Z1,T1,UChi_11)				\
           VMADDSUB(Z2,T1,UChi_12)        \
	   VMOVRDUP(1,%r8,Z3 )					\
	   VMOVRDUP(4,%r8,Z4 )					\
	   VMOVRDUP(7,%r8,Z5 )					\
           /*38*/						\
           VMADDSUB(Z3,Chi_01,UChi_00)    \
           VMADDSUB(Z4,Chi_01,UChi_01)    \
           VMADDSUB(Z5,Chi_01,UChi_02)    \
	   VMOVRDUP(1,%r9,Z3 )					\
	   VMOVRDUP(4,%r9,Z4 )					\
	   VMOVRDUP(7,%r9,Z5 )					\
           VMADDSUB(Z3,Chi_11,UChi_10)				\
           VMADDSUB(Z4,Chi_11,UChi_11)    \
           VMADDSUB(Z5,Chi_11,UChi_12)				\
	   /*48*/						\
	   VSHUF(Chi_02,T0)					\
	   VSHUF(Chi_12,T1)					\
	   VMOVIDUP(2,%r8,Z0 )					\
	   VMOVIDUP(5,%r8,Z1 )					\
	   VMOVIDUP(8,%r8,Z2 )					\
           VMADDSUB(Z0,T0,UChi_00)				\
           VMADDSUB(Z1,T0,UChi_01)			      \
           VMADDSUB(Z2,T0,UChi_02)			      \
	   VMOVIDUP(2,%r9,Z0 )					\
	   VMOVIDUP(5,%r9,Z1 )					\
	   VMOVIDUP(8,%r9,Z2 )					\
           VMADDSUB(Z0,T1,UChi_10)			      \
           VMADDSUB(Z1,T1,UChi_11)			      \
           VMADDSUB(Z2,T1,UChi_12)			      \
	   /*55*/					      \
	   VMOVRDUP(2,%r8,Z3 )		  \
	   VMOVRDUP(5,%r8,Z4 )					\
	   VMOVRDUP(8,%r8,Z5 )				      \
           VMADDSUB(Z3,Chi_02,UChi_00)			      \
           VMADDSUB(Z4,Chi_02,UChi_01)			      \
           VMADDSUB(Z5,Chi_02,UChi_02)			      \
	   VMOVRDUP(2,%r9,Z3 )		  \
	   VMOVRDUP(5,%r9,Z4 )					\
	   VMOVRDUP(8,%r9,Z5 )				      \
           VMADDSUB(Z3,Chi_12,UChi_10)			      \
           VMADDSUB(Z4,Chi_12,UChi_11)			      \
           VMADDSUB(Z5,Chi_12,UChi_12)			      \
	   /*61 insns*/							);

#define MULT_XYZTa(g0,g1)					\
  asm ( "movq %0, %%r8 \n\t"					\
	"movq %1, %%r9 \n\t" :  : "r"(g0), "r"(g1) : "%r8","%r9" ); \
  __asm__ (							  \
   VSHUFMEM(0,%r8,Z00)		   VSHUFMEM(0,%r9,Z10)	  \
   VRDUP(Chi_00,T0)           VIDUP(Chi_00,Chi_00)	          \
   VRDUP(Chi_10,T1)           VIDUP(Chi_10,Chi_10)		  \
   VMUL(Z00,Chi_00,Z1)        VMUL(Z10,Chi_10,Z2)		  \
   VSHUFMEM(3,%r8,Z00)	      VSHUFMEM(3,%r9,Z10)		  \
   VMUL(Z00,Chi_00,Z3)        VMUL(Z10,Chi_10,Z4)		  \
   VSHUFMEM(6,%r8,Z00)	      VSHUFMEM(6,%r9,Z10)		  \
   VMUL(Z00,Chi_00,Z5)        VMUL(Z10,Chi_10,Z6)		  \
   VMULMEM(0,%r8,T0,UChi_00)  VMULMEM(0,%r9,T1,UChi_10)		  \
   VMULMEM(3,%r8,T0,UChi_01)  VMULMEM(3,%r9,T1,UChi_11)		  \
   VMULMEM(6,%r8,T0,UChi_02)  VMULMEM(6,%r9,T1,UChi_12)		  \
   VSHUFMEM(1,%r8,Z00)	      VSHUFMEM(1,%r9,Z10)		  \
   VRDUP(Chi_01,T0)           VIDUP(Chi_01,Chi_01)		  \
   VRDUP(Chi_11,T1)           VIDUP(Chi_11,Chi_11)		  \
   VMADD(Z00,Chi_01,Z1)       VMADD(Z10,Chi_11,Z2)		  \
   VSHUFMEM(4,%r8,Z00)	      VSHUFMEM(4,%r9,Z10)		  \
   VMADD(Z00,Chi_01,Z3)       VMADD(Z10,Chi_11,Z4)		  \
   VSHUFMEM(7,%r8,Z00)	      VSHUFMEM(7,%r9,Z10)		  \
   VMADD(Z00,Chi_01,Z5)       VMADD(Z10,Chi_11,Z6)		  \
   VMADDMEM(1,%r8,T0,UChi_00) VMADDMEM(1,%r9,T1,UChi_10)	  \
   VMADDMEM(4,%r8,T0,UChi_01) VMADDMEM(4,%r9,T1,UChi_11)	  \
   VMADDMEM(7,%r8,T0,UChi_02) VMADDMEM(7,%r9,T1,UChi_12)	  \
   VSHUFMEM(2,%r8,Z00)	      VSHUFMEM(2,%r9,Z10)			\
   VRDUP(Chi_02,T0)           VIDUP(Chi_02,Chi_02)			\
   VRDUP(Chi_12,T1)           VIDUP(Chi_12,Chi_12)			\
   VMADD(Z00,Chi_02,Z1)       VMADD(Z10,Chi_12,Z2)		  \
   VSHUFMEM(5,%r8,Z00)	      VSHUFMEM(5,%r9,Z10)		  \
   VMADD(Z00,Chi_02,Z3)       VMADD(Z10,Chi_12,Z4)		  \
   VSHUFMEM(8,%r8,Z00)	      VSHUFMEM(8,%r9,Z10)		  \
   VMADD(Z00,Chi_02,Z5)       VMADD(Z10,Chi_12,Z6)		  \
   VMADDSUBMEM(2,%r8,T0,Z1)   VMADDSUBMEM(2,%r9,T1,Z2)		  \
   VMADDSUBMEM(5,%r8,T0,Z3)   VMADDSUBMEM(5,%r9,T1,Z4)	          \
   VMADDSUBMEM(8,%r8,T0,Z5)   VMADDSUBMEM(8,%r9,T1,Z6)	       \
   VADD(Z1,UChi_00,UChi_00)   VADD(Z2,UChi_10,UChi_10)	       \
   VADD(Z3,UChi_01,UChi_01)   VADD(Z4,UChi_11,UChi_11)	       \
   VADD(Z5,UChi_02,UChi_02)   VADD(Z6,UChi_12,UChi_12) );


#define LOAD_CHI(a0,a1,a2,a3)						\
  asm (									\
       "movq %0, %%r8 \n\t"						\
       VLOAD(0,%%r8,pChi_00)						\
       VLOAD(1,%%r8,pChi_01)						\
       VLOAD(2,%%r8,pChi_02)						\
       : : "r" (a0) : "%r8" );						\
  asm (									\
       "movq %0, %%r8 \n\t"						\
       VLOAD(0,%%r8,pChi_10)						\
       VLOAD(1,%%r8,pChi_11)						\
       VLOAD(2,%%r8,pChi_12)						\
       : : "r" (a1) : "%r8" );						\
  asm (									\
       "movq %0, %%r8 \n\t"						\
       VLOAD(0,%%r8,pChi_20)						\
       VLOAD(1,%%r8,pChi_21)						\
       VLOAD(2,%%r8,pChi_22)						\
       : : "r" (a2) : "%r8" );						\
  asm (									\
       "movq %0, %%r8 \n\t"						\
       VLOAD(0,%%r8,pChi_30)						\
       VLOAD(1,%%r8,pChi_31)						\
       VLOAD(2,%%r8,pChi_32)						\
       : : "r" (a3) : "%r8" );						

#define LOAD_CHIa(a0,a1)						\
  asm (									\
       "movq %0, %%r8 \n\t"						\
       VLOAD(0,%%r8,pChi_00)						\
       VLOAD(1,%%r8,pChi_01)						\
       VLOAD(2,%%r8,pChi_02)						\
       : : "r" (a0) : "%r8" );						\
  asm (									\
       "movq %0, %%r8 \n\t"						\
       VLOAD(0,%%r8,pChi_10)						\
       VLOAD(1,%%r8,pChi_11)						\
       VLOAD(2,%%r8,pChi_12)						\
       : : "r" (a1) : "%r8" );						

#define PF_CHI(a0)							
#define PF_CHIa(a0)							\
  asm (									\
       "movq %0, %%r8 \n\t"						\
       VPREFETCH1(0,%%r8)						\
       VPREFETCH1(1,%%r8)						\
       VPREFETCH1(2,%%r8)						\
       : : "r" (a0) : "%r8" );						\

#define PF_GAUGE_XYZT(a0)							
#define PF_GAUGE_XYZTa(a0)						\
  asm (									\
       "movq %0, %%r8 \n\t"						\
       VPREFETCH1(0,%%r8)						\
       VPREFETCH1(1,%%r8)						\
       VPREFETCH1(2,%%r8)						\
       VPREFETCH1(3,%%r8)						\
       VPREFETCH1(4,%%r8)						\
       VPREFETCH1(5,%%r8)						\
       VPREFETCH1(6,%%r8)						\
       VPREFETCH1(7,%%r8)						\
       VPREFETCH1(8,%%r8)						\
       : : "r" (a0) : "%r8" );						\

#define PF_GAUGE_LS(a0)							
#define PF_GAUGE_LSa(a0)							\
  asm (									\
       "movq %0, %%r8 \n\t"						\
       VPREFETCH1(0,%%r8)						\
       VPREFETCH1(1,%%r8)						\
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
  asm (								\
       VSTORE(0,%0,pUChi_00)					\
       VSTORE(1,%0,pUChi_01)					\
       VSTORE(2,%0,pUChi_02)					\
       : : "r" (out) : "memory" );

#define nREDUCE(out)							\
  asm (									\
       VADD(UChi_00,UChi_10,UChi_00)					\
       VADD(UChi_01,UChi_11,UChi_01)					\
       VADD(UChi_02,UChi_12,UChi_02)					\
       VADD(UChi_30,UChi_20,UChi_30)					\
       VADD(UChi_31,UChi_21,UChi_31)					\
       VADD(UChi_32,UChi_22,UChi_32)					\
       VADD(UChi_00,UChi_30,UChi_00)					\
       VADD(UChi_01,UChi_31,UChi_01)					\
       VADD(UChi_02,UChi_32,UChi_02)				);	\
  asm (VZERO(Chi_00)							\
       VSUB(UChi_00,Chi_00,UChi_00)					\
       VSUB(UChi_01,Chi_00,UChi_01)					\
       VSUB(UChi_02,Chi_00,UChi_02)				);	\
  asm (								\
       VSTORE(0,%0,pUChi_00)					\
       VSTORE(1,%0,pUChi_01)					\
       VSTORE(2,%0,pUChi_02)					\
       : : "r" (out) : "memory" );

#define REDUCEa(out)					\
  asm (							\
  VADD(UChi_00,UChi_10,UChi_00)				\
  VADD(UChi_01,UChi_11,UChi_01)				\
  VADD(UChi_02,UChi_12,UChi_02)	);			\
  asm (									\
       VSTORE(0,%0,pUChi_00)						\
       VSTORE(1,%0,pUChi_01)						\
       VSTORE(2,%0,pUChi_02)						\
       : : "r" (out) : "memory" );

// FIXME is sign right in the VSUB ?
#define nREDUCEa(out)					\
  asm (							\
  VADD(UChi_00,UChi_10,UChi_00)				\
  VADD(UChi_01,UChi_11,UChi_01)				\
  VADD(UChi_02,UChi_12,UChi_02)	);			\
  asm (VZERO(Chi_00)							\
       VSUB(UChi_00,Chi_00,UChi_00)					\
       VSUB(UChi_01,Chi_00,UChi_01)					\
       VSUB(UChi_02,Chi_00,UChi_02)				);	\
  asm (									\
       VSTORE(0,%0,pUChi_00)				\
       VSTORE(1,%0,pUChi_01)				\
       VSTORE(2,%0,pUChi_02)				\
       : : "r" (out) : "memory" );

#define PERMUTE_DIR(dir)			\
      permute##dir(Chi_0,Chi_0);\
      permute##dir(Chi_1,Chi_1);\
      permute##dir(Chi_2,Chi_2);

namespace Grid {
namespace QCD {

template <class Impl>
void StaggeredKernels<Impl>::DhopSiteAsm(StencilImpl &st, LebesgueOrder &lo, 
					 DoubledGaugeField &U, DoubledGaugeField &UUU,
					 SiteSpinor *buf, int LLs, int sU, 
					 const FermionField &in, FermionField &out,int dag) 
{
  assert(0);
};


//#define CONDITIONAL_MOVE(l,o,out) if ( l ) { out = (uint64_t) &in._odata[o] ; } else { out =(uint64_t) &buf[o]; }

#define CONDITIONAL_MOVE(l,o,out) { const SiteSpinor *ptr = l? in_p : buf; out = (uint64_t) &ptr[o]; }

#define PREPARE_XYZT(X,Y,Z,T,skew,UU)			\
  PREPARE(X,Y,Z,T,skew,UU);				\
  PF_GAUGE_XYZT(gauge0);					\
  PF_GAUGE_XYZT(gauge1);					\
  PF_GAUGE_XYZT(gauge2);					\
  PF_GAUGE_XYZT(gauge3);					

#define PREPARE_LS(X,Y,Z,T,skew,UU)			\
  PREPARE(X,Y,Z,T,skew,UU);				\
  PF_GAUGE_LS(gauge0);					\
  PF_GAUGE_LS(gauge1);					\
  PF_GAUGE_LS(gauge2);					\
  PF_GAUGE_LS(gauge3);					

#define PREPARE(X,Y,Z,T,skew,UU)					\
  SE0=st.GetEntry(ptype,X+skew,sF);					\
  o0 = SE0->_offset;							\
  l0 = SE0->_is_local;							\
  p0 = SE0->_permute;							\
  CONDITIONAL_MOVE(l0,o0,addr0);					\
  PF_CHI(addr0);							\
  									\
  SE1=st.GetEntry(ptype,Y+skew,sF);					\
  o1 = SE1->_offset;							\
  l1 = SE1->_is_local;							\
  p1 = SE1->_permute;							\
  CONDITIONAL_MOVE(l1,o1,addr1);					\
  PF_CHI(addr1);							\
  									\
  SE2=st.GetEntry(ptype,Z+skew,sF);					\
  o2 = SE2->_offset;							\
  l2 = SE2->_is_local;							\
  p2 = SE2->_permute;							\
  CONDITIONAL_MOVE(l2,o2,addr2);					\
  PF_CHI(addr2);							\
  									\
  SE3=st.GetEntry(ptype,T+skew,sF);					\
  o3 = SE3->_offset;							\
  l3 = SE3->_is_local;							\
  p3 = SE3->_permute;							\
  CONDITIONAL_MOVE(l3,o3,addr3);					\
  PF_CHI(addr3);							\
  									\
  gauge0 =(uint64_t)&UU._odata[sU]( X );				\
  gauge1 =(uint64_t)&UU._odata[sU]( Y );				\
  gauge2 =(uint64_t)&UU._odata[sU]( Z );				\
  gauge3 =(uint64_t)&UU._odata[sU]( T ); 
  
  // This is the single precision 5th direction vectorised kernel
#include <simd/Intel512single.h>
template <> void StaggeredKernels<StaggeredVec5dImplF>::DhopSiteAsm(StencilImpl &st, LebesgueOrder &lo, 
								    DoubledGaugeField &U, DoubledGaugeField &UUU,
								    SiteSpinor *buf, int LLs, int sU, 
								    const FermionField &in, FermionField &out,int dag) 
{
#ifdef AVX512
  uint64_t gauge0,gauge1,gauge2,gauge3;
  uint64_t addr0,addr1,addr2,addr3;
  const SiteSpinor *in_p; in_p = &in._odata[0];

  int o0,o1,o2,o3; // offsets
  int l0,l1,l2,l3; // local 
  int p0,p1,p2,p3; // perm
  int ptype;
  StencilEntry *SE0;
  StencilEntry *SE1;
  StencilEntry *SE2;
  StencilEntry *SE3;

   for(int s=0;s<LLs;s++){

    int sF=s+LLs*sU;
    // Xp, Yp, Zp, Tp
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

    addr0 = (uint64_t) &out._odata[sF];
    if ( dag ) {
      nREDUCE(addr0);
    } else { 
      REDUCE(addr0);
    }
   }
#else 
    assert(0);
#endif
   
}

#include <simd/Intel512double.h>
template <> void StaggeredKernels<StaggeredVec5dImplD>::DhopSiteAsm(StencilImpl &st, LebesgueOrder &lo, 
								    DoubledGaugeField &U, DoubledGaugeField &UUU,
								    SiteSpinor *buf, int LLs, int sU, 
								    const FermionField &in, FermionField &out,int dag) 
{
#ifdef AVX512
  uint64_t gauge0,gauge1,gauge2,gauge3;
  uint64_t addr0,addr1,addr2,addr3;
  const SiteSpinor *in_p; in_p = &in._odata[0];

  int o0,o1,o2,o3; // offsets
  int l0,l1,l2,l3; // local 
  int p0,p1,p2,p3; // perm
  int ptype;
  StencilEntry *SE0;
  StencilEntry *SE1;
  StencilEntry *SE2;
  StencilEntry *SE3;

  for(int s=0;s<LLs;s++){
    int sF=s+LLs*sU;
    // Xp, Yp, Zp, Tp
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

    addr0 = (uint64_t) &out._odata[sF];
    if ( dag ) {
      nREDUCE(addr0);
    } else { 
      REDUCE(addr0);
    }
  }
#else 
  assert(0);
#endif
}
   
   


#define PERMUTE_DIR3 __asm__ (	\
  VPERM3(Chi_00,Chi_00)	\
  VPERM3(Chi_01,Chi_01)	\
  VPERM3(Chi_02,Chi_02)	);

#define PERMUTE_DIR2 __asm__ (	\
  VPERM2(Chi_10,Chi_10)	\
  VPERM2(Chi_11,Chi_11)	\
  VPERM2(Chi_12,Chi_12) );

#define PERMUTE_DIR1 __asm__ (	\
  VPERM1(Chi_00,Chi_00)	\
  VPERM1(Chi_01,Chi_01)	\
  VPERM1(Chi_02,Chi_02)	);

#define PERMUTE_DIR0 __asm__ (			\
  VPERM0(Chi_10,Chi_10)	\
  VPERM0(Chi_11,Chi_11)	\
  VPERM0(Chi_12,Chi_12) );

#define PERMUTE01 \
  if ( p0 ) { PERMUTE_DIR3; }\
  if ( p1 ) { PERMUTE_DIR2; }

#define PERMUTE23 \
  if ( p2 ) { PERMUTE_DIR1; }\
  if ( p3 ) { PERMUTE_DIR0; }

  // This is the single precision 5th direction vectorised kernel

#include <simd/Intel512single.h>
template <> void StaggeredKernels<StaggeredImplF>::DhopSiteAsm(StencilImpl &st, LebesgueOrder &lo, 
							       DoubledGaugeField &U, DoubledGaugeField &UUU,
							       SiteSpinor *buf, int LLs, int sU, 
							       const FermionField &in, FermionField &out,int dag) 
{
#ifdef AVX512
  uint64_t gauge0,gauge1,gauge2,gauge3;
  uint64_t addr0,addr1,addr2,addr3;
  const SiteSpinor *in_p; in_p = &in._odata[0];

  int o0,o1,o2,o3; // offsets
  int l0,l1,l2,l3; // local 
  int p0,p1,p2,p3; // perm
  int ptype;
  StencilEntry *SE0;
  StencilEntry *SE1;
  StencilEntry *SE2;
  StencilEntry *SE3;

  for(int s=0;s<LLs;s++){
    
    int sF=s+LLs*sU;
    // Xp, Yp, Zp, Tp
    PREPARE(Xp,Yp,Zp,Tp,0,U);
    LOAD_CHIa(addr0,addr1);
    PERMUTE01;
    MULT_XYZT(gauge0,gauge1);
    LOAD_CHIa(addr2,addr3);
    PERMUTE23;
    MULT_ADD_XYZT(gauge2,gauge3);  

    PREPARE(Xm,Ym,Zm,Tm,0,U);
    LOAD_CHIa(addr0,addr1);
    PERMUTE01;
    MULT_ADD_XYZT(gauge0,gauge1);
    LOAD_CHIa(addr2,addr3);
    PERMUTE23;
    MULT_ADD_XYZT(gauge2,gauge3);  

    PREPARE(Xp,Yp,Zp,Tp,8,UUU);
    LOAD_CHIa(addr0,addr1);
    PERMUTE01;
    MULT_ADD_XYZT(gauge0,gauge1);
    LOAD_CHIa(addr2,addr3);
    PERMUTE23;
    MULT_ADD_XYZT(gauge2,gauge3);  
    
    PREPARE(Xm,Ym,Zm,Tm,8,UUU);
    LOAD_CHIa(addr0,addr1);
    PERMUTE01;
    MULT_ADD_XYZT(gauge0,gauge1);
    LOAD_CHIa(addr2,addr3);
    PERMUTE23;
    MULT_ADD_XYZT(gauge2,gauge3);  

    addr0 = (uint64_t) &out._odata[sF];
    if ( dag ) { 
      nREDUCEa(addr0);
    } else { 
      REDUCEa(addr0);
    }
  }
#else 
  assert(0);
#endif
}

#include <simd/Intel512double.h>
template <> void StaggeredKernels<StaggeredImplD>::DhopSiteAsm(StencilImpl &st, LebesgueOrder &lo, 
							       DoubledGaugeField &U, DoubledGaugeField &UUU,
							       SiteSpinor *buf, int LLs, int sU, 
							       const FermionField &in, FermionField &out,int dag) 
{
#ifdef AVX512
  uint64_t gauge0,gauge1,gauge2,gauge3;
  uint64_t addr0,addr1,addr2,addr3;
  const SiteSpinor *in_p; in_p = &in._odata[0];

  int o0,o1,o2,o3; // offsets
  int l0,l1,l2,l3; // local 
  int p0,p1,p2,p3; // perm
  int ptype;
  StencilEntry *SE0;
  StencilEntry *SE1;
  StencilEntry *SE2;
  StencilEntry *SE3;

  for(int s=0;s<LLs;s++){
    
    int sF=s+LLs*sU;
    // Xp, Yp, Zp, Tp
    PREPARE(Xp,Yp,Zp,Tp,0,U);
    LOAD_CHIa(addr0,addr1);
    PERMUTE01;
    MULT_XYZT(gauge0,gauge1);
    LOAD_CHIa(addr2,addr3);
    PERMUTE23;
    MULT_ADD_XYZT(gauge2,gauge3);  
    
    PREPARE(Xm,Ym,Zm,Tm,0,U);
    LOAD_CHIa(addr0,addr1);
    PERMUTE01;
    MULT_ADD_XYZT(gauge0,gauge1);
    LOAD_CHIa(addr2,addr3);
    PERMUTE23;
    MULT_ADD_XYZT(gauge2,gauge3);  
    
    PREPARE(Xp,Yp,Zp,Tp,8,UUU);
    LOAD_CHIa(addr0,addr1);
    PERMUTE01;
    MULT_ADD_XYZT(gauge0,gauge1);
    LOAD_CHIa(addr2,addr3);
    PERMUTE23;
    MULT_ADD_XYZT(gauge2,gauge3);  
    
    PREPARE(Xm,Ym,Zm,Tm,8,UUU);
    LOAD_CHIa(addr0,addr1);
    PERMUTE01;
    MULT_ADD_XYZT(gauge0,gauge1);
    LOAD_CHIa(addr2,addr3);
    PERMUTE23;
    MULT_ADD_XYZT(gauge2,gauge3);  
    
    addr0 = (uint64_t) &out._odata[sF];
    if ( dag ) {
      nREDUCEa(addr0);
    } else { 
      REDUCEa(addr0);
    }
  }
#else 
  assert(0);
#endif
}

#define KERNEL_INSTANTIATE(CLASS,FUNC,IMPL)			    \
  template void CLASS<IMPL>::FUNC(StencilImpl &st, LebesgueOrder &lo,	\
				  DoubledGaugeField &U,			\
				  DoubledGaugeField &UUU,		\
				  SiteSpinor *buf, int LLs,		\
				  int sU, const FermionField &in, FermionField &out,int dag);

KERNEL_INSTANTIATE(StaggeredKernels,DhopSiteAsm,StaggeredImplD);
KERNEL_INSTANTIATE(StaggeredKernels,DhopSiteAsm,StaggeredImplF);
KERNEL_INSTANTIATE(StaggeredKernels,DhopSiteAsm,StaggeredVec5dImplD);
KERNEL_INSTANTIATE(StaggeredKernels,DhopSiteAsm,StaggeredVec5dImplF);

}}

