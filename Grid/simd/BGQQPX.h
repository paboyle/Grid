   /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/simd/BGQQPX.h

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
#ifndef GRID_ASM_BGQ_QPX_H
#define GRID_ASM_BGQ_QPX_H

#include <stddint.h>

/*********************************************************
 * Architectural macros
 *********************************************************/
#define VLOADf(OFF,PTR,DEST)         "qvlfsux     " #DEST  "," #OFF "," #PTR ") ;\n"
#define VLOADd(OFF,PTR,DEST)         "qvlfdux     " #DEST  "," #OFF "," #PTR ") ;\n"
#define VSTOREf(OFF,PTR,SRC)         "qvstfsux    " #SRC  "," #OFF "," #PTR ") ;\n"
#define VSTOREd(OFF,PTR,SRC)         "qvstfdux    " #SRC  "," #OFF "," #PTR ") ;\n"
#define VSPLATf(A,B,DEST)            "qvlfcdxa    " #A "," #B "," #DEST  ";\n"
#define VSPLATd(A,B,DEST)            "qvlfcsxa    " #A "," #B "," #DEST  ";\n"

#define LOAD64(A,ptr)  
#define VZERO(DEST)                  "qvfclr      " #DEST "; \n"
#define VONE (DEST)                  "qvfset      " #DEST "; \n" 
#define VNEG (SRC,DEST)              "qvfneg      " #DEST "," #SRC "; \n"
#define VMOV(A,DEST)                 "qvfmr       " #DEST, "," #A   ";\n"

#define VADD(A,B,DEST)               "qvfadd      " #DEST "," #A "," #B  ";\n"
#define VSUB(A,B,DEST)               "qvfsub      " #DEST "," #A "," #B  ";\n"
#define VMUL(A,B,DEST)               "qvfmul      " #DEST "," #A "," #B  ";\n"
#define VMUL_RR_RI(A,B,DEST)         "qvfxmul     " #DEST "," #A "," #B  ";\n" 
#define VMADD(A,B,C,DEST)            "qvfmadd     " #DEST "," #A "," #B ","#C ";\n"
#define VMADD_RR_RI(A,B,C,DEST)      "qvfxmadd    " #DEST "," #A "," #B ","#C ";\n" 
#define VMADD_MII_IR(A,B,C,DEST)     "qvfxxnpmadd " #DEST "," #A "," #B ","#C ";\n" 
#define VMADD_II_MIR(A,B,C,DEST)     "qvfmadd     " #DEST "," #A "," #B ","#C ";\n"  

#define CACHE_LOCK  (PTR)    asm (" dcbtls  %%r0, %0 \n" : : "r" (PTR) );
#define CACHE_UNLOCK(PTR)    asm (" dcblc   %%r0, %0 \n" : : "r" (PTR) );
#define CACHE_FLUSH (PTR)    asm (" dcbf    %%r0, %0 \n" : : "r" (PTR) );
#define CACHE_TOUCH (PTR)    asm (" dcbt    %%r0, %0 \n" : : "r" (PTR) );

// Gauge field locking 2 x 9 complex == 18*8 / 16 bytes per link
// This is 144/288 bytes == 4.5; 9 lines
#define MASK_REGS   /*NOOP ON BGQ*/
#define PF_GAUGE(A) /*NOOP ON BGQ*/
#define PREFETCH1_CHIMU(base) /*NOOP ON BGQ*/
#define PREFETCH_CHIMU(base) /*NOOP ON BGQ*/

/*********************************************************
 * Register definitions
 *********************************************************/
#define psi_00 0
#define psi_01 1
#define psi_02 2
  
#define psi_10 3
#define psi_11 4
#define psi_12 5

#define psi_20 6
#define psi_21 7
#define psi_22 8

#define psi_30 9
#define psi_31 10
#define psi_32 11

#define Chi_00 12
#define Chi_01 13
#define Chi_02 14

#define Chi_10 15
#define Chi_11 16
#define Chi_12 17  

#define UChi_00 18 
#define UChi_01 19
#define UChi_02 20

#define UChi_10 21
#define UChi_11 22
#define UChi_12 23 

#define U0 24
#define U1 25
#define U2 26
#define one 27

#define REP  %%r16
#define IMM  %%r17

/*Alias regs*/
#define Chimu_00 Chi_00
#define Chimu_01 Chi_01
#define Chimu_02 Chi_02
#define Chimu_10 Chi_10
#define Chimu_11 Chi_11
#define Chimu_12 Chi_02
#define Chimu_20 UChi_00
#define Chimu_21 UChi_01
#define Chimu_22 UChi_02
#define Chimu_30 UChi_10
#define Chimu_31 UChi_11
#define Chimu_32 UChi_02

/*********************************************************
 * Macro sequences encoding QCD
 *********************************************************/
#define LOCK_GAUGE(dir)							\
  {									\
    uint8_t *byte_addr = (uint8_t *)&U._odata[sU](dir);			\
    for(int i=0;i< 18*2*BYTES_PER_WORD*8;i+=32){			\
      CACHE_LOCK(&byte_addr[i]);					\
    }									\
  }

#define UNLOCK_GAUGE(dir)						\
  {									\
    uint8_t *byte_addr = (uint8_t *)&U._odata[sU](dir);			\
    for(int i=0;i< 18*2*BYTES_PER_WORD*8;i+=32){			\
      CACHE_UNLOCK(&byte_addr[i]);					\
    }									\
  }

#define MAYBEPERM(A,B)

#define PERMUTE_DIR3 
#define PERMUTE_DIR2 
#define PERMUTE_DIR1 
#define PERMUTE_DIR0 

#define MULT_2SPIN_DIR_PFXP(A,p) MULT_2SPIN(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFYP(A,p) MULT_2SPIN(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFZP(A,p) MULT_2SPIN(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFTP(A,p) MULT_2SPIN(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFXM(A,p) MULT_2SPIN(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFYM(A,p) MULT_2SPIN(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFZM(A,p) MULT_2SPIN(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFTM(A,p) MULT_2SPIN(&U._odata[sU](A),p)

#define MULT_SPIN(ptr,p) {					\
    uint64_t ub = ((uint64_t)base);				\
    asm (							\
         VLOAD(%0,%3,U0)					\
         VLOAD(%1,%3,U1)					\
         VLOAD(%2,%3,U2)					\
	 VMUL_RR_RI(U0,Chi_00,UChi_00)					\
	 VMUL_RR_RI(U1,Chi_00,UChi_01)					\
	 VMUL_RR_RI(U2,Chi_00,UChi_02)					\
	 VMUL_RR_RI(U0,Chi_10,UChi_10)					\
	 VMUL_RR_RI(U1,Chi_10,UChi_11)					\
	 VMUL_RR_RI(U2,Chi_10,UChi_12)					\
	 VMADD_MII_IR(U0,Chi_00,UChi_00,UChi_00)			\
	 VMADD_MII_IR(U1,Chi_00,UChi_01,UChi_01)			\
	 VMADD_MII_IR(U2,Chi_00,UChi_02,UChi_02)			\
	 VMADD_MII_IR(U0,Chi_10,UChi_10,UChi_10)			\
	 VMADD_MII_IR(U1,Chi_10,UChi_11,UChi_11)			\
	 VMADD_MII_IR(U2,Chi_10,UChi_12,UChi_12)			\
	 : : "r" (0), "r" (32*3), "r" (32*6), "r" (ub ));		\
    asm (								\
         VLOAD(%0,%3,U0)						\
         VLOAD(%1,%3,U1)						\
         VLOAD(%2,%3,U2)						\
	 VMADD_RR_RI(U0,Chi_01,UChi_00,UChi_00)				\
	 VMADD_RR_RI(U1,Chi_01,UChi_01,UChi_01)				\
	 VMADD_RR_RI(U2,Chi_01,UChi_02,UChi_02)				\
	 VMADD_RR_RI(U0,Chi_11,UChi_10,UChi_10)				\
	 VMADD_RR_RI(U1,Chi_11,UChi_11,UChi_11)				\
	 VMADD_RR_RI(U2,Chi_11,UChi_12,UChi_12)				\
	 VMADD_MII_IR(U0,Chi_01,UChi_00,UChi_00)			\
	 VMADD_MII_IR(U1,Chi_01,UChi_01,UChi_01)			\
	 VMADD_MII_IR(U2,Chi_01,UChi_02,UChi_02)			\
	 VMADD_MII_IR(U0,Chi_11,UChi_10,UChi_10)			\
	 VMADD_MII_IR(U1,Chi_11,UChi_11,UChi_11)			\
	 VMADD_MII_IR(U2,Chi_11,UChi_12,UChi_12)			\
	 : : "r" (32), "r" (32*4), "r" (32*7), "r" (ub ));		\
    asm (								\
         VLOAD(%0,%3,U0)						\
         VLOAD(%1,%3,U1)						\
         VLOAD(%2,%3,U2)						\
	 VMADD_RR_RI(U0,Chi_02,UChi_00,UChi_00)				\
	 VMADD_RR_RI(U1,Chi_02,UChi_01,UChi_01)				\
	 VMADD_RR_RI(U2,Chi_02,UChi_02,UChi_02)				\
	 VMADD_RR_RI(U0,Chi_12,UChi_10,UChi_10)				\
	 VMADD_RR_RI(U1,Chi_12,UChi_11,UChi_11)				\
	 VMADD_RR_RI(U2,Chi_12,UChi_12,UChi_12)				\
	 VMADD_MII_IR(U0,Chi_02,UChi_00,UChi_00)			\
	 VMADD_MII_IR(U1,Chi_02,UChi_01,UChi_01)			\
	 VMADD_MII_IR(U2,Chi_02,UChi_02,UChi_02)			\
	 VMADD_MII_IR(U0,Chi_12,UChi_10,UChi_10)			\
	 VMADD_MII_IR(U1,Chi_12,UChi_11,UChi_11)			\
	 VMADD_MII_IR(U2,Chi_12,UChi_12,UChi_12)			\
	 : : "r" (32*2), "r" (32*5), "r" (32*8), "r" (ub ));		\
  }

#define SAVE_RESULT(base,basep) {\
    uint64_t ub = ((uint64_t)base)  - 32;				\
    asm("mr %0,"REP";\n\t"						\
	"li " IMM ",32;\n\t"						\
	VSTORE(IMM,REP,psi_00)						\
	VSTORE(IMM,REP,psi_01)						\
	VSTORE(IMM,REP,psi_02)						\
	VSTORE(IMM,REP,psi_10)						\
	VSTORE(IMM,REP,psi_11)						\
	VSTORE(IMM,REP,psi_12)						\
	VSTORE(IMM,REP,psi_20)						\
	VSTORE(IMM,REP,psi_21)						\
	VSTORE(IMM,REP,psi_22)						\
	VSTORE(IMM,REP,psi_30)						\
	VSTORE(IMM,REP,psi_31)						\
	VSTORE(IMM,REP,psi_32)						\
	);								\
}

/*
 *Annoying BG/Q loads with no immediat indexing and big performance hit
 *when second miss to a L1 line occurs
 */
#define LOAD_CHI(base) {						\
    uint64_t ub = ((uint64_t)base)  - 64;				\
    asm("mr  %0,"REP";\n\t"						\
	"li " IMM ",64;\n\t"				    		\
	VLOAD(IMM,REP,Chi_00)						\
	VLOAD(IMM,REP,Chi_02)						\
	VLOAD(IMM,REP,Chi_11) : : "r" (ub) 	     );			\
    ub = ((uint64_t)base)  - 32;					\
    asm("mr  %0,"REP";\n\t"						\
	"li IMM,64;\n\t"				    		\
	VLOAD(IMM,REP,Chimu_01)						\
	VLOAD(IMM,REP,Chimu_10)						\
	VLOAD(IMM,REP,Chimu_12)	: : "r" (ub) 	     );			\
  }

#define LOAD_CHIMU(base) {						\
    uint64_t ub = ((uint64_t)base)  - 64;				\
    asm("mr  %0,"REP";\n\t"						\
	"li IMM,64;\n\t"				    		\
	VLOAD(IMM,REP,Chimu_00)						\
	VLOAD(IMM,REP,Chimu_02)						\
	VLOAD(IMM,REP,Chimu_11)						\
	VLOAD(IMM,REP,Chimu_20)						\
	VLOAD(IMM,REP,Chimu_22)						\
	VLOAD(IMM,REP,Chimu_31) : : "r" (ub) 	     );			\
    ub = ((uint64_t)base)  - 32;					\
    asm("mr  %0,"REP";\n\t"						\
	"li IMM,64;\n\t"				    		\
	VLOAD(IMM,REP,Chimu_01)						\
	VLOAD(IMM,REP,Chimu_10)						\
	VLOAD(IMM,REP,Chimu_12)						\
	VLOAD(IMM,REP,Chimu_21)						\
	VLOAD(IMM,REP,Chimu_30)						\
	VLOAD(IMM,REP,Chimu_32)	: : "r" (ub) 	     );			\
  }

//      hspin(0)=fspin(0)+timesI(fspin(3));
//      hspin(1)=fspin(1)+timesI(fspin(2));
#define XP_PROJMEM(base) {					\
    LOAD_CHIMU(base);						\
    asm (							\
         VONE(one)						\
	 VMADD_MII_IR(one,Chimu_30,Chimu_00,Chi_00)		\
	 VMADD_MII_IR(one,Chimu_31,Chimu_01,Chi_01)		\
	 VMADD_MII_IR(one,Chimu_32,Chimu_02,Chi_02)		\
	 VMADD_MII_IR(one,Chimu_20,Chimu_10,Chi_10)		\
	 VMADD_MII_IR(one,Chimu_21,Chimu_11,Chi_11)		\
	 VMADD_MII_IR(one,Chimu_22,Chimu_12,Chi_12)		\
							);	\
  }

#define XM_PROJMEM(base) {				\
    LOAD_CHIMU(base);					\
    asm (						\
         VONE(one)						\
	 VMADD_II_MIR(one,Chimu_30,Chimu_00,Chi_00)	\
	 VMADD_II_MIR(one,Chimu_31,Chimu_01,Chi_01)	\
	 VMADD_II_MIR(one,Chimu_32,Chimu_02,Chi_02)	\
	 VMADD_II_MIR(one,Chimu_20,Chimu_10,Chi_10)	\
	 VMADD_II_MIR(one,Chimu_21,Chimu_11,Chi_11)	\
	 VMADD_II_MIR(one,Chimu_22,Chimu_12,Chi_12)	\
							);	\
  }

//      hspin(0)=fspin(0)-fspin(3);
//      hspin(1)=fspin(1)+fspin(2);
#define YP_PROJMEM(base) {  \
    LOAD_CHIMU(base);						\
    asm (							\
	 VSUB(Chimu_00,Chimu_00,Chi_30)		\
	 VSUB(Chimu_01,Chimu_01,Chi_31)		\
	 VSUB(Chimu_02,Chimu_02,Chi_32)		\
	 VADD(Chimu_10,Chimu_10,Chi_20)		\
	 VADD(Chimu_11,Chimu_11,Chi_21)		\
	 VADD(Chimu_12,Chimu_12,Chi_22)		\
							);	\
  }

#define YM_PROJMEM(base) {			\
    LOAD_CHIMU(base);						\
    asm (							\
	 VADD(Chimu_00,Chimu_00,Chi_30)		\
	 VADD(Chimu_01,Chimu_01,Chi_31)		\
	 VADD(Chimu_02,Chimu_02,Chi_32)		\
	 VSUB(Chimu_10,Chimu_10,Chi_20)		\
	 VSUB(Chimu_11,Chimu_11,Chi_21)		\
	 VSUB(Chimu_12,Chimu_12,Chi_22)		\
							);	\
  }

	    /*Gz
	     *  0 0  i  0   [0]+-i[2]
	     *  0 0  0 -i   [1]-+i[3]
	     * -i 0  0  0
	     *  0 i  0  0
	     */
#define ZP_PROJMEM(base) {  \
    LOAD_CHIMU(base);						\
    asm (							\
         VONE(one)						\
	 VMADD_MII_IR(one,Chimu_20,Chimu_00,Chi_00)		\
	 VMADD_MII_IR(one,Chimu_21,Chimu_01,Chi_01)		\
	 VMADD_MII_IR(one,Chimu_22,Chimu_02,Chi_02)		\
	 VMADD_II_MIR(one,Chimu_30,Chimu_10,Chi_10)		\
	 VMADD_II_MIR(one,Chimu_31,Chimu_11,Chi_11)		\
	 VMADD_II_MIR(one,Chimu_32,Chimu_12,Chi_12)		\
							);	\
  }

#define ZM_PROJMEM(base) {  \
    LOAD_CHIMU(base);						\
    asm (							\
         VONE(one)						\
	 VMADD_II_MIR(one,Chimu_20,Chimu_00,Chi_00)		\
	 VMADD_II_MIR(one,Chimu_21,Chimu_01,Chi_01)		\
	 VMADD_II_MIR(one,Chimu_22,Chimu_02,Chi_02)		\
	 VMADD_MII_IR(one,Chimu_30,Chimu_10,Chi_10)		\
	 VMADD_MII_IR(one,Chimu_31,Chimu_11,Chi_11)		\
	 VMADD_MII_IR(one,Chimu_32,Chimu_12,Chi_12)		\
							);	\
  }
	    /*Gt
	     *  0 0  1  0 [0]+-[2]
	     *  0 0  0  1 [1]+-[3]
	     *  1 0  0  0
	     *  0 1  0  0
	     */
#define TP_PROJMEM(base) {  \
    LOAD_CHIMU(base);						\
    asm (							\
	 VADD(Chimu_00,Chimu_00,Chi_20)		\
	 VADD(Chimu_01,Chimu_01,Chi_21)		\
	 VADD(Chimu_02,Chimu_02,Chi_22)		\
	 VADD(Chimu_10,Chimu_10,Chi_30)		\
	 VADD(Chimu_11,Chimu_11,Chi_31)		\
	 VADD(Chimu_12,Chimu_12,Chi_32)		\
							);	\
  }

#define TM_PROJMEM(base) {  \
    LOAD_CHIMU(base);						\
    asm (							\
	 VSUB(Chimu_00,Chimu_00,Chi_20)		\
	 VSUB(Chimu_01,Chimu_01,Chi_21)		\
	 VSUB(Chimu_02,Chimu_02,Chi_22)		\
	 VSUB(Chimu_10,Chimu_10,Chi_30)		\
	 VSUB(Chimu_11,Chimu_11,Chi_31)		\
	 VSUB(Chimu_12,Chimu_12,Chi_32)		\
							);	\
  }

/*
      fspin(0)=hspin(0);
      fspin(1)=hspin(1);
      fspin(2)=timesMinusI(hspin(1));
      fspin(3)=timesMinusI(hspin(0));

      fspin(0)+=hspin(0);
      fspin(1)+=hspin(1);
      fspin(2)-=timesI(hspin(1));
      fspin(3)-=timesI(hspin(0));
 */
#define XP_RECON {				\
    asm(\
	VONE(one)\
	VMOV(psi_00,UChi_00) 	VMOV(psi_01,UChi_01)	VMOV(psi_02,UChi_02)\
	VMOV(psi_10,UChi_10) 	VMOV(psi_11,UChi_11)	VMOV(psi_12,UChi_12)\
	VZERO(psi_20)	VZERO(psi_21)	VZERO(psi_22) \
	VZERO(psi_30) 	VZERO(psi_31)   VZERO(psi_32) \
	VMADD_II_MIR(one,UChi_10,psi_20,psi_20)	      \
	VMADD_II_MIR(one,UChi_11,psi_21,psi_21)	      \
	VMADD_II_MIR(one,UChi_12,psi_22,psi_22)	      \
	VMADD_II_MIR(one,UChi_00,psi_30,psi_30)	      \
	VMADD_II_MIR(one,UChi_01,psi_31,psi_31)	      \
	VMADD_II_MIR(one,UChi_02,psi_32,psi_32)	      \
	);		     \
  }

#define XM_RECON {				\
    asm(\
	VONE(one)\
	VMOV(psi_00,UChi_00) 	VMOV(psi_01,UChi_01)	VMOV(psi_02,UChi_02)\
	VMOV(psi_10,UChi_10) 	VMOV(psi_11,UChi_11)	VMOV(psi_12,UChi_12)\
	VZERO(psi_20)	VZERO(psi_21)	VZERO(psi_22) \
	VZERO(psi_30) 	VZERO(psi_31)   VZERO(psi_32) \
	VMADD_MII_IR(one,UChi_10,psi_20,psi_20)	      \
	VMADD_MII_IR(one,UChi_11,psi_21,psi_21)	      \
	VMADD_MII_IR(one,UChi_12,psi_22,psi_22)	      \
	VMADD_MII_IR(one,UChi_00,psi_30,psi_30)	      \
	VMADD_MII_IR(one,UChi_01,psi_31,psi_31)	      \
	VMADD_MII_IR(one,UChi_02,psi_32,psi_32)	      \
	);		     \
  }

#define XP_RECON_ACCUM {				\
    asm(\
	VONE(one)\
	VADD(psi_00,UChi_00,psi_00) 	VADD(psi_01,UChi_01,psi_01)	VADD(psi_02,UChi_02,psi_02) \
	VADD(psi_10,UChi_10,psi_10) 	VADD(psi_11,UChi_11,psi_11)	VADD(psi_12,UChi_12,psi_12) \
	VMADD_II_MIR(one,UChi_10,psi_20,psi_20)	      \
	VMADD_II_MIR(one,UChi_11,psi_21,psi_21)	      \
	VMADD_II_MIR(one,UChi_12,psi_22,psi_22)	      \
	VMADD_II_MIR(one,UChi_00,psi_30,psi_30)	      \
	VMADD_II_MIR(one,UChi_01,psi_31,psi_31)	      \
	VMADD_II_MIR(one,UChi_02,psi_32,psi_32)	      \
	);		     \
  }

#define XM_RECON_ACCUM {				\
    asm(\
	VONE(one)\
	VADD(psi_00,UChi_00,psi_00) 	VADD(psi_01,UChi_01,psi_01)	VADD(psi_02,UChi_02,psi_02) \
	VADD(psi_10,UChi_10,psi_10) 	VADD(psi_11,UChi_11,psi_11)	VADD(psi_12,UChi_12,psi_12) \
	VMADD_MII_IR(one,UChi_10,psi_20,psi_20)	      \
	VMADD_MII_IR(one,UChi_11,psi_21,psi_21)	      \
	VMADD_MII_IR(one,UChi_12,psi_22,psi_22)	      \
	VMADD_MII_IR(one,UChi_00,psi_30,psi_30)	      \
	VMADD_MII_IR(one,UChi_01,psi_31,psi_31)	      \
	VMADD_MII_IR(one,UChi_02,psi_32,psi_32)	      \
	);		     \
  }

//      fspin(2)+=hspin(1);
//      fspin(3)-=hspin(0);
#define YP_RECON_ACCUM {\
    asm(\
	VADD(psi_00,UChi_00,psi_00) 	VADD(psi_01,UChi_01,psi_01)	VADD(psi_02,UChi_02,psi_02) \
	VADD(psi_10,UChi_10,psi_10) 	VADD(psi_11,UChi_11,psi_11)	VADD(psi_12,UChi_12,psi_12) \
	VADD(psi_20,UChi_10,psi_20) 	VADD(psi_21,UChi_11,psi_21)	VADD(psi_22,UChi_12,psi_22) \
	VSUB(psi_30,UChi_00,psi_30) 	VSUB(psi_31,UChi_01,psi_31)	VSUB(psi_32,UChi_02,psi_32) \
	);\
 }
#define YM_RECON_ACCUM {\
    asm(\
	VADD(psi_00,UChi_00,psi_00) 	VADD(psi_01,UChi_01,psi_01)	VADD(psi_02,UChi_02,psi_02) \
	VADD(psi_10,UChi_10,psi_10) 	VADD(psi_11,UChi_11,psi_11)	VADD(psi_12,UChi_12,psi_12) \
	VSUB(psi_20,UChi_10,psi_20) 	VSUB(psi_21,UChi_11,psi_21)	VSUB(psi_22,UChi_12,psi_22) \
	VADD(psi_30,UChi_00,psi_30) 	VADD(psi_31,UChi_01,psi_31)	VADD(psi_32,UChi_02,psi_32) \
	);\
 }

//      fspin(2)-=timesI(hspin(0));
//      fspin(3)+=timesI(hspin(1));
#define ZP_RECON_ACCUM {\
    asm(\
	VONE(one)\
	VADD(psi_00,UChi_00,psi_00) 	VADD(psi_01,UChi_01,psi_01)	VADD(psi_02,UChi_02,psi_02) \
	VADD(psi_10,UChi_10,psi_10) 	VADD(psi_11,UChi_11,psi_11)	VADD(psi_12,UChi_12,psi_12) \
	VMADD_II_MIR(one,UChi_00,psi_20,psi_20)	      \
	VMADD_II_MIR(one,UChi_01,psi_21,psi_21)	      \
	VMADD_II_MIR(one,UChi_02,psi_22,psi_22)	      \
	VMADD_MII_IR(one,UChi_10,psi_30,psi_30)	      \
	VMADD_MII_IR(one,UChi_11,psi_31,psi_31)	      \
	VMADD_MII_IR(one,UChi_12,psi_32,psi_32)	      \
	);\
 }

#define ZM_RECON_ACCUM {\
    asm(\
	VONE(one)\
	VADD(psi_00,UChi_00,psi_00) 	VADD(psi_01,UChi_01,psi_01)	VADD(psi_02,UChi_02,psi_02) \
	VADD(psi_10,UChi_10,psi_10) 	VADD(psi_11,UChi_11,psi_11)	VADD(psi_12,UChi_12,psi_12) \
	VMADD_MII_IR(one,UChi_00,psi_20,psi_20)	      \
	VMADD_MII_IR(one,UChi_01,psi_21,psi_21)	      \
	VMADD_MII_IR(one,UChi_02,psi_22,psi_22)	      \
	VMADD_II_MIR(one,UChi_10,psi_30,psi_30)	      \
	VMADD_II_MIR(one,UChi_11,psi_31,psi_31)	      \
	VMADD_II_MIR(one,UChi_12,psi_32,psi_32)	      \
	);\
 }

//      fspin(2)+=hspin(0);
//      fspin(3)+=hspin(1);
#define TP_RECON_ACCUM {\
    asm(\
	VADD(psi_00,UChi_00,psi_00) 	VADD(psi_01,UChi_01,psi_01)	VADD(psi_02,UChi_02,psi_02) \
	VADD(psi_10,UChi_10,psi_10) 	VADD(psi_11,UChi_11,psi_11)	VADD(psi_12,UChi_12,psi_12) \
	VADD(psi_20,UChi_00,psi_20) 	VADD(psi_21,UChi_01,psi_21)	VADD(psi_22,UChi_02,psi_22) \
	VADD(psi_30,UChi_10,psi_30) 	VADD(psi_31,UChi_11,psi_31)	VADD(psi_32,UChi_12,psi_32) \
	);\
 }

#define TM_RECON_ACCUM {\
    asm(\
	VONE(one)\
	VADD(psi_00,UChi_00,psi_00) 	VADD(psi_01,UChi_01,psi_01)	VADD(psi_02,UChi_02,psi_02) \
	VADD(psi_10,UChi_10,psi_10) 	VADD(psi_11,UChi_11,psi_11)	VADD(psi_12,UChi_12,psi_12) \
	VSUB(psi_20,UChi_00,psi_20) 	VSUB(psi_21,UChi_01,psi_21)	VSUB(psi_22,UChi_02,psi_22) \
	VSUB(psi_30,UChi_10,psi_30) 	VSUB(psi_31,UChi_11,psi_31)	VSUB(psi_32,UChi_12,psi_32) \
	);\
 }

uint64_t GetPFInfo(int nent,int plocal);
uint64_t GetInfo(int ptype,int local,int perm,int Xp,int ent,int plocal); 

#define COMPLEX_TYPE int; 
int signs[4];

void testme(int osites,int ssU)
{
  int local,perm, ptype;
  uint64_t base;
  uint64_t basep;
  const uint64_t plocal =(uint64_t) & in._odata[0];

  //  vComplexF isigns[2] = { signs[0], signs[1] };
  //COMPLEX_TYPE is vComplexF of vComplexD depending 
  //on the chosen precision
  COMPLEX_TYPE *isigns = &signs[0];

  MASK_REGS;
  int nmax=osites;
  for(int site=0;site<Ns;site++) {
    int sU =ssU;
  int ssn=ssU+1; 
  if(ssn>=nmax) ssn=0;
  int sUn=ssn;
  for(int s=0;s<Ls;s++) {
  ss =sU*Ls+s;
  ssn=sUn*Ls+s; 
  ////////////////////////////////
  // Xp
  ////////////////////////////////
  int  ent=ss*8;// 2*Ndim
  int nent=ssn*8;

  PF_GAUGE(Xp); 
  base  = GetInfo(ptype,local,perm,Xp,ent,plocal); ent++;
  PREFETCH1_CHIMU(base);

  basep = GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);
#ifdef KERNEL_DAG
    XP_PROJMEM(base);
#else 
    XM_PROJMEM(base);
#endif
    MAYBEPERM(PERMUTE_DIR3,perm);
  } else { 
    LOAD_CHI(base);
  }
  base = GetInfo(ptype,local,perm,Yp,ent,plocal); ent++;
  PREFETCH_CHIMU(base);
  {
    MULT_2SPIN_DIR_PFXP(Xp,basep);
  }
  LOAD64(%r10,isigns);
#ifdef KERNEL_DAG
  XP_RECON;
#else
  XM_RECON;
#endif
  ////////////////////////////////
  // Yp
  ////////////////////////////////
  basep = GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
    YP_PROJMEM(base);
#else
    YM_PROJMEM(base);
#endif
    MAYBEPERM(PERMUTE_DIR2,perm);
  } else { 
    LOAD_CHI(base);
  }
  base  = GetInfo(ptype,local,perm,Zp,ent,plocal); ent++;
  PREFETCH_CHIMU(base);
  {
    MULT_2SPIN_DIR_PFYP(Yp,basep);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
  YP_RECON_ACCUM;
#else
  YM_RECON_ACCUM;
#endif

  ////////////////////////////////
  // Zp
  ////////////////////////////////
  basep = GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
    ZP_PROJMEM(base);
#else
    ZM_PROJMEM(base);
#endif
    MAYBEPERM(PERMUTE_DIR1,perm);
  } else { 
    LOAD_CHI(base);
  }
  base  = GetInfo(ptype,local,perm,Tp,ent,plocal); ent++;
  PREFETCH_CHIMU(base);
  {
    MULT_2SPIN_DIR_PFZP(Zp,basep);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
  ZP_RECON_ACCUM;
#else
  ZM_RECON_ACCUM;
#endif

  ////////////////////////////////
  // Tp
  ////////////////////////////////
  basep = GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
    TP_PROJMEM(base);
#else
    TM_PROJMEM(base);
#endif
    MAYBEPERM(PERMUTE_DIR0,perm);
  } else { 
    LOAD_CHI(base);
  }
  base = GetInfo(ptype,local,perm,Xm,ent,plocal); ent++;
  PREFETCH_CHIMU(base);
  {
    MULT_2SPIN_DIR_PFTP(Tp,basep);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
  TP_RECON_ACCUM;
#else
  TM_RECON_ACCUM;
#endif

  ////////////////////////////////
  // Xm
  ////////////////////////////////
#ifndef STREAM_STORE
  basep= (uint64_t) &out._odata[ss];
#endif
  //  basep= GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
    XM_PROJMEM(base);
#else
    XP_PROJMEM(base);
#endif
    MAYBEPERM(PERMUTE_DIR3,perm);
  } else { 
    LOAD_CHI(base);
  }
  base = GetInfo(ptype,local,perm,Ym,ent,plocal); ent++;
  PREFETCH_CHIMU(base);
  {
    MULT_2SPIN_DIR_PFXM(Xm,basep);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
  XM_RECON_ACCUM;
#else
  XP_RECON_ACCUM;
#endif

  ////////////////////////////////
  // Ym
  ////////////////////////////////
  basep= GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
    YM_PROJMEM(base);
#else
    YP_PROJMEM(base);
#endif
    MAYBEPERM(PERMUTE_DIR2,perm);
  } else { 
    LOAD_CHI(base);
  }
  base = GetInfo(ptype,local,perm,Zm,ent,plocal); ent++;
  PREFETCH_CHIMU(base);
  {
    MULT_2SPIN_DIR_PFYM(Ym,basep);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
  YM_RECON_ACCUM;
#else
  YP_RECON_ACCUM;
#endif

  ////////////////////////////////
  // Zm
  ////////////////////////////////
  basep= GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
    ZM_PROJMEM(base);
#else
    ZP_PROJMEM(base);
#endif
    MAYBEPERM(PERMUTE_DIR1,perm);
  } else { 
    LOAD_CHI(base);
  }
  base = GetInfo(ptype,local,perm,Tm,ent,plocal); ent++;
  PREFETCH_CHIMU(base);
  {
    MULT_2SPIN_DIR_PFZM(Zm,basep);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
  ZM_RECON_ACCUM;
#else
  ZP_RECON_ACCUM;
#endif

  ////////////////////////////////
  // Tm
  ////////////////////////////////
  basep= GetPFInfo(nent,plocal); nent++;
  if ( local ) {
    LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
    TM_PROJMEM(base);
#else
    TP_PROJMEM(base);
#endif
    MAYBEPERM(PERMUTE_DIR0,perm);
  } else { 
    LOAD_CHI(base);
  }
  base= (uint64_t) &out._odata[ss];
#ifndef STREAM_STORE
  PREFETCH_CHIMU(base);
#endif
  {
    MULT_2SPIN_DIR_PFTM(Tm,basep);
  }
  LOAD64(%r10,isigns);  // times i => shuffle and xor the real part sign bit
#ifdef KERNEL_DAG
  TM_RECON_ACCUM;
#else
  TP_RECON_ACCUM;
#endif

  basep= GetPFInfo(nent,plocal); nent++;
  SAVE_RESULT(base,basep);
  
  }
  ssU++;
  }
}


#endif
