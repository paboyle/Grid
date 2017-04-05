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

#include <stdint.h>

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
#define perm_reg 28

#define REP  %%r16
#define IMM  %%r17
#define pREP  %r16
#define pIMM  %r17

#define PPC_INST_DCBTLS 0x7c00014c
#define PPC_INST_DCBLC  0x7c00030c
#define __PPC_CT(t)     (((t) & 0x0f) << 21)
#define ___PPC_RA(a)    (((a) & 0x1f) << 16)
#define ___PPC_RB(b)    (((b) & 0x1f) << 11)

#define LOCK_SET   ".long (" HASH(PPC_INST_DCBTLS) "|"  HASH(___PPC_RB(16)) ")\n"
#define LOCK_CLEAR ".long (" HASH(PPC_INST_DCBLC) "|"  HASH(___PPC_RB(16)) ")\n"

/*Alias regs for incoming fourspinor on neighbour site*/
#define Chi_20 UChi_00
#define Chi_21 UChi_01
#define Chi_22 UChi_02
#define Chi_30 UChi_10
#define Chi_31 UChi_11
#define Chi_32 UChi_12

/*********************************************************
 * Architectural macros
 *********************************************************/
#define HASHit(A)  #A 
#define HASH(A)  HASHit(A)
#define LOAD64(A,ptr)  


#define MASK_REGS             /*NOOP ON BGQ*/
#define PF_GAUGE(A)           /*NOOP ON BGQ*/
#define PREFETCH1_CHIMU(base) /*NOOP ON BGQ*/
#define PREFETCH_CHIMU(base)  /*NOOP ON BGQ*/

#define VLOADf(OFF,PTR,DEST)         "qvlfsx      " #DEST  "," #PTR "," #OFF " ;\n"
#define VLOADuf(OFF,PTR,DEST)        "qvlfsux     " #DEST  "," #PTR "," #OFF " ;\n"
#define VSTOREf(OFF,PTR,SRC)         "qvstfsx     " #SRC  "," #PTR "," #OFF " ;\n"
#define VSTOREuf(OFF,PTR,SRC)        "qvstfsux    " #SRC  "," #PTR "," #OFF " ;\n"
#define VSPLATf(A,B,DEST)            "qvlfcsxa    " #DEST  "," #A "," #B ";\n"
#define VSIZEf (16)

#define VPERMIi(p)                 "qvgpci   " #p ", 1217;\n"
#define VPERMi(A,p)                "qvfperm  " #A "," #A "," #A "," #p ";\n"
#define VPERMI(p)                 VPERMIi(p)                 
#define VPERM(A,p)                VPERMi(A,p)                

#define VLOADd(OFF,PTR,DEST)         "qvlfdx      " #DEST  "," #PTR "," #OFF " ;\n"
#define VLOADud(OFF,PTR,DEST)        "qvlfdux     " #DEST  "," #PTR "," #OFF " ;\n"
#define VSTOREd(OFF,PTR,SRC)         "qvstfdx     " #SRC  "," #PTR "," #OFF " ;\n"
#define VSTOREud(OFF,PTR,SRC)        "qvstfdux    " #SRC  "," #PTR "," #OFF " ;\n"
#define VSPLATd(A,B,DEST)            "qvlfcdxa    " #DEST  "," #A "," #B ";\n"
#define VSIZEd (32)

// QPX manual ordering QRT comes first (dest)
#define VZEROi(DEST)                  "qvfset       " #DEST "; \n qvfsub " #DEST ","  #DEST ","  #DEST ";\n" 
#define VONEi(DEST)                   "qvfset       " #DEST "; \n" 
#define VMOVi(DEST,A)                 "qvfmr        " #DEST "," #A   ";\n"
#define VADDi(DEST,A,B)               "qvfadd       " #DEST "," #A "," #B  ";\n"
#define VSUBi(DEST,A,B)               "qvfsub       " #DEST "," #A "," #B  ";\n"
#define VMULi(DEST,A,B)               "qvfmul       " #DEST "," #A "," #B  ";\n"
#define VMUL_RR_RIi(DEST,A,B)         "qvfxmul      " #DEST "," #A "," #B  ";\n" 
#define VMADDi(DEST,A,B,C)            "qvfmadd      " #DEST "," #A "," #B ","#C ";\n"
#define VMADD_RR_RIi(DEST,A,B,C)      "qvfxmadd     " #DEST "," #A "," #B ","#C ";\n" 
#define VMADD_MII_IRi(DEST,A,B,C)     "qvfxxnpmadd  " #DEST "," #B "," #A ","#C ";\n" 
#define VMADD_II_MIRi(DEST,A,B,C)     "qvfxxcpnmadd " #DEST "," #B "," #A ","#C ";\n"  

#define VZERO(C)                  VZEROi(C)                  
#define VONE(C)                   VONEi(C)                   
#define VMOV(C,A)                 VMOVi(C,A)                 
#define VADD(A,B,C)               VADDi(A,B,C)               
#define VSUB(A,B,C)               VSUBi(A,B,C)               
#define VMUL(A,B,C)               VMULi(A,B,C)               
#define VMUL_RR_RI(A,B,C)         VMUL_RR_RIi(A,B,C)         
#define VMADD(A,B,C,D)            VMADDi(A,B,C,D)            
#define VMADD_RR_RI(A,B,C,D)      VMADD_RR_RIi(A,B,C,D)      
#define VMADD_MII_IR(A,B,C,D)     VMADD_MII_IRi(A,B,C,D)     
#define VMADD_II_MIR(A,B,C,D)     VMADD_II_MIRi(A,B,C,D)     

/*********************************************************
 * Macro sequences encoding QCD
 *********************************************************/
#define LOCK_GAUGE(dir)							\
  {									\
    uint64_t byte_addr = (uint64_t)&U._odata[sU];			\
    int count = (sizeof(U._odata[0])+63)/64;				\
    asm (" mtctr %0 \n"							\
	 " mr " HASH(REP) ", %1\n"					\
	 " li " HASH(IMM) ", 64\n"					\
	 "0:\n"							\
	 LOCK_SET							\
	 "  add " HASH(REP) "," HASH(IMM) "," HASH(REP) "\n"		\
	 "  bdnz 0b\n"						\
	 : : "b" (count), "b" (byte_addr) );					\
  }

#define UNLOCK_GAUGE(dir)						\
  {									\
    uint64_t byte_addr = (uint64_t)&U._odata[sU];			\
    int count = (sizeof(U._odata[0])+63)/64;				\
    asm (" mtctr %0 \n"							\
	 " mr " HASH(REP) ", %1\n"					\
	 " li " HASH(IMM) ", 64\n"					\
	 "0:\n"								\
	 LOCK_CLEAR							\
	 "  add " HASH(REP) "," HASH(IMM) "," HASH(REP) "\n"		\
	 "  bdnz 0b\n"						\
	 : : "b" (count), "b" (byte_addr) );					\
  }

#define ZERO_PSI				\
  VZERO(psi_00)					\
  VZERO(psi_01)					\
  VZERO(psi_02)					\
  VZERO(psi_10)					\
  VZERO(psi_11)					\
  VZERO(psi_12)					\
  VZERO(psi_20)					\
  VZERO(psi_21)					\
  VZERO(psi_22)					\
  VZERO(psi_30)					\
  VZERO(psi_31)					\
  VZERO(psi_32)

#define MULT_2SPIN_QPX_LSd(ptr,p) MULT_2SPIN_QPX_INTERNAL(ptr,p,VSPLAT,16) 
#define MULT_2SPIN_QPX_LSf(ptr,p) MULT_2SPIN_QPX_INTERNAL(ptr,p,VSPLAT,8) 
#define MULT_2SPIN_QPXd(ptr,p)    MULT_2SPIN_QPX_INTERNAL(ptr,p,VLOAD,32) 
#define MULT_2SPIN_QPXf(ptr,p)    MULT_2SPIN_QPX_INTERNAL(ptr,p,VLOAD,16) 

#define MULT_2SPIN_QPX_INTERNAL(ptr,p,ULOAD,USKIP) {			\
    uint64_t ub = ((uint64_t)ptr);				\
    asm (							\
         ULOAD(%0,%3,U0)					\
         ULOAD(%1,%3,U1)					\
         ULOAD(%2,%3,U2)					\
	 VMUL_RR_RI(UChi_00,U0,Chi_00)					\
	 VMUL_RR_RI(UChi_01,U1,Chi_00)					\
	 VMUL_RR_RI(UChi_02,U2,Chi_00)					\
	 VMUL_RR_RI(UChi_10,U0,Chi_10)					\
	 VMUL_RR_RI(UChi_11,U1,Chi_10)					\
	 VMUL_RR_RI(UChi_12,U2,Chi_10)					\
	 VMADD_MII_IR(UChi_00,U0,Chi_00,UChi_00)			\
	 VMADD_MII_IR(UChi_01,U1,Chi_00,UChi_01)			\
	 VMADD_MII_IR(UChi_02,U2,Chi_00,UChi_02)			\
	 VMADD_MII_IR(UChi_10,U0,Chi_10,UChi_10)			\
	 VMADD_MII_IR(UChi_11,U1,Chi_10,UChi_11)			\
	 VMADD_MII_IR(UChi_12,U2,Chi_10,UChi_12)			\
	 : : "b" (0), "b" (USKIP*3), "b" (USKIP*6), "b" (ub ));		\
    asm (								\
         ULOAD(%0,%3,U0)						\
         ULOAD(%1,%3,U1)						\
         ULOAD(%2,%3,U2)						\
	 VMADD_RR_RI(UChi_00,U0,Chi_01,UChi_00)				\
	 VMADD_RR_RI(UChi_01,U1,Chi_01,UChi_01)				\
	 VMADD_RR_RI(UChi_02,U2,Chi_01,UChi_02)				\
	 VMADD_RR_RI(UChi_10,U0,Chi_11,UChi_10)				\
	 VMADD_RR_RI(UChi_11,U1,Chi_11,UChi_11)				\
	 VMADD_RR_RI(UChi_12,U2,Chi_11,UChi_12)				\
	 VMADD_MII_IR(UChi_00,U0,Chi_01,UChi_00)			\
	 VMADD_MII_IR(UChi_01,U1,Chi_01,UChi_01)			\
	 VMADD_MII_IR(UChi_02,U2,Chi_01,UChi_02)			\
	 VMADD_MII_IR(UChi_10,U0,Chi_11,UChi_10)			\
	 VMADD_MII_IR(UChi_11,U1,Chi_11,UChi_11)			\
	 VMADD_MII_IR(UChi_12,U2,Chi_11,UChi_12)			\
	 : : "b" (USKIP*1), "b" (USKIP*4), "b" (USKIP*7), "b" (ub ));		\
    asm (								\
         ULOAD(%0,%3,U0)						\
         ULOAD(%1,%3,U1)						\
         ULOAD(%2,%3,U2)						\
	 VMADD_RR_RI(UChi_00,U0,Chi_02,UChi_00)				\
	 VMADD_RR_RI(UChi_01,U1,Chi_02,UChi_01)				\
	 VMADD_RR_RI(UChi_02,U2,Chi_02,UChi_02)				\
	 VMADD_RR_RI(UChi_10,U0,Chi_12,UChi_10)				\
	 VMADD_RR_RI(UChi_11,U1,Chi_12,UChi_11)				\
	 VMADD_RR_RI(UChi_12,U2,Chi_12,UChi_12)				\
	 VMADD_MII_IR(UChi_00,U0,Chi_02,UChi_00)			\
	 VMADD_MII_IR(UChi_01,U1,Chi_02,UChi_01)			\
	 VMADD_MII_IR(UChi_02,U2,Chi_02,UChi_02)			\
	 VMADD_MII_IR(UChi_10,U0,Chi_12,UChi_10)			\
	 VMADD_MII_IR(UChi_11,U1,Chi_12,UChi_11)			\
	 VMADD_MII_IR(UChi_12,U2,Chi_12,UChi_12)			\
	 : : "b" (USKIP*2), "b" (USKIP*5), "b" (USKIP*8), "b" (ub ));		\
  }


#define MULT_2SPIN_DIR_PF(A,p) MULT_2SPIN_PF(&U._odata[sU](A),p)
#define MULT_2SPIN_PF(ptr,pf) MULT_2SPIN(ptr,pf)

#define SAVE_RESULT(base,basep) {\
    uint64_t ub = ((uint64_t)base)  - (VSIZE);			\
    asm("mr " HASH(REP)  ", %0;\n"					\
	"li " HASH(IMM)      "," HASH(VSIZE)" ;\n"				\
	VSTOREu(IMM,REP,psi_00)						\
	VSTOREu(IMM,REP,psi_01)						\
	VSTOREu(IMM,REP,psi_02)						\
	VSTOREu(IMM,REP,psi_10)						\
	VSTOREu(IMM,REP,psi_11)						\
	VSTOREu(IMM,REP,psi_12)						\
	VSTOREu(IMM,REP,psi_20)						\
	VSTOREu(IMM,REP,psi_21)						\
	VSTOREu(IMM,REP,psi_22)						\
	VSTOREu(IMM,REP,psi_30)						\
	VSTOREu(IMM,REP,psi_31)						\
	VSTOREu(IMM,REP,psi_32)						\
	: : "b" (ub) : HASH(pIMM), HASH(pREP) );				\
  }


/*
 *Annoying BG/Q loads with no immediat indexing and big performance hit
 *when second miss to a L1 line occurs
 */
#define LOAD_CHI(base) {						\
    uint64_t ub = ((uint64_t)base)  - (2*VSIZE);			\
    asm("mr  " HASH(REP) ",%0 ;\n"					\
	"li  " HASH(IMM) ",(2*" HASH(VSIZE) ");\n"			\
	VLOADu(IMM,REP,Chi_00)						\
	VLOADu(IMM,REP,Chi_02)						\
	VLOADu(IMM,REP,Chi_11) : : "b" (ub)  : HASH(pIMM), HASH(pREP)  ); \
    ub = ((uint64_t)base)  - VSIZE;					\
    asm("mr  " HASH(REP) ", %0;\n"					\
	"li  " HASH(IMM) ",(2*" HASH(VSIZE) ");\n"			\
	VLOADu(IMM,REP,Chi_01)						\
	VLOADu(IMM,REP,Chi_10)						\
	VLOADu(IMM,REP,Chi_12)	: : "b" (ub)  : HASH(pIMM), HASH(pREP) );	\
  }

#define LOAD_CHIMU(base) {						\
    uint64_t ub = ((uint64_t)base)  - (2*VSIZE);			\
    asm("mr " HASH(REP) ",%0;\n"					\
	"li " HASH(IMM) ",(2*" HASH(VSIZE) ");\n"			\
	VLOADu(IMM,REP,Chi_00)						\
	VLOADu(IMM,REP,Chi_02)						\
	VLOADu(IMM,REP,Chi_11)						\
	VLOADu(IMM,REP,Chi_20)						\
	VLOADu(IMM,REP,Chi_22)						\
	VLOADu(IMM,REP,Chi_31) : : "b" (ub)  : HASH(pIMM), HASH(pREP) ); \
    ub = ((uint64_t)base)  - VSIZE;					\
    asm("mr " HASH(REP) ", %0;\n"					\
	"li " HASH(IMM) ", (2*" HASH(VSIZE) ");\n"			\
	VLOADu(IMM,REP,Chi_01)						\
	VLOADu(IMM,REP,Chi_10)						\
	VLOADu(IMM,REP,Chi_12)						\
	VLOADu(IMM,REP,Chi_21)						\
	VLOADu(IMM,REP,Chi_30)						\
	VLOADu(IMM,REP,Chi_32)	: : "b" (ub)  : HASH(pIMM), HASH(pREP) );	\
  }

//      hspin(0)=fspin(0)+timesI(fspin(3));
//      hspin(1)=fspin(1)+timesI(fspin(2));
#define XP_PROJMEM(base) {					\
    LOAD_CHIMU(base);						\
    asm (							\
         VONE(one)						\
	 VMADD_MII_IR(Chi_00,one,Chi_30,Chi_00)			\
	 VMADD_MII_IR(Chi_01,one,Chi_31,Chi_01)			\
	 VMADD_MII_IR(Chi_02,one,Chi_32,Chi_02)			\
	 VMADD_MII_IR(Chi_10,one,Chi_20,Chi_10)			\
	 VMADD_MII_IR(Chi_11,one,Chi_21,Chi_11)			\
	 VMADD_MII_IR(Chi_12,one,Chi_22,Chi_12)			\
							);	\
  }

#define XM_PROJMEM(base) {				\
    LOAD_CHIMU(base);					\
    asm (						\
         VONE(one)						\
	 VMADD_II_MIR(Chi_00,one,Chi_30,Chi_00)			\
	 VMADD_II_MIR(Chi_01,one,Chi_31,Chi_01)			\
	 VMADD_II_MIR(Chi_02,one,Chi_32,Chi_02)			\
	 VMADD_II_MIR(Chi_10,one,Chi_20,Chi_10)			\
	 VMADD_II_MIR(Chi_11,one,Chi_21,Chi_11)			\
	 VMADD_II_MIR(Chi_12,one,Chi_22,Chi_12)			\
							);	\
  }

//      hspin(0)=fspin(0)-fspin(3);
//      hspin(1)=fspin(1)+fspin(2);
#define YP_PROJMEM(base) {  \
    LOAD_CHIMU(base);						\
    asm (							\
	 VSUB(Chi_00,Chi_00,Chi_30)				\
	 VSUB(Chi_01,Chi_01,Chi_31)				\
	 VSUB(Chi_02,Chi_02,Chi_32)				\
	 VADD(Chi_10,Chi_10,Chi_20)				\
	 VADD(Chi_11,Chi_11,Chi_21)				\
	 VADD(Chi_12,Chi_12,Chi_22)				\
							);	\
  }

#define YM_PROJMEM(base) {			\
    LOAD_CHIMU(base);						\
    asm (							\
	 VADD(Chi_00,Chi_00,Chi_30)				\
	 VADD(Chi_01,Chi_01,Chi_31)				\
	 VADD(Chi_02,Chi_02,Chi_32)				\
	 VSUB(Chi_10,Chi_10,Chi_20)				\
	 VSUB(Chi_11,Chi_11,Chi_21)				\
	 VSUB(Chi_12,Chi_12,Chi_22)			);	\
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
	 VMADD_MII_IR(Chi_00,one,Chi_20,Chi_00)			\
	 VMADD_MII_IR(Chi_01,one,Chi_21,Chi_01)			\
	 VMADD_MII_IR(Chi_02,one,Chi_22,Chi_02)			\
	 VMADD_II_MIR(Chi_10,one,Chi_30,Chi_10)			\
	 VMADD_II_MIR(Chi_11,one,Chi_31,Chi_11)			\
	 VMADD_II_MIR(Chi_12,one,Chi_32,Chi_12)			\
							);	\
  }

#define ZM_PROJMEM(base) {  \
    LOAD_CHIMU(base);						\
    asm (							\
         VONE(one)						\
	 VMADD_II_MIR(Chi_00,one,Chi_20,Chi_00)			\
	 VMADD_II_MIR(Chi_01,one,Chi_21,Chi_01)			\
	 VMADD_II_MIR(Chi_02,one,Chi_22,Chi_02)			\
	 VMADD_MII_IR(Chi_10,one,Chi_30,Chi_10)			\
	 VMADD_MII_IR(Chi_11,one,Chi_31,Chi_11)			\
	 VMADD_MII_IR(Chi_12,one,Chi_32,Chi_12)			\
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
	 VADD(Chi_00,Chi_00,Chi_20)				\
	 VADD(Chi_01,Chi_01,Chi_21)				\
	 VADD(Chi_02,Chi_02,Chi_22)				\
	 VADD(Chi_10,Chi_10,Chi_30)				\
	 VADD(Chi_11,Chi_11,Chi_31)				\
	 VADD(Chi_12,Chi_12,Chi_32)				\
							);	\
  }

#define TM_PROJMEM(base) {  \
    LOAD_CHIMU(base);						\
    asm (							\
	 VSUB(Chi_00,Chi_00,Chi_20)				\
	 VSUB(Chi_01,Chi_01,Chi_21)				\
	 VSUB(Chi_02,Chi_02,Chi_22)				\
	 VSUB(Chi_10,Chi_10,Chi_30)				\
	 VSUB(Chi_11,Chi_11,Chi_31)				\
	 VSUB(Chi_12,Chi_12,Chi_32)				\
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
	VMADD_II_MIR(psi_20,one,UChi_10,psi_20)	      \
	VMADD_II_MIR(psi_21,one,UChi_11,psi_21)	      \
	VMADD_II_MIR(psi_22,one,UChi_12,psi_22)	      \
	VMADD_II_MIR(psi_30,one,UChi_00,psi_30)	      \
	VMADD_II_MIR(psi_31,one,UChi_01,psi_31)	      \
	VMADD_II_MIR(psi_32,one,UChi_02,psi_32)	      \
	);		     \
  }

#define XM_RECON {				\
    asm(\
	VONE(one)\
	VMOV(psi_00,UChi_00) 	VMOV(psi_01,UChi_01)	VMOV(psi_02,UChi_02)\
	VMOV(psi_10,UChi_10) 	VMOV(psi_11,UChi_11)	VMOV(psi_12,UChi_12)\
	VZERO(psi_20)	VZERO(psi_21)	VZERO(psi_22) \
	VZERO(psi_30) 	VZERO(psi_31)   VZERO(psi_32) \
	VMADD_MII_IR(psi_20,one,UChi_10,psi_20)	      \
	VMADD_MII_IR(psi_21,one,UChi_11,psi_21)	      \
	VMADD_MII_IR(psi_22,one,UChi_12,psi_22)	      \
	VMADD_MII_IR(psi_30,one,UChi_00,psi_30)	      \
	VMADD_MII_IR(psi_31,one,UChi_01,psi_31)	      \
	VMADD_MII_IR(psi_32,one,UChi_02,psi_32)	      \
	);		     \
  }

#define XP_RECON_ACCUM {				\
    asm(\
	VONE(one)\
	VADD(psi_00,psi_00,UChi_00) 	VADD(psi_01,psi_01,UChi_01)	VADD(psi_02,psi_02,UChi_02) \
	VADD(psi_10,psi_10,UChi_10) 	VADD(psi_11,psi_11,UChi_11)	VADD(psi_12,psi_12,UChi_12) \
	VMADD_II_MIR(psi_20,one,UChi_10,psi_20)	      \
	VMADD_II_MIR(psi_21,one,UChi_11,psi_21)	      \
	VMADD_II_MIR(psi_22,one,UChi_12,psi_22)	      \
	VMADD_II_MIR(psi_30,one,UChi_00,psi_30)	      \
	VMADD_II_MIR(psi_31,one,UChi_01,psi_31)	      \
	VMADD_II_MIR(psi_32,one,UChi_02,psi_32)	      \
	);		     \
  }

#define XM_RECON_ACCUM {				\
    asm(\
	VONE(one)\
	VADD(psi_00,psi_00,UChi_00) 	VADD(psi_01,psi_01,UChi_01)	VADD(psi_02,psi_02,UChi_02) \
	VADD(psi_10,psi_10,UChi_10) 	VADD(psi_11,psi_11,UChi_11)	VADD(psi_12,psi_12,UChi_12) \
	VMADD_MII_IR(psi_20,one,UChi_10,psi_20)	      \
	VMADD_MII_IR(psi_21,one,UChi_11,psi_21)	      \
	VMADD_MII_IR(psi_22,one,UChi_12,psi_22)	      \
	VMADD_MII_IR(psi_30,one,UChi_00,psi_30)	      \
	VMADD_MII_IR(psi_31,one,UChi_01,psi_31)	      \
	VMADD_MII_IR(psi_32,one,UChi_02,psi_32)	      \
	);		     \
  }

//      fspin(2)+=hspin(1);
//      fspin(3)-=hspin(0);
#define YP_RECON_ACCUM {\
    asm(\
	VADD(psi_00,psi_00,UChi_00) 	VADD(psi_01,psi_01,UChi_01)	VADD(psi_02,psi_02,UChi_02) \
	VADD(psi_10,psi_10,UChi_10) 	VADD(psi_11,psi_11,UChi_11)	VADD(psi_12,psi_12,UChi_12) \
	VADD(psi_20,psi_20,UChi_10) 	VADD(psi_21,psi_21,UChi_11)	VADD(psi_22,psi_22,UChi_12) \
	VSUB(psi_30,psi_30,UChi_00) 	VSUB(psi_31,psi_31,UChi_01)	VSUB(psi_32,psi_32,UChi_02) \
	);\
 }
#define YM_RECON_ACCUM {\
    asm(\
	VADD(psi_00,psi_00,UChi_00) 	VADD(psi_01,psi_01,UChi_01)	VADD(psi_02,psi_02,UChi_02) \
	VADD(psi_10,psi_10,UChi_10) 	VADD(psi_11,psi_11,UChi_11)	VADD(psi_12,psi_12,UChi_12) \
	VSUB(psi_20,psi_20,UChi_10) 	VSUB(psi_21,psi_21,UChi_11)	VSUB(psi_22,psi_22,UChi_12) \
	VADD(psi_30,psi_30,UChi_00) 	VADD(psi_31,psi_31,UChi_01)	VADD(psi_32,psi_32,UChi_02) \
	);\
 }

//      fspin(2)-=timesI(hspin(0));
//      fspin(3)+=timesI(hspin(1));
#define ZP_RECON_ACCUM {\
    asm(\
	VONE(one)\
	VADD(psi_00,psi_00,UChi_00) 	VADD(psi_01,psi_01,UChi_01)	VADD(psi_02,psi_02,UChi_02) \
	VADD(psi_10,psi_10,UChi_10) 	VADD(psi_11,psi_11,UChi_11)	VADD(psi_12,psi_12,UChi_12) \
	VMADD_II_MIR(psi_20,one,UChi_00,psi_20)				\
	VMADD_II_MIR(psi_21,one,UChi_01,psi_21)				\
	VMADD_II_MIR(psi_22,one,UChi_02,psi_22)				\
	VMADD_MII_IR(psi_30,one,UChi_10,psi_30)				\
	VMADD_MII_IR(psi_31,one,UChi_11,psi_31)				\
	VMADD_MII_IR(psi_32,one,UChi_12,psi_32)				\
	);\
 }

#define ZM_RECON_ACCUM {\
    asm(\
	VONE(one)\
	VADD(psi_00,psi_00,UChi_00) 	VADD(psi_01,psi_01,UChi_01)	VADD(psi_02,psi_02,UChi_02) \
	VADD(psi_10,psi_10,UChi_10) 	VADD(psi_11,psi_11,UChi_11)	VADD(psi_12,psi_12,UChi_12) \
	VMADD_MII_IR(psi_20,one,UChi_00,psi_20)				\
	VMADD_MII_IR(psi_21,one,UChi_01,psi_21)				\
	VMADD_MII_IR(psi_22,one,UChi_02,psi_22)				\
	VMADD_II_MIR(psi_30,one,UChi_10,psi_30)				\
	VMADD_II_MIR(psi_31,one,UChi_11,psi_31)				\
	VMADD_II_MIR(psi_32,one,UChi_12,psi_32)				\
	);\
 }

//      fspin(2)+=hspin(0);
//      fspin(3)+=hspin(1);
#define TP_RECON_ACCUM {\
    asm(\
	VADD(psi_00,psi_00,UChi_00) 	VADD(psi_01,psi_01,UChi_01)	VADD(psi_02,psi_02,UChi_02) \
	VADD(psi_10,psi_10,UChi_10) 	VADD(psi_11,psi_11,UChi_11)	VADD(psi_12,psi_12,UChi_12) \
	VADD(psi_20,psi_20,UChi_00) 	VADD(psi_21,psi_21,UChi_01)	VADD(psi_22,psi_22,UChi_02) \
	VADD(psi_30,psi_30,UChi_10) 	VADD(psi_31,psi_31,UChi_11)	VADD(psi_32,psi_32,UChi_12) \
	);\
 }

#define TM_RECON_ACCUM {\
    asm(\
	VADD(psi_00,psi_00,UChi_00) 	VADD(psi_01,psi_01,UChi_01)	VADD(psi_02,psi_02,UChi_02) \
	VADD(psi_10,psi_10,UChi_10) 	VADD(psi_11,psi_11,UChi_11)	VADD(psi_12,psi_12,UChi_12) \
	VSUB(psi_20,psi_20,UChi_00) 	VSUB(psi_21,psi_21,UChi_01)	VSUB(psi_22,psi_22,UChi_02) \
	VSUB(psi_30,psi_30,UChi_10) 	VSUB(psi_31,psi_31,UChi_11)	VSUB(psi_32,psi_32,UChi_12) \
	);\
 }


#define ADD_RESULTi(PTR,pf)						\
  LOAD_CHIMU(PTR)							\
  asm(									\
  VADD(psi_00,chi_00,psi_00)  VADD(psi_01,chi_01,psi_01)  VADD(psi_02,chi_02,psi_02) \
  VADD(psi_10,chi_10,psi_10)  VADD(psi_11,chi_11,psi_11)  VADD(psi_12,chi_12,psi_12) \
  VADD(psi_20,chi_20,psi_20)  VADD(psi_21,chi_21,psi_21)  VADD(psi_22,chi_22,psi_22) \
  VADD(psi_30,chi_30,psi_30)  VADD(psi_31,chi_31,psi_31)  VADD(psi_32,chi_32,psi_32) ); \
  SAVE_RESULT(PTR,pf);


#define PERMUTE_DIR3
#define PERMUTE_DIR2
#define PERMUTE_DIR1

#define PERMUTE_DIR0 {							\
    asm(								\
	VPERMI(perm_reg)							\
	VPERM(Chi_00,perm_reg)   VPERM(Chi_01,perm_reg)   VPERM(Chi_02,perm_reg)	\
	VPERM(Chi_10,perm_reg)   VPERM(Chi_11,perm_reg)   VPERM(Chi_12,perm_reg) );	\
  }

#endif
