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
#ifndef GRID_ASM_AV512_H
#define GRID_ASM_AV512_H

// Serialisation elimination: 
// i)  ZEND -> ZEND1, ZEND2, 6 fold round robin.
// ii) TimesI -> TimesI_1, TimesI_2, 6 fold round robin
//

//////////////////////////////////////////////////////////////////////////////////////////
// Register allocations for Wilson Kernel are precision and IMCI/AVX512 indept
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

//////////////////////////////////////////////////////////////////////////////////////////
// CONFIG IMCI/AVX512
//////////////////////////////////////////////////////////////////////////////////////////

#define ASM_IMCI
#undef  ASM_AVX512

////////////////////////////////////////////////////////////////////////////////////////////////////
// Opcodes common to AVX512 and IMCI
////////////////////////////////////////////////////////////////////////////////////////////////////
#define MASK_REGS \
  __asm__ ("mov     $0xAAAA, %%eax \n"\
	   "kmov    %%eax, %%k6 \n"\
	   "knot     %%k6, %%k7 \n" : : : "%eax");

#define VZEROf(A)       "vpxorq " #A ","  #A "," #A ";\n"
#define VZEROd(A)       "vpxorq " #A ","  #A "," #A ";\n"

#define VTIMESIf(A,DEST, Z) \
  VTIMESI0f(A,DEST, Z) \
  VTIMESI1f(A,DEST, Z) \
  VTIMESI2f(A,DEST, Z) 

#define VTIMESId(A,DEST, Z) \
  VTIMESI0d(A,DEST, Z) \
  VTIMESI1d(A,DEST, Z) \
  VTIMESI2d(A,DEST, Z) 

#define VTIMESMINUSIf(A,DEST, Z) \
        VTIMESMINUSI0f(A,DEST, Z) \
        VTIMESMINUSI1f(A,DEST, Z) \
        VTIMESMINUSI2f(A,DEST, Z) 

#define VTIMESMINUSId(A,DEST, Z) \
        VTIMESMINUSI0d(A,DEST, Z) \
        VTIMESMINUSI1d(A,DEST, Z) \
        VTIMESMINUSI2d(A,DEST, Z) 

#define VACCTIMESIf(A,ACC,tmp)			\
 VACCTIMESI0f(A,ACC,tmp)			\
 VACCTIMESI1f(A,ACC,tmp)			\
 VACCTIMESI2f(A,ACC,tmp)			

#define  VACCTIMESI1MEMf(A,ACC,O,P)  "vaddps  " #O"*64("#P"),"#A "," #ACC"{%k7}" ";\n"
#define  VACCTIMESI2MEMf(A,ACC,O,P)  "vsubrps  " #O"*64("#P"),"#A "," #ACC"{%k6}" ";\n"
#define  VACCTIMESMINUSI1MEMf(A,ACC,O,P)  "vsubrps  " #O"*64("#P"),"#A "," #ACC"{%k7}" ";\n"
#define  VACCTIMESMINUSI2MEMf(A,ACC,O,P)  "vaddps  " #O"*64("#P"),"#A "," #ACC"{%k6}" ";\n"

#define VACCTIMESId(A,ACC,tmp)			\
 VACCTIMESI0d(A,ACC,tmp)			\
 VACCTIMESI1d(A,ACC,tmp)			\
 VACCTIMESI2d(A,ACC,tmp)			

#define VACCTIMESMINUSIf(A,ACC,tmp)			\
  VACCTIMESMINUSI0f(A,ACC,tmp)				\
  VACCTIMESMINUSI1f(A,ACC,tmp)				\
  VACCTIMESMINUSI2f(A,ACC,tmp)			

#define VACCTIMESMINUSId(A,ACC,tmp)			\
  VACCTIMESMINUSI0d(A,ACC,tmp)				\
  VACCTIMESMINUSI1d(A,ACC,tmp)				\
  VACCTIMESMINUSI2d(A,ACC,tmp)			

#define LOAD64i(A,ptr)  __asm__ ( "movq %0, %" #A :  : "r"(ptr)  : #A  );
#define LOAD64(A,ptr)          LOAD64i(A,ptr)

#define VMOVf(A,DEST)   "vmovaps  " #A ", " #DEST  ";\n"
#define VMOVd(A,DEST)   "vmovapd  " #A ", " #DEST  ";\n"

// Field prefetch
#define VPREFETCHNTA(O,A) "vprefetchnta "#O"*64("#A");\n" "vprefetch1 ("#O"+12)*64("#A");\n"
#define VPREFETCH(O,A)    "vprefetch0 "#O"*64("#A");\n" "vprefetch1 ("#O"+12)*64("#A");\n"
#define VPREFETCHG(O,A) 
#define VPREFETCHW(O,A) 
//"vprefetche0 "#O"*64("#A");\n" "vprefetche1 ("#O"+12)*64("#A");\n"
#define VEVICT(O,A)   
//  "clevict0 "#O"*64("#A");\n" 

#define VLOADf(OFF,PTR,DEST)   "vmovaps  " #OFF "*64(" #PTR "), " #DEST  ";\n"
#define VLOADd(OFF,PTR,DEST)   "vmovapd  " #OFF "*64(" #PTR "), " #DEST  ";\n"

#define VADDf(A,B,DEST)        "vaddps   " #A "," #B "," #DEST  ";\n"
#define VADDd(A,B,DEST)        "vaddpd   " #A "," #B "," #DEST  ";\n"

#define VSUBf(A,B,DEST)        "vsubps   " #A "," #B "," #DEST  ";\n"
#define VSUBd(A,B,DEST)        "vsubpd   " #A "," #B "," #DEST  ";\n"

#define VADDMEMf(O,A,B,DEST)        "vaddps   "#O"*64("#A ")," #B "," #DEST  ";\n"
#define VADDMEMd(O,A,B,DEST)        "vaddpd   "#O"*64("#A ")," #B "," #DEST  ";\n"

#define VSUBMEMf(O,A,B,DEST)        "vsubps   "#O"*64("#A ")," #B "," #DEST  ";\n"
#define VSUBMEMd(O,A,B,DEST)        "vsubpd   "#O"*64("#A ")," #B "," #DEST  ";\n"

#define VMULf(A,B,DEST)        "vmulps   " #A "," #B "," #DEST  ";\n"
#define VMULd(A,B,DEST)        "vmulpd   " #A "," #B "," #DEST  ";\n"

#define VMADDf(A,B,DEST)       "vfmadd231ps   " #A "," #B "," #DEST  ";\n"
#define VMADDd(A,B,DEST)       "vfmadd231pd   " #A "," #B "," #DEST  ";\n"

#define VMULMEMf(O,A,B,DEST)   "vmulps   " #O"*64("#A ")," #B "," #DEST  ";\n"
#define VMULMEMd(O,A,B,DEST)   "vmulpd   " #O"*64("#A ")," #B "," #DEST  ";\n"

#define VMADDMEMf(O,A,B,DEST)       "vfmadd231ps   " #O"*64("#A "),"#B "," #DEST  ";\n"
#define VMADDMEMd(O,A,B,DEST)       "vfmadd231pd   " #O"*64("#A "),"#B "," #DEST  ";\n"

#define ZLOADf(OFF,PTR,ri,ir)  VLOADf(OFF,PTR,ir)  VSHUFf(ir,ri)
#define ZLOADd(OFF,PTR,ri,ir)  VLOADd(OFF,PTR,ir)  VSHUFd(ir,ri)

#define ZMULf(Ari,Air,B,Criir,Ciirr)  VMULf(Ari,B,Criir)  VMULf(Air,B,Ciirr)
#define ZMULd(Ari,Air,B,Criir,Ciirr)  VMULd(Ari,B,Criir)  VMULd(Air,B,Ciirr)

#define ZMADDf(Ari,Air,B,Criir,Ciirr) VMADDf(Ari,B,Criir) VMADDf(Air,B,Ciirr)
#define ZMADDd(Ari,Air,B,Criir,Ciirr) VMADDd(Ari,B,Criir) VMADDd(Air,B,Ciirr)

#define ZENDf(Criir,Ciirr, tmp) ZEND1f(Criir,Ciirr, tmp) ZEND2f(Criir,Ciirr, tmp)
#define ZENDd(Criir,Ciirr, tmp) ZEND1d(Criir,Ciirr, tmp) ZEND2d(Criir,Ciirr, tmp)

// Need VSHUFMULMEMf,d for KNC
// AVX512 friendly
#define ZMULMEM2SPf(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr)\
  VSHUFMEMf(O,P,tmp) \
  VMULMEMf(O,P,B,Biirr) \
  VMULMEMf(O,P,C,Ciirr) \
  VMULf(tmp,B,Briir) \
  VMULf(tmp,C,Criir)

#define ZMULMEM2SPd(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr)\
  VSHUFMEMd(O,P,tmp)  \
  VMULMEMd(O,P,B,Biirr)  \ 
  VMULMEMd(O,P,C,Ciirr)  \
  VMULd(tmp,B,Briir)  \
  VMULd(tmp,C,Criir) 

#define ZMADDMEM2SPf(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr)\
  VSHUFMEMf(O,P,tmp) \
  VMADDMEMf(O,P,B,Biirr) \
  VMADDMEMf(O,P,C,Ciirr) \
  VMADDf(tmp,B,Briir) \
  VMADDf(tmp,C,Criir)

#define TRAP " int3 ;\n"

#define ZMADDMEM2SPd(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr)	\
  VSHUFMEMd(O,P,tmp) \
  VMADDMEMd(O,P,B,Biirr) \
  VMADDMEMd(O,P,C,Ciirr) \
  VMADDd(tmp,B,Briir) \
  VMADDd(tmp,C,Criir)

////////////////////////////////////////////////////////////////////////////////////////////////////
// Lane swizzling changed between AVX512 and IMCI and requires arch dependent complex support
////////////////////////////////////////////////////////////////////////////////////////////////////

// AVX512 special (Knights Landing)
#ifdef ASM_AVX512
#define VSTOREf(OFF,PTR,SRC)   "vmovntps " #SRC "," #OFF "*64(" #PTR ")"  ";\n"
#define VSTOREd(OFF,PTR,SRC)   "vmovntpd " #SRC "," #OFF "*64(" #PTR ")"  ";\n"
// Swaps Re/Im
#define VSHUFd(A,DEST)         "vshufpd  $0x5, " #A "," #A "," #DEST  ";\n"
#define VSHUFf(A,DEST)         "vshufps  $0x55," #A "," #A "," #DEST  ";\n"
// Memops are useful for optimisation
#define VSHUFMEMd(OFF,A,DEST)       "vpshufpd  $0x4e, " #OFF"("#A ")," #DEST  ";\n"
#define VSHUFMEMf(OFF,A,DEST)       "vpshufps  $0xb1, " #OFF"("#A ")," #DEST  ";\n"


// Merges accumulation for complex dot chain
// TODO:  12 operation saving:
//	# could SWIZ op 18{cdab} and eliminate temporary // 12cycles
//	# no use KNL though. Fingour something else there.
//	# All swizzles become perms ops, but gain addsub; subadd must use this
//	# uint32_t (0x7F << 23 )
//	# uint64_t (0x3FF<< 52 ) ; vpbroadcast
#define ZEND1f(Criir,Ciirr, tmp)	\
  "vshufps $0xb1," #Ciirr "," #Criir "," #tmp   ";\n"\
  "vaddps  " #Criir "," #tmp "," #Criir"{%k6}"  ";\n"

#define ZEND2f(Criir,Ciirr, tmp) "vsubps  " #Ciirr "," #tmp "," #Criir"{%k7}"  ";\n"

#define ZEND2d(Criir,Ciirr, tmp)	\
  "vshufpd $0x33," #Ciirr "," #Criir "," #tmp  ";\n"\
  "vaddpd  " #Criir "," #tmp "," #Criir"{%k6}" ";\n"
#define ZEND2d(Criir,Ciirr, tmp) "vsubpd  " #Ciirr "," #tmp "," #Criir"{%k7}" ";\n"

// Further opt possible: KNC -- use swizzle operand ; no addsub.
//                       KNL -- addsub. Saves 6 ops, 12 cycles; KNL cost of loading "1" as only fmaddsub
//                              no swizzle on KNL.
#define VTIMESI0f(A,DEST, Z)   VSHUFf(A,DEST)	 
#define VTIMESI1f(A,DEST, Z)   "vaddps  " #DEST "," #Z "," #DEST"{%k7}"  ";\n"
#define VTIMESI2f(A,DEST, Z)   "vsubps  " #DEST "," #Z "," #DEST"{%k6}"  ";\n"

#define VTIMESI0d(A,DEST, Z)   VSHUFd(A,DEST)	 
#define VTIMESI1d(A,DEST, Z)   "vaddpd  " #DEST "," #Z "," #DEST"{%k7}"  ";\n"
#define VTIMESI2d(A,DEST, Z)   "vsubpd  " #DEST "," #Z "," #DEST"{%k6}"  ";\n"

#define VTIMESMINUSI0f(A,DEST,Z)  VSHUFf(A,DEST)					
#define VTIMESMINUSI1f(A,DEST,Z)  "vsubps  " #DEST "," #Z "," #DEST"{%k7}"  ";\n"
#define VTIMESMINUSI2f(A,DEST,Z)  "vaddps  " #DEST "," #Z "," #DEST"{%k6}"  ";\n"

#define VTIMESMINUSI0d(A,DEST,Z)  VSHUFd(A,DEST)					
#define VTIMESMINUSI1d(A,DEST,Z)  "vsubpd  " #DEST "," #Z "," #DEST"{%k7}"  ";\n"
#define VTIMESMINUSI2d(A,DEST,Z)  "vaddpd  " #DEST "," #Z "," #DEST"{%k6}"  ";\n"

#define VACCTIMESMINUSI0f(A,ACC,tmp)  VSHUFf(A,tmp)					
#define VACCTIMESMINUSI1f(A,ACC,tmp)  "vsubps  " #tmp "," #ACC "," #ACC"{%k7}" ";\n"
#define VACCTIMESMINUSI2f(A,ACC,tmp)  "vaddps  " #tmp "," #ACC "," #ACC"{%k6}" ";\n"

#define VACCTIMESMINUSI0d(A,ACC,tmp)  VSHUFd(A,tmp)					
#define VACCTIMESMINUSI1d(A,ACC,tmp)  "vsubpd  " #tmp "," #ACC "," #ACC"{%k7}" ";\n"
#define VACCTIMESMINUSI2d(A,ACC,tmp)  "vaddpd  " #tmp "," #ACC "," #ACC"{%k6}" ";\n"

#define  VACCTIMESI0f(A,ACC,tmp)	VSHUFf(A,tmp)	
#define  VACCTIMESI1f(A,ACC,tmp)  "vaddps  " #tmp "," #ACC "," #ACC"{%k7}" ";\n"
#define  VACCTIMESI2f(A,ACC,tmp)  "vsubps  " #tmp "," #ACC "," #ACC"{%k6}" ";\n"

#define  VACCTIMESI0d(A,ACC,tmp)	VSHUFd(A,tmp)	
#define  VACCTIMESI1d(A,ACC,tmp)  "vaddpd  " #tmp "," #ACC "," #ACC"{%k7}" ";\n"
#define  VACCTIMESI2d(A,ACC,tmp)  "vsubpd  " #tmp "," #ACC "," #ACC"{%k6}" ";\n"

#define VPERM0f(A,B) "vshuff32x4 " #A "," #B "," "#B" ", " #_MM_SELECT_FOUR_FOUR(1,0,3,2) ";\n"
#define VPERM1f(A,B) "vshuff32x4 " #A "," #B "," "#B" ", " #_MM_SELECT_FOUR_FOUR(2,3,0,1) ";\n"
#define VPERM2f(A,B) "vshufps    " #A "," #B "," "#B" ", " #_MM_SELECT_FOUR_FOUR(1,0,3,2) ";\n"
#define VPERM3f(A,B) "vshufps    " #A "," #B "," "#B" ", " #_MM_SELECT_FOUR_FOUR(2,3,0,1) ";\n"

#define VPERM0d(A,B) "vshuff64x2 " #A "," #B "," "#B" ", " #_MM_SELECT_FOUR_FOUR(1,0,3,2) ";\n"
#define VPERM1d(A,B) "vshuff64x2 " #A "," #B "," "#B" ", " #_MM_SELECT_FOUR_FOUR(2,3,0,1) ";\n"
#define VPERM2d(A,B) "vshufpd    " #A "," #B "," "#B" ", " 0x55 ";\n"
#define VPERM3d(A,B) VMOVd(A,B)

#endif

// Knights Corner specials
#ifdef ASM_IMCI
#define VSTOREf(OFF,PTR,SRC)   "vmovnrngoaps " #SRC "," #OFF "*64(" #PTR ")"  ";\n"
#define VSTOREd(OFF,PTR,SRC)   "vmovnrngoapd " #SRC "," #OFF "*64(" #PTR ")"  ";\n"
  //#define VSTOREf(OFF,PTR,SRC)   "vmovaps " #SRC "," #OFF "*64(" #PTR ")"  ";\n"
  //#define VSTOREd(OFF,PTR,SRC)   "vmovapd " #SRC "," #OFF "*64(" #PTR ")"  ";\n"

#define VSHUFf(A,DEST)         "vmovaps " #A "{cdab} , " #DEST  ";\n"
#define VSHUFd(A,DEST)         "vmovapd " #A "{cdab} , " #DEST  ";\n"

// Memops are useful for optimisation
#define VSHUFMEMd(OFF,A,DEST)  "vpshufd  $0x4e, " #OFF"*64("#A ")," #DEST  ";\n"
#define VSHUFMEMf(OFF,A,DEST)  "vpshufd  $0xb1, " #OFF"*64("#A ")," #DEST  ";\n"

#define ZEND1d(Criir,Ciirr, tmp) "vaddpd  " #Criir "{cdab} ," #Criir "," #Criir"{%k6}"  ";\n"
#define ZEND2d(Criir,Ciirr, tmp) "vsubpd  " #Ciirr "{cdab} ," #Ciirr "," #Criir"{%k7}"  ";\n"

#define ZEND1f(Criir,Ciirr, tmp) "vaddps  " #Criir "{cdab} ," #Criir "," #Criir"{%k6}"  ";\n"
#define ZEND2f(Criir,Ciirr, tmp) "vsubps  " #Ciirr "{cdab} ," #Ciirr "," #Criir"{%k7}"  ";\n"

// Further opt possible: KNC -- use swizzle operand ; no addsub.
//                       KNL -- addsub. Saves 6 ops, 12 cycles; KNL cost of loading "1" as only fmaddsub
//                              no swizzle on KNL.
#define VTIMESI0f(A,DEST, Z)   
#define VTIMESI1f(A,DEST, Z)   "vaddps  " #A "{cdab}," #Z "," #DEST"{%k7}"  ";\n"
#define VTIMESI2f(A,DEST, Z)   "vsubps  " #A "{cdab}," #Z "," #DEST"{%k6}"  ";\n"

#define VTIMESI0d(A,DEST, Z)   
#define VTIMESI1d(A,DEST, Z)   "vaddpd  " #A "{cdab}," #Z "," #DEST"{%k7}"  ";\n"
#define VTIMESI2d(A,DEST, Z)   "vsubpd  " #A "{cdab}," #Z "," #DEST"{%k6}"  ";\n"

#define VTIMESMINUSI0f(A,DEST,Z)  
#define VTIMESMINUSI1f(A,DEST,Z)  "vsubps  " #A "{cdab}," #Z "," #DEST"{%k7}"  ";\n"
#define VTIMESMINUSI2f(A,DEST,Z)  "vaddps  " #A "{cdab}," #Z "," #DEST"{%k6}"  ";\n"

#define VTIMESMINUSI0d(A,DEST,Z)  
#define VTIMESMINUSI1d(A,DEST,Z)  "vsubpd  " #A "{cdab}," #Z "," #DEST"{%k7}"  ";\n"
#define VTIMESMINUSI2d(A,DEST,Z)  "vaddpd  " #A "{cdab}," #Z "," #DEST"{%k6}"  ";\n"

#define  VACCTIMESI0f(A,ACC,tmp)
#define  VACCTIMESI1f(A,ACC,tmp)  "vaddps  " #A "{cdab}," #ACC "," #ACC"{%k7}" ";\n"
#define  VACCTIMESI2f(A,ACC,tmp)  "vsubps  " #A "{cdab}," #ACC "," #ACC"{%k6}" ";\n"

#define  VACCTIMESI0d(A,ACC,tmp)
#define  VACCTIMESI1d(A,ACC,tmp)  "vaddpd  " #A "{cdab}," #ACC "," #ACC"{%k7}" ";\n"
#define  VACCTIMESI2d(A,ACC,tmp)  "vsubpd  " #A "{cdab}," #ACC "," #ACC"{%k6}" ";\n"

#define VACCTIMESMINUSI0f(A,ACC,tmp)  
#define VACCTIMESMINUSI1f(A,ACC,tmp)  "vsubps  " #A "{cdab}," #ACC "," #ACC"{%k7}" ";\n"
#define VACCTIMESMINUSI2f(A,ACC,tmp)  "vaddps  " #A "{cdab}," #ACC "," #ACC"{%k6}" ";\n"

#define VACCTIMESMINUSI0d(A,ACC,tmp)  
#define VACCTIMESMINUSI1d(A,ACC,tmp)  "vsubpd  " #A "{cdab}," #ACC "," #ACC"{%k7}" ";\n"
#define VACCTIMESMINUSI2d(A,ACC,tmp)  "vaddpd  " #A "{cdab}," #ACC "," #ACC"{%k6}" ";\n"

//#define ZENDf(Criir,Ciirr, tmp)	

//((1<<6)|(0<<4)|(3<<2)|(2)) == 0100,1110 = 0x4e
//((2<<6)|(3<<4)|(0<<2)|(1)) == 1011,0001 = 0xb1

#define VPERM0f(A,B) "vpermf32x4  $0x4e," #A "," #B ";\n"
#define VPERM1f(A,B) "vpermf32x4  $0xb1," #A "," #B ";\n"
#define VPERM2f(A,B) "vmovaps     " #A "{badc}," #B ";\n"
#define VPERM3f(A,B) "vmovaps     " #A "{cdab}," #B ";\n"

#define VPERM0d(A,B) "vpermf32x4  $0x4e," #A "," #B ";\n"
#define VPERM1d(A,B) "vmovapd     " #A "{badc}," #B ";\n"
#define VPERM2d(A,B) "vmovapd     " #A "{cdab}," #B ";\n"
#define VPERM3d(A,B) VMOVd(A,B)

#endif

                    
//  const SiteSpinor * ptr = & in._odata[offset];	
#define LOAD_CHIMU(PTR)	 LOAD_CHIMUi(PTR)
#define LOAD_CHI(PTR)	 LOAD_CHIi(PTR) 
#define SAVE_UCHI(PTR)	 SAVE_UCHIi(PTR)
#define SAVE_CHI(PTR)	 SAVE_CHIi(PTR)
#define SAVE_RESULT(PTR) SAVE_RESULTi(PTR)

#define LOAD_CHIMUi(PTR) \
	   LOAD64(%r8,PTR)	\
	   __asm__ (\
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

#define LOAD_CHIi(PTR)				\
  LOAD64(%r8,PTR)					\
  __asm__ (						\
  VLOAD(0,%r8,Chi_00)					\
  VLOAD(1,%r8,Chi_01)					\
  VLOAD(2,%r8,Chi_02)					\
  VLOAD(3,%r8,Chi_10)					\
  VLOAD(4,%r8,Chi_11)					\
  VLOAD(5,%r8,Chi_12)					\
							);

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

//   auto ptr = &U._odata[sU](A);		
// A plan for lifting loads 
//  can use Z2/3/4/5/U/U for U field in first step. 
//  can use Chi_00, Chi_10, U U for U field in second step
//  can use Chi_00, Chi_10, Chi_01,11, U U for U field in third step
// Enables to lift ALL loads earlier by a few cycles and alleviate OoO pressure if needed.
// KNL is DUAL issue for FP, and lifting these loads is potentially important.
// Need detailed profile data to be sure.

#define PREFETCH_U(A) \
  LOAD64(%r8,&U._odata[sU](A)) \
  __asm__ (		       \
  VPREFETCHG(0,%r8)	       \
  VPREFETCHG(1,%r8)	       \
  VPREFETCHG(2,%r8)	       \
  VPREFETCHG(3,%r8)	       \
  VPREFETCHG(4,%r8)	       \
  VPREFETCHG(5,%r8)	       \
  VPREFETCHG(6,%r8)	       \
  VPREFETCHG(7,%r8)	       \
  VPREFETCHG(8,%r8)	       );

#define PREFETCH_R(A)  \
  LOAD64(%r8,&out._odata[ss]) \
  __asm__ (		       \
  VPREFETCHW(0,%r8)	       \
  VPREFETCHW(1,%r8)	       \
  VPREFETCHW(2,%r8)	       \
  VPREFETCHW(3,%r8)	       \
  VPREFETCHW(4,%r8)	       \
  VPREFETCHW(5,%r8)	       \
  VPREFETCHW(6,%r8)	       \
  VPREFETCHW(7,%r8)	       \
  VPREFETCHW(8,%r8)	       \
  VPREFETCHW(9,%r8)	       \
  VPREFETCHW(10,%r8)	       \
  VPREFETCHW(11,%r8)	       );

 
#define MULT_2SPIN_DIR(A) MULT_2SPIN(&U._odata[sU](A))

#define MULT_2SPIN_DIR_PFXP(A,p) MULT_2SPIN_PFXP(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFYP(A,p) MULT_2SPIN_PFYP(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFZP(A,p) MULT_2SPIN_PFZP(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFTP(A,p) MULT_2SPIN_PFTP(&U._odata[sU](A),p)

#define MULT_2SPIN_DIR_PFXM(A,p) MULT_2SPIN_PFXM(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFYM(A,p) MULT_2SPIN_PFYM(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFZM(A,p) MULT_2SPIN_PFZM(&U._odata[sU](A),p)
#define MULT_2SPIN_DIR_PFTM(A,p) MULT_2SPIN_PFTM(&U._odata[sU](A),p)

#define MULT_2SPIN_UNOPT(ptr)				\
	   LOAD64(%r8,ptr)			\
  __asm__ (					\
	   ZLOAD (0,%r8,UChi_01,UChi_11)	\
	   ZLOAD (3,%r8,UChi_02,UChi_12)	\
	   ZLOAD (6,%r8,Uri,Uir)		\
	   ZMUL (UChi_01,UChi_11,Chi_00,UChi_00,Z0)	\
	   ZMUL (UChi_01,UChi_11,Chi_10,UChi_10,Z1)	\
	   ZMUL (UChi_02,UChi_12,Chi_00,UChi_01,Z2)	\
	   ZMUL (UChi_02,UChi_12,Chi_10,UChi_11,Z3)	\
	   ZMUL (Uri,Uir,        Chi_00,UChi_02,Z4)	\
	   ZMUL (Uri,Uir,        Chi_10,UChi_12,Z5)	\
	   						\
	   ZLOAD (1,%r8,Uri,Uir)			\
	   ZLOAD (4,%r8,Chi_00, Chi_10)		     	\
	   ZMADD (Uri,Uir,       Chi_01,UChi_00,Z0)	\
	   ZMADD (Uri,Uir,       Chi_11,UChi_10,Z1)	\
	   ZLOAD (7,%r8,Uri,Uir)			\
	   ZMADD (Chi_00, Chi_10,Chi_01,UChi_01,Z2)	\
	   ZMADD (Chi_00, Chi_10,Chi_11,UChi_11,Z3)	\
	   ZLOAD (2,%r8,Chi_00,Chi_10)			\
	   ZMADD(Uri,Uir,        Chi_01,UChi_02,Z4)	\
	   ZMADD(Uri,Uir,        Chi_11,UChi_12,Z5)	\
							\
	   ZLOAD  (5,%r8,Uri,Uir)			\
	   ZMADD (Chi_00,Chi_10, Chi_02,UChi_00,Z0)	\
	   ZMADD (Chi_00,Chi_10, Chi_12,UChi_10,Z1)	\
	   ZLOAD  (8,%r8,Chi_00,Chi_10)			\
	   ZMADD (Uri,Uir,       Chi_02,UChi_01,Z2)    	\
	   ZMADD (Uri,Uir,       Chi_12,UChi_11,Z3)	\
	   ZMADD(Chi_00,Chi_10,  Chi_02,UChi_02,Z4)	\
	   ZMADD(Chi_00,Chi_10,  Chi_12,UChi_12,Z5)	\
	   						\
	   ZEND1(UChi_00,Z0,Chi_01)			\
	   ZEND1(UChi_10,Z1,Chi_11)			\
	   ZEND1(UChi_01,Z2,Chi_00)			\
	   ZEND1(UChi_11,Z3,Chi_10)			\
	   ZEND1(UChi_02,Z4,Chi_02)			\
	   ZEND1(UChi_12,Z5,Chi_12)			\
	   ZEND2(UChi_00,Z0,Chi_01)			\
	   ZEND2(UChi_10,Z1,Chi_11)			\
	   ZEND2(UChi_01,Z2,Chi_00)			\
	   ZEND2(UChi_11,Z3,Chi_10)			\
	   ZEND2(UChi_02,Z4,Chi_02)			\
	   ZEND2(UChi_12,Z5,Chi_12)	     );

#define MULT_2SPIN(ptr)	MULT_2SPIN_PF(ptr,ptr,VPREFETCHG);

#define MULT_2SPIN_PFXM(ptr,pf) MULT_2SPIN_PF(ptr,pf,VPREFETCHNTA)
#define MULT_2SPIN_PFYM(ptr,pf) MULT_2SPIN_PF(ptr,pf,VPREFETCHNTA)
#define MULT_2SPIN_PFZM(ptr,pf) MULT_2SPIN_PF(ptr,pf,VPREFETCHNTA)
#define MULT_2SPIN_PFTM(ptr,pf) MULT_2SPIN_PF(ptr,pf,VPREFETCHNTA)
#define MULT_2SPIN_PFTP(ptr,pf) MULT_2SPIN_PF(ptr,pf,VPREFETCHNTA)
#define MULT_2SPIN_PFZP(ptr,pf) MULT_2SPIN_PF(ptr,pf,VPREFETCHNTA)
#define MULT_2SPIN_PFYP(ptr,pf) MULT_2SPIN_PF(ptr,pf,VPREFETCH)
#define MULT_2SPIN_PFXP(ptr,pf) MULT_2SPIN_PF(ptr,pf,VPREFETCHNTA)


#define MULT_2SPIN_PF(ptr,pf,VPF)			\
	   LOAD64(%r8,ptr)			\
	   LOAD64(%r9,pf)			\
  __asm__ (					\
	   ZMULMEM2SP(0,%r8,Uri,Chi_00,Chi_10,UChi_00,Z0,UChi_10,Z1)	\
	   VPF(0,%r9)						\
	   ZMULMEM2SP(3,%r8,Uri,Chi_00,Chi_10,UChi_01,Z2,UChi_11,Z3)	\
	   VPF(1,%r9)						\
	   ZMULMEM2SP(6,%r8,Uri,Chi_00,Chi_10,UChi_02,Z4,UChi_12,Z5)	\
	   VPF(2,%r9)						\
									\
	   ZMADDMEM2SP(1,%r8,Uri,Chi_01,Chi_11,UChi_00,Z0,UChi_10,Z1)	\
	   VPF(3,%r9)						\
	   ZMADDMEM2SP(4,%r8,Uri,Chi_01,Chi_11,UChi_01,Z2,UChi_11,Z3)	\
	   VPF(4,%r9)						\
	   ZMADDMEM2SP(7,%r8,Uri,Chi_01,Chi_11,UChi_02,Z4,UChi_12,Z5)	\
	   VPF(5,%r9)						\
									\
	   ZMADDMEM2SP(2,%r8,Uri,Chi_02,Chi_12,UChi_00,Z0,UChi_10,Z1)	\
	   VPF(6,%r9)						\
	   ZMADDMEM2SP(5,%r8,Uri,Chi_02,Chi_12,UChi_01,Z2,UChi_11,Z3)	\
	   VPF(7,%r9)						\
	   ZMADDMEM2SP(8,%r8,Uri,Chi_02,Chi_12,UChi_02,Z4,UChi_12,Z5)	\
	   VPF(8,%r9)						\
	   						\
	   ZEND1(UChi_00,Z0,Chi_01)			\
	   ZEND1(UChi_10,Z1,Chi_11)			\
	   ZEND1(UChi_01,Z2,Chi_00)			\
	   ZEND1(UChi_11,Z3,Chi_10)			\
	   VPF(9,%r9)						\
	   ZEND1(UChi_02,Z4,Chi_02)			\
	   ZEND1(UChi_12,Z5,Chi_12)			\
	   ZEND2(UChi_00,Z0,Chi_01)			\
	   ZEND2(UChi_10,Z1,Chi_11)			\
	   VPF(10,%r9)						\
	   ZEND2(UChi_01,Z2,Chi_00)			\
	   ZEND2(UChi_11,Z3,Chi_10)			\
	   ZEND2(UChi_02,Z4,Chi_02)			\
	   VPF(11,%r9)						\
	   ZEND2(UChi_12,Z5,Chi_12)	     );


#define MULT_2SPIN_PFNONE(ptr,pf,VPF)			\
	   LOAD64(%r8,ptr)			\
	   LOAD64(%r9,pf)			\
  __asm__ (					\
	   VPF(0,%r9)						\
	   VPF(1,%r9)						\
	   VPF(2,%r9)						\
	   							\
	   VPF(3,%r9)						\
	   VPF(4,%r9)						\
	   VPF(5,%r9)						\
	   							\
	   VPF(6,%r9)						\
	   VPF(7,%r9)						\
	   VPF(8,%r9)						\
	   							\
	   VPF(9,%r9)						\
	   VPF(10,%r9)						\
	   VPF(11,%r9)						);


// Pretty much Perfectly Pipelined

//////////////////////////////////////////////////////////////////
// Dirac algebra
//////////////////////////////////////////////////////////////////

//      hspin(0)=fspin(0)+timesI(fspin(3));
//      hspin(1)=fspin(1)+timesI(fspin(2));
//define VTIMESIf(A,DEST, Z)		
// These don't work if DEST==Z. FIXME.	
#define XP_PROJ __asm__ (		\
			 VACCTIMESI(Chimu_30,Chi_00,Z0)	\
			 VACCTIMESI(Chimu_31,Chi_01,Z1)	\
			 VACCTIMESI(Chimu_32,Chi_02,Z2)	\
			 VACCTIMESI(Chimu_20,Chi_10,Z3)	\
			 VACCTIMESI(Chimu_21,Chi_11,Z4)	\
			 VACCTIMESI(Chimu_22,Chi_12,Z5)	);

#define XP_PROJMEM(PTR) \
  LOAD64(%r8,PTR)							\
  __asm__ (								\
	   SHUF_CHIMU23i						\
	   VACCTIMESI1MEM(Chimu_30,Chi_00,0,%r8)			\
	   VACCTIMESI1MEM(Chimu_31,Chi_01,1,%r8)			\
	   VACCTIMESI1MEM(Chimu_32,Chi_02,2,%r8)			\
	   VACCTIMESI1MEM(Chimu_20,Chi_10,3,%r8)			\
	   VACCTIMESI1MEM(Chimu_21,Chi_11,4,%r8)			\
	   VACCTIMESI1MEM(Chimu_22,Chi_12,5,%r8)			\
	   VACCTIMESI2MEM(Chimu_30,Chi_00,0,%r8)			\
	   VACCTIMESI2MEM(Chimu_31,Chi_01,1,%r8)			\
	   VACCTIMESI2MEM(Chimu_32,Chi_02,2,%r8)			\
	   VACCTIMESI2MEM(Chimu_20,Chi_10,3,%r8)			\
	   VACCTIMESI2MEM(Chimu_21,Chi_11,4,%r8)			\
	   VACCTIMESI2MEM(Chimu_22,Chi_12,5,%r8)			);


#define YP_PROJ __asm__ ( 			\
  VSUB(Chimu_30,Chimu_00,Chi_00)\
  VSUB(Chimu_31,Chimu_01,Chi_01)\
  VSUB(Chimu_32,Chimu_02,Chi_02)\
  VADD(Chimu_10,Chimu_20,Chi_10)\
  VADD(Chimu_11,Chimu_21,Chi_11)\
  VADD(Chimu_12,Chimu_22,Chi_12) );

#define EVICT_SPINOR(reg) \
  VEVICT(0,reg)		  \
  VEVICT(1,reg)		  \
  VEVICT(2,reg)		  \
  VEVICT(3,reg)		  \
  VEVICT(4,reg)		  \
  VEVICT(5,reg)		  \
  VEVICT(6,reg)		  \
  VEVICT(7,reg)		  \
  VEVICT(8,reg)		  \
  VEVICT(9,reg)		  \
  VEVICT(9,reg)		  \
  VEVICT(10,reg)	  \
  VEVICT(11,reg)	  


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
  //  EVICT_SPINOR(%r8) 			);
  
#define ZP_PROJ __asm__ ( \
			 VACCTIMESI(Chimu_20,Chi_00,Z0)	\
			 VACCTIMESI(Chimu_21,Chi_01,Z1)	\
			 VACCTIMESI(Chimu_22,Chi_02,Z2)	\
			 VACCTIMESMINUSI(Chimu_30,Chi_10,Z3)	\
			 VACCTIMESMINUSI(Chimu_31,Chi_11,Z4)	\
			 VACCTIMESMINUSI(Chimu_32,Chi_12,Z5) );

#define ZP_PROJMEM(PTR) \
  LOAD64(%r8,PTR)							\
  __asm__ (								\
	   SHUF_CHIMU23i						\
	   VACCTIMESI1MEM(Chimu_20,Chi_00,0,%r8)			\
	   VACCTIMESI1MEM(Chimu_21,Chi_01,1,%r8)			\
	   VACCTIMESI1MEM(Chimu_22,Chi_02,2,%r8)			\
	   VACCTIMESMINUSI1MEM(Chimu_30,Chi_10,3,%r8)			\
	   VACCTIMESMINUSI1MEM(Chimu_31,Chi_11,4,%r8)			\
	   VACCTIMESMINUSI1MEM(Chimu_32,Chi_12,5,%r8)			\
	   VACCTIMESI2MEM(Chimu_20,Chi_00,0,%r8)			\
	   VACCTIMESI2MEM(Chimu_21,Chi_01,1,%r8)			\
	   VACCTIMESI2MEM(Chimu_22,Chi_02,2,%r8)			\
	   VACCTIMESMINUSI2MEM(Chimu_30,Chi_10,3,%r8)			\
	   VACCTIMESMINUSI2MEM(Chimu_31,Chi_11,4,%r8)			\
	   VACCTIMESMINUSI2MEM(Chimu_32,Chi_12,5,%r8)			\
	   EVICT_SPINOR(%r8) 			);



#define TP_PROJ __asm__ (			\
			 VADD(Chimu_00,Chimu_20,Chi_00)	\
			 VADD(Chimu_01,Chimu_21,Chi_01)	\
			 VADD(Chimu_02,Chimu_22,Chi_02)	\
			 VADD(Chimu_10,Chimu_30,Chi_10)	\
			 VADD(Chimu_11,Chimu_31,Chi_11)	\
			 VADD(Chimu_12,Chimu_32,Chi_12) );


#define TP_PROJMEM(ptr)				\
  LOAD64(%r8,ptr)				\
  __asm__ (					\
	   LOAD_CHIMU01i			\
	   VADDMEM(6,%r8 ,Chimu_00,Chi_00)	\
	   VADDMEM(7,%r8,Chimu_01,Chi_01)	\
	   VADDMEM(8,%r8,Chimu_02,Chi_02)	\
	   VADDMEM(9,%r8,Chimu_10,Chi_10)	\
	   VADDMEM(10,%r8,Chimu_11,Chi_11)	\
	   VADDMEM(11,%r8,Chimu_12,Chi_12)	\
	   EVICT_SPINOR(%r8) 			);


//      hspin(0)=fspin(0)-timesI(fspin(3))
//      hspin(1)=fspin(1)-timesI(fspin(2))
#define XM_PROJ __asm__ ( \
 VACCTIMESMINUSI(Chimu_30,Chi_00,Z0)		\
 VACCTIMESMINUSI(Chimu_31,Chi_01,Z1)		\
 VACCTIMESMINUSI(Chimu_32,Chi_02,Z2)		\
 VACCTIMESMINUSI(Chimu_20,Chi_10,Z3)		\
 VACCTIMESMINUSI(Chimu_21,Chi_11,Z4)		\
 VACCTIMESMINUSI(Chimu_22,Chi_12,Z5)		);

#define XM_PROJMEM(PTR) \
  LOAD64(%r8,PTR)							\
  __asm__ (								\
	   SHUF_CHIMU23i						\
	   VACCTIMESMINUSI1MEM(Chimu_30,Chi_00,0,%r8)			\
	   VACCTIMESMINUSI1MEM(Chimu_31,Chi_01,1,%r8)			\
	   VACCTIMESMINUSI1MEM(Chimu_32,Chi_02,2,%r8)			\
	   VACCTIMESMINUSI1MEM(Chimu_20,Chi_10,3,%r8)			\
	   VACCTIMESMINUSI1MEM(Chimu_21,Chi_11,4,%r8)			\
	   VACCTIMESMINUSI1MEM(Chimu_22,Chi_12,5,%r8)			\
	   VACCTIMESMINUSI2MEM(Chimu_30,Chi_00,0,%r8)			\
	   VACCTIMESMINUSI2MEM(Chimu_31,Chi_01,1,%r8)			\
	   VACCTIMESMINUSI2MEM(Chimu_32,Chi_02,2,%r8)			\
	   VACCTIMESMINUSI2MEM(Chimu_20,Chi_10,3,%r8)			\
	   VACCTIMESMINUSI2MEM(Chimu_21,Chi_11,4,%r8)			\
	   VACCTIMESMINUSI2MEM(Chimu_22,Chi_12,5,%r8)			);

#define YM_PROJ __asm__ ( \
  VADD(Chimu_00,Chimu_30,Chi_00)\
  VADD(Chimu_01,Chimu_31,Chi_01)\
  VADD(Chimu_02,Chimu_32,Chi_02)\
  VSUB(Chimu_20,Chimu_10,Chi_10)\
  VSUB(Chimu_21,Chimu_11,Chi_11)\
  VSUB(Chimu_22,Chimu_12,Chi_12) );

#define YM_PROJMEM(ptr)				\
  LOAD64(%r8,ptr)				\
  __asm__ (					\
  LOAD_CHIMU01i					\
  VADDMEM(9,%r8 ,Chimu_00,Chi_00)		\
  VADDMEM(10,%r8,Chimu_01,Chi_01)		\
  VADDMEM(11,%r8,Chimu_02,Chi_02)		\
  VSUBMEM(6,%r8,Chimu_10,Chi_10)		\
  VSUBMEM(7,%r8,Chimu_11,Chi_11)		\
  VSUBMEM(8,%r8,Chimu_12,Chi_12)			\
  EVICT_SPINOR(%r8) 			);


#define ZM_PROJ __asm__ ( \
  VACCTIMESMINUSI(Chimu_20,Chi_00,Z0)\
  VACCTIMESMINUSI(Chimu_21,Chi_01,Z1)\
  VACCTIMESMINUSI(Chimu_22,Chi_02,Z2)\
  VACCTIMESI(Chimu_30,Chi_10,Z3)\
  VACCTIMESI(Chimu_31,Chi_11,Z4)\
  VACCTIMESI(Chimu_32,Chi_12,Z5));

#define ZM_PROJMEM(PTR) \
  LOAD64(%r8,PTR)							\
  __asm__ (								\
	   SHUF_CHIMU23i						\
	   VACCTIMESMINUSI1MEM(Chimu_20,Chi_00,0,%r8)			\
	   VACCTIMESMINUSI1MEM(Chimu_21,Chi_01,1,%r8)			\
	   VACCTIMESMINUSI1MEM(Chimu_22,Chi_02,2,%r8)			\
	   VACCTIMESI1MEM(Chimu_30,Chi_10,3,%r8)			\
	   VACCTIMESI1MEM(Chimu_31,Chi_11,4,%r8)			\
	   VACCTIMESI1MEM(Chimu_32,Chi_12,5,%r8)			\
	   VACCTIMESMINUSI2MEM(Chimu_20,Chi_00,0,%r8)			\
	   VACCTIMESMINUSI2MEM(Chimu_21,Chi_01,1,%r8)			\
	   VACCTIMESMINUSI2MEM(Chimu_22,Chi_02,2,%r8)			\
	   VACCTIMESI2MEM(Chimu_30,Chi_10,3,%r8)			\
	   VACCTIMESI2MEM(Chimu_31,Chi_11,4,%r8)			\
	   VACCTIMESI2MEM(Chimu_32,Chi_12,5,%r8)			\
	   EVICT_SPINOR(%r8) 			);


#define TM_PROJ __asm__ ( \
  VSUB(Chimu_20,Chimu_00,Chi_00)\
  VSUB(Chimu_21,Chimu_01,Chi_01)\
  VSUB(Chimu_22,Chimu_02,Chi_02)\
  VSUB(Chimu_30,Chimu_10,Chi_10)\
  VSUB(Chimu_31,Chimu_11,Chi_11)\
  VSUB(Chimu_32,Chimu_12,Chi_12) );


#define TM_PROJMEM(ptr)				\
  LOAD64(%r8,ptr)				\
  __asm__ (					\
  LOAD_CHIMU01i					\
  VSUBMEM(6,%r8 ,Chimu_00,Chi_00)		\
  VSUBMEM(7,%r8,Chimu_01,Chi_01)		\
  VSUBMEM(8,%r8,Chimu_02,Chi_02)		\
  VSUBMEM(9,%r8,Chimu_10,Chi_10)		\
  VSUBMEM(10,%r8,Chimu_11,Chi_11)		\
  VSUBMEM(11,%r8,Chimu_12,Chi_12)		\
  EVICT_SPINOR(%r8) );

//      fspin(0)=hspin(0)
//      fspin(1)=hspin(1)
//      fspin(2)=timesMinusI(hspin(1))
//      fspin(3)=timesMinusI(hspin(0))
#define XP_RECON __asm__ (			\
			  VZERO(TMP)		\
			  VMOV(UChi_00,result_00)	\
			  VMOV(UChi_01,result_01)	\
			  VMOV(UChi_02,result_02)	\
			  VMOV(UChi_10,result_10)	\
			  VMOV(UChi_11,result_11)	\
			  VMOV(UChi_12,result_12)	\
			  VTIMESMINUSI0(UChi_10,result_20,TMP)	\
			  VTIMESMINUSI0(UChi_11,result_21,TMP)	\
			  VTIMESMINUSI0(UChi_12,result_22,TMP)	\
			  VTIMESMINUSI0(UChi_00,result_30,TMP)	\
			  VTIMESMINUSI0(UChi_01,result_31,TMP)	\
			  VTIMESMINUSI0(UChi_02,result_32,TMP)   \
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
  VADD(UChi_00,result_00,result_00)\
  VADD(UChi_01,result_01,result_01)\
  VADD(UChi_02,result_02,result_02)\
  VADD(UChi_10,result_10,result_10)\
  VADD(UChi_11,result_11,result_11)\
  VADD(UChi_12,result_12,result_12)\
  VACCTIMESMINUSI0(UChi_10,result_20,Z0)\
  VACCTIMESMINUSI0(UChi_11,result_21,Z1)\
  VACCTIMESMINUSI0(UChi_12,result_22,Z2)\
  VACCTIMESMINUSI0(UChi_00,result_30,Z3)\
  VACCTIMESMINUSI0(UChi_01,result_31,Z4)\
  VACCTIMESMINUSI0(UChi_02,result_32,Z5)\
  VACCTIMESMINUSI1(UChi_10,result_20,Z0)\
  VACCTIMESMINUSI1(UChi_11,result_21,Z1)\
  VACCTIMESMINUSI1(UChi_12,result_22,Z2)\
  VACCTIMESMINUSI1(UChi_00,result_30,Z3)\
  VACCTIMESMINUSI1(UChi_01,result_31,Z4)\
  VACCTIMESMINUSI1(UChi_02,result_32,Z5)\
  VACCTIMESMINUSI2(UChi_10,result_20,Z0)\
  VACCTIMESMINUSI2(UChi_11,result_21,Z1)\
  VACCTIMESMINUSI2(UChi_12,result_22,Z2)\
  VACCTIMESMINUSI2(UChi_00,result_30,Z3)\
  VACCTIMESMINUSI2(UChi_01,result_31,Z4)\
  VACCTIMESMINUSI2(UChi_02,result_32,Z5)\
				 );

#define XM_RECON __asm__ ( \
  VZERO(TMP)\
  VMOV(UChi_00,result_00)\
  VMOV(UChi_01,result_01)\
  VMOV(UChi_02,result_02)\
  VMOV(UChi_10,result_10)\
  VMOV(UChi_11,result_11)\
  VMOV(UChi_12,result_12)\
  VTIMESI0(UChi_10,result_20,TMP)\
  VTIMESI0(UChi_11,result_21,TMP)\
  VTIMESI0(UChi_12,result_22,TMP)\
  VTIMESI0(UChi_00,result_30,TMP)\
  VTIMESI0(UChi_01,result_31,TMP)\
  VTIMESI0(UChi_02,result_32,TMP)\
  VTIMESI1(UChi_10,result_20,TMP)\
  VTIMESI1(UChi_11,result_21,TMP)\
  VTIMESI1(UChi_12,result_22,TMP)\
  VTIMESI1(UChi_00,result_30,TMP)\
  VTIMESI1(UChi_01,result_31,TMP)\
  VTIMESI1(UChi_02,result_32,TMP)\
  VTIMESI2(UChi_10,result_20,TMP)\
  VTIMESI2(UChi_11,result_21,TMP)\
  VTIMESI2(UChi_12,result_22,TMP)\
  VTIMESI2(UChi_00,result_30,TMP)\
  VTIMESI2(UChi_01,result_31,TMP)\
  VTIMESI2(UChi_02,result_32,TMP)\
			   );

#define XM_RECON_ACCUM __asm__ ( \
  VADD(UChi_00,result_00,result_00)\
  VADD(UChi_01,result_01,result_01)\
  VADD(UChi_02,result_02,result_02)\
  VADD(UChi_10,result_10,result_10)\
  VADD(UChi_11,result_11,result_11)\
  VADD(UChi_12,result_12,result_12)\
  VACCTIMESI0(UChi_10,result_20,Z0)\
  VACCTIMESI0(UChi_11,result_21,Z1)\
  VACCTIMESI0(UChi_12,result_22,Z2)\
  VACCTIMESI0(UChi_00,result_30,Z3)\
  VACCTIMESI0(UChi_01,result_31,Z4)\
  VACCTIMESI0(UChi_02,result_32,Z5)\
  VACCTIMESI1(UChi_10,result_20,Z0)\
  VACCTIMESI1(UChi_11,result_21,Z1)\
  VACCTIMESI1(UChi_12,result_22,Z2)\
  VACCTIMESI1(UChi_00,result_30,Z3)\
  VACCTIMESI1(UChi_01,result_31,Z4)\
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
  VADD(UChi_01,result_01,result_01)\
  VADD(UChi_02,result_02,result_02)\
  VADD(UChi_10,result_10,result_10)\
  VADD(UChi_11,result_11,result_11)\
  VADD(UChi_12,result_12,result_12)\
  VADD(UChi_10,result_20,result_20)\
  VADD(UChi_11,result_21,result_21)\
  VADD(UChi_12,result_22,result_22)\
  VSUB(UChi_00,result_30,result_30)\
  VSUB(UChi_01,result_31,result_31)\
  VSUB(UChi_02,result_32,result_32) );

#define YM_RECON_ACCUM __asm__ ( \
  VADD(UChi_00,result_00,result_00)\
  VADD(UChi_01,result_01,result_01)\
  VADD(UChi_02,result_02,result_02)\
  VADD(UChi_10,result_10,result_10)\
  VADD(UChi_11,result_11,result_11)\
  VADD(UChi_12,result_12,result_12)\
  VSUB(UChi_10,result_20,result_20)\
  VSUB(UChi_11,result_21,result_21)\
  VSUB(UChi_12,result_22,result_22)\
  VADD(UChi_00,result_30,result_30)\
  VADD(UChi_01,result_31,result_31)\
  VADD(UChi_02,result_32,result_32) );

#define ZP_RECON_ACCUM __asm__ ( \
  VADD(UChi_00,result_00,result_00)\
  VADD(UChi_01,result_01,result_01)\
  VADD(UChi_02,result_02,result_02)\
  VADD(UChi_10,result_10,result_10)\
  VADD(UChi_11,result_11,result_11)\
  VADD(UChi_12,result_12,result_12)\
  VACCTIMESMINUSI0(UChi_00,result_20,Z0)\
  VACCTIMESMINUSI0(UChi_01,result_21,Z1)\
  VACCTIMESMINUSI0(UChi_02,result_22,Z2)\
  VACCTIMESI0(UChi_10,result_30,Z3)\
  VACCTIMESI0(UChi_11,result_31,Z4)\
  VACCTIMESI0(UChi_12,result_32,Z5)\
  VACCTIMESMINUSI1(UChi_00,result_20,Z0)\
  VACCTIMESMINUSI1(UChi_01,result_21,Z1)\
  VACCTIMESMINUSI1(UChi_02,result_22,Z2)\
  VACCTIMESI1(UChi_10,result_30,Z3)\
  VACCTIMESI1(UChi_11,result_31,Z4)\
  VACCTIMESI1(UChi_12,result_32,Z5)\
  VACCTIMESMINUSI2(UChi_00,result_20,Z0)\
  VACCTIMESMINUSI2(UChi_01,result_21,Z1)\
  VACCTIMESMINUSI2(UChi_02,result_22,Z2)\
  VACCTIMESI2(UChi_10,result_30,Z3)\
  VACCTIMESI2(UChi_11,result_31,Z4)\
  VACCTIMESI2(UChi_12,result_32,Z5)\
				 );

#define ZM_RECON_ACCUM __asm__ ( \
  VADD(UChi_00,result_00,result_00)\
  VADD(UChi_01,result_01,result_01)\
  VADD(UChi_02,result_02,result_02)\
  VADD(UChi_10,result_10,result_10)\
  VADD(UChi_11,result_11,result_11)\
  VADD(UChi_12,result_12,result_12)\
  VACCTIMESI0(UChi_00,result_20,Z0)\
  VACCTIMESI0(UChi_01,result_21,Z1)\
  VACCTIMESI0(UChi_02,result_22,Z2)\
  VACCTIMESMINUSI0(UChi_10,result_30,Z3)\
  VACCTIMESMINUSI0(UChi_11,result_31,Z4)\
  VACCTIMESMINUSI0(UChi_12,result_32,Z5)\
  VACCTIMESI1(UChi_00,result_20,Z0)\
  VACCTIMESI1(UChi_01,result_21,Z1)\
  VACCTIMESI1(UChi_02,result_22,Z2)\
  VACCTIMESMINUSI1(UChi_10,result_30,Z3)\
  VACCTIMESMINUSI1(UChi_11,result_31,Z4)\
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
  VADD(UChi_01,result_01,result_01)\
  VADD(UChi_02,result_02,result_02)\
  VADD(UChi_10,result_10,result_10)\
  VADD(UChi_11,result_11,result_11)\
  VADD(UChi_12,result_12,result_12)\
  VADD(UChi_00,result_20,result_20)\
  VADD(UChi_01,result_21,result_21)\
  VADD(UChi_02,result_22,result_22)\
  VADD(UChi_10,result_30,result_30)\
  VADD(UChi_11,result_31,result_31)\
  VADD(UChi_12,result_32,result_32) );

#define TM_RECON_ACCUM __asm__ ( \
  VADD(UChi_00,result_00,result_00)\
  VADD(UChi_01,result_01,result_01)\
  VADD(UChi_02,result_02,result_02)\
  VADD(UChi_10,result_10,result_10)\
  VADD(UChi_11,result_11,result_11)\
  VADD(UChi_12,result_12,result_12)\
  VSUB(UChi_00,result_20,result_20)\
  VSUB(UChi_01,result_21,result_21)\
  VSUB(UChi_02,result_22,result_22)\
  VSUB(UChi_10,result_30,result_30)\
  VSUB(UChi_11,result_31,result_31)\
  VSUB(UChi_12,result_32,result_32) );

#define PREFETCH_CHIMU(A) 

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

#endif
