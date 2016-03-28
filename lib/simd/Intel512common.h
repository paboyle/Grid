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
#ifndef GRID_ASM_INTEL_COMMON_512_H
#define GRID_ASM_INTEL_COMMON_512_H

////////////////////////////////////////////////////////////////////////////////////////////////////
// Opcodes common 
////////////////////////////////////////////////////////////////////////////////////////////////////
#define MASK_REGS \
  __asm__ ("mov     $0xAAAA, %%eax \n"\ 
           "kmovw    %%eax, %%k6 \n"\
           "mov     $0x5555, %%eax \n"\
           "kmovw    %%eax, %%k7 \n" : : : "%eax");

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
#define LOAD64(A,ptr)  LOAD64i(A,ptr)

#define VMOVf(A,DEST)   "vmovaps  " #A ", " #DEST  ";\n"
#define VMOVd(A,DEST)   "vmovapd  " #A ", " #DEST  ";\n"

#define VPREFETCHG(O,A) 
#define VPREFETCHW(O,A) 
#define VEVICT(O,A)   

//"vprefetche0 "#O"*64("#A");\n" "vprefetche1 ("#O"+12)*64("#A");\n"
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

#define VPREFETCHNTA(O,A) 
#define VPREFETCH(O,A)    

#define VSTOREf(OFF,PTR,SRC)   "vmovaps " #SRC "," #OFF "*64(" #PTR ")"  ";\n"
#define VSTOREd(OFF,PTR,SRC)   "vmovapd " #SRC "," #OFF "*64(" #PTR ")"  ";\n"

// Swaps Re/Im ; could unify this with IMCI
#define VSHUFd(A,DEST)         "vpshufd  $0x4e," #A "," #DEST  ";\n"    
#define VSHUFf(A,DEST)         "vpshufd  $0xb1," #A "," #DEST  ";\n"    
#define VSHUFMEMd(OFF,A,DEST)  "vpshufd  $0x4e, " #OFF"*64("#A ")," #DEST  ";\n" // 32 bit level: 1,0,3,2
#define VSHUFMEMf(OFF,A,DEST)  "vpshufd  $0xb1, " #OFF"*64("#A ")," #DEST  ";\n" // 32 bit level: 2,3,0,1

#define TRAP " int3 ;\n"

#endif
