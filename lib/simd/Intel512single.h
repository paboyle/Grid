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
// No guard can be multiply included as undef clearge of macros
#undef VZERO
#undef VMOV
#undef VLOAD
#undef VSTORE
#define VZERO(A)                  VZEROf(A)
#define VMOV(A,B)                 VMOVf(A,B)
#define VLOAD(OFF,PTR,DEST)       VLOADf(OFF,PTR,DEST)
#define VSTORE(OFF,PTR,SRC)       VSTOREf(OFF,PTR,SRC)

#undef VADD
#undef VSUB
#undef VMUL
#undef VMADD
#define VADD(A,B,C)               VADDf(A,B,C)
#define VSUB(A,B,C)               VSUBf(A,B,C)
#define VMUL(Uri,Uir,Chi)         VMULf(Uri,Uir,Chi)
#define VMADD(Uri,Uir,Chi)        VMADDf(Uri,Uir,Chi)


#undef VTIMESI
#undef VTIMESI0 
#undef VTIMESI1
#undef VTIMESI2 
#define VTIMESI(A,B,C)                 VTIMESIf(A,B,C)
#define VTIMESI0(A,B,C)                VTIMESI0f(A,B,C)
#define VTIMESI1(A,B,C)                VTIMESI1f(A,B,C)
#define VTIMESI2(A,B,C)                VTIMESI2f(A,B,C)

#undef VTIMESMINUSI
#undef VTIMESMINUSI0
#undef VTIMESMINUSI1
#undef VTIMESMINUSI2
#define VTIMESMINUSI(A,B,C)            VTIMESMINUSIf(A,B,C)
#define VTIMESMINUSI0(A,B,C)           VTIMESMINUSI0f(A,B,C)
#define VTIMESMINUSI1(A,B,C)           VTIMESMINUSI1f(A,B,C)
#define VTIMESMINUSI2(A,B,C)           VTIMESMINUSI2f(A,B,C)

#undef VACCTIMESI
#undef VACCTIMESI0
#undef VACCTIMESI1
#undef VACCTIMESI2
#define VACCTIMESI(A,B,C)         VACCTIMESIf(A,B,C)
#define VACCTIMESI0(A,B,C)             VACCTIMESI0f(A,B,C)
#define VACCTIMESI1(A,B,C)             VACCTIMESI1f(A,B,C)
#define VACCTIMESI2(A,B,C)             VACCTIMESI2f(A,B,C)

#undef VACCTIMESMINUSI
#undef VACCTIMESMINUSI0
#undef VACCTIMESMINUSI1
#undef VACCTIMESMINUSI2
#define VACCTIMESMINUSI(A,B,C)    VACCTIMESMINUSIf(A,B,C)
#define VACCTIMESMINUSI0(A,B,C)        VACCTIMESMINUSI0f(A,B,C)
#define VACCTIMESMINUSI1(A,B,C)        VACCTIMESMINUSI1f(A,B,C)
#define VACCTIMESMINUSI2(A,B,C)        VACCTIMESMINUSI2f(A,B,C)

#undef VACCTIMESI1MEM
#undef VACCTIMESI2MEM
#define VACCTIMESI1MEM(A,ACC,O,P)      VACCTIMESI1MEMf(A,ACC,O,P)
#define VACCTIMESI2MEM(A,ACC,O,P)      VACCTIMESI2MEMf(A,ACC,O,P)

#undef VACCTIMESMINUSI1MEM
#undef VACCTIMESMINUSI2MEM
#define VACCTIMESMINUSI1MEM(A,ACC,O,P) VACCTIMESMINUSI1MEMf(A,ACC,O,P)
#define VACCTIMESMINUSI2MEM(A,ACC,O,P) VACCTIMESMINUSI2MEMf(A,ACC,O,P)

#undef VPERM0
#undef VPERM1
#undef VPERM2
#undef VPERM3
#define VPERM0(A,B)               VPERM0f(A,B)
#define VPERM1(A,B)               VPERM1f(A,B)
#define VPERM2(A,B)               VPERM2f(A,B)
#define VPERM3(A,B)               VPERM3f(A,B)

#undef VSHUFMEM
#undef VADDMEM
#undef VSUBMEM
#define VSHUFMEM(OFF,A,DEST)      VSHUFMEMf(OFF,A,DEST)
#define VADDMEM(O,A,B,C)                                 VADDMEMf(O,A,B,C)
#define VSUBMEM(O,A,B,C)                                 VSUBMEMf(O,A,B,C)

#undef VMOVIDUP
#undef VMOVRDUP
#undef VMADDSUB
#undef VSHUF
#define VMOVIDUP(A,B,C)                                  VMOVIDUPf(A,B,C)
#define VMOVRDUP(A,B,C)                                  VMOVRDUPf(A,B,C)
#define VMADDSUB(A,B,accum)                              VMADDSUBf(A,B,accum) 
#define VSHUF(A,B)                                       VSHUFf(A,B)


#undef ZEND1
#undef ZEND2
#undef ZLOAD
#undef ZMUL
#undef ZMADD
#undef ZMULMEM2SP
#undef ZMADDMEM2SP

#define ZEND1(A,B,C)                                     ZEND1f(A,B,C)
#define ZEND2(A,B,C)                                     ZEND2f(A,B,C)
#define ZLOAD(A,B,C,D)                                   ZLOADf(A,B,C,D)
#define ZMUL(A,B,C,D,E)                                  ZMULf(A,B,C,D,E)
#define ZMADD(A,B,C,D,E)                                 ZMADDf(A,B,C,D,E)
#define ZMULMEM2SP(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr)  ZMULMEM2SPf(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr) 
#define ZMADDMEM2SP(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr) ZMADDMEM2SPf(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr) 

