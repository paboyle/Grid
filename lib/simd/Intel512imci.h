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

////////////////////////////////////////////////////////////	  
// Knights Corner specials
////////////////////////////////////////////////////////////	  

#define ZLOADf(OFF,PTR,ri,ir)  VLOADf(OFF,PTR,ir)  VSHUFf(ir,ri)
#define ZLOADd(OFF,PTR,ri,ir)  VLOADd(OFF,PTR,ir)  VSHUFd(ir,ri)

#define ZMULf(Ari,Air,B,Criir,Ciirr)  VMULf(Ari,B,Criir)  VMULf(Air,B,Ciirr)
#define ZMULd(Ari,Air,B,Criir,Ciirr)  VMULd(Ari,B,Criir)  VMULd(Air,B,Ciirr)

#define ZMADDf(Ari,Air,B,Criir,Ciirr) VMADDf(Ari,B,Criir) VMADDf(Air,B,Ciirr)
#define ZMADDd(Ari,Air,B,Criir,Ciirr) VMADDd(Ari,B,Criir) VMADDd(Air,B,Ciirr)

#define ZENDf(Criir,Ciirr, tmp) ZEND1f(Criir,Ciirr, tmp) ZEND2f(Criir,Ciirr, tmp)
#define ZENDd(Criir,Ciirr, tmp) ZEND1d(Criir,Ciirr, tmp) ZEND2d(Criir,Ciirr, tmp)

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

#define ZMADDMEM2SPd(O,P,tmp,B,C,Briir,Biirr,Criir,Ciirr)	\
  VSHUFMEMd(O,P,tmp) \
  VMADDMEMd(O,P,B,Biirr) \
  VMADDMEMd(O,P,C,Ciirr) \
  VMADDd(tmp,B,Briir) \
  VMADDd(tmp,C,Criir)

#define ZEND1d(Criir,Ciirr, tmp) "vaddpd  " #Criir "{cdab} ," #Criir "," #Criir"{%k6}"  ";\n"
#define ZEND2d(Criir,Ciirr, tmp) "vsubpd  " #Ciirr "{cdab} ," #Ciirr "," #Criir"{%k7}"  ";\n"

#define ZEND1f(Criir,Ciirr, tmp) "vaddps  " #Criir "{cdab} ," #Criir "," #Criir"{%k6}"  ";\n"
#define ZEND2f(Criir,Ciirr, tmp) "vsubps  " #Ciirr "{cdab} ," #Ciirr "," #Criir"{%k7}"  ";\n"

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

	   // Acc = Acc - i A
#define VACCTIMESMINUSI0d(A,ACC,tmp)  
#define VACCTIMESMINUSI1d(A,ACC,tmp)  "vsubpd  " #A "{cdab}," #ACC "," #ACC"{%k7}" ";\n"
#define VACCTIMESMINUSI2d(A,ACC,tmp)  "vaddpd  " #A "{cdab}," #ACC "," #ACC"{%k6}" ";\n"

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
