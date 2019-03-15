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
// No guard; ok multi-include
#undef VSIZE
#undef VLOAD
#undef VLOADu
#undef VSPLAT
#undef VSTORE
#undef VSTOREu
#undef MULT_2SPIN_QPX_LS
#undef MULT_2SPIN_QPX

#define VSIZE VSIZEd
#define VLOAD(A,B,C)     VLOADd(A,B,C)
#define VLOADu(A,B,C)    VLOADud(A,B,C)
#define VSPLAT(A,B,DEST) VSPLATd(A,B,DEST)
#define VSTORE(A,B,C)    VSTOREd(A,B,C)
#define VSTOREu(A,B,C)   VSTOREud(A,B,C)
#define MULT_2SPIN_QPX_LS(ptr,p) MULT_2SPIN_QPX_LSd(ptr,p)
#define MULT_2SPIN_QPX(ptr,p)    MULT_2SPIN_QPXd(ptr,p)

