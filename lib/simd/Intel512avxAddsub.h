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
#ifndef GRID_ASM_AV512_ADDSUB_H
#define GRID_ASM_AV512_ADDSUB_H


////////////////////////////////////////////////////////////////
// Building blocks for SU3 x 2spinor
// Load columns of U
// 18   U DUP's  rr/ii
//  6 Chi shuffles ir,ri
// 6muls, 30 fmaddsubs
////////////////////////////////////////////////////////////////
#define MULT_ADDSUB_2SPIN(ptr)	\
	   LOAD64(%r8,ptr)			\
  __asm__ (					\
	   VMOVIDUPf(0,%r8,Z0 ) \
	   VMOVIDUPf(3,%r8,Z1 )\
	   VMOVIDUPf(6,%r8,Z2 )\
           VSHUFf(Chi_00,T1)    \
           VSHUFf(Chi_10,T2)    \
				\
	   VMULf(Z0,T1,UChi_00)	           VMOVRDUPf(0,%r8,Z3 ) \
	   VMULf(Z0,T2,UChi_10)	           VMOVRDUPf(3,%r8,Z4 ) \
	   VMULf(Z1,T1,UChi_01)	           VMOVRDUPf(6,%r8,Z5 ) \
	   VMULf(Z1,T2,UChi_11)	           VMOVIDUPf(1,%r8,Z0 ) \
	   VMULf(Z2,T1,UChi_02)            VMOVIDUPf(4,%r8,Z1 ) \
	   VMULf(Z2,T2,UChi_12)	           VMOVIDUPf(7,%r8,Z2 ) \
	   			\
	   VMADDSUBf(Z3,Chi_00,UChi_00)    VSHUFf(Chi_01,T1)    \
	   VMADDSUBf(Z3,Chi_10,UChi_10)    VSHUFf(Chi_11,T2)    \
	   VMADDSUBf(Z4,Chi_00,UChi_01)    VMOVRDUPf(1,%r8,Z3 ) \
	   VMADDSUBf(Z4,Chi_10,UChi_11)\
	   VMADDSUBf(Z5,Chi_00,UChi_02)    VMOVRDUPf(4,%r8,Z4 ) \
	   VMADDSUBf(Z5,Chi_10,UChi_12)\
	   			       \
	   VMADDSUBf(Z0,T1,UChi_00) 	   VMOVRDUPf(7,%r8,Z5 ) \ 
	   VMADDSUBf(Z0,T2,UChi_10)\
	   VMADDSUBf(Z1,T1,UChi_01) 	   VMOVIDUPf(2,%r8,Z0 ) \
	   VMADDSUBf(Z1,T2,UChi_11)\
	   VMADDSUBf(Z2,T1,UChi_02)        VMOVIDUPf(5,%r8,Z1 ) \
	   VMADDSUBf(Z2,T2,UChi_12)        VMOVIDUPf(8,%r8,Z2 ) \
					   			\
	   VMADDSUBf(Z3,Chi_01,UChi_00)    VSHUFf(Chi_02,T1)    \
	   VMADDSUBf(Z3,Chi_11,UChi_10)    VSHUFf(Chi_12,T2)    \
	   VMADDSUBf(Z4,Chi_01,UChi_01)	   VMOVRDUPf(2,%r8,Z3 ) \
	   VMADDSUBf(Z4,Chi_11,UChi_11)\
	   VMADDSUBf(Z5,Chi_01,UChi_02)	   VMOVRDUPf(5,%r8,Z4 ) \
	   VMADDSUBf(Z5,Chi_11,UChi_12)\
	   			\
	   VMADDSUBf(Z0,T1,UChi_00) 	   VMOVRDUPf(8,%r8,Z5 ) \
	   VMADDSUBf(Z0,T2,UChi_10)\
	   VMADDSUBf(Z1,T1,UChi_01)\
	   VMADDSUBf(Z1,T2,UChi_11)\
	   VMADDSUBf(Z2,T1,UChi_02)\
	   VMADDSUBf(Z2,T2,UChi_12)\
		   		   \
	   VMADDSUBf(Z3,Chi_02,UChi_00)\
	   VMADDSUBf(Z3,Chi_12,UChi_10)\
	   VMADDSUBf(Z4,Chi_02,UChi_01)\
	   VMADDSUBf(Z4,Chi_12,UChi_11)\
	   VMADDSUBf(Z5,Chi_02,UChi_02)\
	   VMADDSUBf(Z5,Chi_12,UChi_12)\
						);


#endif
