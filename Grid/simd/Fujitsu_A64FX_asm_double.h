/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: XXX

    Copyright (C) 2020

Author: Nils Meyer <nils.meyer@ur.de>

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
#define LOAD_CHIMU_A64FXd(x)           LOAD_CHIMU_INTERLEAVED_A64FXd(x)  
#define PREFETCH_CHIMU_L1(A)  
#define PREFETCH_GAUGE_L1(A)  
#define PREFETCH_CHIMU_L2(A)  
#define PREFETCH_GAUGE_L2(A)  
#define PF_GAUGE(A)  
#define PREFETCH1_CHIMU(A)  
#define PREFETCH_CHIMU(A)  
#define LOCK_GAUGE(A)  
#define UNLOCK_GAUGE(A)  
#define MASK_REGS                      DECLARATIONS_A64FXd(A)  
#define COMPLEX_SIGNS(A)  
#define LOAD64(A,B)  
#define SAVE_RESULT(A,B)               RESULT_A64FXd(A)  
#define MULT_2SPIN_DIR_PF(A,B)         MULT_2SPIN_A64FXd(A)  
#define MAYBEPERM(A,perm)              if (perm) { A ; }  
#define LOAD_CHI(base)                 LOAD_CHI_A64FXd(base)  
#define ZERO_PSI                       ZERO_PSI_A64FXd  
#define XP_PROJMEM(base)               LOAD_CHIMU_A64FXd(base);   XP_PROJ_A64FXd  
#define YP_PROJMEM(base)               LOAD_CHIMU_A64FXd(base);   YP_PROJ_A64FXd  
#define ZP_PROJMEM(base)               LOAD_CHIMU_A64FXd(base);   ZP_PROJ_A64FXd  
#define TP_PROJMEM(base)               LOAD_CHIMU_A64FXd(base);   TP_PROJ_A64FXd  
#define XM_PROJMEM(base)               LOAD_CHIMU_A64FXd(base);   XM_PROJ_A64FXd  
#define YM_PROJMEM(base)               LOAD_CHIMU_A64FXd(base);   YM_PROJ_A64FXd  
#define ZM_PROJMEM(base)               LOAD_CHIMU_A64FXd(base);   ZM_PROJ_A64FXd  
#define TM_PROJMEM(base)               LOAD_CHIMU_A64FXd(base);   TM_PROJ_A64FXd  
#define XP_RECON                       XP_RECON_A64FXd  
#define XM_RECON                       XM_RECON_A64FXd  
#define YM_RECON_ACCUM                 YM_RECON_ACCUM_A64FXd  
#define ZM_RECON_ACCUM                 ZM_RECON_ACCUM_A64FXd  
#define TM_RECON_ACCUM                 TM_RECON_ACCUM_A64FXd  
#define XP_RECON_ACCUM                 XP_RECON_ACCUM_A64FXd  
#define YP_RECON_ACCUM                 YP_RECON_ACCUM_A64FXd  
#define ZP_RECON_ACCUM                 ZP_RECON_ACCUM_A64FXd  
#define TP_RECON_ACCUM                 TP_RECON_ACCUM_A64FXd  
#define PERMUTE_DIR0                   PERM0_A64FXd  
#define PERMUTE_DIR1                   PERM1_A64FXd  
#define PERMUTE_DIR2                   PERM2_A64FXd  
#define PERMUTE_DIR3                   PERM3_A64FXd  
// DECLARATIONS
#define DECLARATIONS_A64FXd(x)  \
    const uint64_t lut[4][8] = { \
        {4, 5, 6, 7, 0, 1, 2, 3}, \
        {2, 3, 0, 1, 6, 7, 4, 5}, \
        {1, 0, 3, 2, 5, 4, 7, 6}, \
        {0, 1, 2, 4, 5, 6, 7, 8} };\
asm ( \
    "fmov z31.d , 0 \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// RESULT
#define RESULT_A64FXd(base)  \
{ \
asm ( \
    "stnt1d { z0.d }, p5, [%[storeptr], -6, mul vl] \n\t" \
    "stnt1d { z1.d }, p5, [%[storeptr], -5, mul vl] \n\t" \
    "stnt1d { z2.d }, p5, [%[storeptr], -4, mul vl] \n\t" \
    "stnt1d { z3.d }, p5, [%[storeptr], -3, mul vl] \n\t" \
    "stnt1d { z4.d }, p5, [%[storeptr], -2, mul vl] \n\t" \
    "stnt1d { z5.d }, p5, [%[storeptr], -1, mul vl] \n\t" \
    "stnt1d { z6.d }, p5, [%[storeptr], 0, mul vl] \n\t" \
    "stnt1d { z7.d }, p5, [%[storeptr], 1, mul vl] \n\t" \
    "stnt1d { z8.d }, p5, [%[storeptr], 2, mul vl] \n\t" \
    "stnt1d { z9.d }, p5, [%[storeptr], 3, mul vl] \n\t" \
    "stnt1d { z10.d }, p5, [%[storeptr], 4, mul vl] \n\t" \
    "stnt1d { z11.d }, p5, [%[storeptr], 5, mul vl] \n\t" \
    :  \
    : [storeptr] "r" (base + 2 * 3 * 64) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// PREFETCH_CHIMU_L2 (prefetch to L2)
#define PREFETCH_CHIMU_L2_INTERNAL_A64FXd(base)  \
{ \
asm ( \
    "prfd PLDL2STRM, p5, [%[fetchptr], 0, MUL VL] \n\t" \
    "prfd PLDL2STRM, p5, [%[fetchptr], 4, MUL VL] \n\t" \
    "prfd PLDL2STRM, p5, [%[fetchptr], 8, MUL VL] \n\t" \
    :  \
    : [fetchptr] "r" (base) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// PREFETCH_CHIMU_L1 (prefetch to L1)
#define PREFETCH_CHIMU_L1_INTERNAL_A64FXd(base)  \
{ \
asm ( \
    "prfd PLDL1STRM, p5, [%[fetchptr], 0, MUL VL] \n\t" \
    "prfd PLDL1STRM, p5, [%[fetchptr], 4, MUL VL] \n\t" \
    "prfd PLDL1STRM, p5, [%[fetchptr], 8, MUL VL] \n\t" \
    :  \
    : [fetchptr] "r" (base) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// PREFETCH_GAUGE_L2 (prefetch to L2)
#define PREFETCH_GAUGE_L2_INTERNAL_A64FXd(A)  \
{ \
    const auto & ref(U[sUn][A]); uint64_t baseU = (uint64_t)&ref[0][0]; \
asm ( \
    "prfd PLDL2STRM, p5, [%[fetchptr], 0, MUL VL] \n\t" \
    "prfd PLDL2STRM, p5, [%[fetchptr], 4, MUL VL] \n\t" \
    "prfd PLDL2STRM, p5, [%[fetchptr], 8, MUL VL] \n\t" \
    :  \
    : [fetchptr] "r" (baseU) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// PREFETCH_GAUGE_L1 (prefetch to L1)
#define PREFETCH_GAUGE_L1_INTERNAL(A)_A64FXd  \
{ \
    const auto & ref(U[sU][A]); uint64_t baseU = (uint64_t)&ref[0][0]; \
asm ( \
    "prfd PLDL1STRM, p5, [%[fetchptr], 0, MUL VL] \n\t" \
    "prfd PLDL1STRM, p5, [%[fetchptr], 4, MUL VL] \n\t" \
    "prfd PLDL1STRM, p5, [%[fetchptr], 8, MUL VL] \n\t" \
    :  \
    : [fetchptr] "r" (baseU) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// LOAD_CHI
#define LOAD_CHI_A64FXd(base)  \
{ \
asm ( \
    "ld1d { z12.d }, p5/z, [%[fetchptr], -6, mul vl] \n\t" \
    "ld1d { z13.d }, p5/z, [%[fetchptr], -5, mul vl] \n\t" \
    "ld1d { z14.d }, p5/z, [%[fetchptr], -4, mul vl] \n\t" \
    "ld1d { z15.d }, p5/z, [%[fetchptr], -3, mul vl] \n\t" \
    "ld1d { z16.d }, p5/z, [%[fetchptr], -2, mul vl] \n\t" \
    "ld1d { z17.d }, p5/z, [%[fetchptr], -1, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (base + 2 * 3 * 64) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// LOAD_CHIMU
#define LOAD_CHIMU_INTERLEAVED_A64FXd(base)  \
{ \
asm ( \
    "ptrue p5.d \n\t" \
    "ld1d { z12.d }, p5/z, [%[fetchptr], -6, mul vl] \n\t" \
    "ld1d { z27.d }, p5/z, [%[fetchptr], 3, mul vl] \n\t" \
    "ld1d { z15.d }, p5/z, [%[fetchptr], -3, mul vl] \n\t" \
    "ld1d { z24.d }, p5/z, [%[fetchptr], 0, mul vl] \n\t" \
    "ld1d { z13.d }, p5/z, [%[fetchptr], -5, mul vl] \n\t" \
    "ld1d { z28.d }, p5/z, [%[fetchptr], 4, mul vl] \n\t" \
    "ld1d { z16.d }, p5/z, [%[fetchptr], -2, mul vl] \n\t" \
    "ld1d { z25.d }, p5/z, [%[fetchptr], 1, mul vl] \n\t" \
    "ld1d { z14.d }, p5/z, [%[fetchptr], -4, mul vl] \n\t" \
    "ld1d { z29.d }, p5/z, [%[fetchptr], 5, mul vl] \n\t" \
    "ld1d { z17.d }, p5/z, [%[fetchptr], -1, mul vl] \n\t" \
    "ld1d { z26.d }, p5/z, [%[fetchptr], 2, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (base + 2 * 3 * 64) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// LOAD_CHIMU_0213
#define LOAD_CHIMU_0213_A64FXd  \
{ \
    const SiteSpinor & ref(in[offset]); \
asm ( \
    "ptrue p5.d \n\t" \
    "ld1d { z12.d }, p5/z, [%[fetchptr], -6, mul vl] \n\t" \
    "ld1d { z24.d }, p5/z, [%[fetchptr], 0, mul vl] \n\t" \
    "ld1d { z13.d }, p5/z, [%[fetchptr], -5, mul vl] \n\t" \
    "ld1d { z25.d }, p5/z, [%[fetchptr], 1, mul vl] \n\t" \
    "ld1d { z14.d }, p5/z, [%[fetchptr], -4, mul vl] \n\t" \
    "ld1d { z26.d }, p5/z, [%[fetchptr], 2, mul vl] \n\t" \
    "ld1d { z15.d }, p5/z, [%[fetchptr], -3, mul vl] \n\t" \
    "ld1d { z27.d }, p5/z, [%[fetchptr], 3, mul vl] \n\t" \
    "ld1d { z16.d }, p5/z, [%[fetchptr], -2, mul vl] \n\t" \
    "ld1d { z28.d }, p5/z, [%[fetchptr], 4, mul vl] \n\t" \
    "ld1d { z17.d }, p5/z, [%[fetchptr], -1, mul vl] \n\t" \
    "ld1d { z29.d }, p5/z, [%[fetchptr], 5, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (&ref[2][0]) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// LOAD_CHIMU_0312
#define LOAD_CHIMU_0312_A64FXd  \
{ \
    const SiteSpinor & ref(in[offset]); \
asm ( \
    "ptrue p5.d \n\t" \
    "ld1d { z12.d }, p5/z, [%[fetchptr], -6, mul vl] \n\t" \
    "ld1d { z27.d }, p5/z, [%[fetchptr], 3, mul vl] \n\t" \
    "ld1d { z13.d }, p5/z, [%[fetchptr], -5, mul vl] \n\t" \
    "ld1d { z28.d }, p5/z, [%[fetchptr], 4, mul vl] \n\t" \
    "ld1d { z14.d }, p5/z, [%[fetchptr], -4, mul vl] \n\t" \
    "ld1d { z29.d }, p5/z, [%[fetchptr], 5, mul vl] \n\t" \
    "ld1d { z15.d }, p5/z, [%[fetchptr], -3, mul vl] \n\t" \
    "ld1d { z24.d }, p5/z, [%[fetchptr], 0, mul vl] \n\t" \
    "ld1d { z16.d }, p5/z, [%[fetchptr], -2, mul vl] \n\t" \
    "ld1d { z25.d }, p5/z, [%[fetchptr], 1, mul vl] \n\t" \
    "ld1d { z17.d }, p5/z, [%[fetchptr], -1, mul vl] \n\t" \
    "ld1d { z26.d }, p5/z, [%[fetchptr], 2, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (&ref[2][0]) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// PERM0
#define PERM0_A64FXd  \
asm ( \
    "ld1d { z30.d }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "tbl z12.d, { z12.d }, z30.d \n\t"  \
    "tbl z13.d, { z13.d }, z30.d \n\t"  \
    "tbl z14.d, { z14.d }, z30.d \n\t"  \
    "tbl z15.d, { z15.d }, z30.d \n\t"  \
    "tbl z16.d, { z16.d }, z30.d \n\t"  \
    "tbl z17.d, { z17.d }, z30.d \n\t"  \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (0) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// PERM1
#define PERM1_A64FXd  \
asm ( \
    "ld1d { z30.d }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "tbl z12.d, { z12.d }, z30.d \n\t"  \
    "tbl z13.d, { z13.d }, z30.d \n\t"  \
    "tbl z14.d, { z14.d }, z30.d \n\t"  \
    "tbl z15.d, { z15.d }, z30.d \n\t"  \
    "tbl z16.d, { z16.d }, z30.d \n\t"  \
    "tbl z17.d, { z17.d }, z30.d \n\t"  \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (1) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// PERM2
#define PERM2_A64FXd  \
asm ( \
    "ld1d { z30.d }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "tbl z12.d, { z12.d }, z30.d \n\t"  \
    "tbl z13.d, { z13.d }, z30.d \n\t"  \
    "tbl z14.d, { z14.d }, z30.d \n\t"  \
    "tbl z15.d, { z15.d }, z30.d \n\t"  \
    "tbl z16.d, { z16.d }, z30.d \n\t"  \
    "tbl z17.d, { z17.d }, z30.d \n\t"  \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (2) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// PERM3
#define PERM3_A64FXd  

// MULT_2SPIN
#define MULT_2SPIN_A64FXd(A)  \
{ \
    const auto & ref(U[sU][A]); \
asm ( \
    "ld1d { z24.d }, p5/z, [%[fetchptr], -6, mul vl] \n\t" \
    "ld1d { z25.d }, p5/z, [%[fetchptr], -3, mul vl] \n\t" \
    "ld1d { z26.d }, p5/z, [%[fetchptr], 0, mul vl] \n\t" \
    "ld1d { z27.d }, p5/z, [%[fetchptr], -5, mul vl] \n\t" \
    "ld1d { z28.d }, p5/z, [%[fetchptr], -2, mul vl] \n\t" \
    "ld1d { z29.d }, p5/z, [%[fetchptr], 1, mul vl] \n\t" \
    "fmov z18.d , 0 \n\t" \
    "fmov z21.d , 0 \n\t" \
    "fmov z19.d , 0 \n\t" \
    "fmov z22.d , 0 \n\t" \
    "fmov z20.d , 0 \n\t" \
    "fmov z23.d , 0 \n\t" \
    "fcmla z18.d, p5/m, z24.d, z12.d, 0 \n\t" \
    "fcmla z21.d, p5/m, z24.d, z15.d, 0 \n\t" \
    "fcmla z19.d, p5/m, z25.d, z12.d, 0 \n\t" \
    "fcmla z22.d, p5/m, z25.d, z15.d, 0 \n\t" \
    "fcmla z20.d, p5/m, z26.d, z12.d, 0 \n\t" \
    "fcmla z23.d, p5/m, z26.d, z15.d, 0 \n\t" \
    "fcmla z18.d, p5/m, z24.d, z12.d, 90 \n\t" \
    "fcmla z21.d, p5/m, z24.d, z15.d, 90 \n\t" \
    "fcmla z19.d, p5/m, z25.d, z12.d, 90 \n\t" \
    "fcmla z22.d, p5/m, z25.d, z15.d, 90 \n\t" \
    "fcmla z20.d, p5/m, z26.d, z12.d, 90 \n\t" \
    "fcmla z23.d, p5/m, z26.d, z15.d, 90 \n\t" \
    "ld1d { z24.d }, p5/z, [%[fetchptr], -4, mul vl] \n\t" \
    "ld1d { z25.d }, p5/z, [%[fetchptr], -1, mul vl] \n\t" \
    "ld1d { z26.d }, p5/z, [%[fetchptr], 2, mul vl] \n\t" \
    "fcmla z18.d, p5/m, z27.d, z13.d, 0 \n\t" \
    "fcmla z21.d, p5/m, z27.d, z16.d, 0 \n\t" \
    "fcmla z19.d, p5/m, z28.d, z13.d, 0 \n\t" \
    "fcmla z22.d, p5/m, z28.d, z16.d, 0 \n\t" \
    "fcmla z20.d, p5/m, z29.d, z13.d, 0 \n\t" \
    "fcmla z23.d, p5/m, z29.d, z16.d, 0 \n\t" \
    "fcmla z18.d, p5/m, z27.d, z13.d, 90 \n\t" \
    "fcmla z21.d, p5/m, z27.d, z16.d, 90 \n\t" \
    "fcmla z19.d, p5/m, z28.d, z13.d, 90 \n\t" \
    "fcmla z22.d, p5/m, z28.d, z16.d, 90 \n\t" \
    "fcmla z20.d, p5/m, z29.d, z13.d, 90 \n\t" \
    "fcmla z23.d, p5/m, z29.d, z16.d, 90 \n\t" \
    "fcmla z18.d, p5/m, z24.d, z14.d, 0 \n\t" \
    "fcmla z21.d, p5/m, z24.d, z17.d, 0 \n\t" \
    "fcmla z19.d, p5/m, z25.d, z14.d, 0 \n\t" \
    "fcmla z22.d, p5/m, z25.d, z17.d, 0 \n\t" \
    "fcmla z20.d, p5/m, z26.d, z14.d, 0 \n\t" \
    "fcmla z23.d, p5/m, z26.d, z17.d, 0 \n\t" \
    "fcmla z18.d, p5/m, z24.d, z14.d, 90 \n\t" \
    "fcmla z21.d, p5/m, z24.d, z17.d, 90 \n\t" \
    "fcmla z19.d, p5/m, z25.d, z14.d, 90 \n\t" \
    "fcmla z22.d, p5/m, z25.d, z17.d, 90 \n\t" \
    "fcmla z20.d, p5/m, z26.d, z14.d, 90 \n\t" \
    "fcmla z23.d, p5/m, z26.d, z17.d, 90 \n\t" \
    :  \
    : [fetchptr] "r" ((uint64_t)&ref[2][0]) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// XP_PROJ
#define XP_PROJ_A64FXd  \
{ \
asm ( \
    "fcadd z12.d, p5/m, z12.d, z27.d, 90 \n\t" \
    "fcadd z13.d, p5/m, z13.d, z28.d, 90 \n\t" \
    "fcadd z14.d, p5/m, z14.d, z29.d, 90 \n\t" \
    "fcadd z15.d, p5/m, z15.d, z24.d, 90 \n\t" \
    "fcadd z16.d, p5/m, z16.d, z25.d, 90 \n\t" \
    "fcadd z17.d, p5/m, z17.d, z26.d, 90 \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// XP_RECON
#define XP_RECON_A64FXd  \
asm ( \
    "fcadd z6.d, p5/m, z6.d, z21.d, 270 \n\t" \
    "fcadd z7.d, p5/m, z7.d, z22.d, 270 \n\t" \
    "fcadd z8.d, p5/m, z8.d, z23.d, 270 \n\t" \
    "fcadd z9.d, p5/m, z9.d, z18.d, 270 \n\t" \
    "fcadd z10.d, p5/m, z10.d, z19.d, 270 \n\t" \
    "fcadd z11.d, p5/m, z11.d, z20.d, 270 \n\t" \
    "mov z0.d, z18.d \n\t" \
    "mov z1.d, z19.d \n\t" \
    "mov z2.d, z20.d \n\t" \
    "mov z3.d, z21.d \n\t" \
    "mov z4.d, z22.d \n\t" \
    "mov z5.d, z23.d \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// XP_RECON_ACCUM
#define XP_RECON_ACCUM_A64FXd  \
asm ( \
    "fcadd z9.d, p5/m, z9.d, z18.d, 270 \n\t" \
    "fadd z0.d, p5/m, z0.d, z18.d \n\t"  \
    "fcadd z10.d, p5/m, z10.d, z19.d, 270 \n\t" \
    "fadd z1.d, p5/m, z1.d, z19.d \n\t"  \
    "fcadd z11.d, p5/m, z11.d, z20.d, 270 \n\t" \
    "fadd z2.d, p5/m, z2.d, z20.d \n\t"  \
    "fcadd z6.d, p5/m, z6.d, z21.d, 270 \n\t" \
    "fadd z3.d, p5/m, z3.d, z21.d \n\t"  \
    "fcadd z7.d, p5/m, z7.d, z22.d, 270 \n\t" \
    "fadd z4.d, p5/m, z4.d, z22.d \n\t"  \
    "fcadd z8.d, p5/m, z8.d, z23.d, 270 \n\t" \
    "fadd z5.d, p5/m, z5.d, z23.d \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// YP_PROJ
#define YP_PROJ_A64FXd  \
{ \
asm ( \
    "ld1d { z30.d }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "fsub z12.d, p5/m, z12.d, z27.d \n\t" \
    "fsub z13.d, p5/m, z13.d, z28.d \n\t" \
    "fsub z14.d, p5/m, z14.d, z29.d \n\t" \
    "fadd z15.d, p5/m, z15.d, z24.d \n\t"  \
    "fadd z16.d, p5/m, z16.d, z25.d \n\t"  \
    "fadd z17.d, p5/m, z17.d, z26.d \n\t"  \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (2) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// ZP_PROJ
#define ZP_PROJ_A64FXd  \
{ \
asm ( \
    "ld1d { z30.d }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "fcadd z12.d, p5/m, z12.d, z24.d, 90 \n\t" \
    "fcadd z13.d, p5/m, z13.d, z25.d, 90 \n\t" \
    "fcadd z14.d, p5/m, z14.d, z26.d, 90 \n\t" \
    "fcadd z15.d, p5/m, z15.d, z27.d, 270 \n\t" \
    "fcadd z16.d, p5/m, z16.d, z28.d, 270 \n\t" \
    "fcadd z17.d, p5/m, z17.d, z29.d, 270 \n\t" \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (1) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// TP_PROJ
#define TP_PROJ_A64FXd  \
{ \
asm ( \
    "ld1d { z30.d }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "fadd z12.d, p5/m, z12.d, z24.d \n\t"  \
    "fadd z13.d, p5/m, z13.d, z25.d \n\t"  \
    "fadd z14.d, p5/m, z14.d, z26.d \n\t"  \
    "fadd z15.d, p5/m, z15.d, z27.d \n\t"  \
    "fadd z16.d, p5/m, z16.d, z28.d \n\t"  \
    "fadd z17.d, p5/m, z17.d, z29.d \n\t"  \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (0) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// XM_PROJ
#define XM_PROJ_A64FXd  \
{ \
asm ( \
    "fcadd z12.d, p5/m, z12.d, z27.d, 270 \n\t" \
    "fcadd z13.d, p5/m, z13.d, z28.d, 270 \n\t" \
    "fcadd z14.d, p5/m, z14.d, z29.d, 270 \n\t" \
    "fcadd z15.d, p5/m, z15.d, z24.d, 270 \n\t" \
    "fcadd z16.d, p5/m, z16.d, z25.d, 270 \n\t" \
    "fcadd z17.d, p5/m, z17.d, z26.d, 270 \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// XM_RECON
#define XM_RECON_A64FXd  \
asm ( \
    "fcadd z6.d, p5/m, z6.d, z21.d, 90 \n\t" \
    "fcadd z7.d, p5/m, z7.d, z22.d, 90 \n\t" \
    "fcadd z8.d, p5/m, z8.d, z23.d, 90 \n\t" \
    "fcadd z9.d, p5/m, z9.d, z18.d, 90 \n\t" \
    "fcadd z10.d, p5/m, z10.d, z19.d, 90 \n\t" \
    "fcadd z11.d, p5/m, z11.d, z20.d, 90 \n\t" \
    "mov z0.d, z18.d \n\t" \
    "mov z1.d, z19.d \n\t" \
    "mov z2.d, z20.d \n\t" \
    "mov z3.d, z21.d \n\t" \
    "mov z4.d, z22.d \n\t" \
    "mov z5.d, z23.d \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// YM_PROJ
#define YM_PROJ_A64FXd  \
{ \
asm ( \
    "ld1d { z30.d }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "fadd z12.d, p5/m, z12.d, z27.d \n\t"  \
    "fadd z13.d, p5/m, z13.d, z28.d \n\t"  \
    "fadd z14.d, p5/m, z14.d, z29.d \n\t"  \
    "fsub z15.d, p5/m, z15.d, z24.d \n\t" \
    "fsub z16.d, p5/m, z16.d, z25.d \n\t" \
    "fsub z17.d, p5/m, z17.d, z26.d \n\t" \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (2) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// ZM_PROJ
#define ZM_PROJ_A64FXd  \
{ \
asm ( \
    "ld1d { z30.d }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "fcadd z12.d, p5/m, z12.d, z24.d, 270 \n\t" \
    "fcadd z13.d, p5/m, z13.d, z25.d, 270 \n\t" \
    "fcadd z14.d, p5/m, z14.d, z26.d, 270 \n\t" \
    "fcadd z15.d, p5/m, z15.d, z27.d, 90 \n\t" \
    "fcadd z16.d, p5/m, z16.d, z28.d, 90 \n\t" \
    "fcadd z17.d, p5/m, z17.d, z29.d, 90 \n\t" \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (1) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// TM_PROJ
#define TM_PROJ_A64FXd  \
{ \
asm ( \
    "ld1d { z30.d }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "fsub z12.d, p5/m, z12.d, z24.d \n\t" \
    "fsub z13.d, p5/m, z13.d, z25.d \n\t" \
    "fsub z14.d, p5/m, z14.d, z26.d \n\t" \
    "fsub z15.d, p5/m, z15.d, z27.d \n\t" \
    "fsub z16.d, p5/m, z16.d, z28.d \n\t" \
    "fsub z17.d, p5/m, z17.d, z29.d \n\t" \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (0) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// XM_RECON_ACCUM
#define XM_RECON_ACCUM_A64FXd  \
asm ( \
    "fcadd z9.d, p5/m, z9.d, z18.d, 90 \n\t" \
    "fcadd z10.d, p5/m, z10.d, z19.d, 90 \n\t" \
    "fcadd z11.d, p5/m, z11.d, z20.d, 90 \n\t" \
    "fcadd z6.d, p5/m, z6.d, z21.d, 90 \n\t" \
    "fcadd z7.d, p5/m, z7.d, z22.d, 90 \n\t" \
    "fcadd z8.d, p5/m, z8.d, z23.d, 90 \n\t" \
    "mov z0.d, z18.d \n\t" \
    "mov z1.d, z19.d \n\t" \
    "mov z2.d, z20.d \n\t" \
    "mov z3.d, z21.d \n\t" \
    "mov z4.d, z22.d \n\t" \
    "mov z5.d, z23.d \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// YP_RECON_ACCUM
#define YP_RECON_ACCUM_A64FXd  \
asm ( \
    "fadd z0.d, p5/m, z0.d, z18.d \n\t"  \
    "fsub z9.d, p5/m, z9.d, z18.d \n\t" \
    "fadd z1.d, p5/m, z1.d, z19.d \n\t"  \
    "fsub z10.d, p5/m, z10.d, z19.d \n\t" \
    "fadd z2.d, p5/m, z2.d, z20.d \n\t"  \
    "fsub z11.d, p5/m, z11.d, z20.d \n\t" \
    "fadd z3.d, p5/m, z3.d, z21.d \n\t"  \
    "fadd z6.d, p5/m, z6.d, z21.d \n\t"  \
    "fadd z4.d, p5/m, z4.d, z22.d \n\t"  \
    "fadd z7.d, p5/m, z7.d, z22.d \n\t"  \
    "fadd z5.d, p5/m, z5.d, z23.d \n\t"  \
    "fadd z8.d, p5/m, z8.d, z23.d \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// YM_RECON_ACCUM
#define YM_RECON_ACCUM_A64FXd  \
asm ( \
    "fadd z0.d, p5/m, z0.d, z18.d \n\t"  \
    "fadd z9.d, p5/m, z9.d, z18.d \n\t"  \
    "fadd z1.d, p5/m, z1.d, z19.d \n\t"  \
    "fadd z10.d, p5/m, z10.d, z19.d \n\t"  \
    "fadd z2.d, p5/m, z2.d, z20.d \n\t"  \
    "fadd z11.d, p5/m, z11.d, z20.d \n\t"  \
    "fadd z3.d, p5/m, z3.d, z21.d \n\t"  \
    "fsub z6.d, p5/m, z6.d, z21.d \n\t" \
    "fadd z4.d, p5/m, z4.d, z22.d \n\t"  \
    "fsub z7.d, p5/m, z7.d, z22.d \n\t" \
    "fadd z5.d, p5/m, z5.d, z23.d \n\t"  \
    "fsub z8.d, p5/m, z8.d, z23.d \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// ZP_RECON_ACCUM
#define ZP_RECON_ACCUM_A64FXd  \
asm ( \
    "fcadd z6.d, p5/m, z6.d, z18.d, 270 \n\t" \
    "fadd z0.d, p5/m, z0.d, z18.d \n\t"  \
    "fcadd z7.d, p5/m, z7.d, z19.d, 270 \n\t" \
    "fadd z1.d, p5/m, z1.d, z19.d \n\t"  \
    "fcadd z8.d, p5/m, z8.d, z20.d, 270 \n\t" \
    "fadd z2.d, p5/m, z2.d, z20.d \n\t"  \
    "fcadd z9.d, p5/m, z9.d, z21.d, 90 \n\t" \
    "fadd z3.d, p5/m, z3.d, z21.d \n\t"  \
    "fcadd z10.d, p5/m, z10.d, z22.d, 90 \n\t" \
    "fadd z4.d, p5/m, z4.d, z22.d \n\t"  \
    "fcadd z11.d, p5/m, z11.d, z23.d, 90 \n\t" \
    "fadd z5.d, p5/m, z5.d, z23.d \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// ZM_RECON_ACCUM
#define ZM_RECON_ACCUM_A64FXd  \
asm ( \
    "fcadd z6.d, p5/m, z6.d, z18.d, 90 \n\t" \
    "fadd z0.d, p5/m, z0.d, z18.d \n\t"  \
    "fcadd z7.d, p5/m, z7.d, z19.d, 90 \n\t" \
    "fadd z1.d, p5/m, z1.d, z19.d \n\t"  \
    "fcadd z8.d, p5/m, z8.d, z20.d, 90 \n\t" \
    "fadd z2.d, p5/m, z2.d, z20.d \n\t"  \
    "fcadd z9.d, p5/m, z9.d, z21.d, 270 \n\t" \
    "fadd z3.d, p5/m, z3.d, z21.d \n\t"  \
    "fcadd z10.d, p5/m, z10.d, z22.d, 270 \n\t" \
    "fadd z4.d, p5/m, z4.d, z22.d \n\t"  \
    "fcadd z11.d, p5/m, z11.d, z23.d, 270 \n\t" \
    "fadd z5.d, p5/m, z5.d, z23.d \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// TP_RECON_ACCUM
#define TP_RECON_ACCUM_A64FXd  \
asm ( \
    "ptrue p5.d \n\t" \
    "fadd z0.d, p5/m, z0.d, z18.d \n\t"  \
    "fadd z6.d, p5/m, z6.d, z18.d \n\t"  \
    "fadd z1.d, p5/m, z1.d, z19.d \n\t"  \
    "fadd z7.d, p5/m, z7.d, z19.d \n\t"  \
    "fadd z2.d, p5/m, z2.d, z20.d \n\t"  \
    "fadd z8.d, p5/m, z8.d, z20.d \n\t"  \
    "fadd z3.d, p5/m, z3.d, z21.d \n\t"  \
    "fadd z9.d, p5/m, z9.d, z21.d \n\t"  \
    "fadd z4.d, p5/m, z4.d, z22.d \n\t"  \
    "fadd z10.d, p5/m, z10.d, z22.d \n\t"  \
    "fadd z5.d, p5/m, z5.d, z23.d \n\t"  \
    "fadd z11.d, p5/m, z11.d, z23.d \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// TM_RECON_ACCUM
#define TM_RECON_ACCUM_A64FXd  \
asm ( \
    "fadd z0.d, p5/m, z0.d, z18.d \n\t"  \
    "fsub z6.d, p5/m, z6.d, z18.d \n\t" \
    "fadd z1.d, p5/m, z1.d, z19.d \n\t"  \
    "fsub z7.d, p5/m, z7.d, z19.d \n\t" \
    "fadd z2.d, p5/m, z2.d, z20.d \n\t"  \
    "fsub z8.d, p5/m, z8.d, z20.d \n\t" \
    "fadd z3.d, p5/m, z3.d, z21.d \n\t"  \
    "fsub z9.d, p5/m, z9.d, z21.d \n\t" \
    "fadd z4.d, p5/m, z4.d, z22.d \n\t"  \
    "fsub z10.d, p5/m, z10.d, z22.d \n\t" \
    "fadd z5.d, p5/m, z5.d, z23.d \n\t"  \
    "fsub z11.d, p5/m, z11.d, z23.d \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// ZERO_PSI
#define ZERO_PSI_A64FXd  \
asm ( \
    "ptrue p5.d \n\t" \
    "fmov z0.d , 0 \n\t" \
    "fmov z1.d , 0 \n\t" \
    "fmov z2.d , 0 \n\t" \
    "fmov z3.d , 0 \n\t" \
    "fmov z4.d , 0 \n\t" \
    "fmov z5.d , 0 \n\t" \
    "fmov z6.d , 0 \n\t" \
    "fmov z7.d , 0 \n\t" \
    "fmov z8.d , 0 \n\t" \
    "fmov z9.d , 0 \n\t" \
    "fmov z10.d , 0 \n\t" \
    "fmov z11.d , 0 \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

