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
#define LOAD_CHIMU_A64FXf(x)           LOAD_CHIMU_INTERLEAVED_A64FXf(x)  
#define PREFETCH_CHIMU_L1(A)  
#define PREFETCH_GAUGE_L1(A)  
#define PREFETCH_CHIMU_L2(A)  
#define PREFETCH_GAUGE_L2(A)  
#define PF_GAUGE(A)  
#define PREFETCH1_CHIMU(A)  
#define PREFETCH_CHIMU(A)  
#define LOCK_GAUGE(A)  
#define UNLOCK_GAUGE(A)  
#define MASK_REGS                      DECLARATIONS_A64FXf  
#define COMPLEX_SIGNS(A)  
#define LOAD64(A,B)  
#define SAVE_RESULT(A,B)               RESULT_A64FXf(A)  
#define MULT_2SPIN_DIR_PF(A,B)         MULT_2SPIN_A64FXf(A)  
#define MAYBEPERM(A,perm)              if (perm) { A ; }  
#define LOAD_CHI(base)                 LOAD_CHI_A64FXf(base)  
#define ZERO_PSI                       ZERO_PSI_A64FXf  
#define XP_PROJMEM(base)               LOAD_CHIMU_A64FXf(base);   XP_PROJ_A64FXf  
#define YP_PROJMEM(base)               LOAD_CHIMU_A64FXf(base);   YP_PROJ_A64FXf  
#define ZP_PROJMEM(base)               LOAD_CHIMU_A64FXf(base);   ZP_PROJ_A64FXf  
#define TP_PROJMEM(base)               LOAD_CHIMU_A64FXf(base);   TP_PROJ_A64FXf  
#define XM_PROJMEM(base)               LOAD_CHIMU_A64FXf(base);   XM_PROJ_A64FXf  
#define YM_PROJMEM(base)               LOAD_CHIMU_A64FXf(base);   YM_PROJ_A64FXf  
#define ZM_PROJMEM(base)               LOAD_CHIMU_A64FXf(base);   ZM_PROJ_A64FXf  
#define TM_PROJMEM(base)               LOAD_CHIMU_A64FXf(base);   TM_PROJ_A64FXf  
#define XP_RECON                       XP_RECON_A64FXf  
#define XM_RECON                       XM_RECON_A64FXf  
#define YM_RECON_ACCUM                 YM_RECON_ACCUM_A64FXf  
#define ZM_RECON_ACCUM                 ZM_RECON_ACCUM_A64FXf  
#define TM_RECON_ACCUM                 TM_RECON_ACCUM_A64FXf  
#define XP_RECON_ACCUM                 XP_RECON_ACCUM_A64FXf  
#define YP_RECON_ACCUM                 YP_RECON_ACCUM_A64FXf  
#define ZP_RECON_ACCUM                 ZP_RECON_ACCUM_A64FXf  
#define TP_RECON_ACCUM                 TP_RECON_ACCUM_A64FXf  
#define PERMUTE_DIR0                   PERM0_A64FXf  
#define PERMUTE_DIR1                   PERM1_A64FXf  
#define PERMUTE_DIR2                   PERM2_A64FXf  
#define PERMUTE_DIR3                   PERM3_A64FXf  
// DECLARATIONS
#define DECLARATIONS_A64FXf  \
    const uint32_t lut[4][16] = { \
        {8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7}, \
        {4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11}, \
        {2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13}, \
        {1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14} }; \
asm ( \
    "fmov z31.f , 0 \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// RESULT
#define RESULT_A64FXf(base)  \
{ \
asm ( \
    "stnt1d { z0.f }, p5, [%[storeptr], -6, mul vl] \n\t" \
    "stnt1d { z1.f }, p5, [%[storeptr], -5, mul vl] \n\t" \
    "stnt1d { z2.f }, p5, [%[storeptr], -4, mul vl] \n\t" \
    "stnt1d { z3.f }, p5, [%[storeptr], -3, mul vl] \n\t" \
    "stnt1d { z4.f }, p5, [%[storeptr], -2, mul vl] \n\t" \
    "stnt1d { z5.f }, p5, [%[storeptr], -1, mul vl] \n\t" \
    "stnt1d { z6.f }, p5, [%[storeptr], 0, mul vl] \n\t" \
    "stnt1d { z7.f }, p5, [%[storeptr], 1, mul vl] \n\t" \
    "stnt1d { z8.f }, p5, [%[storeptr], 2, mul vl] \n\t" \
    "stnt1d { z9.f }, p5, [%[storeptr], 3, mul vl] \n\t" \
    "stnt1d { z10.f }, p5, [%[storeptr], 4, mul vl] \n\t" \
    "stnt1d { z11.f }, p5, [%[storeptr], 5, mul vl] \n\t" \
    :  \
    : [storeptr] "r" (base + 2 * 3 * 64) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// PREFETCH_CHIMU_L2 (prefetch to L2)
#define PREFETCH_CHIMU_L2_INTERNAL_A64FXf(base)  \
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
#define PREFETCH_CHIMU_L1_INTERNAL_A64FXf(base)  \
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
#define PREFETCH_GAUGE_L2_INTERNAL_A64FXf(A)  \
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
#define PREFETCH_GAUGE_L1_INTERNAL(A)_A64FXf  \
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
#define LOAD_CHI_A64FXf(base)  \
{ \
asm ( \
    "ld1d { z12.f }, p5/z, [%[fetchptr], -6, mul vl] \n\t" \
    "ld1d { z13.f }, p5/z, [%[fetchptr], -5, mul vl] \n\t" \
    "ld1d { z14.f }, p5/z, [%[fetchptr], -4, mul vl] \n\t" \
    "ld1d { z15.f }, p5/z, [%[fetchptr], -3, mul vl] \n\t" \
    "ld1d { z16.f }, p5/z, [%[fetchptr], -2, mul vl] \n\t" \
    "ld1d { z17.f }, p5/z, [%[fetchptr], -1, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (base + 2 * 3 * 64) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// LOAD_CHIMU
#define LOAD_CHIMU_INTERLEAVED_A64FXf(base)  \
{ \
asm ( \
    "ptrue p5.f \n\t" \
    "ld1d { z12.f }, p5/z, [%[fetchptr], -6, mul vl] \n\t" \
    "ld1d { z27.f }, p5/z, [%[fetchptr], 3, mul vl] \n\t" \
    "ld1d { z15.f }, p5/z, [%[fetchptr], -3, mul vl] \n\t" \
    "ld1d { z24.f }, p5/z, [%[fetchptr], 0, mul vl] \n\t" \
    "ld1d { z13.f }, p5/z, [%[fetchptr], -5, mul vl] \n\t" \
    "ld1d { z28.f }, p5/z, [%[fetchptr], 4, mul vl] \n\t" \
    "ld1d { z16.f }, p5/z, [%[fetchptr], -2, mul vl] \n\t" \
    "ld1d { z25.f }, p5/z, [%[fetchptr], 1, mul vl] \n\t" \
    "ld1d { z14.f }, p5/z, [%[fetchptr], -4, mul vl] \n\t" \
    "ld1d { z29.f }, p5/z, [%[fetchptr], 5, mul vl] \n\t" \
    "ld1d { z17.f }, p5/z, [%[fetchptr], -1, mul vl] \n\t" \
    "ld1d { z26.f }, p5/z, [%[fetchptr], 2, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (base + 2 * 3 * 64) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// LOAD_CHIMU_0213
#define LOAD_CHIMU_0213_A64FXf  \
{ \
    const SiteSpinor & ref(in[offset]); \
asm ( \
    "ptrue p5.f \n\t" \
    "ld1d { z12.f }, p5/z, [%[fetchptr], -6, mul vl] \n\t" \
    "ld1d { z24.f }, p5/z, [%[fetchptr], 0, mul vl] \n\t" \
    "ld1d { z13.f }, p5/z, [%[fetchptr], -5, mul vl] \n\t" \
    "ld1d { z25.f }, p5/z, [%[fetchptr], 1, mul vl] \n\t" \
    "ld1d { z14.f }, p5/z, [%[fetchptr], -4, mul vl] \n\t" \
    "ld1d { z26.f }, p5/z, [%[fetchptr], 2, mul vl] \n\t" \
    "ld1d { z15.f }, p5/z, [%[fetchptr], -3, mul vl] \n\t" \
    "ld1d { z27.f }, p5/z, [%[fetchptr], 3, mul vl] \n\t" \
    "ld1d { z16.f }, p5/z, [%[fetchptr], -2, mul vl] \n\t" \
    "ld1d { z28.f }, p5/z, [%[fetchptr], 4, mul vl] \n\t" \
    "ld1d { z17.f }, p5/z, [%[fetchptr], -1, mul vl] \n\t" \
    "ld1d { z29.f }, p5/z, [%[fetchptr], 5, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (&ref[2][0]) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// LOAD_CHIMU_0312
#define LOAD_CHIMU_0312_A64FXf  \
{ \
    const SiteSpinor & ref(in[offset]); \
asm ( \
    "ptrue p5.f \n\t" \
    "ld1d { z12.f }, p5/z, [%[fetchptr], -6, mul vl] \n\t" \
    "ld1d { z27.f }, p5/z, [%[fetchptr], 3, mul vl] \n\t" \
    "ld1d { z13.f }, p5/z, [%[fetchptr], -5, mul vl] \n\t" \
    "ld1d { z28.f }, p5/z, [%[fetchptr], 4, mul vl] \n\t" \
    "ld1d { z14.f }, p5/z, [%[fetchptr], -4, mul vl] \n\t" \
    "ld1d { z29.f }, p5/z, [%[fetchptr], 5, mul vl] \n\t" \
    "ld1d { z15.f }, p5/z, [%[fetchptr], -3, mul vl] \n\t" \
    "ld1d { z24.f }, p5/z, [%[fetchptr], 0, mul vl] \n\t" \
    "ld1d { z16.f }, p5/z, [%[fetchptr], -2, mul vl] \n\t" \
    "ld1d { z25.f }, p5/z, [%[fetchptr], 1, mul vl] \n\t" \
    "ld1d { z17.f }, p5/z, [%[fetchptr], -1, mul vl] \n\t" \
    "ld1d { z26.f }, p5/z, [%[fetchptr], 2, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (&ref[2][0]) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// PERM0
#define PERM0_A64FXf  \
asm ( \
    "ld1d { z30.f }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "tbl z12.f, { z12.f }, z30.f \n\t"  \
    "tbl z13.f, { z13.f }, z30.f \n\t"  \
    "tbl z14.f, { z14.f }, z30.f \n\t"  \
    "tbl z15.f, { z15.f }, z30.f \n\t"  \
    "tbl z16.f, { z16.f }, z30.f \n\t"  \
    "tbl z17.f, { z17.f }, z30.f \n\t"  \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (0) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// PERM1
#define PERM1_A64FXf  \
asm ( \
    "ld1d { z30.f }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "tbl z12.f, { z12.f }, z30.f \n\t"  \
    "tbl z13.f, { z13.f }, z30.f \n\t"  \
    "tbl z14.f, { z14.f }, z30.f \n\t"  \
    "tbl z15.f, { z15.f }, z30.f \n\t"  \
    "tbl z16.f, { z16.f }, z30.f \n\t"  \
    "tbl z17.f, { z17.f }, z30.f \n\t"  \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (1) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// PERM2
#define PERM2_A64FXf  \
asm ( \
    "ld1d { z30.f }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "tbl z12.f, { z12.f }, z30.f \n\t"  \
    "tbl z13.f, { z13.f }, z30.f \n\t"  \
    "tbl z14.f, { z14.f }, z30.f \n\t"  \
    "tbl z15.f, { z15.f }, z30.f \n\t"  \
    "tbl z16.f, { z16.f }, z30.f \n\t"  \
    "tbl z17.f, { z17.f }, z30.f \n\t"  \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (2) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// PERM3
#define PERM3_A64FXf  \
asm ( \
    "ld1d { z30.f }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "tbl z12.f, { z12.f }, z30.f \n\t"  \
    "tbl z13.f, { z13.f }, z30.f \n\t"  \
    "tbl z14.f, { z14.f }, z30.f \n\t"  \
    "tbl z15.f, { z15.f }, z30.f \n\t"  \
    "tbl z16.f, { z16.f }, z30.f \n\t"  \
    "tbl z17.f, { z17.f }, z30.f \n\t"  \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (3) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// MULT_2SPIN
#define MULT_2SPIN_A64FXf(A)  \
{ \
    const auto & ref(U[sU][A]); \
asm ( \
    "ld1d { z24.f }, p5/z, [%[fetchptr], -6, mul vl] \n\t" \
    "ld1d { z25.f }, p5/z, [%[fetchptr], -3, mul vl] \n\t" \
    "ld1d { z26.f }, p5/z, [%[fetchptr], 0, mul vl] \n\t" \
    "ld1d { z27.f }, p5/z, [%[fetchptr], -5, mul vl] \n\t" \
    "ld1d { z28.f }, p5/z, [%[fetchptr], -2, mul vl] \n\t" \
    "ld1d { z29.f }, p5/z, [%[fetchptr], 1, mul vl] \n\t" \
    "fmov z18.f , 0 \n\t" \
    "fmov z21.f , 0 \n\t" \
    "fmov z19.f , 0 \n\t" \
    "fmov z22.f , 0 \n\t" \
    "fmov z20.f , 0 \n\t" \
    "fmov z23.f , 0 \n\t" \
    "fcmla z18.f, p5/m, z24.f, z12.f, 0 \n\t" \
    "fcmla z21.f, p5/m, z24.f, z15.f, 0 \n\t" \
    "fcmla z19.f, p5/m, z25.f, z12.f, 0 \n\t" \
    "fcmla z22.f, p5/m, z25.f, z15.f, 0 \n\t" \
    "fcmla z20.f, p5/m, z26.f, z12.f, 0 \n\t" \
    "fcmla z23.f, p5/m, z26.f, z15.f, 0 \n\t" \
    "fcmla z18.f, p5/m, z24.f, z12.f, 90 \n\t" \
    "fcmla z21.f, p5/m, z24.f, z15.f, 90 \n\t" \
    "fcmla z19.f, p5/m, z25.f, z12.f, 90 \n\t" \
    "fcmla z22.f, p5/m, z25.f, z15.f, 90 \n\t" \
    "fcmla z20.f, p5/m, z26.f, z12.f, 90 \n\t" \
    "fcmla z23.f, p5/m, z26.f, z15.f, 90 \n\t" \
    "ld1d { z24.f }, p5/z, [%[fetchptr], -4, mul vl] \n\t" \
    "ld1d { z25.f }, p5/z, [%[fetchptr], -1, mul vl] \n\t" \
    "ld1d { z26.f }, p5/z, [%[fetchptr], 2, mul vl] \n\t" \
    "fcmla z18.f, p5/m, z27.f, z13.f, 0 \n\t" \
    "fcmla z21.f, p5/m, z27.f, z16.f, 0 \n\t" \
    "fcmla z19.f, p5/m, z28.f, z13.f, 0 \n\t" \
    "fcmla z22.f, p5/m, z28.f, z16.f, 0 \n\t" \
    "fcmla z20.f, p5/m, z29.f, z13.f, 0 \n\t" \
    "fcmla z23.f, p5/m, z29.f, z16.f, 0 \n\t" \
    "fcmla z18.f, p5/m, z27.f, z13.f, 90 \n\t" \
    "fcmla z21.f, p5/m, z27.f, z16.f, 90 \n\t" \
    "fcmla z19.f, p5/m, z28.f, z13.f, 90 \n\t" \
    "fcmla z22.f, p5/m, z28.f, z16.f, 90 \n\t" \
    "fcmla z20.f, p5/m, z29.f, z13.f, 90 \n\t" \
    "fcmla z23.f, p5/m, z29.f, z16.f, 90 \n\t" \
    "fcmla z18.f, p5/m, z24.f, z14.f, 0 \n\t" \
    "fcmla z21.f, p5/m, z24.f, z17.f, 0 \n\t" \
    "fcmla z19.f, p5/m, z25.f, z14.f, 0 \n\t" \
    "fcmla z22.f, p5/m, z25.f, z17.f, 0 \n\t" \
    "fcmla z20.f, p5/m, z26.f, z14.f, 0 \n\t" \
    "fcmla z23.f, p5/m, z26.f, z17.f, 0 \n\t" \
    "fcmla z18.f, p5/m, z24.f, z14.f, 90 \n\t" \
    "fcmla z21.f, p5/m, z24.f, z17.f, 90 \n\t" \
    "fcmla z19.f, p5/m, z25.f, z14.f, 90 \n\t" \
    "fcmla z22.f, p5/m, z25.f, z17.f, 90 \n\t" \
    "fcmla z20.f, p5/m, z26.f, z14.f, 90 \n\t" \
    "fcmla z23.f, p5/m, z26.f, z17.f, 90 \n\t" \
    :  \
    : [fetchptr] "r" ((uint64_t)&ref[2][0]) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// XP_PROJ
#define XP_PROJ_A64FXf  \
{ \
asm ( \
    "ld1d { z30.f }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "fcadd z12.f, p5/m, z12.f, z27.f, 90 \n\t" \
    "fcadd z13.f, p5/m, z13.f, z28.f, 90 \n\t" \
    "fcadd z14.f, p5/m, z14.f, z29.f, 90 \n\t" \
    "fcadd z15.f, p5/m, z15.f, z24.f, 90 \n\t" \
    "fcadd z16.f, p5/m, z16.f, z25.f, 90 \n\t" \
    "fcadd z17.f, p5/m, z17.f, z26.f, 90 \n\t" \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (3) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// XP_RECON
#define XP_RECON_A64FXf  \
asm ( \
    "fcadd z6.f, p5/m, z6.f, z21.f, 270 \n\t" \
    "fcadd z7.f, p5/m, z7.f, z22.f, 270 \n\t" \
    "fcadd z8.f, p5/m, z8.f, z23.f, 270 \n\t" \
    "fcadd z9.f, p5/m, z9.f, z18.f, 270 \n\t" \
    "fcadd z10.f, p5/m, z10.f, z19.f, 270 \n\t" \
    "fcadd z11.f, p5/m, z11.f, z20.f, 270 \n\t" \
    "mov z0.f, z18.f \n\t" \
    "mov z1.f, z19.f \n\t" \
    "mov z2.f, z20.f \n\t" \
    "mov z3.f, z21.f \n\t" \
    "mov z4.f, z22.f \n\t" \
    "mov z5.f, z23.f \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// XP_RECON_ACCUM
#define XP_RECON_ACCUM_A64FXf  \
asm ( \
    "fcadd z9.f, p5/m, z9.f, z18.f, 270 \n\t" \
    "fadd z0.f, p5/m, z0.f, z18.f \n\t"  \
    "fcadd z10.f, p5/m, z10.f, z19.f, 270 \n\t" \
    "fadd z1.f, p5/m, z1.f, z19.f \n\t"  \
    "fcadd z11.f, p5/m, z11.f, z20.f, 270 \n\t" \
    "fadd z2.f, p5/m, z2.f, z20.f \n\t"  \
    "fcadd z6.f, p5/m, z6.f, z21.f, 270 \n\t" \
    "fadd z3.f, p5/m, z3.f, z21.f \n\t"  \
    "fcadd z7.f, p5/m, z7.f, z22.f, 270 \n\t" \
    "fadd z4.f, p5/m, z4.f, z22.f \n\t"  \
    "fcadd z8.f, p5/m, z8.f, z23.f, 270 \n\t" \
    "fadd z5.f, p5/m, z5.f, z23.f \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// YP_PROJ
#define YP_PROJ_A64FXf  \
{ \
asm ( \
    "ld1d { z30.f }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "fsub z12.f, p5/m, z12.f, z27.f \n\t" \
    "fsub z13.f, p5/m, z13.f, z28.f \n\t" \
    "fsub z14.f, p5/m, z14.f, z29.f \n\t" \
    "fadd z15.f, p5/m, z15.f, z24.f \n\t"  \
    "fadd z16.f, p5/m, z16.f, z25.f \n\t"  \
    "fadd z17.f, p5/m, z17.f, z26.f \n\t"  \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (2) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// ZP_PROJ
#define ZP_PROJ_A64FXf  \
{ \
asm ( \
    "ld1d { z30.f }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "fcadd z12.f, p5/m, z12.f, z24.f, 90 \n\t" \
    "fcadd z13.f, p5/m, z13.f, z25.f, 90 \n\t" \
    "fcadd z14.f, p5/m, z14.f, z26.f, 90 \n\t" \
    "fcadd z15.f, p5/m, z15.f, z27.f, 270 \n\t" \
    "fcadd z16.f, p5/m, z16.f, z28.f, 270 \n\t" \
    "fcadd z17.f, p5/m, z17.f, z29.f, 270 \n\t" \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (1) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// TP_PROJ
#define TP_PROJ_A64FXf  \
{ \
asm ( \
    "ld1d { z30.f }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "fadd z12.f, p5/m, z12.f, z24.f \n\t"  \
    "fadd z13.f, p5/m, z13.f, z25.f \n\t"  \
    "fadd z14.f, p5/m, z14.f, z26.f \n\t"  \
    "fadd z15.f, p5/m, z15.f, z27.f \n\t"  \
    "fadd z16.f, p5/m, z16.f, z28.f \n\t"  \
    "fadd z17.f, p5/m, z17.f, z29.f \n\t"  \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (0) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// XM_PROJ
#define XM_PROJ_A64FXf  \
{ \
asm ( \
    "ld1d { z30.f }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "fcadd z12.f, p5/m, z12.f, z27.f, 270 \n\t" \
    "fcadd z13.f, p5/m, z13.f, z28.f, 270 \n\t" \
    "fcadd z14.f, p5/m, z14.f, z29.f, 270 \n\t" \
    "fcadd z15.f, p5/m, z15.f, z24.f, 270 \n\t" \
    "fcadd z16.f, p5/m, z16.f, z25.f, 270 \n\t" \
    "fcadd z17.f, p5/m, z17.f, z26.f, 270 \n\t" \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (3) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// XM_RECON
#define XM_RECON_A64FXf  \
asm ( \
    "fcadd z6.f, p5/m, z6.f, z21.f, 90 \n\t" \
    "fcadd z7.f, p5/m, z7.f, z22.f, 90 \n\t" \
    "fcadd z8.f, p5/m, z8.f, z23.f, 90 \n\t" \
    "fcadd z9.f, p5/m, z9.f, z18.f, 90 \n\t" \
    "fcadd z10.f, p5/m, z10.f, z19.f, 90 \n\t" \
    "fcadd z11.f, p5/m, z11.f, z20.f, 90 \n\t" \
    "mov z0.f, z18.f \n\t" \
    "mov z1.f, z19.f \n\t" \
    "mov z2.f, z20.f \n\t" \
    "mov z3.f, z21.f \n\t" \
    "mov z4.f, z22.f \n\t" \
    "mov z5.f, z23.f \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// YM_PROJ
#define YM_PROJ_A64FXf  \
{ \
asm ( \
    "ld1d { z30.f }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "fadd z12.f, p5/m, z12.f, z27.f \n\t"  \
    "fadd z13.f, p5/m, z13.f, z28.f \n\t"  \
    "fadd z14.f, p5/m, z14.f, z29.f \n\t"  \
    "fsub z15.f, p5/m, z15.f, z24.f \n\t" \
    "fsub z16.f, p5/m, z16.f, z25.f \n\t" \
    "fsub z17.f, p5/m, z17.f, z26.f \n\t" \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (2) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// ZM_PROJ
#define ZM_PROJ_A64FXf  \
{ \
asm ( \
    "ld1d { z30.f }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "fcadd z12.f, p5/m, z12.f, z24.f, 270 \n\t" \
    "fcadd z13.f, p5/m, z13.f, z25.f, 270 \n\t" \
    "fcadd z14.f, p5/m, z14.f, z26.f, 270 \n\t" \
    "fcadd z15.f, p5/m, z15.f, z27.f, 90 \n\t" \
    "fcadd z16.f, p5/m, z16.f, z28.f, 90 \n\t" \
    "fcadd z17.f, p5/m, z17.f, z29.f, 90 \n\t" \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (1) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// TM_PROJ
#define TM_PROJ_A64FXf  \
{ \
asm ( \
    "ld1d { z30.f }, p5/z, [%[tableptr], %[index], mul vl] \n\t" \
    "fsub z12.f, p5/m, z12.f, z24.f \n\t" \
    "fsub z13.f, p5/m, z13.f, z25.f \n\t" \
    "fsub z14.f, p5/m, z14.f, z26.f \n\t" \
    "fsub z15.f, p5/m, z15.f, z27.f \n\t" \
    "fsub z16.f, p5/m, z16.f, z28.f \n\t" \
    "fsub z17.f, p5/m, z17.f, z29.f \n\t" \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (0) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// XM_RECON_ACCUM
#define XM_RECON_ACCUM_A64FXf  \
asm ( \
    "fcadd z9.f, p5/m, z9.f, z18.f, 90 \n\t" \
    "fcadd z10.f, p5/m, z10.f, z19.f, 90 \n\t" \
    "fcadd z11.f, p5/m, z11.f, z20.f, 90 \n\t" \
    "fcadd z6.f, p5/m, z6.f, z21.f, 90 \n\t" \
    "fcadd z7.f, p5/m, z7.f, z22.f, 90 \n\t" \
    "fcadd z8.f, p5/m, z8.f, z23.f, 90 \n\t" \
    "mov z0.f, z18.f \n\t" \
    "mov z1.f, z19.f \n\t" \
    "mov z2.f, z20.f \n\t" \
    "mov z3.f, z21.f \n\t" \
    "mov z4.f, z22.f \n\t" \
    "mov z5.f, z23.f \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// YP_RECON_ACCUM
#define YP_RECON_ACCUM_A64FXf  \
asm ( \
    "fadd z0.f, p5/m, z0.f, z18.f \n\t"  \
    "fsub z9.f, p5/m, z9.f, z18.f \n\t" \
    "fadd z1.f, p5/m, z1.f, z19.f \n\t"  \
    "fsub z10.f, p5/m, z10.f, z19.f \n\t" \
    "fadd z2.f, p5/m, z2.f, z20.f \n\t"  \
    "fsub z11.f, p5/m, z11.f, z20.f \n\t" \
    "fadd z3.f, p5/m, z3.f, z21.f \n\t"  \
    "fadd z6.f, p5/m, z6.f, z21.f \n\t"  \
    "fadd z4.f, p5/m, z4.f, z22.f \n\t"  \
    "fadd z7.f, p5/m, z7.f, z22.f \n\t"  \
    "fadd z5.f, p5/m, z5.f, z23.f \n\t"  \
    "fadd z8.f, p5/m, z8.f, z23.f \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// YM_RECON_ACCUM
#define YM_RECON_ACCUM_A64FXf  \
asm ( \
    "fadd z0.f, p5/m, z0.f, z18.f \n\t"  \
    "fadd z9.f, p5/m, z9.f, z18.f \n\t"  \
    "fadd z1.f, p5/m, z1.f, z19.f \n\t"  \
    "fadd z10.f, p5/m, z10.f, z19.f \n\t"  \
    "fadd z2.f, p5/m, z2.f, z20.f \n\t"  \
    "fadd z11.f, p5/m, z11.f, z20.f \n\t"  \
    "fadd z3.f, p5/m, z3.f, z21.f \n\t"  \
    "fsub z6.f, p5/m, z6.f, z21.f \n\t" \
    "fadd z4.f, p5/m, z4.f, z22.f \n\t"  \
    "fsub z7.f, p5/m, z7.f, z22.f \n\t" \
    "fadd z5.f, p5/m, z5.f, z23.f \n\t"  \
    "fsub z8.f, p5/m, z8.f, z23.f \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// ZP_RECON_ACCUM
#define ZP_RECON_ACCUM_A64FXf  \
asm ( \
    "fcadd z6.f, p5/m, z6.f, z18.f, 270 \n\t" \
    "fadd z0.f, p5/m, z0.f, z18.f \n\t"  \
    "fcadd z7.f, p5/m, z7.f, z19.f, 270 \n\t" \
    "fadd z1.f, p5/m, z1.f, z19.f \n\t"  \
    "fcadd z8.f, p5/m, z8.f, z20.f, 270 \n\t" \
    "fadd z2.f, p5/m, z2.f, z20.f \n\t"  \
    "fcadd z9.f, p5/m, z9.f, z21.f, 90 \n\t" \
    "fadd z3.f, p5/m, z3.f, z21.f \n\t"  \
    "fcadd z10.f, p5/m, z10.f, z22.f, 90 \n\t" \
    "fadd z4.f, p5/m, z4.f, z22.f \n\t"  \
    "fcadd z11.f, p5/m, z11.f, z23.f, 90 \n\t" \
    "fadd z5.f, p5/m, z5.f, z23.f \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// ZM_RECON_ACCUM
#define ZM_RECON_ACCUM_A64FXf  \
asm ( \
    "fcadd z6.f, p5/m, z6.f, z18.f, 90 \n\t" \
    "fadd z0.f, p5/m, z0.f, z18.f \n\t"  \
    "fcadd z7.f, p5/m, z7.f, z19.f, 90 \n\t" \
    "fadd z1.f, p5/m, z1.f, z19.f \n\t"  \
    "fcadd z8.f, p5/m, z8.f, z20.f, 90 \n\t" \
    "fadd z2.f, p5/m, z2.f, z20.f \n\t"  \
    "fcadd z9.f, p5/m, z9.f, z21.f, 270 \n\t" \
    "fadd z3.f, p5/m, z3.f, z21.f \n\t"  \
    "fcadd z10.f, p5/m, z10.f, z22.f, 270 \n\t" \
    "fadd z4.f, p5/m, z4.f, z22.f \n\t"  \
    "fcadd z11.f, p5/m, z11.f, z23.f, 270 \n\t" \
    "fadd z5.f, p5/m, z5.f, z23.f \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// TP_RECON_ACCUM
#define TP_RECON_ACCUM_A64FXf  \
asm ( \
    "ptrue p5.f \n\t" \
    "fadd z0.f, p5/m, z0.f, z18.f \n\t"  \
    "fadd z6.f, p5/m, z6.f, z18.f \n\t"  \
    "fadd z1.f, p5/m, z1.f, z19.f \n\t"  \
    "fadd z7.f, p5/m, z7.f, z19.f \n\t"  \
    "fadd z2.f, p5/m, z2.f, z20.f \n\t"  \
    "fadd z8.f, p5/m, z8.f, z20.f \n\t"  \
    "fadd z3.f, p5/m, z3.f, z21.f \n\t"  \
    "fadd z9.f, p5/m, z9.f, z21.f \n\t"  \
    "fadd z4.f, p5/m, z4.f, z22.f \n\t"  \
    "fadd z10.f, p5/m, z10.f, z22.f \n\t"  \
    "fadd z5.f, p5/m, z5.f, z23.f \n\t"  \
    "fadd z11.f, p5/m, z11.f, z23.f \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// TM_RECON_ACCUM
#define TM_RECON_ACCUM_A64FXf  \
asm ( \
    "fadd z0.f, p5/m, z0.f, z18.f \n\t"  \
    "fsub z6.f, p5/m, z6.f, z18.f \n\t" \
    "fadd z1.f, p5/m, z1.f, z19.f \n\t"  \
    "fsub z7.f, p5/m, z7.f, z19.f \n\t" \
    "fadd z2.f, p5/m, z2.f, z20.f \n\t"  \
    "fsub z8.f, p5/m, z8.f, z20.f \n\t" \
    "fadd z3.f, p5/m, z3.f, z21.f \n\t"  \
    "fsub z9.f, p5/m, z9.f, z21.f \n\t" \
    "fadd z4.f, p5/m, z4.f, z22.f \n\t"  \
    "fsub z10.f, p5/m, z10.f, z22.f \n\t" \
    "fadd z5.f, p5/m, z5.f, z23.f \n\t"  \
    "fsub z11.f, p5/m, z11.f, z23.f \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// ZERO_PSI
#define ZERO_PSI_A64FXf  \
asm ( \
    "ptrue p5.f \n\t" \
    "fmov z0.f , 0 \n\t" \
    "fmov z1.f , 0 \n\t" \
    "fmov z2.f , 0 \n\t" \
    "fmov z3.f , 0 \n\t" \
    "fmov z4.f , 0 \n\t" \
    "fmov z5.f , 0 \n\t" \
    "fmov z6.f , 0 \n\t" \
    "fmov z7.f , 0 \n\t" \
    "fmov z8.f , 0 \n\t" \
    "fmov z9.f , 0 \n\t" \
    "fmov z10.f , 0 \n\t" \
    "fmov z11.f , 0 \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

