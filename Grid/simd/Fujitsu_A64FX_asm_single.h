/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: Fujitsu_A64FX_asm_single.h

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
#define LOAD_CHIMU(base)               LOAD_CHIMU_INTERLEAVED_A64FXf(base)  
#define PREFETCH_CHIMU_L1(A)           PREFETCH_CHIMU_L1_INTERNAL_A64FXf(A)  
#define PREFETCH_GAUGE_L1(A)           PREFETCH_GAUGE_L1_INTERNAL_A64FXf(A)  
#define PREFETCH_CHIMU_L2(A)           PREFETCH_CHIMU_L2_INTERNAL_A64FXf(A)  
#define PREFETCH_GAUGE_L2(A)           PREFETCH_GAUGE_L2_INTERNAL_A64FXf(A)  
#define PF_GAUGE(A)  
#define PREFETCH_RESULT_L2_STORE(A)    PREFETCH_RESULT_L2_STORE_INTERNAL_A64FXf(A)  
#define PREFETCH_RESULT_L1_STORE(A)    PREFETCH_RESULT_L1_STORE_INTERNAL_A64FXf(A)  
#define PREFETCH1_CHIMU(A)             PREFETCH_CHIMU_L1(A)  
#define PREFETCH_CHIMU(A)              PREFETCH_CHIMU_L1(A)  
#define LOCK_GAUGE(A)  
#define UNLOCK_GAUGE(A)  
#define MASK_REGS                      DECLARATIONS_A64FXf  
#define SAVE_RESULT(A,B)               RESULT_A64FXf(A); PREFETCH_RESULT_L2_STORE(B)  
#define MULT_2SPIN_1(Dir)              MULT_2SPIN_1_A64FXf(Dir)  
#define MULT_2SPIN_2                   MULT_2SPIN_2_A64FXf  
#define LOAD_CHI(base)                 LOAD_CHI_A64FXf(base)  
#define ADD_RESULT(base,basep)         LOAD_CHIMU(base); ADD_RESULT_INTERNAL_A64FXf; RESULT_A64FXf(base)  
#define XP_PROJ                        XP_PROJ_A64FXf  
#define YP_PROJ                        YP_PROJ_A64FXf  
#define ZP_PROJ                        ZP_PROJ_A64FXf  
#define TP_PROJ                        TP_PROJ_A64FXf  
#define XM_PROJ                        XM_PROJ_A64FXf  
#define YM_PROJ                        YM_PROJ_A64FXf  
#define ZM_PROJ                        ZM_PROJ_A64FXf  
#define TM_PROJ                        TM_PROJ_A64FXf  
#define XP_RECON                       XP_RECON_A64FXf  
#define XM_RECON                       XM_RECON_A64FXf  
#define XM_RECON_ACCUM                 XM_RECON_ACCUM_A64FXf  
#define YM_RECON_ACCUM                 YM_RECON_ACCUM_A64FXf  
#define ZM_RECON_ACCUM                 ZM_RECON_ACCUM_A64FXf  
#define TM_RECON_ACCUM                 TM_RECON_ACCUM_A64FXf  
#define XP_RECON_ACCUM                 XP_RECON_ACCUM_A64FXf  
#define YP_RECON_ACCUM                 YP_RECON_ACCUM_A64FXf  
#define ZP_RECON_ACCUM                 ZP_RECON_ACCUM_A64FXf  
#define TP_RECON_ACCUM                 TP_RECON_ACCUM_A64FXf  
#define PERMUTE_DIR0                   0  
#define PERMUTE_DIR1                   1  
#define PERMUTE_DIR2                   2  
#define PERMUTE_DIR3                   3  
#define PERMUTE                        PERMUTE_A64FXf;  
#define LOAD_TABLE(Dir)                if (Dir == 0) { LOAD_TABLE0; } else if (Dir == 1) { LOAD_TABLE1 } else if (Dir == 2) { LOAD_TABLE2; } else if (Dir == 3) { LOAD_TABLE3; }  
#define MAYBEPERM(A,perm)              if (perm) { PERMUTE; }  
// DECLARATIONS
#define DECLARATIONS_A64FXf  \
    const uint32_t lut[4][16] = { \
        {8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7}, \
        {4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11}, \
        {2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13}, \
        {1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14} }; \
asm ( \
    "fmov z31.s , 0 \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// RESULT
#define RESULT_A64FXf(base)  \
{ \
asm ( \
    "str z0, [%[storeptr], -6, mul vl] \n\t" \
    "str z1, [%[storeptr], -5, mul vl] \n\t" \
    "str z2, [%[storeptr], -4, mul vl] \n\t" \
    "str z3, [%[storeptr], -3, mul vl] \n\t" \
    "str z4, [%[storeptr], -2, mul vl] \n\t" \
    "str z5, [%[storeptr], -1, mul vl] \n\t" \
    "str z6, [%[storeptr], 0, mul vl] \n\t" \
    "str z7, [%[storeptr], 1, mul vl] \n\t" \
    "str z8, [%[storeptr], 2, mul vl] \n\t" \
    "str z9, [%[storeptr], 3, mul vl] \n\t" \
    "str z10, [%[storeptr], 4, mul vl] \n\t" \
    "str z11, [%[storeptr], 5, mul vl] \n\t" \
    :  \
    : [storeptr] "r" (base + 2 * 3 * 64) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// PREFETCH_CHIMU_L2 (prefetch to L2)
#define PREFETCH_CHIMU_L2_INTERNAL_A64FXf(base)  \
{ \
asm ( \
    "prfd PLDL2STRM, p5, [%[fetchptr], 0, mul vl] \n\t" \
    "prfd PLDL2STRM, p5, [%[fetchptr], 4, mul vl] \n\t" \
    "prfd PLDL2STRM, p5, [%[fetchptr], 8, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (base) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// PREFETCH_CHIMU_L1 (prefetch to L1)
#define PREFETCH_CHIMU_L1_INTERNAL_A64FXf(base)  \
{ \
asm ( \
    "prfd PLDL1STRM, p5, [%[fetchptr], 0, mul vl] \n\t" \
    "prfd PLDL1STRM, p5, [%[fetchptr], 4, mul vl] \n\t" \
    "prfd PLDL1STRM, p5, [%[fetchptr], 8, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (base) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// PREFETCH_GAUGE_L2 (prefetch to L2)
#define PREFETCH_GAUGE_L2_INTERNAL_A64FXf(A)  \
{ \
    const auto & ref(U[sUn](A)); uint64_t baseU = (uint64_t)&ref + 3 * 3 * 64; \
asm ( \
    "prfd PLDL2STRM, p5, [%[fetchptr], -4, mul vl] \n\t" \
    "prfd PLDL2STRM, p5, [%[fetchptr], 0, mul vl] \n\t" \
    "prfd PLDL2STRM, p5, [%[fetchptr], 4, mul vl] \n\t" \
    "prfd PLDL2STRM, p5, [%[fetchptr], 8, mul vl] \n\t" \
    "prfd PLDL2STRM, p5, [%[fetchptr], 12, mul vl] \n\t" \
    "prfd PLDL2STRM, p5, [%[fetchptr], 16, mul vl] \n\t" \
    "prfd PLDL2STRM, p5, [%[fetchptr], 20, mul vl] \n\t" \
    "prfd PLDL2STRM, p5, [%[fetchptr], 24, mul vl] \n\t" \
    "prfd PLDL2STRM, p5, [%[fetchptr], 28, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (baseU) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// PREFETCH_GAUGE_L1 (prefetch to L1)
#define PREFETCH_GAUGE_L1_INTERNAL_A64FXf(A)  \
{ \
    const auto & ref(U[sU](A)); uint64_t baseU = (uint64_t)&ref; \
asm ( \
    "prfd PLDL1STRM, p5, [%[fetchptr], 0, mul vl] \n\t" \
    "prfd PLDL1STRM, p5, [%[fetchptr], 4, mul vl] \n\t" \
    "prfd PLDL1STRM, p5, [%[fetchptr], 8, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (baseU) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// LOAD_CHI
#define LOAD_CHI_A64FXf(base)  \
{ \
asm ( \
    "ldr z12, [%[fetchptr], 0, mul vl] \n\t" \
    "ldr z13, [%[fetchptr], 1, mul vl] \n\t" \
    "ldr z14, [%[fetchptr], 2, mul vl] \n\t" \
    "ldr z15, [%[fetchptr], 3, mul vl] \n\t" \
    "ldr z16, [%[fetchptr], 4, mul vl] \n\t" \
    "ldr z17, [%[fetchptr], 5, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (base) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// LOAD_CHIMU
#define LOAD_CHIMU_INTERLEAVED_A64FXf(base)  \
{ \
asm ( \
    "ptrue p5.s \n\t" \
    "ldr z12, [%[fetchptr], -6, mul vl] \n\t" \
    "ldr z21, [%[fetchptr], 3, mul vl] \n\t" \
    "ldr z15, [%[fetchptr], -3, mul vl] \n\t" \
    "ldr z18, [%[fetchptr], 0, mul vl] \n\t" \
    "ldr z13, [%[fetchptr], -5, mul vl] \n\t" \
    "ldr z22, [%[fetchptr], 4, mul vl] \n\t" \
    "ldr z16, [%[fetchptr], -2, mul vl] \n\t" \
    "ldr z19, [%[fetchptr], 1, mul vl] \n\t" \
    "ldr z14, [%[fetchptr], -4, mul vl] \n\t" \
    "ldr z23, [%[fetchptr], 5, mul vl] \n\t" \
    "ldr z17, [%[fetchptr], -1, mul vl] \n\t" \
    "ldr z20, [%[fetchptr], 2, mul vl] \n\t" \
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
    "ptrue p5.s \n\t" \
    "ldr z12, [%[fetchptr], -6, mul vl] \n\t" \
    "ldr z18, [%[fetchptr], 0, mul vl] \n\t" \
    "ldr z13, [%[fetchptr], -5, mul vl] \n\t" \
    "ldr z19, [%[fetchptr], 1, mul vl] \n\t" \
    "ldr z14, [%[fetchptr], -4, mul vl] \n\t" \
    "ldr z20, [%[fetchptr], 2, mul vl] \n\t" \
    "ldr z15, [%[fetchptr], -3, mul vl] \n\t" \
    "ldr z21, [%[fetchptr], 3, mul vl] \n\t" \
    "ldr z16, [%[fetchptr], -2, mul vl] \n\t" \
    "ldr z22, [%[fetchptr], 4, mul vl] \n\t" \
    "ldr z17, [%[fetchptr], -1, mul vl] \n\t" \
    "ldr z23, [%[fetchptr], 5, mul vl] \n\t" \
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
    "ptrue p5.s \n\t" \
    "ldr z12, [%[fetchptr], -6, mul vl] \n\t" \
    "ldr z21, [%[fetchptr], 3, mul vl] \n\t" \
    "ldr z13, [%[fetchptr], -5, mul vl] \n\t" \
    "ldr z22, [%[fetchptr], 4, mul vl] \n\t" \
    "ldr z14, [%[fetchptr], -4, mul vl] \n\t" \
    "ldr z23, [%[fetchptr], 5, mul vl] \n\t" \
    "ldr z15, [%[fetchptr], -3, mul vl] \n\t" \
    "ldr z18, [%[fetchptr], 0, mul vl] \n\t" \
    "ldr z16, [%[fetchptr], -2, mul vl] \n\t" \
    "ldr z19, [%[fetchptr], 1, mul vl] \n\t" \
    "ldr z17, [%[fetchptr], -1, mul vl] \n\t" \
    "ldr z20, [%[fetchptr], 2, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (&ref[2][0]) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// LOAD_TABLE0
#define LOAD_TABLE0  \
asm ( \
    "ldr z30, [%[tableptr], %[index], mul vl] \n\t" \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (0) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// LOAD_TABLE1
#define LOAD_TABLE1  \
asm ( \
    "ldr z30, [%[tableptr], %[index], mul vl] \n\t" \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (1) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// LOAD_TABLE2
#define LOAD_TABLE2  \
asm ( \
    "ldr z30, [%[tableptr], %[index], mul vl] \n\t" \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (2) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// LOAD_TABLE3
#define LOAD_TABLE3  \
asm ( \
    "ldr z30, [%[tableptr], %[index], mul vl] \n\t" \
    :  \
    : [tableptr] "r" (&lut[0]),[index] "i" (3) \
    : "memory","cc","p5","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// PERMUTE
#define PERMUTE_A64FXf  \
asm ( \
    "tbl z12.s, { z12.s }, z30.s \n\t"  \
    "tbl z13.s, { z13.s }, z30.s \n\t"  \
    "tbl z14.s, { z14.s }, z30.s \n\t"  \
    "tbl z15.s, { z15.s }, z30.s \n\t"  \
    "tbl z16.s, { z16.s }, z30.s \n\t"  \
    "tbl z17.s, { z17.s }, z30.s \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// LOAD_GAUGE
#define LOAD_GAUGE  \
    const auto & ref(U[sU](A)); uint64_t baseU = (uint64_t)&ref; \
{ \
asm ( \
    "ptrue p5.s \n\t" \
    "ldr z24, [%[fetchptr], -6, mul vl] \n\t" \
    "ldr z25, [%[fetchptr], -3, mul vl] \n\t" \
    "ldr z26, [%[fetchptr], 0, mul vl] \n\t" \
    "ldr z27, [%[fetchptr], -5, mul vl] \n\t" \
    "ldr z28, [%[fetchptr], -2, mul vl] \n\t" \
    "ldr z29, [%[fetchptr], 1, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (baseU + 2 * 3 * 64) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// MULT_2SPIN
#define MULT_2SPIN_1_A64FXf(A)  \
{ \
    const auto & ref(U[sU](A)); uint64_t baseU = (uint64_t)&ref; \
asm ( \
    "ldr z24, [%[fetchptr], -6, mul vl] \n\t" \
    "ldr z25, [%[fetchptr], -3, mul vl] \n\t" \
    "ldr z26, [%[fetchptr], 0, mul vl] \n\t" \
    "ldr z27, [%[fetchptr], -5, mul vl] \n\t" \
    "ldr z28, [%[fetchptr], -2, mul vl] \n\t" \
    "ldr z29, [%[fetchptr], 1, mul vl] \n\t" \
    "movprfx z18.s, p5/m, z31.s \n\t" \
    "fcmla z18.s, p5/m, z24.s, z12.s, 0 \n\t" \
    "movprfx z21.s, p5/m, z31.s \n\t" \
    "fcmla z21.s, p5/m, z24.s, z15.s, 0 \n\t" \
    "movprfx z19.s, p5/m, z31.s \n\t" \
    "fcmla z19.s, p5/m, z25.s, z12.s, 0 \n\t" \
    "movprfx z22.s, p5/m, z31.s \n\t" \
    "fcmla z22.s, p5/m, z25.s, z15.s, 0 \n\t" \
    "movprfx z20.s, p5/m, z31.s \n\t" \
    "fcmla z20.s, p5/m, z26.s, z12.s, 0 \n\t" \
    "movprfx z23.s, p5/m, z31.s \n\t" \
    "fcmla z23.s, p5/m, z26.s, z15.s, 0 \n\t" \
    "fcmla z18.s, p5/m, z24.s, z12.s, 90 \n\t" \
    "fcmla z21.s, p5/m, z24.s, z15.s, 90 \n\t" \
    "fcmla z19.s, p5/m, z25.s, z12.s, 90 \n\t" \
    "fcmla z22.s, p5/m, z25.s, z15.s, 90 \n\t" \
    "fcmla z20.s, p5/m, z26.s, z12.s, 90 \n\t" \
    "fcmla z23.s, p5/m, z26.s, z15.s, 90 \n\t" \
    "ldr z24, [%[fetchptr], -4, mul vl] \n\t" \
    "ldr z25, [%[fetchptr], -1, mul vl] \n\t" \
    "ldr z26, [%[fetchptr], 2, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (baseU + 2 * 3 * 64) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// MULT_2SPIN_BACKEND
#define MULT_2SPIN_2_A64FXf  \
{ \
asm ( \
    "fcmla z18.s, p5/m, z27.s, z13.s, 0 \n\t" \
    "fcmla z21.s, p5/m, z27.s, z16.s, 0 \n\t" \
    "fcmla z19.s, p5/m, z28.s, z13.s, 0 \n\t" \
    "fcmla z22.s, p5/m, z28.s, z16.s, 0 \n\t" \
    "fcmla z20.s, p5/m, z29.s, z13.s, 0 \n\t" \
    "fcmla z23.s, p5/m, z29.s, z16.s, 0 \n\t" \
    "fcmla z18.s, p5/m, z27.s, z13.s, 90 \n\t" \
    "fcmla z21.s, p5/m, z27.s, z16.s, 90 \n\t" \
    "fcmla z19.s, p5/m, z28.s, z13.s, 90 \n\t" \
    "fcmla z22.s, p5/m, z28.s, z16.s, 90 \n\t" \
    "fcmla z20.s, p5/m, z29.s, z13.s, 90 \n\t" \
    "fcmla z23.s, p5/m, z29.s, z16.s, 90 \n\t" \
    "fcmla z18.s, p5/m, z24.s, z14.s, 0 \n\t" \
    "fcmla z21.s, p5/m, z24.s, z17.s, 0 \n\t" \
    "fcmla z19.s, p5/m, z25.s, z14.s, 0 \n\t" \
    "fcmla z22.s, p5/m, z25.s, z17.s, 0 \n\t" \
    "fcmla z20.s, p5/m, z26.s, z14.s, 0 \n\t" \
    "fcmla z23.s, p5/m, z26.s, z17.s, 0 \n\t" \
    "fcmla z18.s, p5/m, z24.s, z14.s, 90 \n\t" \
    "fcmla z21.s, p5/m, z24.s, z17.s, 90 \n\t" \
    "fcmla z19.s, p5/m, z25.s, z14.s, 90 \n\t" \
    "fcmla z22.s, p5/m, z25.s, z17.s, 90 \n\t" \
    "fcmla z20.s, p5/m, z26.s, z14.s, 90 \n\t" \
    "fcmla z23.s, p5/m, z26.s, z17.s, 90 \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// XP_PROJ
#define XP_PROJ_A64FXf  \
{ \
asm ( \
    "fcadd z12.s, p5/m, z12.s, z21.s, 90 \n\t" \
    "fcadd z13.s, p5/m, z13.s, z22.s, 90 \n\t" \
    "fcadd z14.s, p5/m, z14.s, z23.s, 90 \n\t" \
    "fcadd z15.s, p5/m, z15.s, z18.s, 90 \n\t" \
    "fcadd z16.s, p5/m, z16.s, z19.s, 90 \n\t" \
    "fcadd z17.s, p5/m, z17.s, z20.s, 90 \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// XP_RECON
#define XP_RECON_A64FXf  \
asm ( \
    "movprfx z6.s, p5/m, z31.s \n\t" \
    "fcadd z6.s, p5/m, z6.s, z21.s, 270 \n\t" \
    "movprfx z7.s, p5/m, z31.s \n\t" \
    "fcadd z7.s, p5/m, z7.s, z22.s, 270 \n\t" \
    "movprfx z8.s, p5/m, z31.s \n\t" \
    "fcadd z8.s, p5/m, z8.s, z23.s, 270 \n\t" \
    "movprfx z9.s, p5/m, z31.s \n\t" \
    "fcadd z9.s, p5/m, z9.s, z18.s, 270 \n\t" \
    "movprfx z10.s, p5/m, z31.s \n\t" \
    "fcadd z10.s, p5/m, z10.s, z19.s, 270 \n\t" \
    "movprfx z11.s, p5/m, z31.s \n\t" \
    "fcadd z11.s, p5/m, z11.s, z20.s, 270 \n\t" \
    "mov z0.s, p5/m, z18.s \n\t" \
    "mov z1.s, p5/m, z19.s \n\t" \
    "mov z2.s, p5/m, z20.s \n\t" \
    "mov z3.s, p5/m, z21.s \n\t" \
    "mov z4.s, p5/m, z22.s \n\t" \
    "mov z5.s, p5/m, z23.s \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// XP_RECON_ACCUM
#define XP_RECON_ACCUM_A64FXf  \
asm ( \
    "fcadd z9.s, p5/m, z9.s, z18.s, 270 \n\t" \
    "fadd z0.s, p5/m, z0.s, z18.s \n\t"  \
    "fcadd z10.s, p5/m, z10.s, z19.s, 270 \n\t" \
    "fadd z1.s, p5/m, z1.s, z19.s \n\t"  \
    "fcadd z11.s, p5/m, z11.s, z20.s, 270 \n\t" \
    "fadd z2.s, p5/m, z2.s, z20.s \n\t"  \
    "fcadd z6.s, p5/m, z6.s, z21.s, 270 \n\t" \
    "fadd z3.s, p5/m, z3.s, z21.s \n\t"  \
    "fcadd z7.s, p5/m, z7.s, z22.s, 270 \n\t" \
    "fadd z4.s, p5/m, z4.s, z22.s \n\t"  \
    "fcadd z8.s, p5/m, z8.s, z23.s, 270 \n\t" \
    "fadd z5.s, p5/m, z5.s, z23.s \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// YP_PROJ
#define YP_PROJ_A64FXf  \
{ \
asm ( \
    "fsub z12.s, p5/m, z12.s, z21.s \n\t" \
    "fsub z13.s, p5/m, z13.s, z22.s \n\t" \
    "fsub z14.s, p5/m, z14.s, z23.s \n\t" \
    "fadd z15.s, p5/m, z15.s, z18.s \n\t"  \
    "fadd z16.s, p5/m, z16.s, z19.s \n\t"  \
    "fadd z17.s, p5/m, z17.s, z20.s \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// ZP_PROJ
#define ZP_PROJ_A64FXf  \
{ \
asm ( \
    "fcadd z12.s, p5/m, z12.s, z18.s, 90 \n\t" \
    "fcadd z13.s, p5/m, z13.s, z19.s, 90 \n\t" \
    "fcadd z14.s, p5/m, z14.s, z20.s, 90 \n\t" \
    "fcadd z15.s, p5/m, z15.s, z21.s, 270 \n\t" \
    "fcadd z16.s, p5/m, z16.s, z22.s, 270 \n\t" \
    "fcadd z17.s, p5/m, z17.s, z23.s, 270 \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// TP_PROJ
#define TP_PROJ_A64FXf  \
{ \
asm ( \
    "fadd z12.s, p5/m, z12.s, z18.s \n\t"  \
    "fadd z13.s, p5/m, z13.s, z19.s \n\t"  \
    "fadd z14.s, p5/m, z14.s, z20.s \n\t"  \
    "fadd z15.s, p5/m, z15.s, z21.s \n\t"  \
    "fadd z16.s, p5/m, z16.s, z22.s \n\t"  \
    "fadd z17.s, p5/m, z17.s, z23.s \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// XM_PROJ
#define XM_PROJ_A64FXf  \
{ \
asm ( \
    "fcadd z12.s, p5/m, z12.s, z21.s, 270 \n\t" \
    "fcadd z13.s, p5/m, z13.s, z22.s, 270 \n\t" \
    "fcadd z14.s, p5/m, z14.s, z23.s, 270 \n\t" \
    "fcadd z15.s, p5/m, z15.s, z18.s, 270 \n\t" \
    "fcadd z16.s, p5/m, z16.s, z19.s, 270 \n\t" \
    "fcadd z17.s, p5/m, z17.s, z20.s, 270 \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// XM_RECON
#define XM_RECON_A64FXf  \
asm ( \
    "movprfx z6.s, p5/m, z31.s \n\t" \
    "fcadd z6.s, p5/m, z6.s, z21.s, 90 \n\t" \
    "movprfx z7.s, p5/m, z31.s \n\t" \
    "fcadd z7.s, p5/m, z7.s, z22.s, 90 \n\t" \
    "movprfx z8.s, p5/m, z31.s \n\t" \
    "fcadd z8.s, p5/m, z8.s, z23.s, 90 \n\t" \
    "movprfx z9.s, p5/m, z31.s \n\t" \
    "fcadd z9.s, p5/m, z9.s, z18.s, 90 \n\t" \
    "movprfx z10.s, p5/m, z31.s \n\t" \
    "fcadd z10.s, p5/m, z10.s, z19.s, 90 \n\t" \
    "movprfx z11.s, p5/m, z31.s \n\t" \
    "fcadd z11.s, p5/m, z11.s, z20.s, 90 \n\t" \
    "mov z0.s, p5/m, z18.s \n\t" \
    "mov z1.s, p5/m, z19.s \n\t" \
    "mov z2.s, p5/m, z20.s \n\t" \
    "mov z3.s, p5/m, z21.s \n\t" \
    "mov z4.s, p5/m, z22.s \n\t" \
    "mov z5.s, p5/m, z23.s \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// YM_PROJ
#define YM_PROJ_A64FXf  \
{ \
asm ( \
    "fadd z12.s, p5/m, z12.s, z21.s \n\t"  \
    "fadd z13.s, p5/m, z13.s, z22.s \n\t"  \
    "fadd z14.s, p5/m, z14.s, z23.s \n\t"  \
    "fsub z15.s, p5/m, z15.s, z18.s \n\t" \
    "fsub z16.s, p5/m, z16.s, z19.s \n\t" \
    "fsub z17.s, p5/m, z17.s, z20.s \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// ZM_PROJ
#define ZM_PROJ_A64FXf  \
{ \
asm ( \
    "fcadd z12.s, p5/m, z12.s, z18.s, 270 \n\t" \
    "fcadd z13.s, p5/m, z13.s, z19.s, 270 \n\t" \
    "fcadd z14.s, p5/m, z14.s, z20.s, 270 \n\t" \
    "fcadd z15.s, p5/m, z15.s, z21.s, 90 \n\t" \
    "fcadd z16.s, p5/m, z16.s, z22.s, 90 \n\t" \
    "fcadd z17.s, p5/m, z17.s, z23.s, 90 \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// TM_PROJ
#define TM_PROJ_A64FXf  \
{ \
asm ( \
    "ptrue p5.s \n\t" \
    "fsub z12.s, p5/m, z12.s, z18.s \n\t" \
    "fsub z13.s, p5/m, z13.s, z19.s \n\t" \
    "fsub z14.s, p5/m, z14.s, z20.s \n\t" \
    "fsub z15.s, p5/m, z15.s, z21.s \n\t" \
    "fsub z16.s, p5/m, z16.s, z22.s \n\t" \
    "fsub z17.s, p5/m, z17.s, z23.s \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); \
}
// XM_RECON_ACCUM
#define XM_RECON_ACCUM_A64FXf  \
asm ( \
    "fcadd z9.s, p5/m, z9.s, z18.s, 90 \n\t" \
    "fcadd z10.s, p5/m, z10.s, z19.s, 90 \n\t" \
    "fcadd z11.s, p5/m, z11.s, z20.s, 90 \n\t" \
    "fcadd z6.s, p5/m, z6.s, z21.s, 90 \n\t" \
    "fcadd z7.s, p5/m, z7.s, z22.s, 90 \n\t" \
    "fcadd z8.s, p5/m, z8.s, z23.s, 90 \n\t" \
    "fadd z0.s, p5/m, z0.s, z18.s \n\t"  \
    "fadd z1.s, p5/m, z1.s, z19.s \n\t"  \
    "fadd z2.s, p5/m, z2.s, z20.s \n\t"  \
    "fadd z3.s, p5/m, z3.s, z21.s \n\t"  \
    "fadd z4.s, p5/m, z4.s, z22.s \n\t"  \
    "fadd z5.s, p5/m, z5.s, z23.s \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// YP_RECON_ACCUM
#define YP_RECON_ACCUM_A64FXf  \
asm ( \
    "fadd z0.s, p5/m, z0.s, z18.s \n\t"  \
    "fsub z9.s, p5/m, z9.s, z18.s \n\t" \
    "fadd z1.s, p5/m, z1.s, z19.s \n\t"  \
    "fsub z10.s, p5/m, z10.s, z19.s \n\t" \
    "fadd z2.s, p5/m, z2.s, z20.s \n\t"  \
    "fsub z11.s, p5/m, z11.s, z20.s \n\t" \
    "fadd z3.s, p5/m, z3.s, z21.s \n\t"  \
    "fadd z6.s, p5/m, z6.s, z21.s \n\t"  \
    "fadd z4.s, p5/m, z4.s, z22.s \n\t"  \
    "fadd z7.s, p5/m, z7.s, z22.s \n\t"  \
    "fadd z5.s, p5/m, z5.s, z23.s \n\t"  \
    "fadd z8.s, p5/m, z8.s, z23.s \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// YM_RECON_ACCUM
#define YM_RECON_ACCUM_A64FXf  \
asm ( \
    "fadd z0.s, p5/m, z0.s, z18.s \n\t"  \
    "fadd z9.s, p5/m, z9.s, z18.s \n\t"  \
    "fadd z1.s, p5/m, z1.s, z19.s \n\t"  \
    "fadd z10.s, p5/m, z10.s, z19.s \n\t"  \
    "fadd z2.s, p5/m, z2.s, z20.s \n\t"  \
    "fadd z11.s, p5/m, z11.s, z20.s \n\t"  \
    "fadd z3.s, p5/m, z3.s, z21.s \n\t"  \
    "fsub z6.s, p5/m, z6.s, z21.s \n\t" \
    "fadd z4.s, p5/m, z4.s, z22.s \n\t"  \
    "fsub z7.s, p5/m, z7.s, z22.s \n\t" \
    "fadd z5.s, p5/m, z5.s, z23.s \n\t"  \
    "fsub z8.s, p5/m, z8.s, z23.s \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// ZP_RECON_ACCUM
#define ZP_RECON_ACCUM_A64FXf  \
asm ( \
    "fcadd z6.s, p5/m, z6.s, z18.s, 270 \n\t" \
    "fadd z0.s, p5/m, z0.s, z18.s \n\t"  \
    "fcadd z7.s, p5/m, z7.s, z19.s, 270 \n\t" \
    "fadd z1.s, p5/m, z1.s, z19.s \n\t"  \
    "fcadd z8.s, p5/m, z8.s, z20.s, 270 \n\t" \
    "fadd z2.s, p5/m, z2.s, z20.s \n\t"  \
    "fcadd z9.s, p5/m, z9.s, z21.s, 90 \n\t" \
    "fadd z3.s, p5/m, z3.s, z21.s \n\t"  \
    "fcadd z10.s, p5/m, z10.s, z22.s, 90 \n\t" \
    "fadd z4.s, p5/m, z4.s, z22.s \n\t"  \
    "fcadd z11.s, p5/m, z11.s, z23.s, 90 \n\t" \
    "fadd z5.s, p5/m, z5.s, z23.s \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// ZM_RECON_ACCUM
#define ZM_RECON_ACCUM_A64FXf  \
asm ( \
    "fcadd z6.s, p5/m, z6.s, z18.s, 90 \n\t" \
    "fadd z0.s, p5/m, z0.s, z18.s \n\t"  \
    "fcadd z7.s, p5/m, z7.s, z19.s, 90 \n\t" \
    "fadd z1.s, p5/m, z1.s, z19.s \n\t"  \
    "fcadd z8.s, p5/m, z8.s, z20.s, 90 \n\t" \
    "fadd z2.s, p5/m, z2.s, z20.s \n\t"  \
    "fcadd z9.s, p5/m, z9.s, z21.s, 270 \n\t" \
    "fadd z3.s, p5/m, z3.s, z21.s \n\t"  \
    "fcadd z10.s, p5/m, z10.s, z22.s, 270 \n\t" \
    "fadd z4.s, p5/m, z4.s, z22.s \n\t"  \
    "fcadd z11.s, p5/m, z11.s, z23.s, 270 \n\t" \
    "fadd z5.s, p5/m, z5.s, z23.s \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// TP_RECON_ACCUM
#define TP_RECON_ACCUM_A64FXf  \
asm ( \
    "fadd z0.s, p5/m, z0.s, z18.s \n\t"  \
    "fadd z6.s, p5/m, z6.s, z18.s \n\t"  \
    "fadd z1.s, p5/m, z1.s, z19.s \n\t"  \
    "fadd z7.s, p5/m, z7.s, z19.s \n\t"  \
    "fadd z2.s, p5/m, z2.s, z20.s \n\t"  \
    "fadd z8.s, p5/m, z8.s, z20.s \n\t"  \
    "fadd z3.s, p5/m, z3.s, z21.s \n\t"  \
    "fadd z9.s, p5/m, z9.s, z21.s \n\t"  \
    "fadd z4.s, p5/m, z4.s, z22.s \n\t"  \
    "fadd z10.s, p5/m, z10.s, z22.s \n\t"  \
    "fadd z5.s, p5/m, z5.s, z23.s \n\t"  \
    "fadd z11.s, p5/m, z11.s, z23.s \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// TM_RECON_ACCUM
#define TM_RECON_ACCUM_A64FXf  \
asm ( \
    "fadd z0.s, p5/m, z0.s, z18.s \n\t"  \
    "fsub z6.s, p5/m, z6.s, z18.s \n\t" \
    "fadd z1.s, p5/m, z1.s, z19.s \n\t"  \
    "fsub z7.s, p5/m, z7.s, z19.s \n\t" \
    "fadd z2.s, p5/m, z2.s, z20.s \n\t"  \
    "fsub z8.s, p5/m, z8.s, z20.s \n\t" \
    "fadd z3.s, p5/m, z3.s, z21.s \n\t"  \
    "fsub z9.s, p5/m, z9.s, z21.s \n\t" \
    "fadd z4.s, p5/m, z4.s, z22.s \n\t"  \
    "fsub z10.s, p5/m, z10.s, z22.s \n\t" \
    "fadd z5.s, p5/m, z5.s, z23.s \n\t"  \
    "fsub z11.s, p5/m, z11.s, z23.s \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// ZERO_PSI
#define ZERO_PSI_A64FXf  \
asm ( \
    "ptrue p5.s \n\t" \
    "fmov z0.s , 0 \n\t" \
    "fmov z1.s , 0 \n\t" \
    "fmov z2.s , 0 \n\t" \
    "fmov z3.s , 0 \n\t" \
    "fmov z4.s , 0 \n\t" \
    "fmov z5.s , 0 \n\t" \
    "fmov z6.s , 0 \n\t" \
    "fmov z7.s , 0 \n\t" \
    "fmov z8.s , 0 \n\t" \
    "fmov z9.s , 0 \n\t" \
    "fmov z10.s , 0 \n\t" \
    "fmov z11.s , 0 \n\t" \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

// PREFETCH_RESULT_L2_STORE (prefetch store to L2)
#define PREFETCH_RESULT_L2_STORE_INTERNAL_A64FXf(base)  \
{ \
asm ( \
    "prfd PSTL2STRM, p5, [%[fetchptr], 0, mul vl] \n\t" \
    "prfd PSTL2STRM, p5, [%[fetchptr], 4, mul vl] \n\t" \
    "prfd PSTL2STRM, p5, [%[fetchptr], 8, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (base) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// PREFETCH_RESULT_L1_STORE (prefetch store to L1)
#define PREFETCH_RESULT_L1_STORE_INTERNAL_A64FXf(base)  \
{ \
asm ( \
    "prfd PSTL1STRM, p5, [%[fetchptr], 0, mul vl] \n\t" \
    "prfd PSTL1STRM, p5, [%[fetchptr], 4, mul vl] \n\t" \
    "prfd PSTL1STRM, p5, [%[fetchptr], 8, mul vl] \n\t" \
    :  \
    : [fetchptr] "r" (base) \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31","memory" \
); \
}
// ADD_RESULT_INTERNAL
#define ADD_RESULT_INTERNAL_A64FXf  \
asm ( \
    "fadd z0.s, p5/m, z0.s, z12.s \n\t"  \
    "fadd z1.s, p5/m, z1.s, z13.s \n\t"  \
    "fadd z2.s, p5/m, z2.s, z14.s \n\t"  \
    "fadd z3.s, p5/m, z3.s, z15.s \n\t"  \
    "fadd z4.s, p5/m, z4.s, z16.s \n\t"  \
    "fadd z5.s, p5/m, z5.s, z17.s \n\t"  \
    "fadd z6.s, p5/m, z6.s, z18.s \n\t"  \
    "fadd z7.s, p5/m, z7.s, z19.s \n\t"  \
    "fadd z8.s, p5/m, z8.s, z20.s \n\t"  \
    "fadd z9.s, p5/m, z9.s, z21.s \n\t"  \
    "fadd z10.s, p5/m, z10.s, z22.s \n\t"  \
    "fadd z11.s, p5/m, z11.s, z23.s \n\t"  \
    :  \
    :  \
    : "p5","cc","z0","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","z19","z20","z21","z22","z23","z24","z25","z26","z27","z28","z29","z30","z31" \
); 

