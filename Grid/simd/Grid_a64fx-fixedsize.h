    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: Grid_a64fx-fixedsize.h

    Copyright (C) 2020

    Author: Nils Meyer         <nils.meyer@ur.de>           Regensburg University

    with support from Arm
            Richard Sandiford  <richard.sandiford@arm.com>

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

/////////////////////////////////////////////////////
// Using SVE ACLE with fixed-size data types
/////////////////////////////////////////////////////


// gcc 10 features
#if __ARM_FEATURE_SVE_BITS==512
/* gcc 10.0.1 and gcc 10.1 bug using ACLE data types  CAS-159553-Y1K4C6
   workaround: use gcc's internal data types, bugfix expected for gcc 10.2
typedef svbool_t    pred __attribute__((arm_sve_vector_bits(512)));
typedef svfloat16_t vech __attribute__((arm_sve_vector_bits(512)));
typedef svfloat32_t vecf __attribute__((arm_sve_vector_bits(512)));
typedef svfloat64_t vecd __attribute__((arm_sve_vector_bits(512)));
typedef svuint32_t  veci __attribute__((arm_sve_vector_bits(512)));
typedef svuint32_t  lutf __attribute__((arm_sve_vector_bits(512))); // LUTs for float
typedef svuint64_t  lutd __attribute__((arm_sve_vector_bits(512))); // LUTs for double
*/
typedef __SVBool_t    pred __attribute__((arm_sve_vector_bits(512)));
typedef __SVFloat16_t vech __attribute__((arm_sve_vector_bits(512)));
typedef __SVFloat32_t vecf __attribute__((arm_sve_vector_bits(512)));
typedef __SVFloat64_t vecd __attribute__((arm_sve_vector_bits(512)));
typedef __SVUint32_t  veci __attribute__((arm_sve_vector_bits(512)));
typedef __SVUint32_t  lutf __attribute__((arm_sve_vector_bits(512))); // LUTs for float
typedef __SVUint64_t  lutd __attribute__((arm_sve_vector_bits(512))); // LUTs for double
#else
#pragma error("Oops. Illegal SVE vector size!?")
#endif /* __ARM_FEATURE_SVE_BITS */

// low-level API
NAMESPACE_BEGIN(Grid);
NAMESPACE_BEGIN(Optimization);

// convenience union types for tables eliminating loads
union ulutf {
  lutf v;
  uint32_t s[16];
};
union ulutd {
  lutd v;
  uint64_t s[8];
};

template <typename T>
struct acle{};

template <>
struct acle<double>{
  static inline lutd tbl_swap(){
    const ulutd t = { .s = {1, 0, 3, 2, 5, 4, 7, 6} };
    return t.v;
  }
  static inline lutd tbl0(){
    const ulutd t = { .s = {4, 5, 6, 7, 0, 1, 2, 3} };
    return t.v;
  }
  static inline lutd tbl1(){
    const ulutd t = { .s = {2, 3, 0, 1, 6, 7, 4, 5} };
    return t.v;
  }
  static inline lutd tbl_exch1a(){ // Exchange1
    const ulutd t = { .s = {0, 1, 4, 5, 2, 3, 6, 7} };
    return t.v;
  }
  static inline lutd tbl_exch1b(){ // Exchange1
    const ulutd t = { .s = {2, 3, 6, 7, 0, 1, 4, 5} };
    return t.v;
  }
  static inline lutd tbl_exch1c(){ // Exchange1
    const ulutd t = { .s = {4, 5, 0, 1, 6, 7, 2, 3} };
    return t.v;
  }
  static inline pred pg1(){return svptrue_b64();}
  static inline pred pg_even(){return svzip1_b64(svptrue_b64(), svpfalse_b());}
  static inline pred pg_odd() {return svzip1_b64(svpfalse_b(), svptrue_b64());}
  static inline vecd zero(){return svdup_f64(0.);}
};

template <>
struct acle<float>{
  // exchange neighboring elements
  static inline lutf tbl_swap(){
    const ulutf t = { .s = {1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14} };
    return t.v;
  }
  static inline lutf tbl0(){
    const ulutf t = { .s = {8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7} };
    return t.v;
  }
  static inline lutf tbl1(){
    const ulutf t = { .s = {4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11} };
    return t.v;
  }
  static inline lutf tbl2(){
    const ulutf t = { .s = {2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13} };
    return t.v;
  }
  static inline lutf tbl_exch1a(){ // Exchange1
    const ulutf t = { .s = {0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7, 12, 13, 14, 15 } };
    return t.v;
  }
  static inline lutf tbl_exch1b(){ // Exchange1
    const ulutf t = { .s = {4, 5, 6, 7, 12, 13, 14, 15, 0, 1, 2, 3, 8, 9, 10, 11 } };
    return t.v;
  }
  static inline lutf tbl_exch1c(){ // Exchange1
    const ulutf t = { .s = {8, 9, 10, 11, 0, 1, 2, 3, 12, 13, 14, 15, 4, 5, 6, 7} };
    return t.v;
  }
  static inline pred pg1(){return svptrue_b32();}
  static inline pred pg_even(){return svzip1_b32(svptrue_b32(), svpfalse_b());}
  static inline pred pg_odd() {return svzip1_b32(svpfalse_b(), svptrue_b32());}
  static inline vecf zero(){return svdup_f32(0.);}
};

template <>
struct acle<uint16_t>{
  static inline pred pg1(){return svptrue_b16();}
  static inline pred pg_even(){return svzip1_b16(svptrue_b16(), svpfalse_b());}
  static inline pred pg_odd() {return svzip1_b16(svpfalse_b(), svptrue_b16());}
  static inline vech zero(){return svdup_f16(0.);}
};

template <>
struct acle<Integer>{
  //static inline svbool_t pg1(){return svptrue_b16();}
  static inline pred pg1(){return svptrue_b32();}
  static inline pred pg_even(){return svzip1_b32(svptrue_b32(), svpfalse_b());}
  static inline pred pg_odd() {return svzip1_b32(svpfalse_b(), svptrue_b32());}
};

// ---------------------------------------------------

struct Vsplat{
  // Complex float
  inline vecf operator()(float a, float b){
    vecf a_v = svdup_f32(a);
    vecf b_v = svdup_f32(b);
    return svzip1(a_v, b_v);
  }
  // Real float
  inline vecf operator()(float a){
    return svdup_f32(a);
  }
  // Complex double
  inline vecd operator()(double a, double b){
    vecd a_v = svdup_f64(a);
    vecd b_v = svdup_f64(b);
    return svzip1(a_v, b_v);
  }
  // Real double
  inline vecd operator()(double a){
    return svdup_f64(a);
  }
  // Integer
  inline veci operator()(Integer a){
    return svdup_u32(a);
  }
};

struct Vstore{
  // Real float
  inline void operator()(vecf a, float *D){
    pred pg1 = acle<float>::pg1();
    svst1(pg1, D, a);
  }
  // Real double
  inline void operator()(vecd a, double *D){
    pred pg1 = acle<double>::pg1();
    svst1(pg1, D, a);
  }
  // Real float
  inline void operator()(veci a, Integer *D){
    pred pg1 = acle<Integer>::pg1();
    svst1(pg1, D, a);
  }
};

struct Vstream{
  // Real float
  inline void operator()(float * a, vecf b){
    pred pg1 = acle<float>::pg1();
    svstnt1(pg1, a, b);
    //svst1(pg1, a, b);
  }
  // Real double
  inline void operator()(double * a, vecd b){
    pred pg1 = acle<double>::pg1();
    svstnt1(pg1, a, b);
    //svst1(pg1, a, b);
  }
};

struct Vset{
  // Complex float
  inline vecf operator()(Grid::ComplexF *a){
    pred pg1 = acle<float>::pg1();
    return svld1(pg1, (float*)a);
  }
  // Complex double
  inline vecd operator()(Grid::ComplexD *a){
    pred pg1 = acle<double>::pg1();
    return svld1(pg1, (double*)a);
  }
  // Real float
  inline vecf operator()(float *a){
    pred pg1 = acle<float>::pg1();
    return svld1(pg1, a);
  }
  // Real double
  inline vecd operator()(double *a){
    pred pg1 = acle<double>::pg1();
    return svld1(pg1, a);
  }
  // Integer
  inline veci operator()(Integer *a){
    pred pg1 = acle<Integer>::pg1();
    return svld1(pg1, a);
  }
};

/////////////////////////////////////////////////////
// Arithmetic operations
/////////////////////////////////////////////////////

struct Sum{
  // Complex/real float
  inline vecf operator()(vecf a, vecf b){
    pred pg1 = acle<float>::pg1();
    return svadd_x(pg1, a, b);
  }
  // Complex/real double
  inline vecd operator()(vecd a, vecd b){
    pred pg1 = acle<double>::pg1();
    return svadd_x(pg1, a, b);
  }
  // Integer
  inline veci operator()(veci a, veci b){
    pred pg1 = acle<Integer>::pg1();
    return svadd_x(pg1, a, b);
  }
};

struct Sub{
  // Complex/real float
  inline vecf operator()(vecf a, vecf b){
    pred pg1 = acle<float>::pg1();
    return svsub_x(pg1, a, b);
  }
  // Complex/real double
  inline vecd operator()(vecd a, vecd b){
    pred pg1 = acle<double>::pg1();
    return svsub_x(pg1, a, b);
  }
  // Integer
  inline veci operator()(veci a, veci b){
    pred pg1 = acle<Integer>::pg1();
    return svsub_x(pg1, a, b);
  }

};

struct Mult{
  // Real float fma
  inline vecf operator()(vecf a, vecf b, vecf c){
    pred pg1 = acle<float>::pg1();
    return svmad_x(pg1, b, c, a);
  }
  // Real double fma
  inline vecd operator()(vecd a, vecd b, vecd c){
    pred pg1 = acle<double>::pg1();
    return svmad_x(pg1, b, c, a);
  }
  // Real float
  inline vecf operator()(vecf a, vecf b){
    pred pg1 = acle<float>::pg1();
    return svmul_x(pg1, a, b);
  }
  // Real double
  inline vecd operator()(vecd a, vecd b){
    pred pg1 = acle<double>::pg1();
    return svmul_x(pg1, a, b);
  }
  // Integer
  inline veci operator()(veci a, veci b){
    pred pg1 = acle<Integer>::pg1();
    return svmul_x(pg1, a, b);
  }
};

struct MultRealPart{
  // Complex float
  inline vecf operator()(vecf a, vecf b){
    pred pg1 = acle<float>::pg1();
    // using FCMLA
    vecf z_v = acle<float>::zero();
    return svcmla_x(pg1, z_v, a, b, 0);
  }
  // Complex double
  inline vecd operator()(vecd a, vecd b){
    pred pg1 = acle<double>::pg1();
    // using FCMLA
    vecd z_v = acle<double>::zero();
    return svcmla_x(pg1, z_v, a, b, 0);
  }
};

struct MaddRealPart{
  // Complex float
  inline vecf operator()(vecf a, vecf b, vecf c){
    pred pg1 = acle<float>::pg1();
    // using FCMLA
    return svcmla_x(pg1, c, a, b, 0);
  }
  // Complex double
  inline vecd operator()(vecd a, vecd b, vecd c){
    pred pg1 = acle<double>::pg1();
    // using FCMLA
    return svcmla_x(pg1, c, a, b, 0);
  }
};

struct MultComplex{
  // Complex a*b
  // Complex float
  inline vecf operator()(vecf a, vecf b){
    pred pg1 = acle<float>::pg1();
    vecf z = acle<float>::zero();
    // using FCMLA
    vecf r_v = svcmla_x(pg1, z, a, b, 0);
    return svcmla_x(pg1, r_v, a, b, 90);
  }
  // Complex double
  inline vecd operator()(vecd a, vecd b){
    pred pg1 = acle<double>::pg1();
    vecd z = acle<double>::zero();
    // using FCMLA
    vecd r_v = svcmla_x(pg1, z, a, b, 0);
    return svcmla_x(pg1, r_v, a, b, 90);
  }
};

struct MultAddComplex{
  // Complex a*b+c
  // Complex float
  inline vecf operator()(vecf a, vecf b, vecf c){
    pred pg1 = acle<float>::pg1();
    // using FCMLA
    vecf r_v = svcmla_x(pg1, c, a, b, 0);
    return svcmla_x(pg1, r_v, a, b, 90);
  }
  // Complex double
  inline vecd operator()(vecd a, vecd b, vecd c){
    pred pg1 = acle<double>::pg1();
    // using FCMLA
    vecd r_v = svcmla_x(pg1, c, a, b, 0);
    return svcmla_x(pg1, r_v, a, b, 90);
  }
};

struct Div{
  // Real float
  inline vecf operator()(vecf a, vecf b){
    pred pg1 = acle<float>::pg1();
    return svdiv_x(pg1, a, b);
  }
  // Real double
  inline vecd operator()(vecd a, vecd b){
    pred pg1 = acle<double>::pg1();
    return svdiv_x(pg1, a, b);
  }
};

struct Conj{
  // Complex float
  inline vecf operator()(vecf a){
    pred pg_odd = acle<float>::pg_odd();
    //return svneg_x(pg_odd, a);  this is unsafe
    return svneg_m(a, pg_odd, a);
  }
  // Complex double
  inline vecd operator()(vecd a){
    pred pg_odd = acle<double>::pg_odd();
    //return svneg_x(pg_odd, a);  this is unsafe
    return svneg_m(a, pg_odd, a);
  }
};

struct TimesMinusI{
  // Complex float
  inline vecf operator()(vecf a, vecf b){
    lutf tbl_swap = acle<float>::tbl_swap();
    pred pg1 = acle<float>::pg1();
    pred pg_odd = acle<float>::pg_odd();

    vecf a_v = svtbl(a, tbl_swap);
    //return svneg_x(pg_odd, a_v);  this is unsafe
    return svneg_m(a_v, pg_odd, a_v);
  }
  // Complex double
  inline vecd operator()(vecd a, vecd b){
    lutd tbl_swap = acle<double>::tbl_swap();
    pred pg1 = acle<double>::pg1();
    pred pg_odd = acle<double>::pg_odd();

    vecd a_v = svtbl(a, tbl_swap);
    //return svneg_x(pg_odd, a_v);  this is unsafe
    return svneg_m(a_v, pg_odd, a_v);
  }
};

struct TimesI{
  // Complex float
  inline vecf operator()(vecf a, vecf b){
    lutf tbl_swap = acle<float>::tbl_swap();
    pred pg1 = acle<float>::pg1();
    pred pg_even = acle<float>::pg_even();

    vecf a_v = svtbl(a, tbl_swap);
    //return svneg_x(pg_even, a_v);  this is unsafe
    return svneg_m(a_v, pg_even, a_v);
  }
  // Complex double
  inline vecd operator()(vecd a, vecd b){
    lutd tbl_swap = acle<double>::tbl_swap();
    pred pg1 = acle<double>::pg1();
    pred pg_even = acle<double>::pg_even();

    vecd a_v = svtbl(a, tbl_swap);
    //return svneg_x(pg_even, a_v);  this is unsafe
    return svneg_m(a_v, pg_even, a_v);
  }
};

struct PrecisionChange {
  static inline vech StoH (vecf sa, vecf sb) {
    pred pg1s = acle<float>::pg1();
    vech ha_v = svcvt_f16_x(pg1s, sa);
    vech hb_v = svcvt_f16_x(pg1s, sb);
    return svuzp1(ha_v, hb_v);
  }
  static inline void HtoS(vech h,vecf &sa,vecf &sb) {
    pred pg1s = acle<float>::pg1();
    vech ha_v = svzip1(h, h);
    vech hb_v = svzip2(h, h);
    sa = svcvt_f32_x(pg1s, ha_v);
    sb = svcvt_f32_x(pg1s, hb_v);
  }
  static inline vecf DtoS (vecd a,vecd b) {
    pred pg1d = acle<double>::pg1();
    vecf sa_v = svcvt_f32_x(pg1d, a);
    vecf sb_v = svcvt_f32_x(pg1d, b);
    return svuzp1(sa_v, sb_v);
  }
  static inline void StoD (vecf s,vecd &a,vecd &b) {
    pred pg1d = acle<double>::pg1();
    vecf sa_v = svzip1(s, s);
    vecf sb_v = svzip2(s, s);
    a = svcvt_f64_x(pg1d, sa_v);
    b = svcvt_f64_x(pg1d, sb_v);
  }
  static inline vech DtoH (vecd a,vecd b,vecd c,vecd d) {
    pred pg1d = acle<double>::pg1();
    pred pg1h = acle<uint16_t>::pg1();
    vech ha_v = svcvt_f16_x(pg1d, a);
    vech hb_v = svcvt_f16_x(pg1d, b);
    vech hc_v = svcvt_f16_x(pg1d, c);
    vech hd_v = svcvt_f16_x(pg1d, d);
    vech hab_v = svuzp1(ha_v, hb_v);
    vech hcd_v = svuzp1(hc_v, hd_v);
    return svuzp1(hab_v, hcd_v);

/*
    vecf sa,sb;
    sa = DtoS(a,b);
    sb = DtoS(c,d);
    return StoH(sa,sb);
*/
  }
  static inline void HtoD(vech h,vecd &a,vecd &b,vecd &c,vecd &d) {
    pred pg1h = acle<uint16_t>::pg1();
    pred pg1d = acle<double>::pg1();
    vech sa_v = svzip1(h, h);
    vech sb_v = svzip2(h, h);
    vech da_v = svzip1(sa_v, sa_v);
    vech db_v = svzip2(sa_v, sa_v);
    vech dc_v = svzip1(sb_v, sb_v);
    vech dd_v = svzip2(sb_v, sb_v);
    a = svcvt_f64_x(pg1d, da_v);
    b = svcvt_f64_x(pg1d, db_v);
    c = svcvt_f64_x(pg1d, dc_v);
    d = svcvt_f64_x(pg1d, dd_v);

/*
    vecf sa,sb;
    HtoS(h,sa,sb);
    StoD(sa,a,b);
    StoD(sb,c,d);
*/
  }
};

struct Exchange{
  // float
  static inline void Exchange0(vecf &out1, vecf &out2, vecf in1, vecf in2){
    vecf r1_v = svext(in1, in1, (uint64_t)8u);
    vecf r2_v = svext(in2, in2, (uint64_t)8u);
    out1 = svext(r1_v, in2, (uint64_t)8u);
    out2 = svext(in1, r2_v, (uint64_t)8u);
  }
  static inline void Exchange1(vecf &out1, vecf &out2, vecf in1, vecf in2){
    // this one is tricky; svtrn2q* from SVE2 fits best, but it is not available in SVE1
    // alternative: use 4-el structure; expect translation into 4x ldp + 4x stp -> SFI
    lutf tbl_exch1a = acle<float>::tbl_exch1a();
    lutf tbl_exch1b = acle<float>::tbl_exch1b();
    lutf tbl_exch1c = acle<float>::tbl_exch1c();

    vecf a1_v = svtbl(in1, tbl_exch1a);
    vecf a2_v = svtbl(in2, tbl_exch1b);
    vecf b1_v  = svext(a2_v, a1_v, (uint64_t)8u);
    vecf b2_v  = svext(a1_v, a2_v, (uint64_t)8u);
    out1 = svtbl(b1_v, tbl_exch1c);
    out2 = svtbl(b2_v, tbl_exch1a);
  }
  static inline void Exchange2(vecf &out1, vecf &out2, vecf in1, vecf in2){
    out1 = (vecf)svtrn1((vecd)in1, (vecd)in2);
    out2 = (vecf)svtrn2((vecd)in1, (vecd)in2);
  }
  static inline void Exchange3(vecf &out1, vecf &out2, vecf in1, vecf in2){
    out1 = svtrn1(in1, in2);
    out2 = svtrn2(in1, in2);
  }

  // double
  static inline void Exchange0(vecd &out1, vecd &out2, vecd in1, vecd in2){
    vecd r1_v = svext(in1, in1, (uint64_t)4u);
    vecd r2_v = svext(in2, in2, (uint64_t)4u);
    out1 = svext(r1_v, in2, (uint64_t)4u);
    out2 = svext(in1, r2_v, (uint64_t)4u);
  }
  static inline void Exchange1(vecd &out1, vecd &out2, vecd in1, vecd in2){
    // this one is tricky; svtrn2q* from SVE2 fits best, but it is not available in SVE1
    // alternative: use 4-el structure; expect translation into 4x ldp + 4x stp -> SFI
    lutd tbl_exch1a = acle<double>::tbl_exch1a();
    lutd tbl_exch1b = acle<double>::tbl_exch1b();
    lutd tbl_exch1c = acle<double>::tbl_exch1c();

    vecd a1_v = svtbl(in1, tbl_exch1a);
    vecd a2_v = svtbl(in2, tbl_exch1b);
    vecd b1_v = svext(a2_v, a1_v, (uint64_t)4u);
    vecd b2_v = svext(a1_v, a2_v, (uint64_t)4u);
    out1 = svtbl(b1_v, tbl_exch1c);
    out2 = svtbl(b2_v, tbl_exch1a);
  }
  static inline void Exchange2(vecd &out1, vecd &out2, vecd in1, vecd in2){
    out1 = svtrn1(in1, in2);
    out2 = svtrn2(in1, in2);
  }
  static inline void Exchange3(vecd &out1, vecd &out2, vecd in1, vecd in2){
    assert(0);
    return;
  }
};

#undef VECTOR_FOR

struct Permute{
  // float
  static inline vecf Permute0(vecf in) {
    return svext(in, in, (uint64_t)8u);
  }
  static inline vecf Permute1(vecf in) {
    lutf tbl_swap = acle<float>::tbl1();
    return svtbl(in, tbl_swap);
  }
  static inline vecf Permute2(vecf in) {
    lutf tbl_swap = acle<float>::tbl2();
    return svtbl(in, tbl_swap);
  }
  static inline vecf Permute3(vecf in) {
    lutf tbl_swap = acle<float>::tbl_swap();
    return svtbl(in, tbl_swap);
  }

  // double
  static inline vecd Permute0(vecd in) {
    return svext(in, in, (uint64_t)(8u / 2u));
  }
  static inline vecd Permute1(vecd in) {
    lutd tbl_swap = acle<double>::tbl1();
    return svtbl(in, tbl_swap);
  }
  static inline vecd Permute2(vecd in) {
    lutd tbl_swap = acle<double>::tbl_swap();
    return svtbl(in, tbl_swap);
  }
  static inline vecd Permute3(vecd in) {
    return in;
  }
};

struct Rotate{

  static inline vecf rotate(vecf in, int n){
    switch(n){
    case 0:  return tRotate<0>(in); break;
    case 1:  return tRotate<1>(in); break;
    case 2:  return tRotate<2>(in); break;
    case 3:  return tRotate<3>(in); break;
    case 4:  return tRotate<4>(in); break;
    case 5:  return tRotate<5>(in); break;
    case 6:  return tRotate<6>(in); break;
    case 7:  return tRotate<7>(in); break;

    case 8:  return tRotate<8>(in); break;
    case 9:  return tRotate<9>(in); break;
    case 10: return tRotate<10>(in); break;
    case 11: return tRotate<11>(in); break;
    case 12: return tRotate<12>(in); break;
    case 13: return tRotate<13>(in); break;
    case 14: return tRotate<14>(in); break;
    case 15: return tRotate<15>(in); break;
    default: assert(0);
    }
  }
  static inline vecd rotate(vecd in, int n){
    switch(n){
    case 0:  return tRotate<0>(in); break;
    case 1:  return tRotate<1>(in); break;
    case 2:  return tRotate<2>(in); break;
    case 3:  return tRotate<3>(in); break;
    case 4:  return tRotate<4>(in); break;
    case 5:  return tRotate<5>(in); break;
    case 6:  return tRotate<6>(in); break;
    case 7:  return tRotate<7>(in); break;
    default: assert(0);
    }
  }

  template <int n> static inline vecf tRotate(vecf in){
    return svext(in, in, (uint64_t)n);
  }
  template <int n> static inline vecd tRotate(vecd in){
    return svext(in, in, (uint64_t)n);
  }
};

// tree-based reduction
#define svred(pg, v)\
svaddv(pg, v);

// left-to-right reduction
// #define svred(pg, v)\
// svadda(pg, 0, v)

template <typename Out_type, typename In_type>
struct Reduce{
  //Need templated class to overload output type
  //General form must generate error if compiled
  inline Out_type operator()(In_type in){
    printf("Error, using wrong Reduce function\n");
    //exit(1);
    return 0;
  }
};
//Complex float Reduce
template <>
inline Grid::ComplexF Reduce<Grid::ComplexF, vecf>::operator()(vecf in){
  pred pg_even = acle<float>::pg_even();
  pred pg_odd  = acle<float>::pg_odd();
  float a = svred(pg_even, in);
  float b = svred(pg_odd, in);
  return Grid::ComplexF(a, b);
}
//Real float Reduce
template <>
inline Grid::RealF Reduce<Grid::RealF, vecf>::operator()(vecf in){
  pred pg1 = acle<float>::pg1();
  return svred(pg1, in);
}
//Complex double Reduce
template <>
inline Grid::ComplexD Reduce<Grid::ComplexD, vecd>::operator()(vecd in){
  pred pg_even = acle<double>::pg_even();
  pred pg_odd  = acle<double>::pg_odd();
  double a = svred(pg_even, in);
  double b = svred(pg_odd, in);
  return Grid::ComplexD(a, b);
}
//Real double Reduce
template <>
inline Grid::RealD Reduce<Grid::RealD, vecd>::operator()(vecd in){
  pred pg1 = acle<double>::pg1();
  return svred(pg1, in);
}
//Integer Reduce
template <>
inline Integer Reduce<Integer, veci>::operator()(veci in){
  pred pg1 = acle<Integer>::pg1();
  return svred(pg1, in);
}

#undef svred

NAMESPACE_END(Optimization);

//////////////////////////////////////////////////////////////////////////////////////
// Here assign types

typedef vech SIMD_Htype; // Reduced precision type
typedef vecf SIMD_Ftype; // Single precision type
typedef vecd SIMD_Dtype; // Double precision type
typedef veci SIMD_Itype; // Integer type

// prefetch utilities
inline void v_prefetch0(int size, const char *ptr){};
inline void prefetch_HINT_T0(const char *ptr){};

// Function name aliases
typedef Optimization::Vsplat   VsplatSIMD;
typedef Optimization::Vstore   VstoreSIMD;
typedef Optimization::Vset     VsetSIMD;
typedef Optimization::Vstream  VstreamSIMD;
template <typename S, typename T> using ReduceSIMD = Optimization::Reduce<S,T>;

// Arithmetic operations
typedef Optimization::Sum            SumSIMD;
typedef Optimization::Sub            SubSIMD;
typedef Optimization::Div            DivSIMD;
typedef Optimization::Mult           MultSIMD;
typedef Optimization::MultComplex    MultComplexSIMD;
typedef Optimization::MultAddComplex MultAddComplexSIMD;
typedef Optimization::MultRealPart   MultRealPartSIMD;
typedef Optimization::MaddRealPart   MaddRealPartSIMD;
typedef Optimization::Conj           ConjSIMD;
typedef Optimization::TimesMinusI    TimesMinusISIMD;
typedef Optimization::TimesI         TimesISIMD;

NAMESPACE_END(Grid);
