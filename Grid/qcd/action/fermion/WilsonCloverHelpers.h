/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/WilsonCloverHelpers.h

    Copyright (C) 2021 - 2022

    Author: Daniel Richtmann <daniel.richtmann@gmail.com>

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

#pragma once

// Helper routines that implement common clover functionality

NAMESPACE_BEGIN(Grid);

template<class Impl> class WilsonCloverHelpers {
public:

  INHERIT_IMPL_TYPES(Impl);
  INHERIT_CLOVER_TYPES(Impl);

  // Computing C_{\mu \nu}(x) as in Eq.(B.39) in Zbigniew Sroczynski's PhD thesis
  static GaugeLinkField Cmunu(std::vector<GaugeLinkField> &U, GaugeLinkField &lambda, int mu, int nu)
  {
    conformable(lambda.Grid(), U[0].Grid());
    GaugeLinkField out(lambda.Grid()), tmp(lambda.Grid());
    // insertion in upper staple
    // please check redundancy of shift operations

    // C1+
    tmp = lambda * U[nu];
    out = Impl::ShiftStaple(Impl::CovShiftForward(tmp, nu, Impl::CovShiftBackward(U[mu], mu, Impl::CovShiftIdentityBackward(U[nu], nu))), mu);

    // C2+
    tmp = U[mu] * Impl::ShiftStaple(adj(lambda), mu);
    out += Impl::ShiftStaple(Impl::CovShiftForward(U[nu], nu, Impl::CovShiftBackward(tmp, mu, Impl::CovShiftIdentityBackward(U[nu], nu))), mu);

    // C3+
    tmp = U[nu] * Impl::ShiftStaple(adj(lambda), nu);
    out += Impl::ShiftStaple(Impl::CovShiftForward(U[nu], nu, Impl::CovShiftBackward(U[mu], mu, Impl::CovShiftIdentityBackward(tmp, nu))), mu);

    // C4+
    out += Impl::ShiftStaple(Impl::CovShiftForward(U[nu], nu, Impl::CovShiftBackward(U[mu], mu, Impl::CovShiftIdentityBackward(U[nu], nu))), mu) * lambda;

    // insertion in lower staple
    // C1-
    out -= Impl::ShiftStaple(lambda, mu) * Impl::ShiftStaple(Impl::CovShiftBackward(U[nu], nu, Impl::CovShiftBackward(U[mu], mu, U[nu])), mu);

    // C2-
    tmp = adj(lambda) * U[nu];
    out -= Impl::ShiftStaple(Impl::CovShiftBackward(tmp, nu, Impl::CovShiftBackward(U[mu], mu, U[nu])), mu);

    // C3-
    tmp = lambda * U[nu];
    out -= Impl::ShiftStaple(Impl::CovShiftBackward(U[nu], nu, Impl::CovShiftBackward(U[mu], mu, tmp)), mu);

    // C4-
    out -= Impl::ShiftStaple(Impl::CovShiftBackward(U[nu], nu, Impl::CovShiftBackward(U[mu], mu, U[nu])), mu) * lambda;

    return out;
  }

  static CloverField fillCloverYZ(const GaugeLinkField &F)
  {
    CloverField T(F.Grid());
    T = Zero();
    autoView(T_v,T,AcceleratorWrite);
    autoView(F_v,F,AcceleratorRead);
    accelerator_for(i, T.Grid()->oSites(),CloverField::vector_type::Nsimd(),
    {
      coalescedWrite(T_v[i]()(0, 1), coalescedRead(timesMinusI(F_v[i]()())));
      coalescedWrite(T_v[i]()(1, 0), coalescedRead(timesMinusI(F_v[i]()())));
      coalescedWrite(T_v[i]()(2, 3), coalescedRead(timesMinusI(F_v[i]()())));
      coalescedWrite(T_v[i]()(3, 2), coalescedRead(timesMinusI(F_v[i]()())));
    });

    return T;
  }

  static CloverField fillCloverXZ(const GaugeLinkField &F)
  {
    CloverField T(F.Grid());
    T = Zero();
    
    autoView(T_v, T,AcceleratorWrite);
    autoView(F_v, F,AcceleratorRead);
    accelerator_for(i, T.Grid()->oSites(),CloverField::vector_type::Nsimd(),
    {
      coalescedWrite(T_v[i]()(0, 1), coalescedRead(-F_v[i]()()));
      coalescedWrite(T_v[i]()(1, 0), coalescedRead(F_v[i]()()));
      coalescedWrite(T_v[i]()(2, 3), coalescedRead(-F_v[i]()()));
      coalescedWrite(T_v[i]()(3, 2), coalescedRead(F_v[i]()()));
    });

    return T;
  }

  static CloverField fillCloverXY(const GaugeLinkField &F)
  {
    CloverField T(F.Grid());
    T = Zero();

    autoView(T_v,T,AcceleratorWrite);
    autoView(F_v,F,AcceleratorRead);
    accelerator_for(i, T.Grid()->oSites(),CloverField::vector_type::Nsimd(),
    {
      coalescedWrite(T_v[i]()(0, 0), coalescedRead(timesMinusI(F_v[i]()())));
      coalescedWrite(T_v[i]()(1, 1), coalescedRead(timesI(F_v[i]()())));
      coalescedWrite(T_v[i]()(2, 2), coalescedRead(timesMinusI(F_v[i]()())));
      coalescedWrite(T_v[i]()(3, 3), coalescedRead(timesI(F_v[i]()())));
    });

    return T;
  }

  static CloverField fillCloverXT(const GaugeLinkField &F)
  {
    CloverField T(F.Grid());
    T = Zero();

    autoView( T_v , T, AcceleratorWrite);
    autoView( F_v , F, AcceleratorRead);
    accelerator_for(i, T.Grid()->oSites(),CloverField::vector_type::Nsimd(),
    {
      coalescedWrite(T_v[i]()(0, 1), coalescedRead(timesI(F_v[i]()())));
      coalescedWrite(T_v[i]()(1, 0), coalescedRead(timesI(F_v[i]()())));
      coalescedWrite(T_v[i]()(2, 3), coalescedRead(timesMinusI(F_v[i]()())));
      coalescedWrite(T_v[i]()(3, 2), coalescedRead(timesMinusI(F_v[i]()())));
    });

    return T;
  }

  static CloverField fillCloverYT(const GaugeLinkField &F)
  {
    CloverField T(F.Grid());
    T = Zero();
    
    autoView( T_v ,T,AcceleratorWrite);
    autoView( F_v ,F,AcceleratorRead);
    accelerator_for(i, T.Grid()->oSites(),CloverField::vector_type::Nsimd(),
    {
      coalescedWrite(T_v[i]()(0, 1), coalescedRead(-(F_v[i]()())));
      coalescedWrite(T_v[i]()(1, 0), coalescedRead((F_v[i]()())));
      coalescedWrite(T_v[i]()(2, 3), coalescedRead((F_v[i]()())));
      coalescedWrite(T_v[i]()(3, 2), coalescedRead(-(F_v[i]()())));
    });

    return T;
  }

  static CloverField fillCloverZT(const GaugeLinkField &F)
  {
    CloverField T(F.Grid());

    T = Zero();

    autoView( T_v , T,AcceleratorWrite);
    autoView( F_v , F,AcceleratorRead);
    accelerator_for(i, T.Grid()->oSites(),CloverField::vector_type::Nsimd(),
    {
      coalescedWrite(T_v[i]()(0, 0), coalescedRead(timesI(F_v[i]()())));
      coalescedWrite(T_v[i]()(1, 1), coalescedRead(timesMinusI(F_v[i]()())));
      coalescedWrite(T_v[i]()(2, 2), coalescedRead(timesMinusI(F_v[i]()())));
      coalescedWrite(T_v[i]()(3, 3), coalescedRead(timesI(F_v[i]()())));
    });

    return T;
  }

  template<class _Spinor>
  static accelerator_inline void multClover(_Spinor& phi, const SiteClover& C, const _Spinor& chi) {
    auto CC = coalescedRead(C);
    mult(&phi, &CC, &chi);
  }

  template<class _SpinorField>
  inline void multCloverField(_SpinorField& out, const CloverField& C, const _SpinorField& phi) {
    const int Nsimd = SiteSpinor::Nsimd();
    autoView(out_v, out, AcceleratorWrite);
    autoView(phi_v, phi, AcceleratorRead);
    autoView(C_v,   C,   AcceleratorRead);
    typedef decltype(coalescedRead(out_v[0])) calcSpinor;
    accelerator_for(sss,out.Grid()->oSites(),Nsimd,{
      calcSpinor tmp;
      multClover(tmp,C_v[sss],phi_v(sss));
      coalescedWrite(out_v[sss],tmp);
    });
  }
};


template<class Impl> class CompactWilsonCloverHelpers {
public:

  INHERIT_COMPACT_CLOVER_SIZES(Impl);

  INHERIT_IMPL_TYPES(Impl);
  INHERIT_CLOVER_TYPES(Impl);
  INHERIT_COMPACT_CLOVER_TYPES(Impl);

  #if 0
  static accelerator_inline typename SiteCloverTriangle::vector_type triangle_elem(const SiteCloverTriangle& triangle, int block, int i, int j) {
    assert(i != j);
    if(i < j) {
      return triangle()(block)(triangle_index(i, j));
    } else { // i > j
      return conjugate(triangle()(block)(triangle_index(i, j)));
    }
  }
  #else
  template<typename vobj>
  static accelerator_inline vobj triangle_elem(const iImplCloverTriangle<vobj>& triangle, int block, int i, int j) {
    assert(i != j);
    if(i < j) {
      return triangle()(block)(triangle_index(i, j));
    } else { // i > j
      return conjugate(triangle()(block)(triangle_index(i, j)));
    }
  }
  #endif

  static accelerator_inline int triangle_index(int i, int j) {
    if(i == j)
      return 0;
    else if(i < j)
      return Nred * (Nred - 1) / 2 - (Nred - i) * (Nred - i - 1) / 2 + j - i - 1;
    else // i > j
      return Nred * (Nred - 1) / 2 - (Nred - j) * (Nred - j - 1) / 2 + i - j - 1;
  }

  static void MooeeKernel_gpu(int                        Nsite,
                              int                        Ls,
                              const FermionField&        in,
                              FermionField&              out,
                              const CloverDiagonalField& diagonal,
                              const CloverTriangleField& triangle) {
    autoView(diagonal_v, diagonal, AcceleratorRead);
    autoView(triangle_v, triangle, AcceleratorRead);
    autoView(in_v,       in,       AcceleratorRead);
    autoView(out_v,      out,      AcceleratorWrite);

    typedef decltype(coalescedRead(out_v[0])) CalcSpinor;

    const uint64_t NN = Nsite * Ls;

    accelerator_for(ss, NN, Simd::Nsimd(), {
      int sF = ss;
      int sU = ss/Ls;
      CalcSpinor res;
      CalcSpinor in_t = in_v(sF);
      auto diagonal_t = diagonal_v(sU);
      auto triangle_t = triangle_v(sU);
      for(int block=0; block<Nhs; block++) {
        int s_start = block*Nhs;
        for(int i=0; i<Nred; i++) {
          int si = s_start + i/Nc, ci = i%Nc;
          res()(si)(ci) = diagonal_t()(block)(i) * in_t()(si)(ci);
          for(int j=0; j<Nred; j++) {
            if (j == i) continue;
            int sj = s_start + j/Nc, cj = j%Nc;
            res()(si)(ci) = res()(si)(ci) + triangle_elem(triangle_t, block, i, j) * in_t()(sj)(cj);
          };
        };
      };
      coalescedWrite(out_v[sF], res);
    });
  }

  static void MooeeKernel_cpu(int                        Nsite,
                              int                        Ls,
                              const FermionField&        in,
                              FermionField&              out,
                              const CloverDiagonalField& diagonal,
                              const CloverTriangleField& triangle) {
    autoView(diagonal_v, diagonal, CpuRead);
    autoView(triangle_v, triangle, CpuRead);
    autoView(in_v,       in,       CpuRead);
    autoView(out_v,      out,      CpuWrite);

    typedef SiteSpinor CalcSpinor;

#if defined(A64FX) || defined(A64FXFIXEDSIZE)
#define PREFETCH_CLOVER(BASE) {                                     \
    uint64_t base;                                                  \
    int pf_dist_L1 = 1;                                             \
    int pf_dist_L2 = -5; /* -> penalty -> disable */                \
                                                                    \
    if ((pf_dist_L1 >= 0) && (sU + pf_dist_L1 < Nsite)) {           \
      base = (uint64_t)&diag_t()(pf_dist_L1+BASE)(0);               \
      svprfd(svptrue_b64(), (int64_t*)(base +    0), SV_PLDL1STRM); \
      svprfd(svptrue_b64(), (int64_t*)(base +  256), SV_PLDL1STRM); \
      svprfd(svptrue_b64(), (int64_t*)(base +  512), SV_PLDL1STRM); \
      svprfd(svptrue_b64(), (int64_t*)(base +  768), SV_PLDL1STRM); \
      svprfd(svptrue_b64(), (int64_t*)(base + 1024), SV_PLDL1STRM); \
      svprfd(svptrue_b64(), (int64_t*)(base + 1280), SV_PLDL1STRM); \
    }                                                               \
                                                                    \
    if ((pf_dist_L2 >= 0) && (sU + pf_dist_L2 < Nsite)) {           \
      base = (uint64_t)&diag_t()(pf_dist_L2+BASE)(0);               \
      svprfd(svptrue_b64(), (int64_t*)(base +    0), SV_PLDL2STRM); \
      svprfd(svptrue_b64(), (int64_t*)(base +  256), SV_PLDL2STRM); \
      svprfd(svptrue_b64(), (int64_t*)(base +  512), SV_PLDL2STRM); \
      svprfd(svptrue_b64(), (int64_t*)(base +  768), SV_PLDL2STRM); \
      svprfd(svptrue_b64(), (int64_t*)(base + 1024), SV_PLDL2STRM); \
      svprfd(svptrue_b64(), (int64_t*)(base + 1280), SV_PLDL2STRM); \
    }                                                               \
  }
// TODO: Implement/generalize this for other architectures
// I played around a bit on KNL (see below) but didn't bring anything
// #elif defined(AVX512)
// #define PREFETCH_CLOVER(BASE) {                              \
//     uint64_t base;                                           \
//     int pf_dist_L1 = 1;                                      \
//     int pf_dist_L2 = +4;                                     \
//                                                              \
//     if ((pf_dist_L1 >= 0) && (sU + pf_dist_L1 < Nsite)) {    \
//       base = (uint64_t)&diag_t()(pf_dist_L1+BASE)(0);        \
//       _mm_prefetch((const char*)(base +    0), _MM_HINT_T0); \
//       _mm_prefetch((const char*)(base +   64), _MM_HINT_T0); \
//       _mm_prefetch((const char*)(base +  128), _MM_HINT_T0); \
//       _mm_prefetch((const char*)(base +  192), _MM_HINT_T0); \
//       _mm_prefetch((const char*)(base +  256), _MM_HINT_T0); \
//       _mm_prefetch((const char*)(base +  320), _MM_HINT_T0); \
//     }                                                        \
//                                                              \
//     if ((pf_dist_L2 >= 0) && (sU + pf_dist_L2 < Nsite)) {    \
//       base = (uint64_t)&diag_t()(pf_dist_L2+BASE)(0);        \
//       _mm_prefetch((const char*)(base +    0), _MM_HINT_T1); \
//       _mm_prefetch((const char*)(base +   64), _MM_HINT_T1); \
//       _mm_prefetch((const char*)(base +  128), _MM_HINT_T1); \
//       _mm_prefetch((const char*)(base +  192), _MM_HINT_T1); \
//       _mm_prefetch((const char*)(base +  256), _MM_HINT_T1); \
//       _mm_prefetch((const char*)(base +  320), _MM_HINT_T1); \
//     }                                                        \
//   }
#else
#define PREFETCH_CLOVER(BASE)
#endif

    const uint64_t NN = Nsite * Ls;

    thread_for(ss, NN, {
      int sF = ss;
      int sU = ss/Ls;
      CalcSpinor res;
      CalcSpinor in_t = in_v[sF];
      auto diag_t     = diagonal_v[sU]; // "diag" instead of "diagonal" here to make code below easier to read
      auto triangle_t = triangle_v[sU];

      // upper half
      PREFETCH_CLOVER(0);

      auto in_cc_0_0 = conjugate(in_t()(0)(0)); // Nils: reduces number
      auto in_cc_0_1 = conjugate(in_t()(0)(1)); // of conjugates from
      auto in_cc_0_2 = conjugate(in_t()(0)(2)); // 30 to 20
      auto in_cc_1_0 = conjugate(in_t()(1)(0));
      auto in_cc_1_1 = conjugate(in_t()(1)(1));

      res()(0)(0) =               diag_t()(0)( 0) * in_t()(0)(0)
                  +           triangle_t()(0)( 0) * in_t()(0)(1)
                  +           triangle_t()(0)( 1) * in_t()(0)(2)
                  +           triangle_t()(0)( 2) * in_t()(1)(0)
                  +           triangle_t()(0)( 3) * in_t()(1)(1)
                  +           triangle_t()(0)( 4) * in_t()(1)(2);

      res()(0)(1) =           triangle_t()(0)( 0) * in_cc_0_0;
      res()(0)(1) =               diag_t()(0)( 1) * in_t()(0)(1)
                  +           triangle_t()(0)( 5) * in_t()(0)(2)
                  +           triangle_t()(0)( 6) * in_t()(1)(0)
                  +           triangle_t()(0)( 7) * in_t()(1)(1)
                  +           triangle_t()(0)( 8) * in_t()(1)(2)
                  + conjugate(       res()(0)( 1));

      res()(0)(2) =           triangle_t()(0)( 1) * in_cc_0_0
                  +           triangle_t()(0)( 5) * in_cc_0_1;
      res()(0)(2) =               diag_t()(0)( 2) * in_t()(0)(2)
                  +           triangle_t()(0)( 9) * in_t()(1)(0)
                  +           triangle_t()(0)(10) * in_t()(1)(1)
                  +           triangle_t()(0)(11) * in_t()(1)(2)
                  + conjugate(       res()(0)( 2));

      res()(1)(0) =           triangle_t()(0)( 2) * in_cc_0_0
                  +           triangle_t()(0)( 6) * in_cc_0_1
                  +           triangle_t()(0)( 9) * in_cc_0_2;
      res()(1)(0) =               diag_t()(0)( 3) * in_t()(1)(0)
                  +           triangle_t()(0)(12) * in_t()(1)(1)
                  +           triangle_t()(0)(13) * in_t()(1)(2)
                  + conjugate(       res()(1)( 0));

      res()(1)(1) =           triangle_t()(0)( 3) * in_cc_0_0
                  +           triangle_t()(0)( 7) * in_cc_0_1
                  +           triangle_t()(0)(10) * in_cc_0_2
                  +           triangle_t()(0)(12) * in_cc_1_0;
      res()(1)(1) =               diag_t()(0)( 4) * in_t()(1)(1)
                  +           triangle_t()(0)(14) * in_t()(1)(2)
                  + conjugate(       res()(1)( 1));

      res()(1)(2) =           triangle_t()(0)( 4) * in_cc_0_0
                  +           triangle_t()(0)( 8) * in_cc_0_1
                  +           triangle_t()(0)(11) * in_cc_0_2
                  +           triangle_t()(0)(13) * in_cc_1_0
                  +           triangle_t()(0)(14) * in_cc_1_1;
      res()(1)(2) =               diag_t()(0)( 5) * in_t()(1)(2)
                  + conjugate(       res()(1)( 2));

      vstream(out_v[sF]()(0)(0), res()(0)(0));
      vstream(out_v[sF]()(0)(1), res()(0)(1));
      vstream(out_v[sF]()(0)(2), res()(0)(2));
      vstream(out_v[sF]()(1)(0), res()(1)(0));
      vstream(out_v[sF]()(1)(1), res()(1)(1));
      vstream(out_v[sF]()(1)(2), res()(1)(2));

      // lower half
      PREFETCH_CLOVER(1);

      auto in_cc_2_0 = conjugate(in_t()(2)(0));
      auto in_cc_2_1 = conjugate(in_t()(2)(1));
      auto in_cc_2_2 = conjugate(in_t()(2)(2));
      auto in_cc_3_0 = conjugate(in_t()(3)(0));
      auto in_cc_3_1 = conjugate(in_t()(3)(1));

      res()(2)(0) =               diag_t()(1)( 0) * in_t()(2)(0)
                  +           triangle_t()(1)( 0) * in_t()(2)(1)
                  +           triangle_t()(1)( 1) * in_t()(2)(2)
                  +           triangle_t()(1)( 2) * in_t()(3)(0)
                  +           triangle_t()(1)( 3) * in_t()(3)(1)
                  +           triangle_t()(1)( 4) * in_t()(3)(2);

      res()(2)(1) =           triangle_t()(1)( 0) * in_cc_2_0;
      res()(2)(1) =               diag_t()(1)( 1) * in_t()(2)(1)
                  +           triangle_t()(1)( 5) * in_t()(2)(2)
                  +           triangle_t()(1)( 6) * in_t()(3)(0)
                  +           triangle_t()(1)( 7) * in_t()(3)(1)
                  +           triangle_t()(1)( 8) * in_t()(3)(2)
                  + conjugate(       res()(2)( 1));

      res()(2)(2) =           triangle_t()(1)( 1) * in_cc_2_0
                  +           triangle_t()(1)( 5) * in_cc_2_1;
      res()(2)(2) =               diag_t()(1)( 2) * in_t()(2)(2)
                  +           triangle_t()(1)( 9) * in_t()(3)(0)
                  +           triangle_t()(1)(10) * in_t()(3)(1)
                  +           triangle_t()(1)(11) * in_t()(3)(2)
                  + conjugate(       res()(2)( 2));

      res()(3)(0) =           triangle_t()(1)( 2) * in_cc_2_0
                  +           triangle_t()(1)( 6) * in_cc_2_1
                  +           triangle_t()(1)( 9) * in_cc_2_2;
      res()(3)(0) =               diag_t()(1)( 3) * in_t()(3)(0)
                  +           triangle_t()(1)(12) * in_t()(3)(1)
                  +           triangle_t()(1)(13) * in_t()(3)(2)
                  + conjugate(       res()(3)( 0));

      res()(3)(1) =           triangle_t()(1)( 3) * in_cc_2_0
                  +           triangle_t()(1)( 7) * in_cc_2_1
                  +           triangle_t()(1)(10) * in_cc_2_2
                  +           triangle_t()(1)(12) * in_cc_3_0;
      res()(3)(1) =               diag_t()(1)( 4) * in_t()(3)(1)
                  +           triangle_t()(1)(14) * in_t()(3)(2)
                  + conjugate(       res()(3)( 1));

      res()(3)(2) =           triangle_t()(1)( 4) * in_cc_2_0
                  +           triangle_t()(1)( 8) * in_cc_2_1
                  +           triangle_t()(1)(11) * in_cc_2_2
                  +           triangle_t()(1)(13) * in_cc_3_0
                  +           triangle_t()(1)(14) * in_cc_3_1;
      res()(3)(2) =               diag_t()(1)( 5) * in_t()(3)(2)
                  + conjugate(       res()(3)( 2));

      vstream(out_v[sF]()(2)(0), res()(2)(0));
      vstream(out_v[sF]()(2)(1), res()(2)(1));
      vstream(out_v[sF]()(2)(2), res()(2)(2));
      vstream(out_v[sF]()(3)(0), res()(3)(0));
      vstream(out_v[sF]()(3)(1), res()(3)(1));
      vstream(out_v[sF]()(3)(2), res()(3)(2));
    });
  }

  static void MooeeKernel(int                        Nsite,
                          int                        Ls,
                          const FermionField&        in,
                          FermionField&              out,
                          const CloverDiagonalField& diagonal,
                          const CloverTriangleField& triangle) {
#if defined(GRID_CUDA) || defined(GRID_HIP)
    MooeeKernel_gpu(Nsite, Ls, in, out, diagonal, triangle);
#else
    MooeeKernel_cpu(Nsite, Ls, in, out, diagonal, triangle);
#endif
  }

  static void Invert(const CloverDiagonalField& diagonal,
                     const CloverTriangleField& triangle,
                     CloverDiagonalField&       diagonalInv,
                     CloverTriangleField&       triangleInv) {
    conformable(diagonal, diagonalInv);
    conformable(triangle, triangleInv);
    conformable(diagonal, triangle);

    diagonalInv.Checkerboard() = diagonal.Checkerboard();
    triangleInv.Checkerboard() = triangle.Checkerboard();

    GridBase* grid = diagonal.Grid();

    long lsites = grid->lSites();

    typedef typename SiteCloverDiagonal::scalar_object scalar_object_diagonal;
    typedef typename SiteCloverTriangle::scalar_object scalar_object_triangle;

    autoView(diagonal_v,  diagonal,  CpuRead);
    autoView(triangle_v,  triangle,  CpuRead);
    autoView(diagonalInv_v, diagonalInv, CpuWrite);
    autoView(triangleInv_v, triangleInv, CpuWrite);

    thread_for(site, lsites, { // NOTE: Not on GPU because of Eigen & (peek/poke)LocalSite
      Eigen::MatrixXcd clover_inv_eigen = Eigen::MatrixXcd::Zero(Ns*Nc, Ns*Nc);
      Eigen::MatrixXcd clover_eigen = Eigen::MatrixXcd::Zero(Ns*Nc, Ns*Nc);

      scalar_object_diagonal diagonal_tmp     = Zero();
      scalar_object_diagonal diagonal_inv_tmp = Zero();
      scalar_object_triangle triangle_tmp     = Zero();
      scalar_object_triangle triangle_inv_tmp = Zero();

      Coordinate lcoor;
      grid->LocalIndexToLocalCoor(site, lcoor);

      peekLocalSite(diagonal_tmp, diagonal_v, lcoor);
      peekLocalSite(triangle_tmp, triangle_v, lcoor);

      // TODO: can we save time here by inverting the two 6x6 hermitian matrices separately?
      for (long s_row=0;s_row<Ns;s_row++) {
        for (long s_col=0;s_col<Ns;s_col++) {
          if(abs(s_row - s_col) > 1 || s_row + s_col == 3) continue;
          int block       = s_row / Nhs;
          int s_row_block = s_row % Nhs;
          int s_col_block = s_col % Nhs;
          for (long c_row=0;c_row<Nc;c_row++) {
            for (long c_col=0;c_col<Nc;c_col++) {
              int i = s_row_block * Nc + c_row;
              int j = s_col_block * Nc + c_col;
              if(i == j)
                clover_eigen(s_row*Nc+c_row, s_col*Nc+c_col) = static_cast<ComplexD>(TensorRemove(diagonal_tmp()(block)(i)));
              else
                clover_eigen(s_row*Nc+c_row, s_col*Nc+c_col) = static_cast<ComplexD>(TensorRemove(triangle_elem(triangle_tmp, block, i, j)));
            }
          }
        }
      }

      clover_inv_eigen = clover_eigen.inverse();

      for (long s_row=0;s_row<Ns;s_row++) {
        for (long s_col=0;s_col<Ns;s_col++) {
          if(abs(s_row - s_col) > 1 || s_row + s_col == 3) continue;
          int block       = s_row / Nhs;
          int s_row_block = s_row % Nhs;
          int s_col_block = s_col % Nhs;
          for (long c_row=0;c_row<Nc;c_row++) {
            for (long c_col=0;c_col<Nc;c_col++) {
              int i = s_row_block * Nc + c_row;
              int j = s_col_block * Nc + c_col;
              if(i == j)
                diagonal_inv_tmp()(block)(i) = clover_inv_eigen(s_row*Nc+c_row, s_col*Nc+c_col);
              else if(i < j)
                triangle_inv_tmp()(block)(triangle_index(i, j)) = clover_inv_eigen(s_row*Nc+c_row, s_col*Nc+c_col);
              else
                continue;
            }
          }
        }
      }

      pokeLocalSite(diagonal_inv_tmp, diagonalInv_v, lcoor);
      pokeLocalSite(triangle_inv_tmp, triangleInv_v, lcoor);
    });
  }

  static void ConvertLayout(const CloverField&   full,
                            CloverDiagonalField& diagonal,
                            CloverTriangleField& triangle) {
    conformable(full, diagonal);
    conformable(full, triangle);

    diagonal.Checkerboard() = full.Checkerboard();
    triangle.Checkerboard() = full.Checkerboard();

    autoView(full_v,     full,     AcceleratorRead);
    autoView(diagonal_v, diagonal, AcceleratorWrite);
    autoView(triangle_v, triangle, AcceleratorWrite);

    // NOTE: this function cannot be 'private' since nvcc forbids this for kernels
    accelerator_for(ss, full.Grid()->oSites(), 1, {
      for(int s_row = 0; s_row < Ns; s_row++) {
        for(int s_col = 0; s_col < Ns; s_col++) {
          if(abs(s_row - s_col) > 1 || s_row + s_col == 3) continue;
          int block       = s_row / Nhs;
          int s_row_block = s_row % Nhs;
          int s_col_block = s_col % Nhs;
          for(int c_row = 0; c_row < Nc; c_row++) {
            for(int c_col = 0; c_col < Nc; c_col++) {
              int i = s_row_block * Nc + c_row;
              int j = s_col_block * Nc + c_col;
              if(i == j)
                diagonal_v[ss]()(block)(i) = full_v[ss]()(s_row, s_col)(c_row, c_col);
              else if(i < j)
                triangle_v[ss]()(block)(triangle_index(i, j)) = full_v[ss]()(s_row, s_col)(c_row, c_col);
              else
                continue;
            }
          }
        }
      }
    });
  }


  static void ConvertLayout(const CloverDiagonalField& diagonal,
                            const CloverTriangleField& triangle,
                            CloverField&               full) {
    conformable(full, diagonal);
    conformable(full, triangle);

    full.Checkerboard() = diagonal.Checkerboard();

    full = Zero();

    autoView(diagonal_v, diagonal, AcceleratorRead);
    autoView(triangle_v, triangle, AcceleratorRead);
    autoView(full_v,     full,     AcceleratorWrite);

    // NOTE: this function cannot be 'private' since nvcc forbids this for kernels
    accelerator_for(ss, full.Grid()->oSites(), 1, {
      for(int s_row = 0; s_row < Ns; s_row++) {
        for(int s_col = 0; s_col < Ns; s_col++) {
          if(abs(s_row - s_col) > 1 || s_row + s_col == 3) continue;
          int block       = s_row / Nhs;
          int s_row_block = s_row % Nhs;
          int s_col_block = s_col % Nhs;
          for(int c_row = 0; c_row < Nc; c_row++) {
            for(int c_col = 0; c_col < Nc; c_col++) {
              int i = s_row_block * Nc + c_row;
              int j = s_col_block * Nc + c_col;
              if(i == j)
                full_v[ss]()(s_row, s_col)(c_row, c_col) = diagonal_v[ss]()(block)(i);
              else
                full_v[ss]()(s_row, s_col)(c_row, c_col) = triangle_elem(triangle_v[ss], block, i, j);
            }
          }
        }
      }
    });
  }

  static void ModifyBoundaries(CloverDiagonalField& diagonal, CloverTriangleField& triangle, RealD csw_t, RealD cF, RealD diag_mass) {
    // Checks/grid
    double t0 = usecond();
    conformable(diagonal, triangle);
    GridBase* grid = diagonal.Grid();

    // Determine the boundary coordinates/sites
    double t1 = usecond();
    int t_dir = Nd - 1;
    Lattice<iScalar<vInteger>> t_coor(grid);
    LatticeCoordinate(t_coor, t_dir);
    int T = grid->GlobalDimensions()[t_dir];

    // Set off-diagonal parts at boundary to zero -- OK
    double t2 = usecond();
    CloverTriangleField zeroTriangle(grid);
    zeroTriangle.Checkerboard() = triangle.Checkerboard();
    zeroTriangle = Zero();
    triangle = where(t_coor == 0,   zeroTriangle, triangle);
    triangle = where(t_coor == T-1, zeroTriangle, triangle);

    // Set diagonal to unity (scaled correctly) -- OK
    double t3 = usecond();
    CloverDiagonalField tmp(grid);
    tmp.Checkerboard() = diagonal.Checkerboard();
    tmp                = -1.0 * csw_t + diag_mass;
    diagonal           = where(t_coor == 0,   tmp, diagonal);
    diagonal           = where(t_coor == T-1, tmp, diagonal);

    // Correct values next to boundary
    double t4 = usecond();
    if(cF != 1.0) {
      tmp = cF - 1.0;
      tmp += diagonal;
      diagonal = where(t_coor == 1,   tmp, diagonal);
      diagonal = where(t_coor == T-2, tmp, diagonal);
    }

    // Report timings
    double t5 = usecond();
#if 0
    std::cout << GridLogMessage << "CompactWilsonCloverHelpers::ModifyBoundaries timings:"
              << " checks = "          << (t1 - t0) / 1e6
              << ", coordinate = "     << (t2 - t1) / 1e6
              << ", off-diag zero = "  << (t3 - t2) / 1e6
              << ", diagonal unity = " << (t4 - t3) / 1e6
              << ", near-boundary = "  << (t5 - t4) / 1e6
              << ", total = "          << (t5 - t0) / 1e6
              << std::endl;
#endif
  }

  template<class Field, class Mask>
  static strong_inline void ApplyBoundaryMask(Field& f, const Mask& m) {
    conformable(f, m);
    auto grid  = f.Grid();
    const int Nsite = grid->oSites();
    const int Nsimd = grid->Nsimd();
    autoView(f_v, f, AcceleratorWrite);
    autoView(m_v, m, AcceleratorRead);
    // NOTE: this function cannot be 'private' since nvcc forbids this for kernels
    accelerator_for(ss, Nsite, Nsimd, {
      coalescedWrite(f_v[ss], m_v(ss) * f_v(ss));
    });
  }

  template<class MaskField>
  static void SetupMasks(MaskField& full, MaskField& even, MaskField& odd) {
    assert(even.Grid()->_isCheckerBoarded && even.Checkerboard() == Even);
    assert(odd.Grid()->_isCheckerBoarded  && odd.Checkerboard()  == Odd);
    assert(!full.Grid()->_isCheckerBoarded);

    GridBase* grid = full.Grid();
    int t_dir = Nd-1;
    Lattice<iScalar<vInteger>> t_coor(grid);
    LatticeCoordinate(t_coor, t_dir);
    int T = grid->GlobalDimensions()[t_dir];

    MaskField zeroMask(grid); zeroMask = Zero();
    full = 1.0;
    full = where(t_coor == 0,   zeroMask, full);
    full = where(t_coor == T-1, zeroMask, full);

    pickCheckerboard(Even, even, full);
    pickCheckerboard(Odd,  odd,  full);
  }
};

NAMESPACE_END(Grid);
