/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: tests/Test_extended_meson_field.cc

Copyright (C) 2015-2025

Author: Peter Boyle <pboyle@bnl.gov>
Author: Masaaki Tomii <masaaki.tomii@uconn.edu> (original Hadrons kernels)

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
#include "disable_tests_without_instantiations.h"
#ifdef ENABLE_FERMION_INSTANTIATIONS

#include <Grid/Grid.h>
#include <Grid/qcd/utils/A2Autils.h>

using namespace Grid;

typedef WilsonImplD FImpl;
typedef typename FImpl::FermionField FermionField;
typedef typename FImpl::SiteSpinor   vobj;
typedef typename vobj::scalar_type   scalar_type;
typedef typename vobj::vector_type   vector_type;
typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
typedef iSpinColourVector<vector_type> SpinColourVector_v;
typedef iSpinMatrix<vector_type>       SpinMatrix_v;
typedef iSinglet<vector_type>          Scalar_v;
typedef iSinglet<scalar_type>          Scalar_s;
typedef Lattice<SpinColourMatrix_v>    PropagatorField;

// CPU reference + optionally batched GEMM spatial sum, ported from
// Hadrons/Modules/MContraction/A2AExtendedMesonField.hpp
// (M. Tomii, mtomii/Hadrons:local-2025-edits). Hadrons infrastructure removed.
// thread_for / CpuRead / orthogdim=3 preserved.
class A2AExtendedMesonFieldRef
{
public:
  // result is indexed [nt][N_i][N_j].
  // use_blas=true replaces the scalar spatial accumulation with A2ASpatialSum.
  static void compute(
      Eigen::Tensor<ComplexD, 3> &result,
      const std::vector<FermionField> &left,
      const std::vector<FermionField> &right,
      const std::vector<FermionField> &loop1,
      const std::vector<FermionField> &loop2,
      const std::vector<Gamma::Algebra> &gamma1,
      const std::vector<Gamma::Algebra> &gamma2,
      int type,
      bool use_blas = false)
  {
    GridBase *grid = left[0].Grid();

    const int orthogdim = 3;
    int rd     = grid->_rdimensions[orthogdim];
    int ld     = grid->_ldimensions[orthogdim];
    int Nd     = grid->_ndimension;
    int Nsimd  = grid->Nsimd();

    int nt  = result.dimension(0);
    int N_i = (int)left.size();
    int N_j = (int)right.size();

    std::string tag = std::string(use_blas ? "[blas" : "[ref ") + " type=" + std::to_string(type) + "]";
    auto Tms = [](double us) { return us * 1e-3; };
    double t0;

    // ------------------------------------------------------------------
    // Loop propagator: sum_k outerProduct(loop1[k], loop2[k])
    // ------------------------------------------------------------------
    t0 = usecond();
    PropagatorField loop(grid);
    loop = Zero();
    for (unsigned int k = 0; k < loop1.size(); ++k)
      loop += outerProduct(loop1[k], loop2[k]);
    std::cout << GridLogMessage << tag << " loop_build:      " << Tms(usecond()-t0) << " ms\n";

    // ------------------------------------------------------------------
    // Pack conjugated left vectors
    // ------------------------------------------------------------------
    t0 = usecond();
    std::vector<FermionField> leftv(N_i, grid);
    for (int i = 0; i < N_i; i++)
      leftv[i] = conjugate(left[i]);
    std::cout << GridLogMessage << tag << " pack_left:       " << Tms(usecond()-t0) << " ms\n";

    // ------------------------------------------------------------------
    // Per-site loop contraction into PropagatorField tloop (type-dependent)
    // ------------------------------------------------------------------
    t0 = usecond();
    PropagatorField tloop(grid);
    tloop = Zero();
    {
      autoView(tloopv, tloop, CpuWrite);
      autoView(loopv,  loop,  CpuRead);
      if (type == 0) {
        thread_for(ss, grid->oSites(), {
          for (int s1 = 0; s1 < Ns; ++s1)
          for (int s2 = 0; s2 < Ns; ++s2)
            tloopv[ss]()(s1,s2)(0,0) = loopv[ss]()(s1,s2)(0,0)
                                     + loopv[ss]()(s1,s2)(1,1)
                                     + loopv[ss]()(s1,s2)(2,2);
        });
      }
      if (type == 1) {
        thread_for(ss, grid->oSites(), {
          tloopv[ss] = Zero();
          for (int mu = 0; mu < (int)gamma1.size(); ++mu)
            tloopv[ss] = tloopv[ss] + Gamma(gamma1[mu]) * loopv[ss] * Gamma(gamma2[mu]);
        });
      }
      if (type == 2) {
        thread_for(ss, grid->oSites(), {
          tloopv[ss] = Zero();
          for (int mu = 0; mu < (int)gamma2.size(); ++mu) {
            SpinColourMatrix_v tmp = Gamma(gamma2[mu]) * loopv[ss];
            int s1 = mu / Ns;
            int s2 = mu % Ns;
            for (int c1 = 0; c1 < Nc; ++c1)
            for (int c2 = 0; c2 < Nc; ++c2)
              tloopv[ss]()(s1,s2)(c1,c2) = tmp()(0,0)(c1,c2) + tmp()(1,1)(c1,c2)
                                          + tmp()(2,2)(c1,c2) + tmp()(3,3)(c1,c2);
          }
        });
      }
      if (type == 3) {
        thread_for(ss, grid->oSites(), {
          SpinMatrix_v spinLoop = Zero();
          for (int s1 = 0; s1 < Ns; ++s1)
          for (int s2 = 0; s2 < Ns; ++s2)
            spinLoop()(s1,s2)() = loopv[ss]()(s1,s2)(0,0)
                                + loopv[ss]()(s1,s2)(1,1)
                                + loopv[ss]()(s1,s2)(2,2);
          tloopv[ss] = Zero();
          for (int mu = 0; mu < (int)gamma1.size(); ++mu) {
            SpinMatrix_v tmp2 = Gamma(gamma1[mu]) * spinLoop * Gamma(gamma2[mu]);
            for (int s1 = 0; s1 < Ns; ++s1)
            for (int s2 = 0; s2 < Ns; ++s2)
              tloopv[ss]()(s1,s2)(0,0) = tloopv[ss]()(s1,s2)(0,0) + tmp2()(s1,s2)();
          }
        });
      }
    }
    std::cout << GridLogMessage << tag << " tloop:           " << Tms(usecond()-t0) << " ms\n";

    // Select addLoopRight kernel for this type
    std::function<void(SpinColourVector_v &,
                       const SpinColourMatrix_v &,
                       const SpinColourVector_v &,
                       const std::vector<Gamma::Algebra> &,
                       const std::vector<Gamma::Algebra> &)> addLoopRight;

    if (type == 0) {
      addLoopRight = [](SpinColourVector_v &lR,
                        const SpinColourMatrix_v &loopm,
                        const SpinColourVector_v &rightv,
                        const std::vector<Gamma::Algebra> &g1,
                        const std::vector<Gamma::Algebra> &g2) {
        SpinMatrix_v spinLoop = Zero();
        for (int s1 = 0; s1 < Ns; ++s1)
        for (int s2 = 0; s2 < Ns; ++s2)
          spinLoop()(s1,s2)() = loopm()(s1,s2)(0,0);
        for (int mu = 0; mu < (int)g1.size(); ++mu) {
          SpinMatrix_v GLoop = Gamma(g2[mu]) * spinLoop;
          auto trGLoop = GLoop()(0,0)() + GLoop()(1,1)() + GLoop()(2,2)() + GLoop()(3,3)();
          SpinColourVector_v Grightv = Gamma(g1[mu]) * rightv;
          for (int s = 0; s < Ns; ++s)
          for (int c = 0; c < Nc; ++c)
            lR()(s)(c) += Grightv()(s)(c) * trGLoop;
        }
      };
    }
    if (type == 1) {
      addLoopRight = [](SpinColourVector_v &lR,
                        const SpinColourMatrix_v &loopm,
                        const SpinColourVector_v &rightv,
                        const std::vector<Gamma::Algebra> &g1,
                        const std::vector<Gamma::Algebra> &g2) {
        for (int s = 0; s < Ns; ++s)
        for (int c = 0; c < Nc; ++c) {
          lR()(s)(c)
            += loopm()(s,0)(c,0) * rightv()(0)(0)
            +  loopm()(s,0)(c,1) * rightv()(0)(1)
            +  loopm()(s,0)(c,2) * rightv()(0)(2)
            +  loopm()(s,1)(c,0) * rightv()(1)(0)
            +  loopm()(s,1)(c,1) * rightv()(1)(1)
            +  loopm()(s,1)(c,2) * rightv()(1)(2)
            +  loopm()(s,2)(c,0) * rightv()(2)(0)
            +  loopm()(s,2)(c,1) * rightv()(2)(1)
            +  loopm()(s,2)(c,2) * rightv()(2)(2)
            +  loopm()(s,3)(c,0) * rightv()(3)(0)
            +  loopm()(s,3)(c,1) * rightv()(3)(1)
            +  loopm()(s,3)(c,2) * rightv()(3)(2);
        }
      };
    }
    if (type == 2) {
      addLoopRight = [](SpinColourVector_v &lR,
                        const SpinColourMatrix_v &loopm,
                        const SpinColourVector_v &rightv,
                        const std::vector<Gamma::Algebra> &g1,
                        const std::vector<Gamma::Algebra> &g2) {
        for (int mu = 0; mu < (int)g1.size(); ++mu) {
          int s1 = mu / Ns;
          int s2 = mu % Ns;
          SpinColourVector_v Grightv = Gamma(g1[mu]) * rightv;
          for (int s = 0; s < Ns; ++s)
          for (int c = 0; c < Nc; ++c)
            lR()(s)(c) += loopm()(s1,s2)(c,0) * Grightv()(s)(0)
                       +  loopm()(s1,s2)(c,1) * Grightv()(s)(1)
                       +  loopm()(s1,s2)(c,2) * Grightv()(s)(2);
        }
      };
    }
    if (type == 3) {
      addLoopRight = [](SpinColourVector_v &lR,
                        const SpinColourMatrix_v &loopm,
                        const SpinColourVector_v &rightv,
                        const std::vector<Gamma::Algebra> &g1,
                        const std::vector<Gamma::Algebra> &g2) {
        for (int s = 0; s < Ns; ++s)
        for (int c = 0; c < Nc; ++c)
          lR()(s)(c) += loopm()(s,0)(0,0) * rightv()(0)(c)
                     +  loopm()(s,1)(0,0) * rightv()(1)(c)
                     +  loopm()(s,2)(0,0) * rightv()(2)(c)
                     +  loopm()(s,3)(0,0) * rightv()(3)(c);
      };
    }

    // ------------------------------------------------------------------
    // Pack loopRight[j] = type-kernel(tloop, right[j]) per site
    // ------------------------------------------------------------------
    t0 = usecond();
    std::vector<FermionField> loopRight(N_j, grid);
    {
      autoView(tlv, tloop, CpuRead);
      for (int j = 0; j < N_j; j++) {
        loopRight[j] = Zero();
        autoView(lRv, loopRight[j], CpuWrite);
        autoView(rv,  right[j],     CpuRead);
        thread_for(ss, grid->oSites(), {
          addLoopRight(lRv[ss], tlv[ss], rv[ss], gamma1, gamma2);
        });
      }
    }
    std::cout << GridLogMessage << tag << " pack_loopright:  " << Tms(usecond()-t0) << " ms\n";

    if (use_blas) {
      // ------------------------------------------------------------------
      // BLAS path: A2ASpatialSum (Allocate + PackLeft + PackRight + Sum)
      // ------------------------------------------------------------------
      A2ASpatialSum<SpinColourVector_v> spatial_sum;
      double t_blas_start = usecond();

      t0 = usecond();
      spatial_sum.Allocate(N_i, N_j, grid);
      std::cout << GridLogMessage << tag << " Allocate:        " << Tms(usecond()-t0) << " ms\n";

      t0 = usecond();
      spatial_sum.PackLeft(leftv);
      std::cout << GridLogMessage << tag << " PackLeft:        " << Tms(usecond()-t0) << " ms\n";

      t0 = usecond();
      spatial_sum.PackRight(loopRight);
      std::cout << GridLogMessage << tag << " PackRight:       " << Tms(usecond()-t0) << " ms\n";

      t0 = usecond();
      spatial_sum.Sum(result);
      std::cout << GridLogMessage << tag << " Sum (GEMM+MPI):  " << Tms(usecond()-t0) << " ms\n";

      std::cout << GridLogMessage << tag << " A2ASpatialSum:   " << Tms(usecond()-t_blas_start) << " ms  [TOTAL]\n";

    } else {
      // ------------------------------------------------------------------
      // Reference path: SIMD spatial sum + scalar extraction
      // ------------------------------------------------------------------
      int MFrvol = rd * N_i * N_j;
      int MFlvol = ld * N_i * N_j;

      Vector<Scalar_v> lvSum(MFrvol);
      thread_for(r, MFrvol, { lvSum[r] = Zero(); });

      t0 = usecond();
      {
        int e1     = grid->_slice_nblock[orthogdim];
        int e2     = grid->_slice_block [orthogdim];
        int stride = grid->_slice_stride[orthogdim];
        using LView = decltype(leftv[0].View(CpuRead));
        using RView = decltype(loopRight[0].View(CpuRead));
        std::vector<LView> lv_views;
        std::vector<RView> lr_views;
        lv_views.reserve(N_i);
        lr_views.reserve(N_j);
        for (int i = 0; i < N_i; i++) lv_views.push_back(leftv[i].View(CpuRead));
        for (int j = 0; j < N_j; j++) lr_views.push_back(loopRight[j].View(CpuRead));
        thread_for(r, rd, {
          int so   = r * grid->_ostride[orthogdim];
          int base = N_i * N_j * r;
          for (int n = 0; n < e1; n++)
          for (int b = 0; b < e2; b++) {
            int ss = so + n * stride + b;
            for (int ii = 0; ii < N_i; ii++)
            for (int jj = 0; jj < N_j; jj++) {
              int idx = jj + N_j * ii + base;
              for (int s = 0; s < Ns; ++s)
              for (int c = 0; c < Nc; ++c)
                lvSum[idx]()()() += lv_views[ii][ss]()(s)(c) * lr_views[jj][ss]()(s)(c);
            }
          }
        });
        for (auto &v : lv_views) v.ViewClose();
        for (auto &v : lr_views) v.ViewClose();
      }
      std::cout << GridLogMessage << tag << " spatial_sum:     " << Tms(usecond()-t0) << " ms\n";

      Vector<Scalar_s> lsSum(MFlvol);
      thread_for(r, MFlvol, { lsSum[r] = scalar_type(0.0); });

      t0 = usecond();
      thread_for(rt, rd, {
        Coordinate icoor(Nd);
        ExtractBuffer<Scalar_s> extracted(Nsimd);
        for (int ii = 0; ii < N_i; ii++)
        for (int jj = 0; jj < N_j; jj++) {
          int ij_rdx = jj + N_j * (ii + N_i * rt);
          extract(lvSum[ij_rdx], extracted);
          for (int idx = 0; idx < Nsimd; idx++) {
            grid->iCoorFromIindex(icoor, idx);
            int ldx    = rt + icoor[orthogdim] * rd;
            int ij_ldx = jj + N_j * (ii + N_i * ldx);
            lsSum[ij_ldx] = lsSum[ij_ldx] + extracted[idx];
          }
        }
      });
      std::cout << GridLogMessage << tag << " simd_extract:    " << Tms(usecond()-t0) << " ms\n";

      int pd = grid->_processors[orthogdim];
      int pc = grid->_processor_coor[orthogdim];

      t0 = usecond();
      Vector<ComplexD> cache(nt * N_i * N_j, ComplexD(0.0));
      for (int lt = 0; lt < ld; lt++)
      for (int pt = 0; pt < pd; pt++) {
        int t = lt + pt * ld;
        for (int ii = 0; ii < N_i; ii++)
        for (int jj = 0; jj < N_j; jj++) {
          if (pt == pc) {
            int ij_ldx = jj + N_j * (ii + N_i * lt);
            cache[jj + N_j * (ii + N_i * t)] = lsSum[ij_ldx]()()();
          }
        }
      }
      grid->GlobalSumVector(cache.data(), nt * N_i * N_j);
      std::cout << GridLogMessage << tag << " globalsum:       " << Tms(usecond()-t0) << " ms\n";

      for (int t  = 0; t  < nt;  t++)
      for (int ii = 0; ii < N_i; ii++)
      for (int jj = 0; jj < N_j; jj++)
        result(t, ii, jj) = cache[jj + N_j * (ii + N_i * t)];
    }
  }
};

// ================================================================
// Free-function GPU kernels — accelerator_for, v(ss) reads,
// coalescedWrite writes, vobj-level arithmetic throughout.
// Gamma arrays passed as Vector<Gamma::Algebra> (UVM).
// ================================================================

void A2ALoopPropagator(PropagatorField &loop,
                       const std::vector<FermionField> &loop1,
                       const std::vector<FermionField> &loop2)
{
  loop = Zero();
  for (unsigned int k = 0; k < loop1.size(); ++k)
    loop += outerProduct(loop1[k], loop2[k]);
}

void A2APackLeftConjugated(FermionField &out, const FermionField &in)
{
  autoView(outv, out, AcceleratorWrite);
  autoView(inv,  in,  AcceleratorRead);
  uint64_t Osites = in.Grid()->oSites();
  int Nsimd = SpinColourVector_v::Nsimd();
  accelerator_for(ss, Osites, Nsimd, {
    coalescedWrite(outv[ss], conjugate(inv(ss)));
  });
}

// Type 0: colour-trace stored in (s1,s2)(0,0)
void A2ALoopLeftContractionType0(PropagatorField &tloop, const PropagatorField &loop)
{
  autoView(tloopv, tloop, AcceleratorWrite);
  autoView(loopv,  loop,  AcceleratorRead);
  uint64_t Osites = loop.Grid()->oSites();
  int Nsimd = SpinColourMatrix_v::Nsimd();
  accelerator_for(ss, Osites, Nsimd, {
    auto l = loopv(ss);
    auto tmp = l; tmp = Zero();
    for (int s1 = 0; s1 < Ns; ++s1)
    for (int s2 = 0; s2 < Ns; ++s2)
      tmp()(s1,s2)(0,0) = l()(s1,s2)(0,0) + l()(s1,s2)(1,1) + l()(s1,s2)(2,2);
    coalescedWrite(tloopv[ss], tmp);
  });
}

// Type 1: tloop = sum_mu Gamma(g1[mu]) * loop * Gamma(g2[mu])
void A2ALoopLeftContractionType1(PropagatorField &tloop, const PropagatorField &loop,
                                  const Vector<Gamma::Algebra> &gamma1,
                                  const Vector<Gamma::Algebra> &gamma2)
{
  int ng = (int)gamma1.size();
  const Gamma::Algebra *g1 = gamma1.data();
  const Gamma::Algebra *g2 = gamma2.data();
  autoView(tloopv, tloop, AcceleratorWrite);
  autoView(loopv,  loop,  AcceleratorRead);
  uint64_t Osites = loop.Grid()->oSites();
  int Nsimd = SpinColourMatrix_v::Nsimd();
  accelerator_for(ss, Osites, Nsimd, {
    auto l = loopv(ss);
    auto tmp = l; tmp = Zero();
    for (int mu = 0; mu < ng; ++mu)
      tmp = tmp + Gamma(g1[mu]) * l * Gamma(g2[mu]);
    coalescedWrite(tloopv[ss], tmp);
  });
}

// Type 2: for mu=[0..ng), s1=mu/Ns, s2=mu%Ns;
//         tloop(s1,s2)(c1,c2) = Tr_spin( Gamma(g2[mu]) * loop )(c1,c2)
void A2ALoopLeftContractionType2(PropagatorField &tloop, const PropagatorField &loop,
                                  const Vector<Gamma::Algebra> &gamma2)
{
  int ng = (int)gamma2.size();
  const Gamma::Algebra *g2 = gamma2.data();
  autoView(tloopv, tloop, AcceleratorWrite);
  autoView(loopv,  loop,  AcceleratorRead);
  uint64_t Osites = loop.Grid()->oSites();
  int Nsimd = SpinColourMatrix_v::Nsimd();
  accelerator_for(ss, Osites, Nsimd, {
    auto l = loopv(ss);
    auto tmp = l; tmp = Zero();
    for (int mu = 0; mu < ng; ++mu) {
      auto gtmp = Gamma(g2[mu]) * l;
      int s1 = mu / Ns;
      int s2 = mu % Ns;
      for (int c1 = 0; c1 < Nc; ++c1)
      for (int c2 = 0; c2 < Nc; ++c2)
        tmp()(s1,s2)(c1,c2) = gtmp()(0,0)(c1,c2) + gtmp()(1,1)(c1,c2)
                            + gtmp()(2,2)(c1,c2) + gtmp()(3,3)(c1,c2);
    }
    coalescedWrite(tloopv[ss], tmp);
  });
}

// Type 3: colour-trace → spin matrix → sum_mu G1*spinLoop*G2 stored in (s1,s2)(0,0)
void A2ALoopLeftContractionType3(PropagatorField &tloop, const PropagatorField &loop,
                                  const Vector<Gamma::Algebra> &gamma1,
                                  const Vector<Gamma::Algebra> &gamma2)
{
  int ng = (int)gamma1.size();
  const Gamma::Algebra *g1 = gamma1.data();
  const Gamma::Algebra *g2 = gamma2.data();
  autoView(tloopv, tloop, AcceleratorWrite);
  autoView(loopv,  loop,  AcceleratorRead);
  uint64_t Osites = loop.Grid()->oSites();
  int Nsimd = SpinColourMatrix_v::Nsimd();
  accelerator_for(ss, Osites, Nsimd, {
    typedef decltype(coalescedRead(loopv[0])) calcSCMatrix;
    typedef iSpinMatrix<typename calcSCMatrix::vector_type> calcSpinMatrix;
    auto l = loopv(ss);
    calcSpinMatrix spinLoop; spinLoop = Zero();
    for (int s1 = 0; s1 < Ns; ++s1)
    for (int s2 = 0; s2 < Ns; ++s2)
      spinLoop()(s1,s2)() = l()(s1,s2)(0,0) + l()(s1,s2)(1,1) + l()(s1,s2)(2,2);
    auto tmp = l; tmp = Zero();
    for (int mu = 0; mu < ng; ++mu) {
      calcSpinMatrix tmp2 = Gamma(g1[mu]) * spinLoop * Gamma(g2[mu]);
      for (int s1 = 0; s1 < Ns; ++s1)
      for (int s2 = 0; s2 < Ns; ++s2)
        tmp()(s1,s2)(0,0) = tmp()(s1,s2)(0,0) + tmp2()(s1,s2)();
    }
    coalescedWrite(tloopv[ss], tmp);
  });
}

// Type 0: loopRight = sum_mu Tr(G2*spinLoop) * G1*right
//   where spinLoop(s1,s2) = tloop(s1,s2)(0,0)
void A2ALoopRightContractionType0(FermionField &loopRight,
                                   const PropagatorField &tloop,
                                   const FermionField &right,
                                   const Vector<Gamma::Algebra> &gamma1,
                                   const Vector<Gamma::Algebra> &gamma2)
{
  int ng = (int)gamma1.size();
  const Gamma::Algebra *g1 = gamma1.data();
  const Gamma::Algebra *g2 = gamma2.data();
  autoView(lRv, loopRight, AcceleratorWrite);
  autoView(tlv, tloop,     AcceleratorRead);
  autoView(rv,  right,     AcceleratorRead);
  uint64_t Osites = right.Grid()->oSites();
  int Nsimd = SpinColourVector_v::Nsimd();
  accelerator_for(ss, Osites, Nsimd, {
    typedef decltype(coalescedRead(rv[0]))  calcSCVector;
    typedef decltype(coalescedRead(tlv[0])) calcSCMatrix;
    typedef iSpinMatrix<typename calcSCMatrix::vector_type> calcSpinMatrix;
    auto loopm  = tlv(ss);
    auto rightv = rv(ss);
    calcSpinMatrix spinLoop; spinLoop = Zero();
    for (int s1 = 0; s1 < Ns; ++s1)
    for (int s2 = 0; s2 < Ns; ++s2)
      spinLoop()(s1,s2)() = loopm()(s1,s2)(0,0);
    calcSCVector lR; lR = Zero();
    for (int mu = 0; mu < ng; ++mu) {
      auto GLoop   = Gamma(g2[mu]) * spinLoop;
      auto trGLoop = GLoop()(0,0)() + GLoop()(1,1)() + GLoop()(2,2)() + GLoop()(3,3)();
      auto Grightv = Gamma(g1[mu]) * rightv;
      for (int s = 0; s < Ns; ++s)
      for (int c = 0; c < Nc; ++c)
        lR()(s)(c) = lR()(s)(c) + Grightv()(s)(c) * trGLoop;
    }
    coalescedWrite(lRv[ss], lR);
  });
}

// Type 1: loopRight = tloop * right  (SpinColourMatrix * SpinColourVector)
void A2ALoopRightContractionType1(FermionField &loopRight,
                                   const PropagatorField &tloop,
                                   const FermionField &right)
{
  autoView(lRv, loopRight, AcceleratorWrite);
  autoView(tlv, tloop,     AcceleratorRead);
  autoView(rv,  right,     AcceleratorRead);
  uint64_t Osites = right.Grid()->oSites();
  int Nsimd = SpinColourVector_v::Nsimd();
  accelerator_for(ss, Osites, Nsimd, {
    coalescedWrite(lRv[ss], tlv(ss) * rv(ss));
  });
}

// Type 2: loopRight(s)(c) = sum_{mu,c'} tloop(s1,s2)(c,c') * (G(g1[mu])*right)(s)(c')
void A2ALoopRightContractionType2(FermionField &loopRight,
                                   const PropagatorField &tloop,
                                   const FermionField &right,
                                   const Vector<Gamma::Algebra> &gamma1)
{
  int ng = (int)gamma1.size();
  const Gamma::Algebra *g1 = gamma1.data();
  autoView(lRv, loopRight, AcceleratorWrite);
  autoView(tlv, tloop,     AcceleratorRead);
  autoView(rv,  right,     AcceleratorRead);
  uint64_t Osites = right.Grid()->oSites();
  int Nsimd = SpinColourVector_v::Nsimd();
  accelerator_for(ss, Osites, Nsimd, {
    typedef decltype(coalescedRead(rv[0])) calcSCVector;
    auto loopm  = tlv(ss);
    auto rightv = rv(ss);
    calcSCVector lR; lR = Zero();
    for (int mu = 0; mu < ng; ++mu) {
      int s1 = mu / Ns;
      int s2 = mu % Ns;
      auto Grightv = Gamma(g1[mu]) * rightv;
      for (int s = 0; s < Ns; ++s)
      for (int c = 0; c < Nc; ++c)
        lR()(s)(c) = lR()(s)(c)
                   + loopm()(s1,s2)(c,0) * Grightv()(s)(0)
                   + loopm()(s1,s2)(c,1) * Grightv()(s)(1)
                   + loopm()(s1,s2)(c,2) * Grightv()(s)(2);
    }
    coalescedWrite(lRv[ss], lR);
  });
}

// Type 3: loopRight(s)(c) = sum_{s'} tloop(s,s')(0,0) * right(s')(c)
void A2ALoopRightContractionType3(FermionField &loopRight,
                                   const PropagatorField &tloop,
                                   const FermionField &right)
{
  autoView(lRv, loopRight, AcceleratorWrite);
  autoView(tlv, tloop,     AcceleratorRead);
  autoView(rv,  right,     AcceleratorRead);
  uint64_t Osites = right.Grid()->oSites();
  int Nsimd = SpinColourVector_v::Nsimd();
  accelerator_for(ss, Osites, Nsimd, {
    typedef decltype(coalescedRead(rv[0])) calcSCVector;
    auto loopm  = tlv(ss);
    auto rightv = rv(ss);
    calcSCVector lR; lR = Zero();
    for (int s = 0; s < Ns; ++s)
    for (int c = 0; c < Nc; ++c)
      lR()(s)(c) = loopm()(s,0)(0,0) * rightv()(0)(c)
                 + loopm()(s,1)(0,0) * rightv()(1)(c)
                 + loopm()(s,2)(0,0) * rightv()(2)(c)
                 + loopm()(s,3)(0,0) * rightv()(3)(c);
    coalescedWrite(lRv[ss], lR);
  });
}

// ================================================================
// GPU-offloaded extended meson field: accelerator_for contractions
// + A2ASpatialSum GEMM spatial reduction.
// ================================================================
class A2AExtendedMesonFieldGPU
{
public:
  static void compute(
      Eigen::Tensor<ComplexD, 3> &result,
      const std::vector<FermionField> &left,
      const std::vector<FermionField> &right,
      const std::vector<FermionField> &loop1,
      const std::vector<FermionField> &loop2,
      const std::vector<Gamma::Algebra> &gamma1_in,
      const std::vector<Gamma::Algebra> &gamma2_in,
      int type)
  {
    GridBase *grid = left[0].Grid();
    int N_i = (int)left.size();
    int N_j = (int)right.size();

    std::string tag = std::string("[gpu  type=") + std::to_string(type) + "]";
    auto Tms = [](double us) { return us * 1e-3; };
    double t0;

    Vector<Gamma::Algebra> gamma1(gamma1_in.begin(), gamma1_in.end());
    Vector<Gamma::Algebra> gamma2(gamma2_in.begin(), gamma2_in.end());

    t0 = usecond();
    PropagatorField loop(grid);
    A2ALoopPropagator(loop, loop1, loop2);
    std::cout << GridLogMessage << tag << " loop_build:      " << Tms(usecond()-t0) << " ms\n";

    t0 = usecond();
    std::vector<FermionField> leftv(N_i, grid);
    for (int i = 0; i < N_i; i++)
      A2APackLeftConjugated(leftv[i], left[i]);
    std::cout << GridLogMessage << tag << " pack_left:       " << Tms(usecond()-t0) << " ms\n";

    t0 = usecond();
    PropagatorField tloop(grid);
    tloop = Zero();
    switch (type) {
    case 0: A2ALoopLeftContractionType0(tloop, loop);                break;
    case 1: A2ALoopLeftContractionType1(tloop, loop, gamma1, gamma2); break;
    case 2: A2ALoopLeftContractionType2(tloop, loop, gamma2);         break;
    case 3: A2ALoopLeftContractionType3(tloop, loop, gamma1, gamma2); break;
    }
    std::cout << GridLogMessage << tag << " tloop:           " << Tms(usecond()-t0) << " ms\n";

    t0 = usecond();
    std::vector<FermionField> loopRight(N_j, grid);
    for (int j = 0; j < N_j; j++) {
      switch (type) {
      case 0: A2ALoopRightContractionType0(loopRight[j], tloop, right[j], gamma1, gamma2); break;
      case 1: A2ALoopRightContractionType1(loopRight[j], tloop, right[j]);                 break;
      case 2: A2ALoopRightContractionType2(loopRight[j], tloop, right[j], gamma1);         break;
      case 3: A2ALoopRightContractionType3(loopRight[j], tloop, right[j]);                 break;
      }
    }
    std::cout << GridLogMessage << tag << " pack_loopright:  " << Tms(usecond()-t0) << " ms\n";

    A2ASpatialSum<SpinColourVector_v> spatial_sum;
    double t_blas = usecond();

    t0 = usecond();
    spatial_sum.Allocate(N_i, N_j, grid);
    std::cout << GridLogMessage << tag << " Allocate:        " << Tms(usecond()-t0) << " ms\n";

    t0 = usecond();
    spatial_sum.PackLeft(leftv);
    std::cout << GridLogMessage << tag << " PackLeft:        " << Tms(usecond()-t0) << " ms\n";

    t0 = usecond();
    spatial_sum.PackRight(loopRight);
    std::cout << GridLogMessage << tag << " PackRight:       " << Tms(usecond()-t0) << " ms\n";

    t0 = usecond();
    spatial_sum.Sum(result);
    std::cout << GridLogMessage << tag << " Sum (GEMM+MPI):  " << Tms(usecond()-t0) << " ms\n";

    std::cout << GridLogMessage << tag << " A2ASpatialSum:   " << Tms(usecond()-t_blas) << " ms  [TOTAL]\n";
  }
};

int main(int argc, char *argv[])
{
  Grid_init(&argc, &argv);

  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(4, vComplexD::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();
  GridCartesian grid(latt_size, simd_layout, mpi_layout);

  int Nt    = latt_size[Tp];
  int N_i   = 8;
  int N_j   = 8;
  int Nloop = 4;

  GridParallelRNG pRNG(&grid);
  pRNG.SeedFixedIntegers({1, 2, 3, 4});

  std::vector<FermionField> left(N_i,   &grid);
  std::vector<FermionField> right(N_j,  &grid);
  std::vector<FermionField> loop1(Nloop, &grid);
  std::vector<FermionField> loop2(Nloop, &grid);

  for (auto &f : left)  random(pRNG, f);
  for (auto &f : right) random(pRNG, f);
  for (auto &f : loop1) random(pRNG, f);
  for (auto &f : loop2) random(pRNG, f);

  std::vector<Gamma::Algebra> GammaMU = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT
  };

  Eigen::Tensor<ComplexD, 3> result_ref(Nt, N_i, N_j);
  Eigen::Tensor<ComplexD, 3> result_blas(Nt, N_i, N_j);
  Eigen::Tensor<ComplexD, 3> result_gpu(Nt, N_i, N_j);
  double t_ref = 0, t_blas = 0, t_gpu = 0, start, stop;

  for (int type = 0; type < 4; type++) {

    result_ref.setZero();
    start = usecond();
    A2AExtendedMesonFieldRef::compute(result_ref,  left, right, loop1, loop2,
                                      GammaMU, GammaMU, type, false);
    stop = usecond(); t_ref = stop - start;

    result_blas.setZero();
    start = usecond();
    A2AExtendedMesonFieldRef::compute(result_blas, left, right, loop1, loop2,
                                      GammaMU, GammaMU, type, true);
    stop = usecond(); t_blas = stop - start;

    result_gpu.setZero();
    start = usecond();
    A2AExtendedMesonFieldGPU::compute(result_gpu, left, right, loop1, loop2,
                                      GammaMU, GammaMU, type);
    stop = usecond(); t_gpu = stop - start;

    double norm2_ref = 0.0, norm2_blas = 0.0, norm2_gpu = 0.0;
    double norm2_diff_blas = 0.0, norm2_diff_gpu = 0.0;
    for (int t  = 0; t  < Nt;  t++)
    for (int ii = 0; ii < N_i; ii++)
    for (int jj = 0; jj < N_j; jj++) {
      norm2_ref  += norm2(result_ref(t, ii, jj));
      norm2_blas += norm2(result_blas(t, ii, jj));
      norm2_gpu  += norm2(result_gpu(t, ii, jj));
      ComplexD diff_blas = result_ref(t, ii, jj) - result_blas(t, ii, jj);
      ComplexD diff_gpu  = result_ref(t, ii, jj) - result_gpu(t, ii, jj);
      norm2_diff_blas += norm2(diff_blas);
      norm2_diff_gpu  += norm2(diff_gpu);
    }

    double rel_blas = (norm2_ref > 0) ? std::sqrt(norm2_diff_blas / norm2_ref) : 0.0;
    double rel_gpu  = (norm2_ref > 0) ? std::sqrt(norm2_diff_gpu  / norm2_ref) : 0.0;

    std::cout << GridLogMessage
              << "type=" << type
              << "  norm2_ref="  << norm2_ref
              << "  norm2_blas=" << norm2_blas
              << "  norm2_gpu="  << norm2_gpu
              << "  rel_blas="   << rel_blas
              << "  rel_gpu="    << rel_gpu
              << "  t_ref="  << t_ref  * 1e-6 << "s"
              << "  t_blas=" << t_blas * 1e-6 << "s"
              << "  t_gpu="  << t_gpu  * 1e-6 << "s"
              << std::endl;

    GRID_ASSERT(rel_blas < 1e-10);
    GRID_ASSERT(rel_gpu  < 1e-10);
  }

  std::cout << GridLogMessage << "All types passed A2ASpatialSum and GPU regression." << std::endl;

  Grid_finalize();
  return EXIT_SUCCESS;
}

#endif
