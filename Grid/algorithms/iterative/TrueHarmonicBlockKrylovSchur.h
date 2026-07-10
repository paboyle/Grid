/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/TrueHarmonicBlockKrylovSchur.h

Copyright (C) 2015

Author: Chulwoo Jung <chulwoo@bnl.gov>

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
#ifndef GRID_TRUE_HARMONIC_BLOCKED_KRYLOV_SCHUR_H
#define GRID_TRUE_HARMONIC_BLOCKED_KRYLOV_SCHUR_H

#include <Grid/algorithms/iterative/HarmonicBlockKrylovSchur.h>

NAMESPACE_BEGIN(Grid);

/**
 * True harmonic-Ritz block Krylov-Schur eigensolver.
 *
 * Difference from the base class
 * ------------------------------
 * HarmonicBlockKrylovSchur is a *shift-targeted* Krylov-Schur: it Schur-
 * decomposes (H - sigma I) and sorts the standard Ritz values by |lambda -
 * sigma|.  For interior eigenvalues of Hermitian indefinite operators the
 * standard Ritz values can be spurious ("ghosts"), so the sort can select
 * junk.  This variant performs genuine harmonic Ritz extraction instead.
 *
 * Harmonic extraction
 * -------------------
 * From the exact block Arnoldi relation
 *
 *   A V = V H + F B^dag                                            (1)
 *
 * harmonic Ritz pairs (theta, g) w.r.t. shift sigma satisfy the Petrov-
 * Galerkin condition (A - sigma)Vg - (theta - sigma)Vg  perp  (A - sigma)V,
 * which reduces (with Hs = H - sigma I) to the small eigenproblem
 *
 *   Hhat g = (theta - sigma) g,   Hhat = Hs + Hs^{-H} B B^dag      (2)
 *
 * Exact thick restart
 * -------------------
 * From (2), H g = theta g - Hs^{-H} B (B^dag g), so the residual of every
 * harmonic pair lies in the Nblock-dimensional column space of the single
 * block
 *
 *   W = F - V Hs^{-H} B :
 *   A (V g) - theta (V g) = W (B^dag g)                            (3)
 *
 * Hence {V g_1 .. g_Nk} + cols(W) admits an EXACT truncated Krylov-Schur
 * relation  A Y = Y H_k + What Bc^dag  with dense H_k, and block Arnoldi
 * resumes cleanly from What.  Krylov-Schur (Stewart, SIMAX 23(3), 2001)
 * does not require H to be triangular after restart -- only that the
 * decomposition holds exactly, which (3) guarantees.
 *
 * This is the block analogue of the GMRES-DR restart (Morgan, SISC 24(1),
 * 2002; block version: Morgan, Appl. Numer. Math. 54, 2005).  The exact
 * composition "block harmonic Krylov-Schur" is a derivation from these
 * ingredients rather than a published algorithm; the relation (1) is
 * checked explicitly at every restart when doVerify=true, and identity (3)
 * was validated to machine precision in a standalone Eigen prototype.
 *
 * Caveats
 * -------
 * - If sigma is (numerically) an eigenvalue of H, Hs^{-H} is singular;
 *   an assert fires -- perturb the shift slightly.
 * - Reported eigenvalues come from the eigensystem of the retained H_k
 *   (the true projected Rayleigh quotient Y^dag A Y), as in the base
 *   class, so convergence tests and reporting are inherited unchanged.
 *
 * Usage: identical to HarmonicBlockKrylovSchur.
 */
template<class Field>
class TrueHarmonicBlockKrylovSchur : public HarmonicBlockKrylovSchur<Field> {

  typedef HarmonicBlockKrylovSchur<Field> Base;
  typedef typename Base::CMat CMat;
  typedef typename Base::CVec CVec;

  // Members of the dependent base class need explicit scope.
  using Base::Nblock;
  using Base::Nm;
  using Base::Nk;
  using Base::Nstop;
  using Base::MaxIter;
  using Base::Tolerance;
  using Base::shift;
  using Base::Linop;
  using Base::Grid_;
  using Base::ritzFilter;
  using Base::basis;
  using Base::H;
  using Base::F;
  using Base::B;
  using Base::beta_k;
  using Base::rtol;
  using Base::evals;
  using Base::littleEvecs;
  using Base::ritzEstimates;

  using Base::expandStartBlock;
  using Base::blockArnoldiIteration;
  using Base::blockQR;
  using Base::computeEigensystem;
  using Base::converged;
  using Base::approxMaxEval;

public:
  using Base::evecs;
  using Base::doEvalCheck;
  using Base::useParityFlip;
  using Base::useGamma5;
  using Base::gamma5Func;
  using Base::verify;

  TrueHarmonicBlockKrylovSchur(LinearOperatorBase<Field>& _Linop, GridBase* _Grid,
                               RealD _Tolerance, ComplexD _shift = 0.0,
                               RitzFilter _rf = EvalNormSmall)
    : Base(_Linop, _Grid, _Tolerance, _shift, _rf)
  {}

  //--------------------------------------------------------------------
  // Main entry point (same signature and outer structure as the base;
  // the Schur-sort restart is replaced by harmonicRestart()).
  //--------------------------------------------------------------------
  virtual void operator()(const std::vector<Field>& v0, int _maxIter, int _Nm, int _Nk,
                          int _Nstop, int _Nblock = 1, bool doubleOrthog = true,
                          bool doVerify = false) override
  {
    MaxIter = _maxIter;
    Nm      = _Nm;
    Nk      = _Nk;
    Nstop   = _Nstop;
    Nblock  = _Nblock;

    {
      int divisor = 1;
      if (useParityFlip) divisor *= 2;
      if (useGamma5)     divisor *= 2;
      assert(Nblock % divisor == 0 && (int)v0.size() >= Nblock / divisor);
    }
    assert(Nm % Nblock == 0);
    assert(Nk % Nblock == 0);
    assert(Nk < Nm);
    if (useGamma5) assert(gamma5Func && "useGamma5: gamma5Func must be set");

    int N = Nm;

    RealD approxLambdaMax = approxMaxEval(v0[0]);
    rtol = Tolerance * approxLambdaMax;
    std::cout << GridLogMessage
              << "TrueHarmonicBlockKrylovSchur: approx max eval = " << approxLambdaMax
              << ", rtol = " << rtol
              << ", shift = " << shift << std::endl;

    H = CMat::Zero(N, N);
    B = CMat::Zero(N, Nblock);

    int start = 0;
    std::vector<Field> startBlock = expandStartBlock(v0);

    for (int iter = 0; iter < MaxIter; iter++) {
      std::cout << GridLogMessage
                << "TrueHarmonicBlockKrylovSchur: restart iteration " << iter << std::endl;

      // ---- Block Arnoldi: extend from block 'start' to block Nm/Nblock ----
      blockArnoldiIteration(startBlock, Nm/Nblock, start, doubleOrthog);
      std::cout << GridLogMessage << "blockArnoldiIteration done " << std::endl;
      start = Nk/Nblock;

      if (doVerify) {
        std::string lbl = "iter " + std::to_string(iter) + " after Arnoldi";
        verify(lbl);
      }

      // ---- Harmonic extraction + exact thick restart ----
      harmonicRestart();
      std::cout << GridLogMessage << "harmonicRestart done " << std::endl;

      // Restart from the residual block W (exact by identity (3))
      startBlock = F;

      if (doVerify) {
        std::string lbl = "iter " + std::to_string(iter) + " after harmonic restart";
        verify(lbl);
      }

      // ---- Eigensystem of retained H_k = Y^dag A Y for convergence ----
      CMat Hk = H(Eigen::seqN(0, Nk), Eigen::seqN(0, Nk));
      computeEigensystem(Hk, Nk);

      int Nconv = converged(Nk);
      std::cout << GridLogMessage
                << "TrueHarmonicBlockKrylovSchur: converged " << Nconv
                << " / " << Nstop << std::endl;

      if (Nconv >= Nstop || iter == MaxIter - 1) {
        std::cout << GridLogMessage
                  << "TrueHarmonicBlockKrylovSchur: done after " << iter
                  << " restarts, " << Nconv << " converged." << std::endl;
        std::cout << GridLogMessage << "Eigenvalues: " << evals.transpose() << std::endl;

        if (doEvalCheck) {
          Field w(Grid_);
          for (int k = 0; k < (int)evecs.size(); k++) {
            Linop.Op(evecs[k], w);
            ComplexD eval_est = toStdCmplx(innerProduct(evecs[k], w));
            w -= eval_est * evecs[k];
            RealD res = std::sqrt(norm2(w));
            std::cout << GridLogMessage << "TrueHarmonicBlockKrylovSchur: evec[" << k << "]"
                      << "  eval_reported = " << evals[k]
                      << "  eval_est = " << eval_est
                      << "  || A v - eval_est * v || = " << res << std::endl;
          }
        }

        return;
      }
    }
  }

protected:

  //--------------------------------------------------------------------
  // Harmonic thick restart
  //--------------------------------------------------------------------
  // On entry:  A V = V H + F B^dag exact, V = basis (Nm fields), F Nblock
  //            residual fields, H/B the Nm-sized small matrices.
  // On exit:   basis = Y (Nk orthonormal fields spanning the harmonic
  //            Ritz space), F = What (Nblock orthonormal fields, perp Y),
  //            H top-left Nk x Nk = Y^dag A Y (dense), B top Nk rows the
  //            new coupling; relation A Y = Y H_k + F B^dag again exact.
  void harmonicRestart()
  {
    int N = Nm;

    // ---- Small harmonic eigenproblem:  Hhat g = (theta - sigma) g ----
    CMat Hs = H - shift * CMat::Identity(N, N);
    Eigen::FullPivLU<CMat> lu(Hs.adjoint());
    assert(lu.isInvertible() &&
           "TrueHarmonicBlockKrylovSchur: H - shift*I singular; perturb the shift");
    CMat M    = lu.solve(B);              // N x Nblock :  Hs^{-H} B
    CMat Hhat = Hs + M * B.adjoint();     // N x N

    Eigen::ComplexEigenSolver<CMat> es(Hhat);

    // Sort harmonic values (theta - sigma) by the Ritz filter
    // (EvalNormSmall -> closest to the shift).
    ComplexComparator cComp(ritzFilter);
    std::vector<int> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b){
      return cComp(toStdCmplx(es.eigenvalues()(a)), toStdCmplx(es.eigenvalues()(b)));
    });

    CMat G(N, Nk);
    CVec theta(Nk);
    for (int k = 0; k < Nk; k++) {
      G.col(k) = es.eigenvectors().col(idx[k]);
      theta(k) = es.eigenvalues()(idx[k]) + toStdCmplx(shift);
    }

    std::cout << GridLogMessage
              << "TrueHarmonicBlockKrylovSchur: harmonic Ritz values nearest shift:" << std::endl;
    for (int k = 0; k < Nk; k++)
      std::cout << GridLogMessage << "  [" << k << "] " << theta(k) << std::endl;

    // ---- Orthonormalise the little vectors:  G = Gt Rg ----
    Eigen::HouseholderQR<CMat> qr(G);
    CMat Gt = qr.householderQ() * CMat::Identity(N, Nk);
    CMat Rg = qr.matrixQR().topLeftCorner(Nk, Nk).template triangularView<Eigen::Upper>();
    for (int k = 0; k < Nk; k++)
      assert(std::abs(Rg(k, k)) > 1e-14 &&
             "TrueHarmonicBlockKrylovSchur: harmonic vectors linearly dependent");
    CMat RgInv = Rg.inverse();

    // A (V G) = (V G) diag(theta) + W (B^dag G),  W = F - V M
    //   =>  A (V Gt) = (V Gt) T + W Ct
    CMat T  = Rg * theta.asDiagonal() * RgInv;   // Nk x Nk
    CMat Ct = (B.adjoint() * G) * RgInv;         // Nblock x Nk

    // ---- Large-field work ----
    // Y = V Gt : new orthonormal basis (Nk fields)
    std::vector<Field> Y;
    constructThin(Y, basis, Gt, Nk);

    // W = F - V M : Nblock residual directions containing every harmonic
    // residual (identity (3) in the class comment)
    std::vector<Field> W = F;
    for (int t = 0; t < Nblock; t++)
      for (int j = 0; j < N; j++)
        W[t] -= basis[j] * M(j, t);

    // Orthogonalise W against Y (two passes), recording P = Y^dag W
    CMat P = CMat::Zero(Nk, Nblock);
    for (int pass = 0; pass < 2; pass++) {
      for (int j = 0; j < Nk; j++) {
        for (int t = 0; t < Nblock; t++) {
          ComplexD c = innerProduct(Y[j], W[t]);
          P(j, t) += toStdCmplx(c);
          W[t] -= c * Y[j];
        }
      }
    }

    // Block QR of the projected residual: W <- What, returns Rw
    CMat Rw = blockQR(W);

    // ---- Reassemble the exact Krylov-Schur relation ----
    //   A Y = Y (T + P Ct) + What (Rw Ct)
    CMat Hk = T + P * Ct;                        // Nk x Nk, dense
    CMat Bc = Rw * Ct;                           // Nblock x Nk

    H = CMat::Zero(N, N);
    H(Eigen::seqN(0, Nk), Eigen::seqN(0, Nk)) = Hk;
    B = CMat::Zero(N, Nblock);
    for (int j = 0; j < Nk; j++)
      for (int t = 0; t < Nblock; t++)
        B(j, t) = std::conj(Bc(t, j));

    basis = Y;
    F     = W;
    beta_k = Bc.norm();

    std::cout << GridLogMessage
              << "TrueHarmonicBlockKrylovSchur: beta_k = " << beta_k << std::endl;
  }

  //--------------------------------------------------------------------
  // Thin basis combination: out[i] = sum_j U[j] * R(j,i), i < ncol
  //--------------------------------------------------------------------
  void constructThin(std::vector<Field>& out, const std::vector<Field>& U,
                     const CMat& R, int ncol)
  {
    out.clear();
    Field tmp(Grid_);
    for (int i = 0; i < ncol; i++) {
      tmp = Zero();
      for (int j = 0; j < (int)U.size(); j++)
        tmp += U[j] * R(j, i);
      out.push_back(tmp);
    }
  }

};

NAMESPACE_END(Grid);

#endif  // GRID_TRUE_HARMONIC_BLOCKED_KRYLOV_SCHUR_H
