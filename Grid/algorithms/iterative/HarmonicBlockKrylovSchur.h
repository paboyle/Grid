/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/HarmonicBlockKrylovSchur.h

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#ifndef GRID_HARMONIC_BLOCKED_KRYLOV_SCHUR_H
#define GRID_HARMONIC_BLOCKED_KRYLOV_SCHUR_H

#include <iomanip>

NAMESPACE_BEGIN(Grid);

/**
 * Block shift-targeted Krylov-Schur eigensolver.
 *
 * Algorithm
 * ---------
 * Uses a block Arnoldi factorisation:
 *
 *   A V = V H + F B^dag                                         (1)
 *
 * with V orthonormal (Nm columns), H the Nm×Nm block upper-Hessenberg
 * Rayleigh quotient, F the Nblock residual vectors and B the Nm×Nblock
 * coupling matrix.
 *
 * Thick restart
 * -------------
 * To target eigenvalues near shift σ, the Schur decomposition is computed
 * for the shifted Rayleigh quotient:
 *
 *   (H - σI) = Q^dag S Q                                        (2)
 *
 * Sorting the Schur values of (H - σI) by smallest |S(i,i)| = |λ - σ| and
 * retaining the leading Nk is equivalent to selecting the Ritz values of H
 * closest to σ.  Since Q diagonalises H - σI (and hence H itself), the
 * rotated Rayleigh quotient is exactly upper triangular:
 *
 *   H_new = Q H Q^dag = S + σI  (upper triangular)              (3)
 *
 * Truncation to Nk is therefore exact: the off-diagonal coupling block
 * H_new[Nk:, :Nk] = 0 by triangularity.  The Krylov-Schur relation after
 * restart is exact and block Arnoldi resumes cleanly from F.
 *
 * Convergence
 * -----------
 * Convergence is declared when || B^H y_k || < Tolerance * approxLambdaMax
 * for each Ritz pair (λ_k, y_k) of the truncated H.
 *
 * Parameters
 * ----------
 * shift    : target shift σ (default 0.0); Schur values sorted by |λ - σ|
 * Nblock   : block size p
 * Nm       : total Krylov dimension (must be divisible by Nblock)
 * Nk       : total vectors to retain after each restart (must be divisible by Nblock, Nk < Nm)
 * Nstop    : stop when this many eigenpairs converge
 * MaxIter  : maximum outer (restart) iterations
 * Tolerance: relative convergence tolerance
 *
 * Usage
 * -----
 *   HarmonicBlockKrylovSchur<Field> hbks(LinOp, Grid, tol, shift, EvalNormSmall);
 *   std::vector<Field> v0(Nblock, Field(Grid));
 *   // fill v0 with random starting vectors
 *   hbks(v0, maxIter, Nm, Nk, Nstop, Nblock);
 *   auto evals = hbks.getEvals();
 *   auto evecs = hbks.getEvecs();
 */
template<class Field>
class HarmonicBlockKrylovSchur {

  typedef Eigen::MatrixXcd CMat;
  typedef Eigen::VectorXcd CVec;

  //--------------------------------------------------------------------
  // Parameters
  //--------------------------------------------------------------------
  int    Nblock;
  int    Nm;        // total Krylov dimension (multiple of Nblock)
  int    Nk;        // total vectors retained after restart (multiple of Nblock)
  int    Nstop;
  int    MaxIter;
  RealD  Tolerance;
  ComplexD shift;       // target shift σ

  //--------------------------------------------------------------------
  // Internal state
  //--------------------------------------------------------------------
  LinearOperatorBase<Field>& Linop;
  GridBase*                  Grid_;
  RitzFilter                 ritzFilter;

  std::vector<Field> basis;   // Nm*Nblock flat basis
  CMat               H;       // (Nm*Nblock)² block-Hessenberg Rayleigh quotient
  std::vector<Field> F;       // Nblock residual vectors
  CMat               B;       // (Nm*Nblock) × Nblock coupling matrix
  RealD              beta_k;
  RealD              rtol;

  CVec               evals;
  CMat               littleEvecs;
  std::vector<RealD> ritzEstimates;

public:
  std::vector<Field>  evecs;
  bool doEvalCheck = false;

  //--------------------------------------------------------------------
  // Constructor
  //--------------------------------------------------------------------
  HarmonicBlockKrylovSchur(LinearOperatorBase<Field>& _Linop, GridBase* _Grid,
                             RealD _Tolerance, ComplexD _shift = 0.0,
                             RitzFilter _rf = EvalNormSmall)
    : Linop(_Linop), Grid_(_Grid), Tolerance(_Tolerance), shift(_shift),
      ritzFilter(_rf),
      Nblock(-1), Nm(-1), Nk(-1), Nstop(-1), MaxIter(-1),
      beta_k(0.0), rtol(0.0)
  {}

  //--------------------------------------------------------------------
  // Main entry point
  //--------------------------------------------------------------------
  void operator()(const std::vector<Field>& v0, int _maxIter, int _Nm, int _Nk,
                  int _Nstop, int _Nblock = 1, bool doubleOrthog = true,
                  bool doVerify = false)
  {
    MaxIter = _maxIter;
    Nm      = _Nm;
    Nk      = _Nk;
    Nstop   = _Nstop;
    Nblock  = _Nblock;

    assert((int)v0.size() >= Nblock);
    assert(Nm % Nblock == 0);
    assert(Nk % Nblock == 0);
    assert(Nk < Nm);

    int N = Nm;

    RealD approxLambdaMax = approxMaxEval(v0[0]);
    rtol = Tolerance * approxLambdaMax;
    std::cout << GridLogMessage
              << "HarmonicBlockKrylovSchur: approx max eval = " << approxLambdaMax
              << ", rtol = " << rtol
              << ", shift = " << shift << std::endl;

    H = CMat::Zero(N, N);
    B = CMat::Zero(N, Nblock);

    int start = 0;
    std::vector<Field> startBlock(v0.begin(), v0.begin() + Nblock);

    for (int iter = 0; iter < MaxIter; iter++) {
      std::cout << GridLogMessage
                << "HarmonicBlockKrylovSchur: restart iteration " << iter << std::endl;

      // ---- Block Arnoldi: extend from block 'start' to block Nm/Nblock ----
      blockArnoldiIteration(startBlock, Nm/Nblock, start, doubleOrthog);
      start = Nk/Nblock;

      if (doVerify) {
        std::string lbl = "iter " + std::to_string(iter) + " after Arnoldi";
        verify(lbl);
      }

      // ---- Schur decompose (H - σI) to select Schur vectors closest to σ ----
      // Sorting the Schur values of (H - σI) by |S(i,i)| = |λ - σ| gives the
      // Ritz values of H nearest the target shift without any matrix inversion.
      // Because Q diagonalises (H - σI), it also diagonalises H:
      //   H_new = Q H Q^dag = S + σI  (upper triangular)
      // Truncation is therefore exact.
      CMat Hshift = H - shift * CMat::Identity(N, N);
      ComplexSchurDecomposition schur(Hshift, false, ritzFilter);
      schur.schurReorder(Nk);

      std::cout << GridLogMessage
                << "HarmonicBlockKrylovSchur: Ritz values nearest shift (first Nk):" << std::endl;
      CMat S = schur.getMatrixS();
      for (int i = 0; i < Nk; i++)
        std::cout << GridLogMessage << "  [" << i << "] " << S(i, i) + shift << std::endl;

      CMat Q  = schur.getMatrixQ();
      CMat Qt = Q.adjoint();

      // ---- Rotate Krylov basis ----
      std::vector<Field> basis2;
      constructUR(basis2, basis, Qt, N);
      basis = basis2;

      // ---- Update H and B ----
      // H_new = S + σI is upper triangular; off-diagonal block H_new[Nk:,:Nk] = 0
      H = S + shift * CMat::Identity(N, N);
      B = Q * B;

      // ---- Truncate to Nk (exact: H upper triangular) ----
      int Nkeep = Nk;

      CMat Htmp = H(Eigen::seqN(0, Nkeep), Eigen::seqN(0, Nkeep));
      H = CMat::Zero(N, N);
      H(Eigen::seqN(0, Nkeep), Eigen::seqN(0, Nkeep)) = Htmp;

      std::vector<Field> basisTmp(basis.begin(), basis.begin() + Nkeep);
      basis = basisTmp;

      CMat Btmp = B(Eigen::seqN(0, Nkeep), Eigen::all);
      B = CMat::Zero(N, Nblock);
      B(Eigen::seqN(0, Nkeep), Eigen::all) = Btmp;

      beta_k = Btmp.norm();
      std::cout << GridLogMessage
                << "HarmonicBlockKrylovSchur: beta_k = " << beta_k << std::endl;

      // Restart from F (exact: no discarded-basis correction needed)
      startBlock = F;

      if (doVerify) {
        std::string lbl = "iter " + std::to_string(iter) + " after restart+truncation";
        verify(lbl);
      }

      // ---- Eigensystem of truncated H for convergence ----
      CMat Hk = H(Eigen::seqN(0, Nkeep), Eigen::seqN(0, Nkeep));
      computeEigensystem(Hk, Nkeep);

      int Nconv = converged(Nkeep);
      std::cout << GridLogMessage
                << "HarmonicBlockKrylovSchur: converged " << Nconv
                << " / " << Nstop << std::endl;

      if (Nconv >= Nstop || iter == MaxIter - 1) {
        std::cout << GridLogMessage
                  << "HarmonicBlockKrylovSchur: done after " << iter
                  << " restarts, " << Nconv << " converged." << std::endl;
        std::cout << GridLogMessage << "Eigenvalues: " << evals.transpose() << std::endl;

        if (doEvalCheck) {
          Field w(Grid_);
          for (int k = 0; k < (int)evecs.size(); k++) {
            Linop.Op(evecs[k], w);
            ComplexD eval_est = toStdCmplx(innerProduct(evecs[k], w));
            w -= eval_est * evecs[k];
            RealD res = std::sqrt(norm2(w));
            std::cout << GridLogMessage << "HarmonicBlockKrylovSchur: evec[" << k << "]"
                      << "  eval_reported = " << evals[k]
                      << "  eval_est = " << eval_est
                      << "  || A v - eval_est * v || = " << res << std::endl;
          }
        }

        return;
      }
    }
  }

  // Accessors
  std::vector<Field>  getEvecs()         { return evecs;         }
  CVec                getEvals()         { return evals;         }
  std::vector<RealD>  getRitzEstimates() { return ritzEstimates; }

  //--------------------------------------------------------------------
  // Verification: check A V = V H + F B^dag explicitly
  //--------------------------------------------------------------------
  /**
   * Checks the block Arnoldi / Krylov-Schur decomposition
   *
   *   A V = V H + F B^dag                                      (KS)
   *
   * by explicit operator applications.  H here is the standard Rayleigh
   * quotient, so the KS relation is the same as for BlockKrylovSchur.
   *
   * Prints:
   *   - H (current Rayleigh quotient, nBasis × nBasis)
   *   - B (coupling matrix, nBasis × Nblock)
   *   - M (explicit inner product matrix <V | A V>)
   *   - max |H[i,j] - M[i,j]|  (should be O(machine epsilon))
   *   - max |<V_i|V_j> - delta_ij|  (orthonormality check)
   *   - for each basis column j: || A v_j - V H[:,j] - F B[j,:]^* ||
   *   - max |<V_i | F_t>|  (F orthogonal to basis)
   */
  void verify(const std::string& label = "")
  {
    int nBasis = (int)basis.size();
    int nF     = (int)F.size();

    if (nBasis == 0) {
      std::cout << GridLogMessage
                << "HarmonicBlockKrylovSchur::verify [" << label
                << "]: basis is empty." << std::endl;
      return;
    }

    std::cout << GridLogMessage
              << "======== HarmonicBlockKrylovSchur::verify [" << label << "] ========" << std::endl;
    std::cout << GridLogMessage
              << "  nBasis = " << nBasis << "  Nblock = " << Nblock
              << "  nF = " << nF << std::endl;

    // ---- Print H ----
    std::cout << GridLogMessage << "H (" << nBasis << " x " << nBasis << "):" << std::endl;
    for (int i = 0; i < nBasis; i++) {
      for (int j = 0; j < nBasis; j++)
        std::cout << "  " << std::setprecision(4) << std::setw(14) << H(i, j);
      std::cout << std::endl;
    }

    // ---- Print B ----
    std::cout << GridLogMessage << "B (" << nBasis << " x " << nF << "):" << std::endl;
    for (int i = 0; i < nBasis; i++) {
      for (int t = 0; t < nF; t++)
        std::cout << "  " << std::setprecision(4) << std::setw(14) << B(i, t);
      std::cout << std::endl;
    }

    // ---- Compute M[i,j] = <basis[i] | A basis[j]> ----
    CMat M = CMat::Zero(nBasis, nBasis);
    Field w(Grid_);
    for (int j = 0; j < nBasis; j++) {
      Linop.Op(basis[j], w);
      for (int i = 0; i < nBasis; i++)
        M(i, j) = toStdCmplx(innerProduct(basis[i], w));
    }

    std::cout << GridLogMessage << "M = <V|AV> (" << nBasis << " x " << nBasis << "):" << std::endl;
    for (int i = 0; i < nBasis; i++) {
      for (int j = 0; j < nBasis; j++)
        std::cout << "  " << std::setprecision(4) << std::setw(14) << M(i, j);
      std::cout << std::endl;
    }

    // ---- max |H - M| ----
    RealD maxHM = 0.0;
    for (int i = 0; i < nBasis; i++)
      for (int j = 0; j < nBasis; j++)
        maxHM = std::max(maxHM, std::abs(H(i,j) - M(i,j)));
    std::cout << GridLogMessage
              << "  max |H[i,j] - M[i,j]| = " << maxHM << std::endl;

    // ---- Check orthonormality of basis ----
    CMat G = CMat::Zero(nBasis, nBasis);
    for (int i = 0; i < nBasis; i++)
      for (int j = 0; j < nBasis; j++)
        G(i, j) = toStdCmplx(innerProduct(basis[i], basis[j]));
    CMat Gerr = G - CMat::Identity(nBasis, nBasis);
    std::cout << GridLogMessage
              << "  max |<V_i|V_j> - delta_ij| = " << Gerr.cwiseAbs().maxCoeff() << std::endl;

    // ---- Per-column residual: || A v_j - V H[:,j] - F B[j,:]^* || ----
    RealD maxColDev = 0.0;
    for (int j = 0; j < nBasis; j++) {
      Linop.Op(basis[j], w);

      // subtract V H[:,j]
      for (int i = 0; i < nBasis; i++)
        w -= basis[i] * H(i, j);

      // subtract F B[j,:]^*  (F[t] * conj(B[j,t]))
      for (int t = 0; t < nF; t++)
        w -= F[t] * std::conj(B(j, t));

      RealD dev = std::sqrt(norm2(w));
      std::cout << GridLogMessage
                << "  || A v[" << j << "] - V H[:,j] - F B[j,:]* || = " << dev << std::endl;
      maxColDev = std::max(maxColDev, dev);
    }
    std::cout << GridLogMessage
              << "  max column deviation = " << maxColDev << std::endl;

    // ---- Check F block orthogonality against basis ----
    if (nF > 0) {
      RealD maxFV = 0.0;
      for (int t = 0; t < nF; t++)
        for (int i = 0; i < nBasis; i++) {
          RealD ip = std::abs(toStdCmplx(innerProduct(basis[i], F[t])));
          maxFV = std::max(maxFV, ip);
        }
      std::cout << GridLogMessage
                << "  max |<V_i | F_t>| (should be ~0) = " << maxFV << std::endl;
    }

    std::cout << GridLogMessage
              << "======== end verify ========" << std::endl;
  }

private:

  //--------------------------------------------------------------------
  // Block Arnoldi iteration
  //--------------------------------------------------------------------
  void blockArnoldiIteration(std::vector<Field>& startBlock, int endBlock,
                             int startIdx, bool doubleOrthog)
  {
    int N = Nm;

    if (startIdx == 0) {
      basis.clear();
      F.clear();
      H = CMat::Zero(N, N);
      B = CMat::Zero(N, Nblock);

      std::vector<Field> V0 = startBlock;
      blockOrthonormalise(V0);
      for (auto& v : V0) basis.push_back(v);
    } else {
      // Append the new starting block to the retained basis.
      //
      // The exact truncated KS relation is  A V_k = V_k H_k + F B_k^dag,
      // so the coupling rows are  H[Nkeep+t, j] = conj(B_k[j,t]).
      // Since H_new = S + σI is upper triangular, the off-diagonal block
      // H_new[Nkeep:, :Nkeep] = 0 and the restart from F is exact.
      int Nkeep = startIdx * Nblock;
      for (auto& v : startBlock) basis.push_back(v);

      // Fill restart coupling rows into H (exact: H_new is upper triangular,
      // so B encodes the only non-zero coupling to the new block).
      for (int t = 0; t < Nblock; t++)
        for (int j = 0; j < Nkeep; j++)
          H(Nkeep + t, j) = std::conj(B(j, t));

      // Zero out B for the retained columns now that the coupling is in H
      for (int j = 0; j < Nkeep; j++)
        for (int t = 0; t < Nblock; t++)
          B(j, t) = 0.0;
    }

    for (int k = startIdx; k < endBlock; k++)
      blockArnoldiStep(k, doubleOrthog);
  }

  //--------------------------------------------------------------------
  // One block Arnoldi step
  //--------------------------------------------------------------------
  void blockArnoldiStep(int k, bool doubleOrthog)
  {
    int kBase  = k * Nblock;
    int prevN  = kBase + Nblock;
    int N      = Nm;

    std::vector<Field> W(Nblock, Field(Grid_));
    for (int t = 0; t < Nblock; t++)
      Linop.Op(basis[kBase + t], W[t]);

    // Full reorthogonalisation against all current basis vectors
    for (int pass = 0; pass < (doubleOrthog ? 2 : 1); pass++) {
      for (int i = 0; i < prevN; i++) {
        for (int t = 0; t < Nblock; t++) {
          ComplexD coeff = innerProduct(basis[i], W[t]);
          if (pass == 0)
            H(i, kBase + t) = toStdCmplx(coeff);
          else
            H(i, kBase + t) += toStdCmplx(coeff);
          W[t] -= coeff * basis[i];
        }
      }
    }

    F = W;

    if (k == Nm/Nblock - 1) {
      // Last block: record coupling in B as R^H (Hermitian conjugate of QR factor)
      // KS relation requires B[kBase+t, s] = conj(R[s,t])
      CMat R = blockQR(F);
      for (int t = 0; t < Nblock; t++)
        for (int s = 0; s < Nblock; s++)
          B(kBase + t, s) = std::conj(R(s, t));    // B_block = R^H
      beta_k = R.norm();
      return;
    }

    // Not last: QR the residual, extend basis
    CMat R = blockQR(F);

    int nextBase = (k + 1) * Nblock;
    for (int i = 0; i < Nblock; i++)
      for (int j = 0; j < Nblock; j++)
        H(nextBase + i, kBase + j) = R(i, j);

    for (int t = 0; t < Nblock; t++)
      basis.push_back(F[t]);

  }

  //--------------------------------------------------------------------
  // Block QR (modified Gram-Schmidt within the block)
  //--------------------------------------------------------------------
  CMat blockQR(std::vector<Field>& W)
  {
    CMat R = CMat::Zero(Nblock, Nblock);
    const RealD deflThresh = 1e-14;

    for (int j = 0; j < Nblock; j++) {
      for (int i = 0; i < j; i++) {
        ComplexD coeff = innerProduct(W[i], W[j]);
        R(i, j) = toStdCmplx(coeff);
        W[j] -= coeff * W[i];
      }
      RealD nrm = std::sqrt(norm2(W[j]));
      R(j, j)   = nrm;
      if (nrm > deflThresh) {
        W[j] *= (1.0 / nrm);
      } else {
        W[j] = Zero();
        std::cout << GridLogMessage
                  << "HarmonicBlockKrylovSchur: deflation at block column " << j
                  << " (norm = " << nrm << ")" << std::endl;
      }
    }
    return R;
  }

  //--------------------------------------------------------------------
  // Orthonormalise starting block
  //--------------------------------------------------------------------
  void blockOrthonormalise(std::vector<Field>& V)
  {
    for (int j = 0; j < (int)V.size(); j++) {
      for (int i = 0; i < j; i++) {
        ComplexD c = innerProduct(V[i], V[j]);
        V[j] -= c * V[i];
      }
      RealD nrm = std::sqrt(norm2(V[j]));
      assert(nrm > 1e-14);
      V[j] *= (1.0 / nrm);
    }
  }

  //--------------------------------------------------------------------
  // Basis rotation: UR[i] = sum_j U[j] * R(j,i)
  //--------------------------------------------------------------------
  void constructUR(std::vector<Field>& UR, std::vector<Field>& U,
                   CMat& R, int N)
  {
    UR.clear();
    Field tmp(Grid_);
    for (int i = 0; i < N; i++) {
      tmp = Zero();
      for (int j = 0; j < N; j++)
        tmp += U[j] * R(j, i);
      UR.push_back(tmp);
    }
  }

  //--------------------------------------------------------------------
  // Eigensystem of the truncated H (not Hhat)
  //--------------------------------------------------------------------
  /**
   * Eigenvalues of H_k are the standard Ritz values in the retained
   * subspace.  After convergence has been declared via harmonic estimates,
   * the final reported eigenvalues and vectors come from H_k (not Hhat_k),
   * since H_k contains the true projected operator.
   */
  void computeEigensystem(CMat& Hk, int Nkeep)
  {
    Eigen::ComplexEigenSolver<CMat> es;
    es.compute(Hk);
    evals       = es.eigenvalues();
    littleEvecs = es.eigenvectors();

    evecs.clear();
    for (int k = 0; k < Nkeep; k++) {
      CVec vec = littleEvecs.col(k);
      Field tmp(Grid_);
      tmp = Zero();
      for (int j = 0; j < (int)basis.size() && j < Nkeep; j++)
        tmp += vec[j] * basis[j];
      evecs.push_back(tmp);
    }
  }

  //--------------------------------------------------------------------
  // Convergence check
  //--------------------------------------------------------------------
  /**
   * Ritz estimate for eigenpair k: || B^H y_k ||
   * where y_k is the k-th eigenvector of the truncated H.
   * The same bound applies whether using Ritz or harmonic Ritz restart.
   */
  int converged(int Nkeep)
  {
    ritzEstimates.clear();
    int Nconv = 0;

    CMat Bk = B(Eigen::seqN(0, Nkeep), Eigen::all);

    for (int k = 0; k < Nkeep; k++) {
      CVec yk   = littleEvecs.col(k);
      CVec Bty  = Bk.adjoint() * yk;
      RealD res = Bty.norm();
      ritzEstimates.push_back(res);
      std::cout << GridLogMessage
                << "HarmonicBlockKrylovSchur: Ritz estimate[" << k
                << "] = " << res << "  eval = " << evals[k] << std::endl;
      if (res < rtol) Nconv++;
    }
    return Nconv;
  }

  //--------------------------------------------------------------------
  // Approximate maximum eigenvalue (power iteration)
  //--------------------------------------------------------------------
  RealD approxMaxEval(const Field& v0, int MAX_ITER = 50)
  {
    assert(norm2(v0) > 1e-8);
    RealD lam = 0.0, denom = std::sqrt(norm2(v0));
    Field vcur(Grid_), vtmp(Grid_);
    vcur = v0;
    for (int i = 0; i < MAX_ITER; i++) {
      Linop.Op(vcur, vtmp);
      vcur  = vtmp;
      RealD num = std::sqrt(norm2(vcur));
      lam   = num / denom;
      denom = num;
    }
    return lam;
  }

};

NAMESPACE_END(Grid);

#endif  // GRID_HARMONIC_BLOCKED_KRYLOV_SCHUR_H
