/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./Grid/algorithms/iterative/RestartedLanczosBidiagonalization.h

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
#ifndef GRID_RESTARTED_LANCZOS_BIDIAGONALIZATION_H
#define GRID_RESTARTED_LANCZOS_BIDIAGONALIZATION_H

NAMESPACE_BEGIN(Grid);

/**
 * Implicitly Restarted Lanczos Bidiagonalization (IRLBA)
 *
 * Computes the p largest (or p smallest) singular triplets of a linear
 * operator A using the Golub-Kahan-Lanczos bidiagonalization with implicit
 * restart via thick-restart / QR shifts.
 *
 * Algorithm (Baglama & Reichel, SIAM J. Sci. Comput. 27(1):19-42, 2005):
 *
 *   Outer loop:
 *     1. Extend the p-step (or seed) bidiagonalization to k steps:
 *           A  V_k = U_k B_k
 *           A^dag U_k = V_k B_k^T + beta_{k+1} v_{k+1} e_k^T
 *     2. Compute SVD:  B_k = X Sigma Y^T
 *     3. Check convergence of the p desired singular values via
 *           |beta_{k+1} * y_{k,i}|  <  tol * sigma_i
 *        where y_{k,i} is the last component of the i-th right singular vector.
 *     4. Apply k-p implicit QR shifts to implicitly compress the basis
 *        to p steps (Sorensen-Lehoucq thick restart):
 *           B_p^+ = X_p^T B_k Y_p   (upper bidiagonal, p x p)
 *        and update the lattice vectors:
 *           V_p^+ = V_k Y_p
 *           U_p^+ = U_k X_p
 *        The new residual coupling is
 *           beta_p^+ v_{p+1}^+ = beta_{k+1} v_{k+1} * (e_k^T Y_p)_p
 *                               + B_k(p,p+1) * (orthogonal tail from QR)
 *     5. Go to step 1.
 *
 * Template parameter
 * ------------------
 *   Field : lattice field type (must support Grid algebra operations)
 *
 * Usage
 * -----
 *   RestartedLanczosBidiagonalization<Field> irlba(Linop, grid, p, k, tol, maxIter);
 *   irlba.run(src);
 *   // Results available via getters.
 */
template <class Field>
class RestartedLanczosBidiagonalization {

public:
  LinearOperatorBase<Field> &Linop;
  GridBase *Grid;

  int    Nk;       // number of desired singular triplets
  int    Nm;       // Lanczos basis size (Nm > Nk)
  RealD  Tolerance;
  int    MaxIter;
  bool   largest; // if true, target largest singular values; otherwise smallest

  // Converged singular triplets (filled after run())
  std::vector<RealD>  singularValues;   // sigma_0 >= sigma_1 >= ...
  std::vector<Field>  leftVectors;      // approximate left singular vectors
  std::vector<Field>  rightVectors;     // approximate right singular vectors

private:
  // Working bases (size up to Nm+1)
  std::vector<Field>  V;    // right Lanczos vectors
  std::vector<Field>  U;    // left  Lanczos vectors
  std::vector<RealD>  alpha;
  std::vector<RealD>  beta;

  // After a thick restart, the column at index restart_col of U^dag A V
  // has extra non-zero entries (rows 0..restart_col-2) beyond what the
  // upper bidiagonal captures.  fvec[j] = <U[j] | A V[restart_col]> for
  // j = 0..restart_col-1.  (fvec[restart_col-1] == beta[restart_col-1].)
  // reset_col == -1 means no restart has occurred yet (pure bidiagonal).
  std::vector<RealD>  fvec;
  int                 restart_col;

public:

  RestartedLanczosBidiagonalization(LinearOperatorBase<Field> &_Linop,
                                    GridBase *_Grid,
                                    int _Nk, int _Nm,
                                    RealD _tol   = 1.0e-8,
                                    int   _maxIt = 300,
                                    bool  _largest = true)
    : Linop(_Linop), Grid(_Grid),
      Nk(_Nk), Nm(_Nm),
      Tolerance(_tol), MaxIter(_maxIt),
      largest(_largest)
  {
    assert(Nm > Nk);
  }

  /**
   * Run IRLBA starting from src.
   * On exit, singularValues, leftVectors, rightVectors are filled with
   * the Nk converged singular triplets.
   */
  void run(const Field &src)
  {
    assert(norm2(src) > 0.0);

    singularValues.clear();
    leftVectors.clear();
    rightVectors.clear();

    // Allocate working bases
    V.clear(); U.clear();
    alpha.clear(); beta.clear();
    fvec.clear(); restart_col = -1;
    V.reserve(Nm + 1);
    U.reserve(Nm);

    // Seed: v_0 = src / ||src||
    Field vtmp(Grid);
    vtmp = src;
    RealD nrm = std::sqrt(norm2(vtmp));
    vtmp = (1.0 / nrm) * vtmp;
    V.push_back(vtmp);

    int pStart = 0;  // current basis size at start of extension
    RealD betaRestart = 0.0; // coupling from previous restart

    for (int iter = 0; iter < MaxIter; ++iter) {

      // ----------------------------------------------------------------
      // Step 1: extend from pStart steps to Nm steps
      // ----------------------------------------------------------------
      extendBasis(pStart, Nm, betaRestart);
      verify();

      // ----------------------------------------------------------------
      // Step 2: SVD of the Nm x Nm matrix B (non-bidiagonal after restart)
      // ----------------------------------------------------------------
      Eigen::MatrixXd B = buildFullB(Nm);
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(B,
          Eigen::ComputeThinU | Eigen::ComputeThinV);

      Eigen::VectorXd sigma = svd.singularValues();  // descending
      Eigen::MatrixXd X     = svd.matrixU();          // Nm x Nm left SVecs of B
      Eigen::MatrixXd Y     = svd.matrixV();          // Nm x Nm right SVecs of B

      // If targeting smallest, reorder so desired ones come first
      Eigen::VectorXi order = sortOrder(sigma);

      // ----------------------------------------------------------------
      // Step 3: check convergence of the Nk desired singular values
      // ----------------------------------------------------------------
      RealD betaK = beta.back();  // beta_{k+1}
      // In our convention A V = U B (exact), the residual is in the A^dag
      // direction: A^dag u_j - sigma_j v_j = betaK * X[Nm-1,j] * V[Nm].
      // Convergence criterion: |betaK * X[Nm-1, idx]| < tol * sigma_idx.
      int nconv = 0;
      for (int i = 0; i < Nk; ++i) {
        int idx = order(i);
        RealD res = std::abs(betaK * X(Nm - 1, idx));
        RealD thr = Tolerance * std::max(sigma(idx), 1.0e-14);
        std::cout << GridLogMessage
                  << "IRLBA iter " << iter
                  << "  sigma[" << i << "] = " << sigma(idx)
                  << "  res = " << res
                  << "  thr = " << thr << std::endl;
        if (res < thr) ++nconv;
        else break;  // residuals not strictly ordered but break is conservative
      }

      if (nconv >= Nk) {
        std::cout << GridLogMessage
                  << "IRLBA converged: " << nconv << " singular values after "
                  << iter + 1 << " iterations." << std::endl;
        // Collect converged triplets
        extractTriplets(Nm, sigma, X, Y, order, Nk);
        return;
      }

      // ----------------------------------------------------------------
      // Step 4: implicit restart — compress to Nk steps
      // ----------------------------------------------------------------
      implicitRestart(Nm, Nk, sigma, X, Y, order, betaK, betaRestart);
      verify();

      // Lucky breakdown: exact invariant subspace found; convergence is exact.
      // B_p^+ = diag(alpha[0..Nk-1]); extract directly from restart basis.
      if (betaRestart < 1.0e-14) {
        std::cout << GridLogMessage
                  << "IRLBA: lucky breakdown after restart (betaRestart = 0)."
                  << " Extracting " << Nk << " exact Ritz triplets." << std::endl;
        // Re-run SVD on the p-step diagonal B^+ to get sorted Ritz triplets.
        Eigen::MatrixXd Bp = buildBidiagonal(Nk);
        Eigen::JacobiSVD<Eigen::MatrixXd> svdp(Bp,
            Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::VectorXi ordp = sortOrder(svdp.singularValues());
        extractTriplets(Nk, svdp.singularValues(), svdp.matrixU(),
                        svdp.matrixV(), ordp, Nk);
        return;
      }

      pStart = Nk;
    }

    std::cout << GridLogMessage
              << "IRLBA: did not converge in " << MaxIter
              << " iterations. Returning best approximations." << std::endl;

    // Return best available approximations
    Eigen::MatrixXd B = buildFullB((int)alpha.size());
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(B,
        Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd sigma = svd.singularValues();
    Eigen::MatrixXd X     = svd.matrixU();
    Eigen::MatrixXd Y     = svd.matrixV();
    Eigen::VectorXi order = sortOrder(sigma);
    int nout = std::min(Nk, (int)alpha.size());
    extractTriplets((int)alpha.size(), sigma, X, Y, order, nout);
  }

  /* Getters */
  int getNk() const { return (int)singularValues.size(); }
  const std::vector<RealD>&  getSingularValues() const { return singularValues; }
  const std::vector<Field>&  getLeftVectors()    const { return leftVectors; }
  const std::vector<Field>&  getRightVectors()   const { return rightVectors; }

  /**
   * Print B_k and U^dag A V to verify the bidiagonalization relation
   *   A V_m = U_m B_m   (exact in our GK convention)
   * On the first call (pStart=0), max|B - U^dag A V| should be ~machine epsilon.
   * After a restart and extension, the column p of U^dag A V deviates from B
   * by O(betaK): this is expected because the thick restart breaks the Krylov
   * structure at column p, introducing off-diagonal terms proportional to betaK.
   * These terms vanish as betaK -> 0 (convergence), so the algorithm is correct.
   */
  void verify()
  {
    int m = (int)alpha.size();
    if (m == 0) { std::cout << GridLogMessage << "IRLBA verify: empty basis" << std::endl; return; }

    // Print B_k
    Eigen::MatrixXd B = buildBidiagonal(m);
    std::cout << GridLogMessage << "IRLBA verify: B_k (" << m << "x" << m << "):" << std::endl;
    for (int i = 0; i < m; ++i) {
      std::cout << GridLogMessage << "  row " << i << ": ";
      for (int j = 0; j < m; ++j) std::cout << B(i,j) << " ";
      std::cout << std::endl;
    }

    // Compute M[i,j] = <U[i] | A | V[j]>  (should equal B[i,j])
    std::cout << GridLogMessage << "IRLBA verify: U^dag A V (" << m << "x" << m << "):" << std::endl;
    Field Avj(Grid);
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(m, m);
    for (int j = 0; j < m && j < (int)V.size(); ++j) {
      Linop.Op(V[j], Avj);
      for (int i = 0; i < m && i < (int)U.size(); ++i) {
        ComplexD ip = innerProduct(U[i], Avj);
        M(i, j) = ip.real();
      }
    }
    for (int i = 0; i < m; ++i) {
      std::cout << GridLogMessage << "  row " << i << ": ";
      for (int j = 0; j < m; ++j) std::cout << M(i,j) << " ";
      std::cout << std::endl;
    }

    // Print max deviation |B - U^dag A V|
    RealD maxdev = (B - M).cwiseAbs().maxCoeff();
    std::cout << GridLogMessage << "IRLBA verify: max|B - U^dag A V| = " << maxdev << std::endl;

    // Also print beta (residual couplings)
    std::cout << GridLogMessage << "IRLBA verify: beta[0.." << (int)beta.size()-1 << "] = ";
    for (auto b : beta) std::cout << b << " ";
    std::cout << std::endl;
  }

private:

  // ------------------------------------------------------------------
  // Build the m x m upper-bidiagonal matrix from alpha[0..m-1], beta[0..m-2]
  // ------------------------------------------------------------------
  Eigen::MatrixXd buildBidiagonal(int m) const
  {
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m, m);
    for (int k = 0; k < m; ++k) {
      B(k, k) = alpha[k];
      if (k + 1 < m && k < (int)beta.size())
        B(k, k + 1) = beta[k];
    }
    return B;
  }

  // ------------------------------------------------------------------
  // Build the full m x m B matrix, including the non-bidiagonal column
  // at restart_col that arises after a thick restart.
  //
  // After restart, A V[restart_col] has projections onto all U[0..restart_col-1]
  // (not just U[restart_col-1]).  These are stored in fvec[0..restart_col-1]
  // and make column restart_col of U^dag A V non-bidiagonal.
  // ------------------------------------------------------------------
  Eigen::MatrixXd buildFullB(int m) const
  {
    Eigen::MatrixXd B = buildBidiagonal(m);
    if (restart_col >= 0 && restart_col < m && (int)fvec.size() > 0) {
      for (int j = 0; j < restart_col && j < (int)fvec.size(); ++j)
        B(j, restart_col) = fvec[j];
    }
    return B;
  }

  // ------------------------------------------------------------------
  // Return a permutation vector that puts the desired Nk singular values
  // first (largest first if largest==true, smallest first otherwise).
  // Eigen's JacobiSVD already returns sigma in descending order, so for
  // largest we just return 0,1,...,m-1; for smallest we reverse.
  // ------------------------------------------------------------------
  Eigen::VectorXi sortOrder(const Eigen::VectorXd &sigma) const
  {
    int m = (int)sigma.size();
    Eigen::VectorXi ord(m);
    if (largest) {
      for (int i = 0; i < m; ++i) ord(i) = i;
    } else {
      for (int i = 0; i < m; ++i) ord(i) = m - 1 - i;
    }
    return ord;
  }

  // ------------------------------------------------------------------
  // Extend the Lanczos bidiagonalization from pStart to kEnd steps.
  // On first call pStart==0 (V[0] already set).
  // On restart calls V[0..pStart], U[0..pStart-1], alpha[0..pStart-1],
  // beta[0..pStart-1] are already set; betaRestart is the coupling
  // beta_{pStart} that drives the first new U step.
  // ------------------------------------------------------------------
  void extendBasis(int pStart, int kEnd, RealD betaRestart)
  {
    // Truncate containers to pStart (Lattice has no default constructor)
    if ((int)V.size() > pStart + 1) V.erase(V.begin() + pStart + 1, V.end());
    if ((int)U.size() > pStart)     U.erase(U.begin() + pStart,     U.end());
    alpha.resize(pStart);
    beta.resize(pStart);

    Field p(Grid), r(Grid);

    for (int k = pStart; k < kEnd; ++k) {

      // p = A v_k
      Linop.Op(V[k], p);

      // Remove previous left vector coupling
      if (k > 0) {
        p = p - beta[k - 1] * U[k - 1];
      }
      // On the first step after a restart, beta[pStart-1] was already set;
      // but V[pStart] was already constructed including the beta correction,
      // so no extra subtraction needed here beyond the standard recurrence.

      // Reorthogonalize p against U, then alpha_k = ||p||, u_k = p/alpha_k
      reorthogonalize(p, U);
      RealD ak = std::sqrt(norm2(p));
      if (ak < 1.0e-14) {
        std::cout << GridLogMessage
                  << "IRLBA extendBasis: lucky breakdown at step " << k
                  << " (alpha = " << ak << ")" << std::endl;
        alpha.push_back(ak);
        Field zero(Grid); zero = Zero();
        U.push_back(zero);
        beta.push_back(0.0);
        V.push_back(zero);
        break;
      }
      alpha.push_back(ak);

      Field u(Grid);
      u = (1.0 / ak) * p;
      U.push_back(u);

      // r = A^dag u_k - alpha_k v_k, reorthogonalize, then beta_{k+1} = ||r||
      Linop.AdjOp(U[k], r);
      r = r - ak * V[k];
      reorthogonalize(r, V);

      RealD bk = std::sqrt(norm2(r));
      beta.push_back(bk);

      std::cout << GridLogMessage
                << "IRLBA extend step " << k
                << "  alpha = " << ak
                << "  beta  = " << bk << std::endl;

      // Always push v_{k+1} (needed as residual direction for restart)
      if (bk < 1.0e-14) {
        std::cout << GridLogMessage
                  << "IRLBA extendBasis: lucky breakdown (beta = 0) at step "
                  << k << std::endl;
        Field zero(Grid); zero = Zero();
        V.push_back(zero);
        break;
      }
      Field vnext(Grid);
      vnext = (1.0 / bk) * r;
      V.push_back(vnext);

      if (k == kEnd - 1) break;  // v_{k+1} pushed above; stop here
    }
  }

  // ------------------------------------------------------------------
  // Full reorthogonalization of vec against the vectors in basis.
  // Subtracts projections only — does NOT normalize.
  // ------------------------------------------------------------------
  void reorthogonalize(Field &vec, const std::vector<Field> &basis)
  {
    for (int j = 0; j < (int)basis.size(); ++j) {
      ComplexD ip = innerProduct(basis[j], vec);
      vec = vec - ip * basis[j];
    }
    // Second pass for numerical stability
    for (int j = 0; j < (int)basis.size(); ++j) {
      ComplexD ip = innerProduct(basis[j], vec);
      vec = vec - ip * basis[j];
    }
  }

  // ------------------------------------------------------------------
  // Implicit restart: given the Nm-step bidiagonalization and its SVD,
  // compress to Nk steps via implicit QR shifts applied to B_k.
  //
  // The "shifts" are the Nm - Nk singular values we want to deflate
  // (those NOT in the desired set).  We apply them as implicit QR steps
  // to the bidiagonal matrix, then update the lattice bases accordingly.
  //
  // After this call:
  //   V[0..Nk],  U[0..Nk-1],  alpha[0..Nk-1],  beta[0..Nk-1]  are updated.
  //   betaRestart  ← new beta_Nk coupling for the next extension.
  // ------------------------------------------------------------------
  void implicitRestart(int k, int p,
                       const Eigen::VectorXd &sigma,
                       const Eigen::MatrixXd &X,
                       const Eigen::MatrixXd &Y,
                       const Eigen::VectorXi &order,
                       RealD betaK,
                       RealD &betaRestart)
  {
    // Thick restart (Baglama & Reichel, Sec. 2.2):
    //
    // Given B_k = X Sigma Y^T, define the new p-step basis by:
    //   V^+_i = V_k * y_{order(i)}      (right sing. vec. of B_k)
    //   U^+_i = U_k * x_{order(i)}      (left  sing. vec. of B_k)
    //
    // Then A V^+_i = A V_k y_{order(i)} = U_k B_k y_{order(i)}
    //             = sigma_{order(i)} U_k x_{order(i)} = sigma_{order(i)} U^+_i
    //
    // So B_p^+ = diag(sigma_{order(0)}, ..., sigma_{order(p-1)}) — DIAGONAL,
    // all internal betas are zero.
    //
    // The residual coupling comes from A^dag U_k = V_k B_k^T + betaK V[k] e_{k-1}^T:
    //   A^dag U^+_{p-1} - sigma_{order(p-1)} V^+_{p-1}
    //     = V_k (B_k^T x_{order(p-1)} - sigma_{order(p-1)} y_{order(p-1)})
    //       + betaK * X(k-1, order(p-1)) * V[k]
    //     = betaK * X(k-1, order(p-1)) * V[k]   (since B_k^T x_j = sigma_j y_j)
    //
    // Therefore: betaRestart = |betaK * X(k-1, order(p-1))|
    //            V[p] = sign(X(k-1, order(p-1))) * V[k]

    // ---- Build new lattice vectors ----
    std::vector<Field> Vnew, Unew;
    Vnew.reserve(p + 1);
    Unew.reserve(p);

    for (int i = 0; i < p; ++i) {
      int idx = order(i);
      Field vi(Grid); vi = Zero();
      for (int j = 0; j < k; ++j)
        vi = vi + Y(j, idx) * V[j];
      Vnew.push_back(vi);
    }

    for (int i = 0; i < p; ++i) {
      int idx = order(i);
      Field ui(Grid); ui = Zero();
      for (int j = 0; j < k; ++j)
        ui = ui + X(j, idx) * U[j];
      Unew.push_back(ui);
    }

    // New v_{p} (0-indexed: V[p]) = sign * V[k]
    // From A^dag U_k = V_k B_k^T + betaK V[k] e_{k-1}^T:
    //   A^dag U^+_j - sigma_j V^+_j = betaK * X(k-1, order(j)) * V[k]
    // The last Ritz pair (j=p-1) defines betaRestart and the sign of V[p].
    // All p couplings (j=0..p-1) are stored in fvec so that buildFullB can
    // reconstruct the exact column p of U^dag A V after the next extension.
    RealD coeff = betaK * X(k - 1, order(p - 1));
    betaRestart  = std::abs(coeff);
    RealD sgn = (coeff >= 0.0) ? 1.0 : -1.0;

    fvec.resize(p);
    for (int j = 0; j < p; ++j)
      fvec[j] = betaK * X(k - 1, order(j)) * sgn;
    // fvec[p-1] == betaRestart by construction
    restart_col = p;

    Field vp(Grid);
    if (betaRestart > 1.0e-14) {
      vp = sgn * V[k];
    } else {
      betaRestart = 0.0;
      vp = Zero();
    }
    Vnew.push_back(vp);  // V[p]

    // ---- New alpha, beta ----
    // B_p^+ is diagonal: alpha^+_i = sigma_{order(i)}, all internal beta = 0
    std::vector<RealD> alpha_new(p), beta_new(p);
    for (int i = 0; i < p; ++i) alpha_new[i] = sigma(order(i));
    for (int i = 0; i < p - 1; ++i) beta_new[i] = 0.0;
    beta_new[p - 1] = betaRestart;

    // ---- Commit new state ----
    V = Vnew;
    U = Unew;
    alpha = alpha_new;
    beta  = beta_new;

    std::cout << GridLogMessage
              << "IRLBA restart: compressed to " << p << " steps,"
              << "  new beta_p = " << betaRestart << std::endl;
  }

  // ------------------------------------------------------------------
  // Extract the desired singular triplets into the public output vectors.
  // ------------------------------------------------------------------
  void extractTriplets(int m,
                       const Eigen::VectorXd &sigma,
                       const Eigen::MatrixXd &X,
                       const Eigen::MatrixXd &Y,
                       const Eigen::VectorXi &order,
                       int nout)
  {
    singularValues.resize(nout);
    leftVectors.clear();   leftVectors.reserve(nout);
    rightVectors.clear();  rightVectors.reserve(nout);

    for (int i = 0; i < nout; ++i) {
      int idx = order(i);
      singularValues[i] = sigma(idx);

      // Left singular vector of A:  svec_L = U_m * x_i
      Field svL(Grid); svL = Zero();
      for (int j = 0; j < m && j < (int)U.size(); ++j)
        svL = svL + X(j, idx) * U[j];
      leftVectors.push_back(svL);

      // Right singular vector of A:  svec_R = V_m * y_i
      Field svR(Grid); svR = Zero();
      for (int j = 0; j < m && j < (int)V.size(); ++j)
        svR = svR + Y(j, idx) * V[j];
      rightVectors.push_back(svR);
    }
  }
};

NAMESPACE_END(Grid);
#endif
