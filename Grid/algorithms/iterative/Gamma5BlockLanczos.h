/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    γ5-Block Lanczos algorithm for γ5-Hermitian operators (Wilson Dirac).

    Based on: S. Yamamoto, "γ5-Block Krylov (Block Lanczos) Methods
    for the Wilson Dirac Operator", April 2026.

    The Wilson Dirac operator D_W satisfies D_W† = γ5 D_W γ5, so it is
    self-adjoint in the indefinite γ5-inner product (u,v) ≡ u†γ5v.

    Block size is fixed at s=2.  Each starting block is Q_1 = [v, γ5v].
    The three-term block recurrence

        Q_{k+1} B_{k+1} = D_W Q_k − Q_k A_k − Q_{k-1} C_k

    produces a block tridiagonal projected matrix T_m whose eigenvalues
    (via a general non-Hermitian solver) approximate eigenvalues of D_W
    directly — not those of H_W = γ5 D_W.

*************************************************************************************/
#ifndef GRID_GAMMA5_BLOCK_LANCZOS_H
#define GRID_GAMMA5_BLOCK_LANCZOS_H

#include <functional>
#include <numeric>
#include <iomanip>

NAMESPACE_BEGIN(Grid);

template<class Field>
class Gamma5BlockLanczos {
public:
  using Gamma5Func = std::function<void(const Field&, Field&)>;

private:
  typedef Eigen::Matrix2cd   CMat2;
  typedef Eigen::MatrixXcd   CMat;
  typedef Eigen::VectorXcd   CVec;
  typedef Eigen::SelfAdjointEigenSolver<CMat2> SAEigen2;

  LinearOperatorBase<Field>& Linop;
  GridBase*                  Grid_;
  Gamma5Func                 applyGamma5;
  RealD                      Tolerance;

  // Per-step 2×2 coefficient blocks.
  // Index k here corresponds to paper's step k+1.
  // A_blocks[k] = A_{k+1}, B_blocks[k] = B_{k+2}, C_blocks[k] = C_{k+1}.
  std::vector<CMat2> A_blocks;
  std::vector<CMat2> B_blocks;   // B_blocks[k] = normalization block after step k
  std::vector<CMat2> C_blocks;
  std::vector<CMat2> G_blocks;   // G_blocks[k] = G_{k+1}

  // Krylov basis: basis[2k], basis[2k+1] are the two columns of Q_{k+1}.
  std::vector<Field> basis;

  int nSteps;   // number of completed steps

  // Krylov-Schur implicit restart state
  CMat  Hmat_;           // full projected matrix (2*nSteps x 2*nSteps)
  bool  useFullH_;       // true while in Krylov-Schur extension mode
  int   Ncompressed_;    // number of compressed column vectors kept after last KS step
  CMat2 Blast_;          // last normalization block from lanczosStepFull (for Ritz estimate)

  // Output
  CVec               evals_;
  std::vector<Field> evecs_;
  std::vector<RealD> residuals_;

public:
  bool doEvalCheck = false;

  Gamma5BlockLanczos(LinearOperatorBase<Field>& op, GridBase* grid,
                     Gamma5Func g5, RealD tol = 1e-8)
    : Linop(op), Grid_(grid), applyGamma5(g5), Tolerance(tol), nSteps(0),
      useFullH_(false), Ncompressed_(0)
  {
    Blast_ = CMat2::Zero();
  }

  CVec               getEvals()     { return evals_;     }
  std::vector<Field> getEvecs()     { return evecs_;     }
  std::vector<RealD> getResiduals() { return residuals_; }

  /**
   * Run γ5-Block Lanczos.
   *
   * v0       : starting vector (must not be a chiral eigenstate γ5 v = ±v)
   * maxSteps : maximum Lanczos steps (each adds 2 basis vectors and 2 Ritz values)
   * Nstop    : target converged pairs (informational; all pairs are always returned)
   * reorthog : full γ5-reorthogonalisation at each step (fixes finite-precision drift)
   */
  void operator()(const Field& v0, int maxSteps, int Nstop, bool reorthog = false)
  {
    basis.clear();
    A_blocks.clear();
    B_blocks.clear();
    C_blocks.clear();
    G_blocks.clear();
    nSteps = 0;

    // Initialise Q_1 = [v, γ5v]
    Field v(Grid_), g5v(Grid_);
    v = v0;
    RealD nrm = std::sqrt(norm2(v));
    assert(nrm > 1e-14);
    v *= (1.0 / nrm);
    applyGamma5(v, g5v);

    CMat2 G1 = gramMatrix(v, g5v);
    ComplexD detG1 = G1(0,0)*G1(1,1) - G1(0,1)*G1(1,0);
    RealD absdetG1 = std::sqrt(detG1.real()*detG1.real() + detG1.imag()*detG1.imag());
    std::cout << GridLogMessage
              << "Gamma5BlockLanczos: G1 = \n" << G1
              << "\n  det G1 = " << detG1 << std::endl;
    if (absdetG1 < 1e-13) {
      std::cout << GridLogMessage
                << "Gamma5BlockLanczos: abort — degenerate start "
                << "(v is a chiral eigenstate or |det G1| < 1e-13)" << std::endl;
      return;
    }

    G_blocks.push_back(G1);
    basis.push_back(v);
    basis.push_back(g5v);

    for (int step = 0; step < maxSteps; step++) {
      bool ok = lanczosStep(step, reorthog);
      if (!ok) break;
      nSteps = step + 1;

      RealD beta = B_blocks[step].norm();
      std::cout << GridLogMessage << "Gamma5BlockLanczos: step " << step
                << "  beta = " << beta << std::endl;
      if (beta < Tolerance) {
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos: beta < tol, converged at step "
                  << step << std::endl;
        break;
      }
    }

    if (nSteps == 0) return;
    computeRitzPairs(nSteps, Nstop);
  }

  /**
   * Restarted γ5-Block Lanczos (explicit restart).
   *
   * Runs operator() for Nstep steps per pass.  After each pass the 2*Nstep
   * Ritz pairs are sorted by the chosen RitzFilter and the top Nk are the
   * "wanted" set.  The seed for the next pass is the normalised equal-weight
   * sum of the Nstop best Ritz vectors from the wanted set; this keeps the
   * restart in the span of the most-wanted approximate eigenspace.
   *
   * Convergence is declared when at least Nstop of the top-Nk pairs have
   * residual < tolerance.
   *
   * v0          : initial starting vector
   * maxRestarts : maximum number of Lanczos passes
   * Nstep       : Lanczos steps per pass  (produces 2*Nstep Ritz values)
   * Nk          : size of the "wanted" set to track (must be >= Nstop)
   * Nstop       : target converged pairs; also the number of vectors summed
   *               to form the restart seed
   * reorthog    : full γ5-reorthogonalisation within each pass
   * filter      : EvalNormSmall   → sort by |λ| ascending
   *               EvalImNormSmall → sort by |Im(λ)| ascending (default)
   */
  void restart(const Field& v0, int maxRestarts, int Nstep, int Nk, int Nstop,
               bool reorthog = false, RitzFilter filter = EvalImNormSmall)
  {
    assert(Nk >= Nstop);
    Field src(Grid_);
    src = v0;

    for (int iter = 0; iter < maxRestarts; iter++) {
      std::cout << GridLogMessage
                << "Gamma5BlockLanczos: ---- restart " << iter << " ----" << std::endl;

      (*this)(src, Nstep, Nstop, reorthog);

      int nRitz = (int)residuals_.size();
      if (nRitz == 0) {
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos: restart — no Ritz pairs, stopping." << std::endl;
        return;
      }

      // Sort all Ritz indices by the chosen filter criterion.
      std::vector<int> idx(nRitz);
      std::iota(idx.begin(), idx.end(), 0);
      switch (filter) {
        case EvalNormSmall:
          std::sort(idx.begin(), idx.end(), [&](int a, int b){
            return std::abs(evals_(a)) < std::abs(evals_(b));
          });
          break;
        case EvalImNormSmall:
        default:
          std::sort(idx.begin(), idx.end(), [&](int a, int b){
            return std::abs(evals_(a).imag()) < std::abs(evals_(b).imag());
          });
          break;
      }

      // Count converged pairs within the top-Nk wanted set.
      int nKeep = std::min(Nk, nRitz);
      int nconv = 0;
      for (int i = 0; i < nKeep; i++)
        if (residuals_[idx[i]] < Tolerance) nconv++;

      std::cout << GridLogMessage
                << "Gamma5BlockLanczos: restart " << iter
                << "  Ritz = " << nRitz
                << "  wanted = " << nKeep
                << "  converged = " << nconv << " / " << Nstop << std::endl;
      for (int i = 0; i < nKeep; i++)
        std::cout << GridLogMessage
                  << "    wanted[" << i << "]  lambda = " << evals_(idx[i])
                  << "  |res| = " << residuals_[idx[i]]
                  << (residuals_[idx[i]] < Tolerance ? "  *" : "") << std::endl;

      if (nconv >= Nstop) {
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos: converged after " << iter + 1
                  << " restart(s)." << std::endl;
        reorderOutput(idx, nKeep);
        return;
      }

      // Build restart seed: equal-weight sum of the top Nstop Ritz vectors
      // (sorted by filter criterion).  Spans the best part of the wanted
      // eigenspace and avoids locking onto a single approximate eigenvalue.
      int nSeed = std::min(Nstop, nKeep);
      src = Zero();
      for (int i = 0; i < nSeed; i++)
        src += evecs_[idx[i]];

      RealD nrm = std::sqrt(norm2(src));
      assert(nrm > 1e-14);
      src *= (1.0 / nrm);

      std::cout << GridLogMessage
                << "Gamma5BlockLanczos: seed = sum of top " << nSeed
                << " Ritz vectors  ||seed|| = " << nrm << std::endl;
    }

    std::cout << GridLogMessage
              << "Gamma5BlockLanczos: max restarts (" << maxRestarts
              << ") reached without full convergence." << std::endl;
  }

  /**
   * Implicitly Restarted Block Lanczos (Krylov-Schur variant).
   *
   * Each outer cycle:
   *   1. Runs Nmax block Lanczos steps (first cycle from v0; subsequently
   *      extends the Nk-step compressed state).
   *   2. Applies a Krylov-Schur restart: Schur-decomposes T_{Nmax}, reorders
   *      by the chosen filter, and compresses the basis from 2*Nmax to Nk
   *      columns.  The leading Nk×Nk upper-triangular Schur block S_k replaces
   *      the old block-tridiagonal T.
   *   3. Extends from Nk/2 to Nmax steps using full γ5-projection (block
   *      Arnoldi with the complete Nk-vector compressed basis).
   *
   * Convergence is declared when Nstop pairs have residual < tolerance.
   *
   * v0       : initial starting vector
   * maxIter  : maximum restart cycles (excluding initial run)
   * Nmax     : Lanczos steps per cycle  (builds 2*Nmax Ritz pairs)
   * Nk       : steps kept after compression (even, 2 ≤ Nk < Nmax, Nk ≥ Nstop)
   * Nstop    : target converged pairs
   * reorthog : γ5-reorthogonalisation in the initial Nmax-step run
   * filter   : eigenvalue selection criterion (default: EvalImNormSmall)
   */
  void implicitRestart(const Field& v0, int maxIter, int Nmax, int Nk, int Nstop,
                       bool reorthog = false, RitzFilter filter = EvalImNormSmall)
  {
    assert(Nk % 2 == 0 && Nk >= 2 && Nk < Nmax);
    assert(Nk >= Nstop);

    // Initial full block-Lanczos run
    (*this)(v0, Nmax, Nstop, reorthog);

    for (int iter = 0; iter < maxIter; iter++) {
      std::cout << GridLogMessage
                << "Gamma5BlockLanczos::implicitRestart ---- cycle " << iter << " ----" << std::endl;

      // Sort current Ritz pairs by filter, count converged
      int nRitz = (int)residuals_.size();
      std::vector<int> idx = sortedIdx(nRitz, filter);
      int nKeep = std::min(Nk, nRitz);
      int nconv = 0;
      for (int i = 0; i < nKeep; i++)
        if (residuals_[idx[i]] < Tolerance) nconv++;

      std::cout << GridLogMessage
                << "  nRitz=" << nRitz
                << "  nconv=" << nconv << "/" << Nstop << std::endl;
      for (int i = 0; i < nKeep; i++)
        std::cout << GridLogMessage
                  << "  [" << i << "]  lambda=" << evals_(idx[i])
                  << "  |res|=" << residuals_[idx[i]]
                  << (residuals_[idx[i]] < Tolerance ? "  *" : "") << std::endl;

      if (nconv >= Nstop) {
        reorderOutput(idx, nKeep);
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos::implicitRestart: converged after "
                  << iter + 1 << " cycle(s)." << std::endl;
        useFullH_ = false;
        return;
      }

      // Krylov-Schur: compress T_{Nmax} to S_{Nk} (upper-triangular Nk×Nk)
      krylovSchurCompress(Nk, filter);

      // Extend from Nk/2 to Nmax steps via full-projection block Arnoldi
      for (int step = nSteps; step < Nmax; step++) {
        bool ok = lanczosStepFull(step);
        if (!ok) break;
        nSteps = step + 1;
        if (Blast_.norm() < Tolerance) {
          std::cout << GridLogMessage
                    << "Gamma5BlockLanczos::implicitRestart: beta < tol at full step "
                    << step << ", stopping extension." << std::endl;
          break;
        }
      }

      // Ritz pairs from the full projected matrix
      computeRitzPairsFull(nSteps, Nstop);
    }

    // maxIter exhausted
    std::cout << GridLogMessage
              << "Gamma5BlockLanczos::implicitRestart: maxIter=" << maxIter
              << " reached without full convergence." << std::endl;
    int nRitz = (int)residuals_.size();
    if (nRitz > 0) {
      std::vector<int> idx = sortedIdx(nRitz, filter);
      reorderOutput(idx, std::min(Nk, nRitz));
    }
    useFullH_ = false;
  }

  /**
   * Verify the block Lanczos decomposition after operator() has run.
   *
   * Checks the three-term recurrence
   *
   *   D_W Q_k = Q_{k-1} C_k + Q_k A_k + Q_{k+1} B_{k+1}
   *
   * for each completed step k=1..m, plus γ5-orthonormality and
   * off-diagonal γ5-orthogonality between all Krylov blocks.
   *
   * Prints:
   *   - T_m (block tridiagonal projected matrix, 2m × 2m)
   *   - A_k, B_k, C_k, G_k for each step
   *   - max |G_computed[k] - G_stored[k]|  (γ5-Gram diagonal consistency)
   *   - max |Q_i†γ5 Q_j| for i≠j           (inter-block γ5-orthogonality)
   *   - per-column recurrence residual || D_W q - Q_{k-1}C - Q_k A - Q_{k+1}B ||
   */
  void verify(const std::string& label = "")
  {
    int m = nSteps;
    if (m == 0) {
      std::cout << GridLogMessage
                << "Gamma5BlockLanczos::verify: no steps completed." << std::endl;
      return;
    }

    std::cout << GridLogMessage
              << "======== Gamma5BlockLanczos::verify [" << label << "] ========" << std::endl;
    std::cout << GridLogMessage
              << "  m = " << m << " completed steps,  basis vectors = " << basis.size() << std::endl;

    // ---- Assemble and print T_m ----
    int dim = 2 * m;
    CMat Tm = CMat::Zero(dim, dim);
    for (int k = 0; k < m; k++) {
      Tm.block(2*k, 2*k, 2, 2) = A_blocks[k];
      if (k < m - 1) {
        Tm.block(2*k+2, 2*k,   2, 2) = B_blocks[k];
        Tm.block(2*k,   2*k+2, 2, 2) = C_blocks[k+1];
      }
    }
    std::cout << GridLogMessage << "T_m (" << dim << " x " << dim << "):" << std::endl;
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++)
        std::cout << "  " << std::setprecision(6) << std::setw(16) << Tm(i, j);
      std::cout << std::endl;
    }

    // ---- Print per-step coefficient blocks ----
    for (int k = 0; k < m; k++) {
      std::cout << GridLogMessage << "  A[" << k << "] =\n" << A_blocks[k] << std::endl;
      std::cout << GridLogMessage << "  B[" << k << "] =\n" << B_blocks[k] << std::endl;
      std::cout << GridLogMessage << "  C[" << k << "] =\n" << C_blocks[k] << std::endl;
      std::cout << GridLogMessage << "  G[" << k << "] =\n" << G_blocks[k] << std::endl;
    }
    std::cout << GridLogMessage << "  G[" << m << "] =\n" << G_blocks[m] << std::endl;

    // ---- Check γ5-Gram consistency: compare stored G_blocks[k] with recomputed Q_k†γ5Q_k ----
    // basis[2k..2k+1] = code block k (= paper Q_{k+1}); G_blocks[k] = paper G_{k+1}.
    RealD maxGramErr = 0.0;
    for (int k = 0; k <= m; k++) {
      CMat2 Gcomp = gramMatrix(basis[2*k], basis[2*k + 1]);
      CMat2 Gerr  = Gcomp - G_blocks[k];
      RealD err   = Gerr.cwiseAbs().maxCoeff();
      maxGramErr  = std::max(maxGramErr, err);
      std::cout << GridLogMessage
                << "  |G_computed[" << k << "] - G_stored[" << k << "]|_max = " << err << std::endl;
    }
    std::cout << GridLogMessage
              << "  max γ5-Gram diagonal error = " << maxGramErr << std::endl;

    // ---- Check γ5-orthogonality between different blocks: Q_i†γ5 Q_j ≈ 0 for i≠j ----
    RealD maxOffDiag = 0.0;
    for (int i = 0; i <= m; i++) {
      for (int j = 0; j <= m; j++) {
        if (i == j) continue;
        CMat2 Mij = g5InnerBlock(basis[2*i], basis[2*i+1], basis[2*j], basis[2*j+1]);
        RealD err = Mij.cwiseAbs().maxCoeff();
        maxOffDiag = std::max(maxOffDiag, err);
      }
    }
    std::cout << GridLogMessage
              << "  max |Q_i†γ5 Q_j| (i≠j, should be ~0) = " << maxOffDiag << std::endl;

    // ---- Check three-term recurrence for each code block k=0..m-1 ----
    // Recurrence (code indices):
    //   D_W basis[2k..2k+1] = basis[2k..2k+1] * A_blocks[k]
    //                       + basis[2(k-1)..] * C_blocks[k]   (k > 0; C_blocks[0] = 0)
    //                       + basis[2(k+1)..] * B_blocks[k]
    RealD maxRecErr = 0.0;
    Field p1(Grid_), p2(Grid_), r1(Grid_), r2(Grid_);
    for (int k = 0; k < m; k++) {
      const Field& q1  = basis[2*k];
      const Field& q2  = basis[2*k + 1];
      const Field& qn1 = basis[2*(k+1)];
      const Field& qn2 = basis[2*(k+1) + 1];

      Linop.Op(q1, p1);
      Linop.Op(q2, p2);

      // subtract Q_k A_k
      r1 = p1 - q1 * A_blocks[k](0,0) - q2 * A_blocks[k](1,0);
      r2 = p2 - q1 * A_blocks[k](0,1) - q2 * A_blocks[k](1,1);

      // subtract Q_{k-1} C_k  (C_blocks[0] = 0, so k=0 is automatically fine)
      if (k > 0) {
        const Field& qp1 = basis[2*(k-1)];
        const Field& qp2 = basis[2*(k-1) + 1];
        r1 -= qp1 * C_blocks[k](0,0) + qp2 * C_blocks[k](1,0);
        r2 -= qp1 * C_blocks[k](0,1) + qp2 * C_blocks[k](1,1);
      }

      // subtract Q_{k+1} B_{k+1}  (= basis[2(k+1)..] * B_blocks[k])
      r1 -= qn1 * B_blocks[k](0,0) + qn2 * B_blocks[k](1,0);
      r2 -= qn1 * B_blocks[k](0,1) + qn2 * B_blocks[k](1,1);

      RealD dev1 = std::sqrt(norm2(r1));
      RealD dev2 = std::sqrt(norm2(r2));
      std::cout << GridLogMessage
                << "  recurrence k=" << k
                << ":  || D_W q[2k]   - ... || = " << dev1
                << "   || D_W q[2k+1] - ... || = " << dev2 << std::endl;
      maxRecErr = std::max({maxRecErr, dev1, dev2});
    }
    std::cout << GridLogMessage
              << "  max recurrence deviation = " << maxRecErr << std::endl;

    std::cout << GridLogMessage
              << "======== end Gamma5BlockLanczos::verify ========" << std::endl;
  }

private:
  // Return a permutation of [0,n) sorted by the chosen RitzFilter criterion.
  std::vector<int> sortedIdx(int n, RitzFilter filter)
  {
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    switch (filter) {
      case EvalNormSmall:
        std::sort(idx.begin(), idx.end(), [&](int a, int b){
          return std::abs(evals_(a)) < std::abs(evals_(b));
        });
        break;
      case EvalImNormSmall:
      default:
        std::sort(idx.begin(), idx.end(), [&](int a, int b){
          return std::abs(evals_(a).imag()) < std::abs(evals_(b).imag());
        });
        break;
    }
    return idx;
  }

  // Krylov-Schur: compress the current nSteps-step block-Lanczos basis to Nk
  // steps by Schur-decomposing T_{nSteps}, reordering "wanted" eigenvalues to
  // the top, and rotating the field basis accordingly.
  //
  // After this call:
  //   basis[0..Nk-1]  = Nk compressed Schur vectors (ṽ_j = V_m U[:,j])
  //   basis[Nk..Nk+1] = original outer-residual block Q_{m+1} (unchanged)
  //   G_blocks[0..Nk/2] recomputed from the new pairs
  //   Hmat_  = leading Nk×Nk block of the Schur form S (upper triangular)
  //   A/B/C blocks cleared; nSteps = Nk/2; useFullH_ = true
  void krylovSchurCompress(int Nk, RitzFilter filter)
  {
    int m   = nSteps;
    int dim = 2 * m;
    assert(Nk > 0 && Nk % 2 == 0 && Nk < dim);

    // Assemble T_m (block tridiagonal, dim×dim)
    CMat Tm = CMat::Zero(dim, dim);
    for (int k = 0; k < m; k++) {
      Tm.block(2*k, 2*k, 2, 2) = A_blocks[k];
      if (k < m - 1) {
        Tm.block(2*k+2, 2*k,   2, 2) = B_blocks[k];
        Tm.block(2*k,   2*k+2, 2, 2) = C_blocks[k+1];
      }
    }

    // Complex Schur decomposition with wanted eigenvalues first.
    // ComplexSchurDecomposition convention: A = Q† S Q,
    //   getMatrixQ() = U† (so U = Q†  →  U.adjoint() = Q)
    //   getMatrixS() = S  (upper triangular)
    // New basis: ṽ_j = V_m U[:,j]  with U = getMatrixQ().adjoint()
    ComplexSchurDecomposition schur(Tm, false, filter);
    schur.schurReorder(Nk);
    CMat U = schur.getMatrixQ().adjoint();   // rotation matrix (dim×dim unitary)
    CMat S = schur.getMatrixS();             // upper-triangular Schur form

    // Build Nk compressed field vectors + keep Q_{m+1} as "next block"
    std::vector<Field> new_basis;
    new_basis.reserve(Nk + 2);
    for (int j = 0; j < Nk; j++) {
      Field col(Grid_);
      col = Zero();
      for (int k = 0; k < dim; k++)
        col += basis[k] * U(k, j);
      new_basis.push_back(col);
    }
    // Outer residual block Q_{m+1} unchanged (extension starts from here)
    new_basis.push_back(basis[dim]);
    new_basis.push_back(basis[dim + 1]);

    basis = new_basis;   // Nk + 2 field vectors

    // Recompute G_blocks for each pair of new basis vectors
    G_blocks.clear();
    for (int k = 0; k <= Nk / 2; k++)
      G_blocks.push_back(gramMatrix(basis[2*k], basis[2*k+1]));

    // Compressed projected matrix = leading Nk×Nk block of S
    Hmat_ = S.block(0, 0, Nk, Nk);

    // Clear three-term block storage (not valid after rotation)
    A_blocks.clear();
    B_blocks.clear();
    C_blocks.clear();

    useFullH_    = true;
    Ncompressed_ = Nk;
    nSteps       = Nk / 2;

    std::cout << GridLogMessage
              << "Gamma5BlockLanczos::krylovSchurCompress: compressed to Nk=" << Nk
              << "  Hmat=" << Nk << "x" << Nk << std::endl;
  }

  // One block Arnoldi step with full γ5-projection against all previous blocks.
  // Used after krylovSchurCompress to extend from Nk/2 to Nmax steps.
  //
  // Computes all coupling entries T_{j,step} = G_j^{-1} Q_j† γ5 D_W Q_step
  // (j=0..step) via classical Gram-Schmidt, stores them in a new column of
  // Hmat_, subtracts the projections to form the residual, then LDL†-normalises
  // to produce Q_{step+1}.  Appends to basis and G_blocks; sets Blast_.
  bool lanczosStepFull(int step)
  {
    const Field& q1 = basis[2*step];
    const Field& q2 = basis[2*step + 1];
    CMat2 Gk = G_blocks[step];

    int old_dim = 2 * step;
    assert(Hmat_.rows() == old_dim && Hmat_.cols() == old_dim);

    // Apply D_W to current block
    Field p1(Grid_), p2(Grid_);
    Linop.Op(q1, p1);
    Linop.Op(q2, p2);

    // Grow Hmat_ by 2 (new column; lower rows left zero by Arnoldi projection)
    int new_dim = old_dim + 2;
    CMat Hmat_new = CMat::Zero(new_dim, new_dim);
    Hmat_new.block(0, 0, old_dim, old_dim) = Hmat_;

    // Coupling to all previous blocks (classical GS: use original p1,p2)
    Field r1(Grid_), r2(Grid_);
    r1 = p1; r2 = p2;
    for (int j = 0; j < step; j++) {
      CMat2 Mj = g5InnerBlock(basis[2*j], basis[2*j+1], p1, p2);
      CMat2 Hj = invert2x2(G_blocks[j]) * Mj;
      Hmat_new.block(2*j, old_dim, 2, 2) = Hj;
      r1 -= basis[2*j]   * Hj(0,0) + basis[2*j+1] * Hj(1,0);
      r2 -= basis[2*j]   * Hj(0,1) + basis[2*j+1] * Hj(1,1);
    }

    // On-diagonal block A_{step}
    CMat2 Mk = g5InnerBlock(q1, q2, p1, p2);
    CMat2 Ak = invert2x2(Gk) * Mk;
    Hmat_new.block(old_dim, old_dim, 2, 2) = Ak;
    r1 -= q1 * Ak(0,0) + q2 * Ak(1,0);
    r2 -= q1 * Ak(0,1) + q2 * Ak(1,1);

    // LDL† normalisation of the residual block
    CMat2       Gamma_k = gramMatrix(r1, r2);
    SAEigen2    es(Gamma_k);
    Eigen::Vector2d D  = es.eigenvalues();
    CMat2           U2 = es.eigenvectors();

    for (int j = 0; j < 2; j++) {
      if (std::abs(D(j)) < 1e-28) {
        Field rj = r1 * U2(0,j) + r2 * U2(1,j);
        if (std::sqrt(norm2(rj)) < Tolerance) {
          std::cout << GridLogMessage
                    << "Gamma5BlockLanczos: happy breakdown (full step " << step << ")" << std::endl;
        } else {
          std::cout << GridLogMessage
                    << "Gamma5BlockLanczos: serious breakdown (full step " << step
                    << ") — stopping." << std::endl;
          return false;
        }
      }
    }

    CMat2 Gkp1 = CMat2::Zero();
    Gkp1(0,0) = ComplexD((D(0) > 0.0) ? 1.0 : -1.0, 0.0);
    Gkp1(1,1) = ComplexD((D(1) > 0.0) ? 1.0 : -1.0, 0.0);
    double sqd0 = std::sqrt(std::abs(D(0)));
    double sqd1 = std::sqrt(std::abs(D(1)));

    CMat2 Bkp1;
    Bkp1.row(0) = U2.col(0).adjoint() * sqd0;
    Bkp1.row(1) = U2.col(1).adjoint() * sqd1;

    Field qnew1 = (r1 * U2(0,0) + r2 * U2(1,0)) * (1.0 / sqd0);
    Field qnew2 = (r1 * U2(0,1) + r2 * U2(1,1)) * (1.0 / sqd1);

    Hmat_   = Hmat_new;
    Blast_  = Bkp1;
    G_blocks.push_back(Gkp1);
    basis.push_back(qnew1);
    basis.push_back(qnew2);

    RealD beta = Bkp1.norm();
    std::cout << GridLogMessage
              << "Gamma5BlockLanczos: full step " << step << "  beta=" << beta << std::endl;
    return true;
  }

  // Compute Ritz pairs from the full projected matrix Hmat_ (set after
  // krylovSchurCompress + lanczosStepFull calls).  Ritz residual estimated
  // cheaply as ||Blast_ τ_j|| where τ_j = last 2 entries of eigenvector y_j.
  void computeRitzPairsFull(int m, int Nstop)
  {
    int dim = 2 * m;
    assert(Hmat_.rows() == dim && Hmat_.cols() == dim);

    Eigen::ComplexEigenSolver<CMat> ces(Hmat_);
    CVec lambdas = ces.eigenvalues();
    CMat Y       = ces.eigenvectors();

    // Sort by |Im(λ)| ascending (near-real = physical modes first)
    std::vector<int> idx(dim);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b){
      return std::abs(lambdas(a).imag()) < std::abs(lambdas(b).imag());
    });

    evals_.resize(dim);
    evecs_.clear();
    residuals_.clear();

    for (int ji = 0; ji < dim; ji++) {
      int  j  = idx[ji];
      evals_(ji) = lambdas(j);
      CVec yj = Y.col(j);

      // Ritz vector: ũ_j = sum_k basis[k] * yj[k]  (k = 0..dim-1)
      Field uj(Grid_);
      uj = Zero();
      for (int k = 0; k < dim; k++)
        uj += basis[k] * yj(k);
      evecs_.push_back(uj);

      // Ritz estimate: || Blast_ τ_j ||  (τ_j = last 2 entries of y_j)
      Eigen::Vector2cd tau(yj(dim - 2), yj(dim - 1));
      RealD res = (Blast_ * tau).norm();
      residuals_.push_back(res);

      std::cout << GridLogMessage
                << "Gamma5BlockLanczos (full): Ritz[" << ji << "]"
                << "  lambda=" << evals_(ji)
                << "  |res|=" << res << std::endl;
    }

    if (doEvalCheck) {
      Field w(Grid_);
      int nCheck = std::min((int)evecs_.size(), 2 * Nstop);
      for (int k = 0; k < nCheck; k++) {
        Linop.Op(evecs_[k], w);
        ComplexD eval_est = toStdCmplx(innerProduct(evecs_[k], w));
        w -= eval_est * evecs_[k];
        RealD res = std::sqrt(norm2(w));
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos: evec[" << k << "]"
                  << "  eval_reported=" << evals_(k)
                  << "  eval_est=" << eval_est
                  << "  ||Av-eval*v||=" << res << std::endl;
      }
    }
  }

  // Reorder evals_/evecs_/residuals_ so the first nKeep entries follow idx[].
  void reorderOutput(const std::vector<int>& idx, int nKeep)
  {
    CVec               evals_new(evals_.size());
    std::vector<Field> evecs_new;
    std::vector<RealD> res_new;
    evecs_new.reserve(evecs_.size());
    res_new.reserve(residuals_.size());

    for (int i = 0; i < nKeep; i++) {
      evals_new(i) = evals_(idx[i]);
      evecs_new.push_back(evecs_[idx[i]]);
      res_new.push_back(residuals_[idx[i]]);
    }
    // append remaining entries not in idx[0..nKeep-1]
    std::vector<bool> used(evals_.size(), false);
    for (int i = 0; i < nKeep; i++) used[idx[i]] = true;
    int j = nKeep;
    for (int i = 0; i < (int)evals_.size(); i++) {
      if (!used[i]) {
        evals_new(j++) = evals_(i);
        evecs_new.push_back(evecs_[i]);
        res_new.push_back(residuals_[i]);
      }
    }
    evals_     = evals_new;
    evecs_     = evecs_new;
    residuals_ = res_new;
  }

  // One Lanczos step.  On success pushes Q_{step+2} and returns true.
  bool lanczosStep(int step, bool reorthog)
  {
    const Field& q1 = basis[2*step];
    const Field& q2 = basis[2*step + 1];
    CMat2 Gk = G_blocks[step];

    // (i) Two matvecs: P_k = D_W Q_k
    Field p1(Grid_), p2(Grid_);
    Linop.Op(q1, p1);
    Linop.Op(q2, p2);

    // (ii) M_k = Q_k† γ5 P_k
    CMat2 Mk = g5InnerBlock(q1, q2, p1, p2);

    // (iii) A_k = G_k^{-1} M_k
    CMat2 Ak = invert2x2(Gk) * Mk;
    A_blocks.push_back(Ak);

    // (iv) C_k = G_{k-1}^{-1} B_k† G_k  (zero at step 0)
    CMat2 Ck = CMat2::Zero();
    if (step > 0) {
      CMat2 Gkm1 = G_blocks[step - 1];
      CMat2 Bk   = B_blocks[step - 1];
      Ck = invert2x2(Gkm1) * Bk.adjoint() * Gk;
    }
    C_blocks.push_back(Ck);

    // (v) Residual R̂_k = P_k - Q_k A_k - Q_{k-1} C_k
    // Column j: r̂^(j) = p^(j) - sum_i q^(i) A[i,j] - sum_i q_prev^(i) C[i,j]
    Field r1(Grid_), r2(Grid_);
    r1 = p1;
    r2 = p2;
    r1 -= q1 * Ak(0,0) + q2 * Ak(1,0);
    r2 -= q1 * Ak(0,1) + q2 * Ak(1,1);
    if (step > 0) {
      const Field& qp1 = basis[2*(step-1)];
      const Field& qp2 = basis[2*(step-1) + 1];
      r1 -= qp1 * Ck(0,0) + qp2 * Ck(1,0);
      r2 -= qp1 * Ck(0,1) + qp2 * Ck(1,1);
    }

    // Optional full γ5-reorthogonalisation against all previous blocks
    if (reorthog) {
      for (int j = 0; j <= step; j++) {
        CMat2 Mj = g5InnerBlock(basis[2*j], basis[2*j+1], r1, r2);
        CMat2 Hj = invert2x2(G_blocks[j]) * Mj;
        r1 -= basis[2*j] * Hj(0,0) + basis[2*j+1] * Hj(1,0);
        r2 -= basis[2*j] * Hj(0,1) + basis[2*j+1] * Hj(1,1);
      }
    }

    // (vi) Γ_k = R̂_k† γ5 R̂_k  (2×2 Hermitian)
    CMat2 Gamma_k = gramMatrix(r1, r2);

    // Eigendecompose Γ_k = U diag(d) U†
    SAEigen2 es(Gamma_k);
    Eigen::Vector2d D = es.eigenvalues();
    CMat2           U = es.eigenvectors();

    // Breakdown check
    for (int j = 0; j < 2; j++) {
      if (std::abs(D(j)) < 1e-28) {
        Field rj(Grid_);
        rj = r1 * U(0,j) + r2 * U(1,j);
        RealD rjnrm = std::sqrt(norm2(rj));
        if (rjnrm < Tolerance) {
          std::cout << GridLogMessage
                    << "Gamma5BlockLanczos: happy breakdown at step " << step
                    << " direction " << j << std::endl;
        } else {
          std::cout << GridLogMessage
                    << "Gamma5BlockLanczos: SERIOUS breakdown at step " << step
                    << " direction " << j
                    << " (look-ahead not implemented; stopping)" << std::endl;
          return false;
        }
      }
    }

    // (vii) Indefinite LDL† factorisation: G_{k+1}, B_{k+1}, Q_{k+1}
    CMat2 Gkp1 = CMat2::Zero();
    Gkp1(0,0) = ComplexD((D(0) > 0.0) ?  1.0 : -1.0, 0.0);
    Gkp1(1,1) = ComplexD((D(1) > 0.0) ?  1.0 : -1.0, 0.0);

    double sqd0 = std::sqrt(std::abs(D(0)));
    double sqd1 = std::sqrt(std::abs(D(1)));

    // B_{k+1} = diag(|d_0|^{1/2}, |d_1|^{1/2}) U†
    CMat2 Bkp1;
    Bkp1.row(0) = U.col(0).adjoint() * sqd0;
    Bkp1.row(1) = U.col(1).adjoint() * sqd1;

    // Q_{k+1} columns: q^(j) = R̂_k u_j / |d_j|^{1/2}
    Field qnew1(Grid_), qnew2(Grid_);
    qnew1 = (r1 * U(0,0) + r2 * U(1,0)) * (1.0 / sqd0);
    qnew2 = (r1 * U(0,1) + r2 * U(1,1)) * (1.0 / sqd1);

    G_blocks.push_back(Gkp1);
    B_blocks.push_back(Bkp1);
    basis.push_back(qnew1);
    basis.push_back(qnew2);

    return true;
  }

  // Assemble T_m and extract Ritz pairs (Algorithm 2 of the paper).
  void computeRitzPairs(int m, int Nstop)
  {
    int dim = 2 * m;

    // Assemble block tridiagonal T_m (2m × 2m).
    // T[k,k]   = A_{k+1}     = A_blocks[k]
    // T[k+1,k] = B_{k+2}     = B_blocks[k]
    // T[k,k+1] = C_{k+2}     = C_blocks[k+1]
    CMat Tm = CMat::Zero(dim, dim);
    for (int k = 0; k < m; k++) {
      Tm.block(2*k, 2*k, 2, 2) = A_blocks[k];
      if (k < m - 1) {
        Tm.block(2*k+2, 2*k,   2, 2) = B_blocks[k];
        Tm.block(2*k,   2*k+2, 2, 2) = C_blocks[k+1];
      }
    }

    std::cout << GridLogMessage << "Gamma5BlockLanczos: assembled T_m ("
              << dim << " x " << dim << ")" << std::endl;

    // General non-Hermitian eigensolver (T_m is γ5-symmetric, not Hermitian)
    Eigen::ComplexEigenSolver<CMat> ces(Tm);
    CVec lambdas = ces.eigenvalues();
    CMat Y       = ces.eigenvectors();

    // Sort by |Im(λ)| ascending (near-real = physical modes first)
    std::vector<int> idx(dim);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b){
      return std::abs(lambdas(a).imag()) < std::abs(lambdas(b).imag());
    });

    // B_{m+1} = B_blocks[m-1];  Q_{m+1} = basis[2m], basis[2m+1]
    const CMat2& Bm1 = B_blocks[m - 1];

    evals_.resize(dim);
    evecs_.clear();
    residuals_.clear();

    for (int ji = 0; ji < dim; ji++) {
      int  j  = idx[ji];
      evals_(ji) = lambdas(j);
      CVec yj = Y.col(j);

      // Ritz vector ũ_j = V_m y_j = sum_k Q_{k+1} * y_j[2k:2k+2]
      Field uj(Grid_);
      uj = Zero();
      for (int k = 0; k < m; k++) {
        uj += basis[2*k]     * yj(2*k);
        uj += basis[2*k + 1] * yj(2*k + 1);
      }
      evecs_.push_back(uj);

      // True residual: r_j = Q_{m+1} B_{m+1} τ_j,  τ_j = last 2 entries of y_j
      Eigen::Vector2cd tau(yj(dim-2), yj(dim-1));
      Eigen::Vector2cd Btau = Bm1 * tau;
      Field rj(Grid_);
      rj = basis[2*m] * Btau(0) + basis[2*m + 1] * Btau(1);
      RealD res = std::sqrt(norm2(rj));
      residuals_.push_back(res);

      std::cout << GridLogMessage << "Gamma5BlockLanczos: Ritz[" << ji << "]"
                << "  lambda = " << evals_(ji)
                << "  |res| = " << res << std::endl;
    }

    if (doEvalCheck) {
      Field w(Grid_);
      int nCheck = std::min((int)evecs_.size(), 2*Nstop);
      for (int k = 0; k < nCheck; k++) {
        Linop.Op(evecs_[k], w);
        ComplexD eval_est = toStdCmplx(innerProduct(evecs_[k], w));
        w -= eval_est * evecs_[k];
        RealD res = std::sqrt(norm2(w));
        std::cout << GridLogMessage << "Gamma5BlockLanczos: evec[" << k << "]"
                  << "  eval_reported = " << evals_(k)
                  << "  eval_est = " << eval_est
                  << "  || A v - eval_est * v || = " << res << std::endl;
      }
    }
  }

  // γ5-Gram matrix of block [q1, q2]:  G[i,j] = q^(i)† γ5 q^(j)
  CMat2 gramMatrix(const Field& q1, const Field& q2)
  {
    Field g5q1(Grid_), g5q2(Grid_);
    applyGamma5(q1, g5q1);
    applyGamma5(q2, g5q2);
    CMat2 G;
    G(0,0) = toStdCmplx(innerProduct(q1, g5q1));
    G(0,1) = toStdCmplx(innerProduct(q1, g5q2));
    G(1,0) = toStdCmplx(innerProduct(q2, g5q1));
    G(1,1) = toStdCmplx(innerProduct(q2, g5q2));
    return G;
  }

  // M = Q†γ5 P for Q=[q1,q2], P=[p1,p2]:  M[i,j] = q^(i)† γ5 p^(j)
  CMat2 g5InnerBlock(const Field& q1, const Field& q2,
                     const Field& p1, const Field& p2)
  {
    Field g5p1(Grid_), g5p2(Grid_);
    applyGamma5(p1, g5p1);
    applyGamma5(p2, g5p2);
    CMat2 M;
    M(0,0) = toStdCmplx(innerProduct(q1, g5p1));
    M(0,1) = toStdCmplx(innerProduct(q1, g5p2));
    M(1,0) = toStdCmplx(innerProduct(q2, g5p1));
    M(1,1) = toStdCmplx(innerProduct(q2, g5p2));
    return M;
  }

  CMat2 invert2x2(const CMat2& G)
  {
    return G.inverse();
  }
};

NAMESPACE_END(Grid);

#endif // GRID_GAMMA5_BLOCK_LANCZOS_H
