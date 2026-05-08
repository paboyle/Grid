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
  CMat  Blink_;          // linking block B̃ = B_{m+1} * U[2m-2:2m-1, 0:Nk-1] from krylovSchurCompress
  CMat  GramCompressed_; // Nk×Nk γ5-Gram matrix of compressed Schur basis: G_Ṽ = U_Nk† G_full U_Nk

  // Output
  CVec               evals_;
  std::vector<Field> evecs_;
  std::vector<RealD> residuals_;

public:
  bool doEvalCheck = false;
  bool doVerify = false;

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
  void operator()(const Field& v0, const Field& v1, int maxSteps, int Nstop,
                  bool reorthog = false, RitzFilter filter = EvalImNormSmall)
  {
    basis.clear();
    A_blocks.clear();
    B_blocks.clear();
    C_blocks.clear();
    G_blocks.clear();
    nSteps = 0;

    // Initialise Q_1 = [u0, u1] via L2 Gram-Schmidt on (v0, v1)
    Field u0(Grid_), u1(Grid_);
    u0 = v0;
    RealD nrm = std::sqrt(norm2(u0));
    assert(nrm > 1e-14 && "first starting vector is zero");
    u0 *= (1.0 / nrm);

#if 0
//HACK
    applyGamma5(u0,u1);
#else
    u1 = v1;
    ComplexD proj = innerProduct(u0, u1);
    u1 -= u0 * proj;
    nrm = std::sqrt(norm2(u1));
    assert(nrm > 1e-14 && "second starting vector is linearly dependent on first");
    u1 *= (1.0 / nrm);
#endif

    CMat2 G1 = gramMatrix(u0, u1);
    {
      Eigen::SelfAdjointEigenSolver<CMat2> es(G1);
      auto evals = es.eigenvalues();
      std::cout << GridLogMessage
                << "Gamma5BlockLanczos: G1 = " <<G1 <<std::endl;
      std::cout << GridLogMessage
                << "Gamma5BlockLanczos: eigenvalues = "
                << evals(0) << "  " << evals(1) << std::endl;
      if (std::abs(evals(0)) < 1e-13 || std::abs(evals(1)) < 1e-13) {
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos: abort — degenerate start "
                  << "(G1 has eigenvalue with |λ| < 1e-13)" << std::endl;
        return;
      }
    }

    G_blocks.push_back(G1);
    basis.push_back(u0);
    basis.push_back(u1);

    for (int step = 0; step < maxSteps; step++) {
      bool ok = lanczosStep(step, reorthog);
      if (!ok) break;
      nSteps = step + 1;

      {
        Eigen::ComplexEigenSolver<CMat2> esB(B_blocks[step]);
        auto evB = esB.eigenvalues();
        std::cout << GridLogMessage << "Gamma5BlockLanczos: step " << step
                  << "  B eigenvalues = " << evB(0) << "  " << evB(1);
        if (step < (int)C_blocks.size()) {
          Eigen::ComplexEigenSolver<CMat2> esC(C_blocks[step]);
          auto evC = esC.eigenvalues();
          std::cout << "  C eigenvalues = " << evC(0) << "  " << evC(1);
        }
        std::cout << std::endl;
        if (step < (int)G_blocks.size()) {
          Eigen::SelfAdjointEigenSolver<CMat2> esG(G_blocks[step]);
          auto evG = esG.eigenvalues();
          std::cout << GridLogMessage << "Gamma5BlockLanczos: step " << step
                    << "  G eigenvalues = " << evG(0) << "  " << evG(1) << std::endl;
        }
      }

      RealD beta = B_blocks[step].norm();
      if (beta < Tolerance) {
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos: beta < tol, converged at step "
                  << step << std::endl;
        break;
      }
    }

    if (nSteps == 0) return;
    computeRitzPairs(nSteps, Nstop, filter);
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
  void restart(const Field& v0, const Field& v1, int maxRestarts, int Nstep, int Nk, int Nstop,
               bool reorthog = false, RitzFilter filter = EvalImNormSmall)
  {
    assert(Nk >= Nstop);
    Field src(Grid_), src2(Grid_);
    src  = v0;
    src2 = v1;

    for (int iter = 0; iter < maxRestarts; iter++) {
      std::cout << GridLogMessage
                << "Gamma5BlockLanczos: ---- restart " << iter << " ----" << std::endl;

      // Run Lanczos with two starting vectors; L2 GS is applied inside operator().
      (*this)(src, src2, Nstep, Nstop, reorthog, filter);
      if(this->doVerify){ verify("iter= "+std::to_string(iter)); exit(-42);}

      int nRitz = (int)residuals_.size();
      if (nRitz == 0) {
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos: restart — no Ritz pairs, stopping." << std::endl;
        return;
      }

      // evals_/evecs_/residuals_ are already sorted by filter from computeRitzPairs.
      int nKeep = std::min(Nk, nRitz);
      int nconv = 0;
      for (int i = 0; i < nKeep; i++) {
        if (residuals_[i] < Tolerance) nconv++;
        else break;
      }

      std::cout << GridLogMessage
                << "Gamma5BlockLanczos: restart " << iter
                << "  Ritz = " << nRitz
                << "  wanted = " << nKeep
                << "  converged = " << nconv << " / " << Nstop << std::endl;
      for (int i = 0; i < nKeep; i++)
        std::cout << GridLogMessage
                  << "    wanted[" << i << "]  lambda = " << evals_(i)
                  << "  |res| = " << residuals_[i]
                  << (residuals_[i] < Tolerance ? "  *" : "") << std::endl;

      if (nconv >= Nstop) {
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos: converged after " << iter + 1
                  << " restart(s)." << std::endl;
        return;
      }

      // Build two restart seeds from the Ritz vectors.
      // Split the top min(Nstop,nKeep) Ritz vectors in half:
      //   src  = normalised sum of the first half
      //   src2 = normalised sum of the second half
      // L2 GS inside operator() orthogonalises them.
      int nSeed = std::min(Nstop, nKeep) / 2;
      if (nSeed < 1) nSeed = 1;
      std::cout << GridLogMessage << "Gamma5BlockLanczos: Nstop nKeep nSeed "
                << Nstop << " " << nKeep << " " << nSeed << std::endl;

      src = Zero(); src2 = Zero();
      for (int i = 0; i < nSeed; i++) {
        RealD enorm = std::sqrt(norm2(evecs_[i]));
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos: seed1 evecs_[" << i << "]"
                  << "  ||u||_L2 = " << enorm
                  << "  lambda = " << evals_(i) << std::endl;
        if (enorm > 1e-14) {
          src += evecs_[i] * (1. / enorm);
        }
      }
      for (int i = nSeed; i < std::min(2*nSeed, nKeep); i++) {
        RealD enorm = std::sqrt(norm2(evecs_[i]));
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos: seed2 evecs_[" << i << "]"
                  << "  ||u||_L2 = " << enorm
                  << "  lambda = " << evals_(i) << std::endl;
        if (enorm > 1e-14) {
          src2 += evecs_[i] * (1. / enorm);
        }
      }

      RealD nrm = std::sqrt(norm2(src));
      assert(nrm > 1e-14);
      src *= (1.0 / nrm);

      RealD nrm2 = std::sqrt(norm2(src2));
      if (nrm2 < 1e-14) {
        src2 = evecs_[0];
        nrm2 = std::sqrt(norm2(src2));
      }
      src2 *= (1.0 / nrm2);

      std::cout << GridLogMessage
                << "Gamma5BlockLanczos: seeds built from top " << 2*nSeed
                << " Ritz vectors" << std::endl;
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
  /**
   * Implicitly Restarted Block Lanczos (Krylov-Schur + L2-Arnoldi extension).
   *
   * Initial cycle: run γ5-block Lanczos for Nmax steps → block-tridiagonal T_{Nmax}.
   * Each restart cycle:
   *   1. Schur-compress T (or previous Hessenberg) to Nk modes.
   *   2. L2-QR factorize the Nk compressed field vectors → L2-orthonormal basis W.
   *      Transform Hmat: H_L2 = R·S_Nk·R⁻¹  (preserves eigenvalues).
   *   3. L2-orthogonalize the residual F against W to get a fresh starting vector.
   *   4. Extend with scalar L2-Arnoldi for Np = Nmax-Nk steps from F.
   *      Each step orthogonalises against ALL previous W+extension vectors (always
   *      well-conditioned — no indefinite Gram matrix).
   *   5. Eigensolve the resulting upper-Hessenberg H_comb → Ritz pairs.
   *      Residual estimate: beta_last * |y_j[dim-1]|.
   *
   * This avoids the ill-conditioned G_Ṽ inversion of the split approach while
   * retaining the γ5-block Lanczos efficiency for the initial run.
   */
  void implicitRestart(const Field& v0, const Field& v1, int maxIter, int Nmax, int Nk, int Nstop,
                       bool reorthog = false, RitzFilter filter = EvalImNormSmall)
  {
    assert(Nk >= 2 && Nk < Nmax && Nk >= Nstop);

    // ── Initial full block-Lanczos run ────────────────────────────────────
    (*this)(v0, v1, Nmax, Nstop, reorthog);

    // Persistent state across cycles
    CMat  H_hess;             // current Hessenberg (Ncur × Ncur)
    int   Ncur = 0;           // dimension of H_hess (= Nmax after first fill)
    Field F_vec(Grid_);       // Arnoldi residual vector (seed for next extension)
    bool  first_iter = true;

    for (int iter = 0; iter < maxIter; iter++) {
      std::cout << GridLogMessage
                << "Gamma5BlockLanczos::implicitRestart ---- cycle " << iter << " ----" << std::endl;

      // ── Convergence check ─────────────────────────────────────────────────
      int nRitz = (int)residuals_.size();
      std::vector<int> idx = sortedIdx(nRitz, filter);
      int nKeep = std::min(Nk, nRitz);
      int nconv = 0;
      for (int i = 0; i < nKeep; i++)
        if (residuals_[idx[i]] < Tolerance) nconv++;

      std::cout << GridLogMessage
                << "  nRitz=" << nRitz << "  nconv=" << nconv << "/" << Nstop << std::endl;
      for (int i = 0; i < std::min(nKeep, Nstop + 2); i++)
        std::cout << GridLogMessage
                  << "  [" << i << "]  lambda=" << evals_(idx[i])
                  << "  |res|=" << residuals_[idx[i]]
                  << (residuals_[idx[i]] < Tolerance ? "  *" : "") << std::endl;

      if (nconv >= Nstop) {
        reorderOutput(idx, nKeep);
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos::implicitRestart: converged after "
                  << iter + 1 << " cycle(s)." << std::endl;
        return;
      }

      // ── Krylov-Schur + L2-Arnoldi restart ────────────────────────────────
      if (first_iter) {
        // ── First restart: start from the γ5-block Lanczos result ──────────
        // Schur-compress T_{Nmax} → S_Nk, Ṽ_Nk, Q_{m+1}, Blink_
        krylovSchurCompress(Nk, filter);
        // basis[0..Nk-1] = Ṽ_Nk (NOT L2-orthonormal)
        // basis[Nk..Nk+1] = Q_{m+1}

        // L2-QR factorize Ṽ_Nk to get W (L2-orthonormal) and R (upper triangular).
        // V_old = W · R  →  D_W W = W (R·S_Nk·R⁻¹) + Q_{m+1} (Blink_·R⁻¹)
        CMat R_mat = l2QRFactor(0, Nk);
        CMat R_inv = R_mat.triangularView<Eigen::Upper>().solve(CMat::Identity(Nk, Nk));
        CMat H_L2  = R_mat * Hmat_ * R_inv;  // Nk×Nk; eigenvalues preserved

        // L2-orthogonalize Q_{m+1} columns against W to obtain F_vec
        // (prefer the column with larger residual after projection)
        Field q0 = basis[Nk], q1 = basis[Nk + 1];
        for (int i = 0; i < Nk; i++) {
          q0 -= basis[i] * toStdCmplx(innerProduct(basis[i], q0));
          q1 -= basis[i] * toStdCmplx(innerProduct(basis[i], q1));
        }
        // also L2-orthogonalize q1 against q0 (to break near-parallel)
        RealD n0 = std::sqrt(norm2(q0));
        RealD n1 = std::sqrt(norm2(q1));
        Field F_candidate = (n0 >= n1) ? q0 * (1.0/std::max(n0, 1e-30))
                                       : q1 * (1.0/std::max(n1, 1e-30));
        // re-orthogonalize once more for safety
        for (int i = 0; i < Nk; i++)
          F_candidate -= basis[i] * toStdCmplx(innerProduct(basis[i], F_candidate));
        RealD fn = std::sqrt(norm2(F_candidate));
        assert(fn > 1e-14 && "Q_{m+1} collapsed into Ṽ_Nk — try a different seed");
        F_candidate *= (1.0 / fn);
        F_vec = F_candidate;

        // Set up H_hess: Nmax × Nmax with H_L2 in top-left
        Ncur   = Nmax;
        H_hess = CMat::Zero(Ncur, Ncur);
        H_hess.block(0, 0, Nk, Nk) = H_L2;

        // Trim basis to Nk (Q_{m+1} replaced by F_vec below)
        // Use erase instead of resize: resize instantiates _M_default_append
        // which requires Lattice's default ctor (nonexistent).
        basis.erase(basis.begin() + Nk, basis.end());
        first_iter = false;

      } else {
        // ── Subsequent restarts: Schur-compress previous H_hess ────────────
        ComplexSchurDecomposition schur(H_hess.block(0, 0, Ncur, Ncur), false, filter);
        schur.schurReorder(Nk);
        CMat Qt = schur.getMatrixQ().adjoint();  // Ncur×Ncur unitary

        // Rotate field basis: new_basis[j] = Σ_k basis[k] * Qt(k,j)
        std::vector<Field> new_basis;
        new_basis.reserve(Nk);
        for (int j = 0; j < Nk; j++) {
          Field col(Grid_); col = Zero();
          int Nb = std::min((int)basis.size(), Ncur);
          for (int k = 0; k < Nb; k++)
            col += basis[k] * Qt(k, j);
          new_basis.push_back(col);
        }
        basis = new_basis;

        // Re-orthogonalize F_vec against the new (rotated) basis
        for (int i = 0; i < Nk; i++)
          F_vec -= basis[i] * toStdCmplx(innerProduct(basis[i], F_vec));
        RealD fn = std::sqrt(norm2(F_vec));
        if (fn < 1e-14) {
          std::cout << GridLogMessage
                    << "Gamma5BlockLanczos::implicitRestart: F_vec collapsed; using random." << std::endl;
          // Gram-Schmidt will fix this in the next extension
          fn = 1.0;
        }
        F_vec *= (1.0 / fn);

        // Reset H_hess: place S_Nk in top-left
        Ncur   = Nmax;
        H_hess = CMat::Zero(Ncur, Ncur);
        H_hess.block(0, 0, Nk, Nk) = schur.getMatrixS().block(0, 0, Nk, Nk);

        // Trim basis back to Nk (extension will re-grow it)
        basis.erase(basis.begin() + Nk, basis.end());
      }

      // ── L2-Arnoldi extension: add Np = Nmax-Nk vectors from F_vec ────────
      int Np = Nmax - Nk;
      basis.push_back(F_vec);

      RealD beta_last = 0.0;
      int Nsteps_done = 0;
      for (int step = 0; step < Np; step++) {
        int j = Nk + step;   // index of current vector in basis (0-based)

        Field p(Grid_);
        Linop.Op(basis[j], p);

        // L2-orthogonalize against all previous vectors (Schur + extension)
        for (int i = 0; i < j; i++) {
          ComplexD h = toStdCmplx(innerProduct(basis[i], p));
          H_hess(i, j) = h;
          p -= basis[i] * h;
        }

        beta_last  = std::sqrt(norm2(p));
        Nsteps_done = step + 1;

        if (step < Np - 1) {
          H_hess(j + 1, j) = ComplexD(beta_last, 0.0);
          if (beta_last < Tolerance) {
            std::cout << GridLogMessage
                      << "Gamma5BlockLanczos::implicitRestart: Arnoldi happy breakdown at step "
                      << step << "  beta=" << beta_last << std::endl;
            break;
          }
          basis.push_back(p * (1.0 / beta_last));
        } else {
          // Last step: save normalised residual for next cycle
          if (beta_last > 1e-14) F_vec = p * (1.0 / beta_last);
        }
      }

      int Ncur_used = Nk + Nsteps_done;
      std::cout << GridLogMessage
                << "Gamma5BlockLanczos::implicitRestart: Arnoldi extended to dim="
                << Ncur_used << "  beta_last=" << beta_last << std::endl;

      // ── Ritz pairs from H_hess ────────────────────────────────────────────
      computeRitzPairsHessenberg(H_hess, Ncur_used, beta_last, filter, Nstop);
    }

    // maxIter exhausted
    std::cout << GridLogMessage
              << "Gamma5BlockLanczos::implicitRestart: maxIter=" << maxIter
              << " reached without full convergence." << std::endl;
    {
      int nRitz = (int)residuals_.size();
      if (nRitz > 0) {
        std::vector<int> idx = sortedIdx(nRitz, filter);
        reorderOutput(idx, std::min(Nk, nRitz));
      }
    }
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

    // ---- Compare T_m (from A/B/C blocks) against directly constructed V†γ5DV ----
    // Apply D_W to every basis vector.
    std::vector<Field> Dv(dim, Field(Grid_));
    for (int j = 0; j < dim; j++)
      Linop.Op(basis[j], Dv[j]);

    // T_proj[i,j] = basis[i]† γ5 D_W basis[j]  (γ5-inner product with D_W ket)
    // G_full [i,j] = basis[i]† γ5 basis[j]
    // Relation: T_proj = G_full * Tm  (exact if γ5-orthogonality holds)
    CMat T_proj  = CMat::Zero(dim, dim);
    CMat G_full  = CMat::Zero(dim, dim);
    for (int ib = 0; ib < m; ib++) {
      for (int kb = 0; kb < m; kb++) {
        CMat2 tp = g5InnerBlock(basis[2*ib], basis[2*ib+1], Dv[2*kb],    Dv[2*kb+1]);
        CMat2 gf = g5InnerBlock(basis[2*ib], basis[2*ib+1], basis[2*kb], basis[2*kb+1]);
        T_proj.block(2*ib, 2*kb, 2, 2) = tp;
        G_full.block(2*ib, 2*kb, 2, 2) = gf;
      }
    }

    CMat Terr     = T_proj - G_full * Tm;
    RealD maxTerr = Terr.cwiseAbs().maxCoeff();
    RealD maxTnrm = (G_full * Tm).cwiseAbs().maxCoeff();
    std::cout << GridLogMessage
              << "  max |T_proj - G_full*T_blocks| = " << maxTerr
              << "  (rel: " << (maxTnrm > 0 ? maxTerr/maxTnrm : 0.0) << ")" << std::endl;

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

    // Compute GramCompressed_ = U_Nk† G_full U_Nk from the ORIGINAL G_blocks
    // (before they are cleared below).  G_full = block-diag(G_0,...,G_{m-1}).
    // Since each G_k = diag(±1) and U is unitary, G_Ṽ² = I → G_Ṽ^{-1} = G_Ṽ.
    {
      CMat G_full_mat = CMat::Zero(dim, dim);
      for (int k = 0; k < m; k++)
        G_full_mat.block(2*k, 2*k, 2, 2) = G_blocks[k];   // G_blocks still has ORIGINAL values here
      CMat U_Nk = U.leftCols(Nk);
      GramCompressed_ = U_Nk.adjoint() * G_full_mat * U_Nk;
    }

    basis = new_basis;   // Nk + 2 field vectors

    // Recompute G_blocks for each pair of new basis vectors
    G_blocks.clear();
    for (int k = 0; k <= Nk / 2; k++)
      G_blocks.push_back(gramMatrix(basis[2*k], basis[2*k+1]));

    // Compressed projected matrix = leading Nk×Nk block of S
    Hmat_ = S.block(0, 0, Nk, Nk);

    // Linking block: D_W Ṽ_{Nk} = Ṽ_{Nk} S_{Nk} + Q_{m+1} B̃_link
    // where B̃_link = B_{m+1} * (last 2 rows of U for first Nk cols).
    // Needed to populate the (Nk+1)-th block-row of Hmat_ on the first extension step.
    Blink_ = B_blocks[m - 1] * U.bottomRows(2).leftCols(Nk);

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
    // On the first extension step after krylovSchurCompress, populate the linking
    // block row: D_W Ṽ_{Nk} = Ṽ_{Nk} S_{Nk} + Q_{m+1} B̃_link  (Nk = Ncompressed_)
    if (useFullH_ && old_dim == Ncompressed_ && Blink_.cols() == old_dim)
      Hmat_new.block(old_dim, 0, 2, old_dim) = Blink_;

    // Coupling to all previous blocks.
    // After krylovSchurCompress the first Ncompressed_ basis vectors are the
    // compressed Schur vectors.  They are NOT mutually γ5-orthogonal, so we
    // must use the full Gram system G_Ṽ (stored in GramCompressed_) for their
    // block.  Since G_Ṽ² = I (G_Ṽ is an involution), G_Ṽ^{-1} = G_Ṽ exactly,
    // so the system is solved by a matrix-vector product (no ill-conditioning).
    //
    // Extension vectors (j ≥ Ncompressed_/2) are built by this same routine
    // with full γ5-projection against all predecessors, so they ARE mutually
    // γ5-orthogonal → the diagonal-block formula suffices for them.
    Field r1(Grid_), r2(Grid_);
    r1 = p1; r2 = p2;

    int Nc2 = (useFullH_ ? Ncompressed_ / 2 : 0);  // number of compressed blocks

    if (Nc2 > 0 && step >= Nc2) {
      // Compressed blocks: collect coupling, solve with full Gram G_Ṽ.
      CMat Mcol = CMat::Zero(Ncompressed_, 2);
      for (int j = 0; j < Nc2; j++) {
        CMat2 Mj = g5InnerBlock(basis[2*j], basis[2*j+1], p1, p2);
        Mcol.block(2*j, 0, 2, 2) = Mj;
      }
      // G_Ṽ^{-1} = G_Ṽ (since G_Ṽ² = I); use LU for numerical safety
      CMat Hcol = GramCompressed_.lu().solve(Mcol);
      Hmat_new.block(0, old_dim, Ncompressed_, 2) = Hcol;
      for (int k = 0; k < Ncompressed_; k++) {
        r1 -= basis[k] * Hcol(k, 0);
        r2 -= basis[k] * Hcol(k, 1);
      }
    }

    // Extension blocks (j ≥ Nc2): mutually γ5-orthogonal, diagonal formula.
    for (int j = Nc2; j < step; j++) {
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

    // Breakdown check per Eq (53)/(54) of Yamamoto 2026.
    // Eq (53) happy breakdown: rjnrm < Tolerance  →  Krylov space invariant; stop step.
    // Eq (54) serious breakdown: |d_j| << rjnrm²  →  neutral vector; stop step.
    // Both cases return false so the caller terminates the Lanczos loop cleanly.
    // Relative threshold 1e-14 (≈ machine ε) avoids spurious triggers.
    const RealD breakdownEps = 1e-14;
    for (int j = 0; j < 2; j++) {
      Field rj = r1 * U2(0,j) + r2 * U2(1,j);
      RealD rjnrm2 = norm2(rj);
      RealD rjnrm  = std::sqrt(rjnrm2);
      if (rjnrm < Tolerance) {
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos: happy breakdown (full step " << step
                  << " direction " << j
                  << ")  ||R̂_k u_j||=" << rjnrm << std::endl;
        return false;
      } else if (std::abs(D(j)) < breakdownEps * rjnrm2) {
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos: SERIOUS breakdown (full step " << step
                  << " direction " << j
                  << ")  ||R̂_k u_j||=" << rjnrm
                  << "  d_j=" << D(j)
                  << "  |d_j|/||r_j||^2=" << std::abs(D(j)) / rjnrm2
                  << "  (look-ahead not implemented; stopping)" << std::endl;
        return false;
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
    {
    CMat2 gamma1 = gramMatrix(r1, r2);
    Eigen::ComplexEigenSolver<CMat2> esB(gamma1);
    auto evB = esB.eigenvalues();
    std::cout << GridLogMessage << "Gamma5BlockLanczos: step " << step
    << "  gamma eigenvalues = " << evB(0) << "  " << evB(1);
    }

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
      // Guard against NaN from degenerate eigenvectors in the non-symmetric eigensolver
      if (!std::isfinite(res)) res = std::numeric_limits<RealD>::infinity();
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

  // Ritz pairs from the block-lower-triangular combined matrix assembled by
  // implicitRestart (split Krylov-Schur strategy).
  //
  // Hmat_comb  = [ S_{Nk}  |  0        ]   (dim_total × dim_total)
  //              [ Blink_  |  T_fresh   ]
  //
  // basis[] contains:
  //   [0..Nk-1]          → compressed Schur vectors Ṽ_{Nk}
  //   [Nk..Nk+2*Nf-1]    → fresh Lanczos blocks (Q_{m+1} through Q_{m+Nf})
  //   [Nk+2*Nf..Nk+2*Nf+1] → outer residual Q_{m+Nf+1} (not in Hmat_comb)
  //
  // Ritz residual for eigenpair (λ_j, y_j) of Hmat_comb:
  //   ||D_W u_j - λ_j u_j|| ≈ ||Blink_save * y_j.head(Nk)|| + ||Blast_ * τ||
  //   where τ = y_j.tail(2)  (last two entries, last fresh block contribution)
  void computeRitzPairsCombined(const CMat& Hmat_comb, const CMat& Blink_save,
                                 int Nk, int Nstop)
  {
    int dim_total = Hmat_comb.rows();
    int dim_fresh = dim_total - Nk;

    Eigen::ComplexEigenSolver<CMat> ces(Hmat_comb);
    CVec lambdas = ces.eigenvalues();
    CMat Y       = ces.eigenvectors();

    // Sort by |Im(λ)| ascending (near-real physical modes first)
    std::vector<int> idx(dim_total);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b){
      return std::abs(lambdas(a).imag()) < std::abs(lambdas(b).imag());
    });

    evals_.resize(dim_total);
    evecs_.clear();
    residuals_.clear();

    for (int ji = 0; ji < dim_total; ji++) {
      int  j  = idx[ji];
      evals_(ji) = lambdas(j);
      CVec yj = Y.col(j);

      // Ritz vector: u_j = sum_{k=0}^{dim_total-1} basis[k] * y_j(k)
      // basis[0..Nk-1]       = Schur vectors
      // basis[Nk..dim_total-1] = fresh Lanczos vectors (Q_{m+1}..Q_{m+Nf})
      Field uj(Grid_);
      uj = Zero();
      for (int k = 0; k < dim_total; k++)
        uj += basis[k] * yj(k);
      evecs_.push_back(uj);

      // Ritz residual: two contributions
      //   1) Schur leak:  Q_{m+1} Blink_ y_schur  (linking row)
      //   2) Fresh tail:  Q_{m+Nf+1} Blast_ τ     (fresh outer residual)
      CVec y_schur = yj.head(Nk);
      Eigen::Vector2cd Blink_y = Blink_save * y_schur;
      RealD res_schur = Blink_y.norm();

      Eigen::Vector2cd tau(yj(dim_total - 2), yj(dim_total - 1));
      RealD res_fresh = (Blast_ * tau).norm();

      RealD res = res_schur + res_fresh;
      if (!std::isfinite(res)) res = std::numeric_limits<RealD>::infinity();
      residuals_.push_back(res);

      std::cout << GridLogMessage
                << "Gamma5BlockLanczos (combined): Ritz[" << ji << "]"
                << "  lambda=" << evals_(ji)
                << "  |res_schur|=" << res_schur
                << "  |res_fresh|=" << res_fresh << std::endl;
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

  // L2 modified Gram-Schmidt: orthonormalizes basis[start..start+count-1] in L2.
  // Returns upper-triangular R (count×count) such that V_old = W·R.
  CMat l2QRFactor(int start, int count)
  {
    CMat R = CMat::Zero(count, count);
    for (int j = 0; j < count; j++) {
      for (int i = 0; i < j; i++) {
        ComplexD h = toStdCmplx(innerProduct(basis[start + i], basis[start + j]));
        R(i, j) = h;
        basis[start + j] -= basis[start + i] * h;
      }
      RealD nrm = std::sqrt(norm2(basis[start + j]));
      R(j, j) = ComplexD(nrm, 0.0);
      if (nrm > 1e-14) basis[start + j] *= (1.0 / nrm);
    }
    return R;
  }

  // Ritz pairs from the upper-Hessenberg H[0:dim,0:dim] built by L2-Arnoldi.
  // Residual estimate: beta_last * |y_j[dim-1]|  (standard Arnoldi formula).
  void computeRitzPairsHessenberg(const CMat& H, int dim, RealD beta_last,
                                   RitzFilter filter, int Nstop)
  {
    Eigen::ComplexEigenSolver<CMat> ces(H.block(0, 0, dim, dim));
    CVec lambdas = ces.eigenvalues();
    CMat Y       = ces.eigenvectors();

    ComplexComparator cComp(filter);
    std::vector<int> idx(dim);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b){
      return cComp(toStdCmplx(lambdas(a)), toStdCmplx(lambdas(b)));
    });

    evals_.resize(dim);
    evecs_.clear();
    residuals_.clear();

    for (int ji = 0; ji < dim; ji++) {
      int  j  = idx[ji];
      evals_(ji) = lambdas(j);
      CVec yj = Y.col(j);

      Field uj(Grid_);
      uj = Zero();
      for (int k = 0; k < dim && k < (int)basis.size(); k++)
        uj += basis[k] * yj(k);
      evecs_.push_back(uj);

      RealD res = beta_last * std::abs(yj(dim - 1));
      if (!std::isfinite(res)) res = std::numeric_limits<RealD>::infinity();
      residuals_.push_back(res);

      std::cout << GridLogMessage
                << "Gamma5BlockLanczos (Hess): Ritz[" << ji << "]"
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
        RealD res_check = std::sqrt(norm2(w));
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos: evec[" << k << "]"
                  << "  eval_reported=" << evals_(k)
                  << "  eval_est=" << eval_est
                  << "  ||Av-eval*v||=" << res_check << std::endl;
      }
    }
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

    ComplexD detG1 = Ck(0,0)*Ck(1,1) - Ck(0,1)*Ck(1,0);
    std::cout << GridLogMessage
              << "Gamma5BlockLanczos: C "<<step<<"= \n" << Ck 
              << "\n  det = " << detG1 << std::endl;
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

    // Breakdown check per Eq (53)/(54) of Yamamoto 2026.
    // Eq (53) happy breakdown: rjnrm < Tolerance  →  Krylov space invariant; stop step.
    // Eq (54) serious breakdown: |d_j| << rjnrm²  →  neutral vector; stop step.
    // Both cases return false so the outer loop terminates the Lanczos run cleanly.
    // Relative threshold 1e-14 (≈ machine ε) avoids spurious triggers.
    const RealD breakdownEps = 1e-14;
    for (int j = 0; j < 2; j++) {
      Field rj(Grid_);
      rj = r1 * U(0,j) + r2 * U(1,j);
      RealD rjnrm2 = norm2(rj);
      RealD rjnrm  = std::sqrt(rjnrm2);
      if (rjnrm < Tolerance) {
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos: happy breakdown at step " << step
                  << " direction " << j
                  << "  ||R̂_k u_j||=" << rjnrm << std::endl;
        return false;
      } else if (std::abs(D(j)) < breakdownEps * rjnrm2) {
        std::cout << GridLogMessage
                  << "Gamma5BlockLanczos: SERIOUS breakdown at step " << step
                  << " direction " << j
                  << "  ||R̂_k u_j||=" << rjnrm
                  << "  d_j=" << D(j)
                  << "  |d_j|/||r_j||^2=" << std::abs(D(j)) / rjnrm2
                  << "  (look-ahead not implemented; stopping)" << std::endl;
        return false;
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

    {
    CMat2 gamma1 = gramMatrix(r1, r2);
    Eigen::ComplexEigenSolver<CMat2> esB(gamma1);
    auto evB = esB.eigenvalues();
    std::cout << GridLogMessage
              << "Gamma5BlockLanczos: gamma "<<step<<"= \n" << gamma1 <<std::endl
    << "  gamma eigenvalues = " << evB(0) << "  " << evB(1)<<std::endl;
    }

    detG1 = Gkp1(0,0)*Gkp1(1,1) - Gkp1(0,1)*Gkp1(1,0);
    std::cout << GridLogMessage
              << "Gamma5BlockLanczos: G "<<step<<"= \n" << Gkp1
              << "\n  det G = " << detG1 << std::endl;

    G_blocks.push_back(Gkp1);
    detG1 = Bkp1(0,0)*Bkp1(1,1) - Bkp1(0,1)*Bkp1(1,0);
    std::cout << GridLogMessage
              << "Gamma5BlockLanczos: B "<<step<<"= \n" << Bkp1 
              << "\n  det B = " << detG1 << std::endl;
    B_blocks.push_back(Bkp1);

    CMat2 G1 = gramMatrix(basis[0], basis[1],qnew1,qnew2);
    detG1 = G1(0,0)*G1(1,1) - G1(0,1)*G1(1,0);
    RealD absdetG1 = std::sqrt(detG1.real()*detG1.real() + detG1.imag()*detG1.imag());
    std::cout << GridLogMessage
              << "Gamma5BlockLanczos: eta "<<step<<" = \n" << G1
              << "\n  det = " << detG1 << std::endl;
    basis.push_back(qnew1);
    basis.push_back(qnew2);

    return true;
  }

  // Assemble T_m and extract Ritz pairs (Algorithm 2 of the paper).
  void computeRitzPairs(int m, int Nstop, RitzFilter filter = EvalImNormSmall)
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

    // Sort by filter criterion (ComplexComparator applies the same penalty as schurReorder)
    ComplexComparator cComp(filter);
    std::vector<int> idx(dim);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b){
      return cComp(toStdCmplx(lambdas(a)), toStdCmplx(lambdas(b)));
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

  // γ5-Gram matrix of block [q1, q2]:  G[i,j] = q^(i)† γ5 q^(j)
  CMat2 gramMatrix(const Field& q11, const Field& q12, const Field& q21, const Field& q22)
  {
    Field g5q21(Grid_), g5q22(Grid_);
    applyGamma5(q21, g5q21);
    applyGamma5(q22, g5q22);
    CMat2 G;
    G(0,0) = toStdCmplx(innerProduct(q11, g5q21));
    G(0,1) = toStdCmplx(innerProduct(q11, g5q22));
    G(1,0) = toStdCmplx(innerProduct(q12, g5q21));
    G(1,1) = toStdCmplx(innerProduct(q12, g5q22));
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
