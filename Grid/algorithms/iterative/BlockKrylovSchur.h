/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/BlockKrylovSchur.h

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
#ifndef GRID_BLOCKED_KRYLOV_SCHUR_H
#define GRID_BLOCKED_KRYLOV_SCHUR_H

#include <iomanip>
#include <numeric>

NAMESPACE_BEGIN(Grid);

/**
 * Block (block-Arnoldi) restarted Krylov-Schur eigensolver for
 * general non-Hermitian operators.
 *
 * Algorithm
 * ---------
 * Uses a block Arnoldi factorisation of block size Nblock:
 *
 *   A V_k = V_k H_k + F_k B_k^dag
 *
 * where
 *   V_k   = Nm orthonormal basis vectors (stored flat in basis[])
 *   H_k   = Nm x Nm upper block-Hessenberg Rayleigh quotient
 *   F_k   = Nblock residual vectors (the next block beyond V_k)
 *   B_k   = Nm x Nblock coupling matrix (non-zero only in last Nblock rows)
 *
 * Each block Arnoldi step applies A to each of the Nblock vectors in the
 * current block, orthogonalises against all previous basis vectors, and
 * reduces the residual block to upper-triangular form via Householder QR
 * (implemented here as modified Gram-Schmidt within the block).
 *
 * The restart is a thick restart via the Schur decomposition of H_k:
 *   H_k = Q^dag S Q
 * The leading Nk Schur vectors (chosen by RitzFilter) are retained,
 * the basis and Rayleigh quotient are truncated, and block Arnoldi continues
 * from the (Nk/Nblock)-th block.
 *
 * Parameters
 * ----------
 * Nblock  : block size p
 * Nm      : total Krylov dimension (must be divisible by Nblock)
 * Nk      : total vectors to keep after each restart (must be divisible by Nblock, Nk < Nm)
 * Nstop   : declare convergence when this many eigenpairs have converged
 * MaxIter : maximum number of outer (restart) iterations
 * Tolerance : relative convergence tolerance (||r|| < Tolerance * |lambda_max|)
 */
template<class Field>
class BlockKrylovSchur {

protected:
  //--------------------------------------------------------------------
  // Types
  //--------------------------------------------------------------------
  typedef Eigen::MatrixXcd CMat;
  typedef Eigen::VectorXcd CVec;

  //--------------------------------------------------------------------
  // Parameters (set by operator())
  //--------------------------------------------------------------------
  int   Nblock;    // block size
  int   Nm;        // total Krylov dimension (multiple of Nblock)
  int   Nk;        // total vectors retained after restart (multiple of Nblock)
  int   Nstop;
  int   MaxIter;
  RealD Tolerance;

  //--------------------------------------------------------------------
  // Internal state
  //--------------------------------------------------------------------
  LinearOperatorBase<Field>& Linop;
  GridBase*                  Grid_;
  RitzFilter                 ritzFilter;

  // Prefix used in log messages; derived classes override in their constructor
  // (e.g. "HarmonicBlockKrylovSchur") so the shared code reports the right name.
  std::string className = "BlockKrylovSchur";

  // Flat storage: basis[s*Nblock + t] is the t-th vector of block s
  // After construction: basis has Nm entries
  std::vector<Field> basis;

  // Rayleigh quotient  Nm x Nm
  CMat H;

  // Residual block: Nblock vectors (the (Nm/Nblock+1)-th block, unnormalised before
  // QR; normalised and orthogonalised as part of block Arnoldi)
  std::vector<Field> F;

  // Coupling matrix B: Nm x Nblock.
  // In exact arithmetic only the last Nblock rows are non-zero:
  //   B(Nm - Nblock + t, s) = H_{Nm/Nblock+1, Nm/Nblock}(t, s)   (the subdiagonal block)
  // We keep it as a full matrix for generality after restarts.
  CMat B;

  RealD beta_k;          // Frobenius norm of the last subdiagonal block
  RealD rtol;            // absolute tolerance = Tolerance * approxLambdaMax

  // Output
  CVec               evals;
  CMat               littleEvecs;   // Nm columns
  std::vector<RealD> ritzEstimates;

public:
  std::vector<Field> evecs;
  bool doEvalCheck  = false;
  // When true (and Nblock even), only Nblock/2 starting vectors are required.
  // The remaining Nblock/2 slots are filled by pairing each supplied vector
  // with its parity-flipped partner (sign negated on odd checkerboard sites).
  bool useParityFlip = false;
  // Same pairing behaviour but the partner is produced by gamma5Func(v).
  // gamma5Func must be set by the caller before operator() is called.
  // It is a std::function so that BlockKrylovSchur.h does not need to know
  // about the Grid QCD Gamma class (which is defined at a later include layer).
  bool useGamma5     = false;
  std::function<void(const Field&, Field&)> gamma5Func;

  //--------------------------------------------------------------------
  // Constructor
  //--------------------------------------------------------------------
  BlockKrylovSchur(LinearOperatorBase<Field>& _Linop, GridBase* _Grid,
                     RealD _Tolerance, RitzFilter _rf = EvalReSmall)
    : Linop(_Linop), Grid_(_Grid), Tolerance(_Tolerance), ritzFilter(_rf),
      Nblock(-1), Nm(-1), Nk(-1), Nstop(-1), MaxIter(-1),
      beta_k(0.0), rtol(0.0)
  {}

  virtual ~BlockKrylovSchur() = default;

  //--------------------------------------------------------------------
  // Main entry point
  //--------------------------------------------------------------------
  /**
   * Run the blocked Krylov-Schur algorithm.
   *
   * Parameters
   * ----------
   * v0       : block of Nblock starting vectors (size >= Nblock)
   * _maxIter : maximum outer (restart) iterations
   * _Nm      : total Krylov dimension (must be divisible by _Nblock)
   * _Nk      : total vectors to keep after restart (must be divisible by _Nblock, _Nk < _Nm)
   * _Nstop   : stop after _Nstop eigenvalues converged
   * _Nblock  : block size
   */
  virtual void operator()(const std::vector<Field>& v0, int _maxIter, int _Nm, int _Nk,
                          int _Nstop, int _Nblock = 1, bool doubleOrthog = true,
                          bool doVerify = false)
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
      std::cout << GridLogMessage << className << ": divisor = " << divisor << std::endl;
    }
    assert(Nm % Nblock == 0);
    assert(Nk % Nblock == 0);
    assert(Nk < Nm);
    if (useGamma5) assert(gamma5Func && "useGamma5: gamma5Func must be set");
    preRun();  // hook: derived classes add parameter assertions here

    int N = Nm;   // total Krylov dimension

    // Approximate largest eigenvalue for tolerance normalisation
    RealD approxLambdaMax = approxMaxEval(v0[0]);
    rtol = Tolerance * approxLambdaMax;
    std::cout << GridLogMessage << className << ": approx max eval = "
              << approxLambdaMax << ", rtol = " << rtol << std::endl;

    // Initialise
    H = CMat::Zero(N, N);
    B = CMat::Zero(N, Nblock);

    int start = 0;
    std::vector<Field> startBlock = expandStartBlock(v0);

    for (int iter = 0; iter < MaxIter; iter++) {
      std::cout << GridLogMessage << className << ": restart iteration " << iter << std::endl;

      // ---- Block Arnoldi: extend from block start to block Nm/Nblock ----
      blockArnoldiIteration(startBlock, Nm/Nblock, start, doubleOrthog);

      // After first full cycle start from block Nk/Nblock
      start = Nk/Nblock;

      if (doVerify) {
        std::string lbl = "iter " + std::to_string(iter) + " after Arnoldi";
        verify(lbl);
      }

      // ---- Restart rotation: compute Q and the rotated+reordered Rayleigh
      //      quotient Hnew.  Base class uses the plain Schur decomposition of
      //      H; derived classes (e.g. HarmonicBlockKrylovSchur) override this
      //      to target a shift.
      CMat Q, Hnew;
      restartRotation(Q, Hnew);
      CMat Qt = Q.adjoint();

      // Rotate Krylov basis: basis_new[i] = sum_j basis[j] * Qt(j,i)
      std::vector<Field> basis2;
      constructUR(basis2, basis, Qt, N);
      basis = basis2;

      // Update B and H
      B = Q * B;
      H = Hnew;

      // ---- Truncate to Nk ----
      int Nkeep = Nk;

      CMat Htmp = H(Eigen::seqN(0, Nkeep), Eigen::seqN(0, Nkeep));
      H = CMat::Zero(N, N);
      H(Eigen::seqN(0, Nkeep), Eigen::seqN(0, Nkeep)) = Htmp;

      std::vector<Field> basisTmp(basis.begin(), basis.begin() + Nkeep);
      basis = basisTmp;

      CMat Btmp = B(Eigen::seqN(0, Nkeep), Eigen::all);
      B = CMat::Zero(N, Nblock);
      B(Eigen::seqN(0, Nkeep), Eigen::all) = Btmp;

      // beta_k = Frobenius norm of the effective coupling
      beta_k = Btmp.norm();
      std::cout << GridLogMessage << className << ": beta_k = " << beta_k << std::endl;

      // Restart: the new starting block is F (the residual block from Arnoldi)
      startBlock = F;

      if (doVerify) {
        std::string lbl = "iter " + std::to_string(iter) + " after restart+truncation";
        verify(lbl);
      }

      // ---- Compute eigensystem of truncated H for convergence check ----
      CMat Hk = H(Eigen::seqN(0, Nkeep), Eigen::seqN(0, Nkeep));
      computeEigensystem(Hk, Nkeep);

      int Nconv = converged(Nkeep);
      std::cout << GridLogMessage << className << ": converged " << Nconv
                << " / " << Nstop << std::endl;

      if (Nconv >= Nstop || iter == MaxIter - 1) {
        std::cout << GridLogMessage << className << ": done after " << iter
                  << " restarts, " << Nconv << " converged." << std::endl;
        std::cout << GridLogMessage << "Eigenvalues: " << evals << std::endl;

        if (doEvalCheck) {
          Field w(Grid_);
          for (int k = 0; k < (int)evecs.size(); k++) {
            Linop.Op(evecs[k], w);
            ComplexD eval_est = toStdCmplx(innerProduct(evecs[k], w));
            w -= eval_est * evecs[k];
            RealD res = std::sqrt(norm2(w));
            std::cout << GridLogMessage << className << ": evec[" << k << "]"
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

  // Negate the field on all odd-parity (checkerboard) sites.
  // Uses pickCheckerboard/setCheckerboard so the SIMD layout is handled
  // correctly and the negation runs in parallel.
  static void parityFlippedField(const Field& v, Field& out) {
    GridCartesian* cgrid = dynamic_cast<GridCartesian*>(v.Grid());
    assert(cgrid != nullptr);
    GridRedBlackCartesian rbgrid(cgrid);
    Field veven(&rbgrid), vodd(&rbgrid);
    pickCheckerboard(Even, veven, v);
    pickCheckerboard(Odd,  vodd,  v);
    vodd = -vodd;
    setCheckerboard(out, veven);
    setCheckerboard(out, vodd);
  }

  //--------------------------------------------------------------------
  // Verification: print H and B, check A V = V H + F B^dag explicitly
  //--------------------------------------------------------------------
  /**
   * Checks the block Arnoldi / Krylov-Schur decomposition
   *
   *   A V = V H + F B^dag                                      (KS)
   *
   * by explicit operator applications.  For each basis vector j:
   *
   *   w_j = A basis[j]
   *
   * The nBasis × nBasis matrix M of inner products is computed:
   *
   *   M[i, j] = <basis[i] | A basis[j]>
   *
   * and compared column-by-column against H.  Separately, the nBasis × Nblock
   * residual coupling matrix R is computed:
   *
   *   R[j, t] = <basis[j] | F[t]> * ||F[t]||   (scaled by F-block norms)
   *
   * but since F is already normalised, R[j,t] = <basis[j] | F[t]>.
   *
   * The KS relation for column j reads:
   *   w_j = sum_i basis[i] H[i,j]  +  sum_t F[t] B[j,t]*
   * so the deviation in column j is
   *   dev_j = w_j - sum_i basis[i] M[i,j]   (should be zero for exact arithmetic)
   * augmented by the F B^dag term in the last block.
   *
   * Prints:
   *   - H (current Rayleigh quotient, nBasis × nBasis)
   *   - B (coupling matrix, nBasis × Nblock)
   *   - M (explicit inner product matrix <V | A V>)
   *   - max |H[i,j] - M[i,j]|  (should be O(machine epsilon))
   *   - for each basis column j: || A v_j - V H[:,j] - F B[j,:]^* ||
   *
   * Parameters
   * ----------
   * label : string printed at the start (e.g. "after restart 2")
   */
  void verify(const std::string& label = "")
  {
    int nBasis = (int)basis.size();
    int nF     = (int)F.size();

    if (nBasis == 0) {
      std::cout << GridLogMessage
                << className << "::verify [" << label
                << "]: basis is empty." << std::endl;
      return;
    }

    // Full allocated dimension of H (= Nm); may be larger than nBasis after restart+truncation
    int Nfull = (int)H.rows();

    std::cout << GridLogMessage
              << "======== " << className << "::verify [" << label << "] ========" << std::endl;
    std::cout << GridLogMessage
              << "  nBasis = " << nBasis << "  Nfull = " << Nfull
              << "  Nblock = " << Nblock << "  nF = " << nF << std::endl;

    // ---- Print H (full Nm x Nm matrix) ----
    // After restart+truncation nBasis < Nfull; entries outside [0:nBasis,0:nBasis]
    // should be ~0 — printing the full matrix makes this visible.
    std::cout << GridLogMessage << "H (" << Nfull << " x " << Nfull << "):" << std::endl;
    for (int i = 0; i < Nfull; i++) {
      for (int j = 0; j < Nfull; j++)
        std::cout << "  " << std::setprecision(4) << std::setw(14) << H(i, j);
      std::cout << std::endl;
    }

    // ---- Print B (full Nm x Nblock matrix) ----
    std::cout << GridLogMessage << "B (" << Nfull << " x " << nF << "):" << std::endl;
    for (int i = 0; i < Nfull; i++) {
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

    // ---- max |H - M| (for the nBasis x nBasis Rayleigh-quotient block) ----
    RealD maxHM = 0.0;
    for (int i = 0; i < nBasis; i++)
      for (int j = 0; j < nBasis; j++)
        maxHM = std::max(maxHM, std::abs(H(i,j) - M(i,j)));
    std::cout << GridLogMessage
              << "  max |H[i,j] - M[i,j]| (i,j < nBasis) = " << maxHM << std::endl;

    // ---- Check entries of H OUTSIDE the nBasis x nBasis block ----
    // These live in rows i >= nBasis or columns j >= nBasis.
    // They must be ~0 after restart+truncation (before the coupling rows are
    // re-filled by the next blockArnoldiIteration call).
    // During a full Arnoldi run (nBasis == Nfull) this loop is a no-op.
    RealD maxOutside = 0.0;
    for (int i = 0; i < Nfull; i++)
      for (int j = 0; j < Nfull; j++)
        if (i >= nBasis || j >= nBasis)
          maxOutside = std::max(maxOutside, std::abs(H(i,j)));
    if (Nfull > nBasis)
      std::cout << GridLogMessage
                << "  max |H[i,j]| outside [0:" << nBasis << ",0:" << nBasis
                << "] block = " << maxOutside
                << " (should be ~0 after restart+truncation)" << std::endl;

    // ---- Check orthonormality of basis ----
    CMat G = CMat::Zero(nBasis, nBasis);
    for (int i = 0; i < nBasis; i++)
      for (int j = 0; j < nBasis; j++)
        G(i, j) = toStdCmplx(innerProduct(basis[i], basis[j]));
    CMat Gerr = G - CMat::Identity(nBasis, nBasis);
    std::cout << GridLogMessage
              << "  max |<V_i|V_j> - delta_ij| = " << Gerr.cwiseAbs().maxCoeff() << std::endl;

    // ---- Per-column residual: || A v_j - V H[:,j] - F B[j,:]^* || ----
    // For each basis vector j, compute A v_j then subtract V H[:,j] and F B[j,:]^*
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

  //--------------------------------------------------------------------
  // Extension hooks for derived classes
  //--------------------------------------------------------------------

  /**
   * Apply the operator (or poly-filtered operator) to the current Arnoldi
   * block and return the results in W.
   *
   * Base implementation: W[t] = Linop.Op(basis[kBase + t]).
   * Derived classes (e.g. SplitGridBlockKrylovSchur) override this to
   * apply poly(A) via Grid_split / Grid_unsplit instead.
   */
  virtual void applyBlock(std::vector<Field>& W, int kBase)
  {
    for (int t = 0; t < Nblock; t++)
      Linop.Op(basis[kBase + t], W[t]);
  }

  /**
   * Called once at the start of operator(), after parameters are set.
   * Base implementation is a no-op; derived classes override to add
   * parameter assertions (e.g. Nblock % mrhs == 0).
   */
  virtual void preRun() {}

  /**
   * Compute the restart rotation for the current Rayleigh quotient H.
   *
   * Returns Q (unitary) and Hnew = Q H Q^dag reordered so that the wanted
   * Nk eigenvalues occupy the leading block, ready for truncation.
   *
   * Base implementation: the plain Schur decomposition of H, reordered by
   * ritzFilter.  HarmonicBlockKrylovSchur overrides this to Schur-decompose
   * the shifted quotient (H - sigma I) so that eigenvalues nearest sigma are
   * selected.
   */
  virtual void restartRotation(CMat& Q, CMat& Hnew)
  {
    ComplexSchurDecomposition schur(H, false, ritzFilter);
    std::cout << GridLogMessage << className << ": Schur decomposed." << std::endl;
    schur.schurReorder(Nk);
    std::cout << GridLogMessage << className << ": Schur reordered." << std::endl;
    Q    = schur.getMatrixQ();
    Hnew = schur.getMatrixS();
  }

protected:

  //--------------------------------------------------------------------
  // Starting-block expansion (parity flip / gamma5 partners)
  //--------------------------------------------------------------------
  /**
   * Expand Nblock/divisor seed vectors into the full Nblock starting block by
   * pairing each seed with its parity-flipped and/or gamma5 partners.  Order:
   * {v, flip(v)} then doubled by {., g5(.)} giving {v, flip(v), g5(v),
   * g5(flip(v))} when both flags are set.  Requires Nblock, Grid_ and the
   * flags to be set.
   */
  std::vector<Field> expandStartBlock(const std::vector<Field>& v0)
  {
    int divisor = (useParityFlip ? 2 : 1) * (useGamma5 ? 2 : 1);
    std::vector<Field> startBlock;
    startBlock.reserve(Nblock);
    for (int i = 0; i < Nblock / divisor; i++) {
      std::vector<Field> group;
      group.push_back(v0[i]);
      if (useParityFlip) {
        int n = (int)group.size();
        for (int j = 0; j < n; j++) {
          Field fp(Grid_);
          parityFlippedField(group[j], fp);
          group.push_back(std::move(fp));
        }
      }
      if (useGamma5) {
        int n = (int)group.size();
        for (int j = 0; j < n; j++) {
          Field g5v(Grid_);
          gamma5Func(group[j], g5v);
          group.push_back(std::move(g5v));
        }
      }
      for (auto& f : group) startBlock.push_back(std::move(f));
    }
    return startBlock;
  }

protected:

  //--------------------------------------------------------------------
  // Block Arnoldi iteration
  //--------------------------------------------------------------------
  /**
   * Extends the block Arnoldi factorisation from block index 'start' to
   * block index 'Nm'.
   *
   * On entry (start > 0): basis[0..start*Nblock-1] already set,
   *   H[0..start*Nblock-1, 0..start*Nblock-1] already set,
   *   B[start*Nblock-1, :] set (coupling from prior residual block).
   *   startBlock = the normalised residual block F from the previous cycle.
   *
   * On entry (start == 0): initialises everything from startBlock.
   */
  void blockArnoldiIteration(std::vector<Field>& startBlock, int endBlock,
                             int startIdx, bool doubleOrthog)
  {
    int N = Nm;

    if (startIdx == 0) {
      basis.clear();
      F.clear();
      H = CMat::Zero(N, N);
      B = CMat::Zero(N, Nblock);

      // Orthonormalise starting block via modified Gram-Schmidt
      std::vector<Field> V0 = startBlock;
      blockOrthonormalise(V0);
      for (auto& v : V0) basis.push_back(v);
    } else {
      // Append residual block (startBlock = F_old) to basis.
      // The truncated KS relation after restart is:
      //
      //   A V_k = V_k S_k + F_old B_old^dag          (*)
      //
      // where V_k = basis[0:Nkeep], S_k is stored in H[0:Nkeep,0:Nkeep],
      // B_old = B[0:Nkeep,:], F_old = startBlock.
      //
      // Once F_old is appended as basis[Nkeep:Nkeep+Nblock], (*) becomes
      // a statement about the extended H matrix:
      //
      //   H[Nkeep+t, j] = (B_old^dag)[t,j] = conj(B_old[j,t])
      //                   for t=0..Nblock-1, j=0..Nkeep-1
      //
      // These entries are the "restart coupling rows" that connect the new
      // block to all retained Schur vectors and must be set before Arnoldi
      // continues, otherwise A V_k = V_k H[:,0:Nkeep] would be missing the
      // F_old B_old^dag term for those columns.

      int Nkeep = startIdx * Nblock;
      for (auto& v : startBlock) basis.push_back(v);

      // Fill restart coupling rows into H
      for (int t = 0; t < Nblock; t++)
        for (int j = 0; j < Nkeep; j++)
          H(Nkeep + t, j) = std::conj(B(j, t));

      // Zero out B for the retained columns now that the coupling is in H
      for (int j = 0; j < Nkeep; j++)
        for (int t = 0; t < Nblock; t++)
          B(j, t) = 0.0;
    }

    // Main block Arnoldi loop
    for (int k = startIdx; k < endBlock; k++) {
      blockArnoldiStep(k, doubleOrthog);
    }
  }

  //--------------------------------------------------------------------
  /**
   * One block Arnoldi step: extends by one block (Nblock vectors).
   *
   * Computes block column k of H and the next basis block V_{k+1}.
   *
   * Layout of basis (flat):
   *   basis[j*Nblock + t]  = t-th vector of j-th block,  j = 0..k
   *
   * After this call:
   *   H[i, k*Nblock : (k+1)*Nblock] filled for i = 0..(k+1)*Nblock - 1
   *   basis[k*Nblock .. (k+1)*Nblock - 1] normalised (already set on entry)
   *   F = residual block (to become V_{k+1} after this step if k < Nm-1)
   *
   * If k < Nm-1, also:
   *   H[(k+1)*Nblock : (k+2)*Nblock, k*Nblock : (k+1)*Nblock] = subdiag block (from QR of residual)
   *   basis extended by Nblock (the normalised residual vectors)
   */
  void blockArnoldiStep(int k, bool doubleOrthog)
  {
    int kBase  = k * Nblock;          // first flat index of current block
    int prevN  = kBase + Nblock;      // number of basis vectors so far after this step
    int N      = Nm;

    // W[t] = op(basis[kBase + t])  — dispatches to applyBlock() virtual
    std::vector<Field> W(Nblock, Field(Grid_));
    applyBlock(W, kBase);

    // Orthogonalise W against all current basis vectors (full reorthogonalisation)
    // H[i, kBase + t] = <basis[i] | W[t]>
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

    // Store residual block F
    F = W;

    if (k == Nm/Nblock - 1) {
      // Last block: compute coupling matrix B for KS decomp.
      //
      // blockQR modifies F in-place (F → Q orthonormal) and returns R
      // such that  W_orig[t] = sum_s F[s] * R[s,t]  (W_orig = F_after * R).
      //
      // The KS relation for column j = kBase+t requires the coefficient of F[s]
      // to be (B†)[s,j] = conj(B[j,s]).  Matching with R[s,t]:
      //   conj(B[kBase+t, s]) = R[s,t]  →  B[kBase+t, s] = conj(R[s,t])
      //
      // Equivalently the last Nblock rows of B are R^H (Hermitian conjugate of R).
      // Note: for Nblock=1, R is scalar real positive, so this reduces to B = R. ✓
      CMat R = blockQR(F);   // F is modified in-place to become Q; returns R
      for (int t = 0; t < Nblock; t++)
        for (int s = 0; s < Nblock; s++)
          B(kBase + t, s) = std::conj(R(s, t));    // B_block = R^H

      beta_k = R.norm();
      return;
    }

    // Not last block: QR-decompose residual to get V_{k+1}
    CMat R = blockQR(F);   // F orthonormalised in-place, R is upper triangular

    // Subdiagonal block of H: H[(k+1)*Nblock : (k+2)*Nblock, kBase : kBase+Nblock] = R
    int nextBase = (k + 1) * Nblock;
    for (int i = 0; i < Nblock; i++)
      for (int j = 0; j < Nblock; j++)
        H(nextBase + i, kBase + j) = R(i, j);

    // Append normalised residual block to basis
    for (int t = 0; t < Nblock; t++)
      basis.push_back(F[t]);
  }

  //--------------------------------------------------------------------
  // Block QR via modified Gram-Schmidt within the block
  //--------------------------------------------------------------------
  /**
   * Given a block of Nblock vectors W (not necessarily orthonormal),
   * orthonormalises them in-place and returns the upper-triangular R
   * such that W_in = W_out * R.
   *
   * Handles (near-)linear dependence by zeroing vectors below threshold.
   */
  CMat blockQR(std::vector<Field>& W)
  {
    CMat R = CMat::Zero(Nblock, Nblock);
    const RealD deflThresh = 1e-14;

    for (int j = 0; j < Nblock; j++) {
      // Orthogonalise W[j] against W[0..j-1]
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
        // deflation: zero this vector
        W[j] = Zero();
        std::cout << GridLogMessage
                  << className << ": deflation at block column " << j
                  << " (norm = " << nrm << ")" << std::endl;
      }
    }
    return R;
  }

  //--------------------------------------------------------------------
  // Orthonormalise a block against itself (no prior basis)
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
  // Basis rotation: UR[i] = sum_j U[j] * R(j, i)
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
  // Eigensystem of the truncated Rayleigh quotient
  //--------------------------------------------------------------------
  void computeEigensystem(CMat& Hk, int Nkeep)
  {
    Eigen::ComplexEigenSolver<CMat> es;
    es.compute(Hk);

    // Sort to match schurReorder ordering.
    int n = es.eigenvalues().size();
    ComplexComparator cComp(ritzFilter);
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b){
      return cComp(toStdCmplx(es.eigenvalues()(a)), toStdCmplx(es.eigenvalues()(b)));
    });

    evals.resize(n);
    littleEvecs.resize(n, n);
    for (int k = 0; k < n; k++) {
      evals(k)           = es.eigenvalues()(idx[k]);
      littleEvecs.col(k) = es.eigenvectors().col(idx[k]);
    }

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
   * An eigenpair (lambda, y) is converged if the Ritz estimate
   *   r = || B^dag y ||
   * satisfies r < rtol.  Here B is the (Nkeep x Nblock) coupling matrix
   * and y is the little eigenvector (Nkeep-vector) of H.
   */
  int converged(int Nkeep)
  {
    ritzEstimates.clear();
    int Nconv = 0;

    CMat Bk = B(Eigen::seqN(0, Nkeep), Eigen::all);  // Nkeep x Nblock

    for (int k = 0; k < Nkeep; k++) {
      CVec yk   = littleEvecs.col(k);                // Nkeep-vector
      CVec Bty  = Bk.adjoint() * yk;                // Nblock-vector
      RealD res = Bty.norm();
      ritzEstimates.push_back(res);
      std::cout << GridLogMessage << className << ": Ritz estimate[" << k
                << "] = " << res << "  eval = " << evals[k] << std::endl;
      if (res < rtol) Nconv++;
      else break;
    }
    return Nconv;
  }

  //--------------------------------------------------------------------
  // Approximate maximum eigenvalue via power iteration
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

#endif  // GRID_BLOCKED_KRYLOV_SCHUR_H
