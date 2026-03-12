/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./Grid/algorithms/iterative/LanczosBidiagonalization.h

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
#ifndef GRID_LANCZOS_BIDIAGONALIZATION_H
#define GRID_LANCZOS_BIDIAGONALIZATION_H

NAMESPACE_BEGIN(Grid);

/**
 * Lanczos Bidiagonalization (Golub-Kahan)
 *
 * For a linear operator A with adjoint A^dag, constructs the bidiagonal
 * decomposition:
 *
 *   A  V_m = U_m B_m
 *   A^dag U_m = V_m B_m^T + beta_{m+1} v_{m+1} e_m^T
 *
 * where:
 *   V_m = [v_1, ..., v_m]  right Lanczos vectors (orthonormal)
 *   U_m = [u_1, ..., u_m]  left  Lanczos vectors (orthonormal)
 *   B_m is upper bidiagonal with diag(alpha_1,...,alpha_m) and
 *       superdiag(beta_2,...,beta_m)
 *
 * The singular values of A are approximated by those of B_m.
 * The singular values of B_m are the square roots of the eigenvalues of
 * the symmetric tridiagonal matrix B_m^T B_m.
 *
 * Usage:
 *   LanczosBidiagonalization<Field> lb(Linop, grid);
 *   lb.run(src, Nm, tol);
 *   // Access results via getters.
 */
template <class Field>
class LanczosBidiagonalization {

  public: 
  LinearOperatorBase<Field> &Linop;
  GridBase *Grid;

  int Nm;           // number of Lanczos steps taken
  RealD Tolerance;  // convergence threshold on beta_{k+1} / alpha_k

  std::vector<Field>  V;       // right Lanczos vectors v_1 ... v_m
  std::vector<Field>  U;       // left  Lanczos vectors u_1 ... u_m
  std::vector<RealD>  alpha;   // diagonal of bidiagonal matrix
  std::vector<RealD>  beta;    // super-diagonal (beta[k] couples u_k and v_{k+1})

  // SVD of the bidiagonal matrix (filled after computeSVD())
  Eigen::VectorXd  singularValues;
  Eigen::MatrixXd  leftSVecs;   // columns are left  singular vectors of B
  Eigen::MatrixXd  rightSVecs;  // columns are right singular vectors of B

public:

  LanczosBidiagonalization(LinearOperatorBase<Field> &_Linop, GridBase *_Grid,
                           RealD _tol = 1.0e-8)
    : Linop(_Linop), Grid(_Grid), Tolerance(_tol), Nm(0)
  {}

  /**
   * Run the Golub-Kahan Lanczos bidiagonalization.
   *
   * Parameters
   * ----------
   * src  : starting vector (need not be normalised)
   * Nmax : maximum number of Lanczos steps
   * reorth : if true, full reorthogonalisation of both V and U bases
   */
  void run(const Field &src, int Nmax, bool reorth = true)
  {
    assert(norm2(src) > 0.0);

    V.clear(); U.clear();
    alpha.clear(); beta.clear();
    Nm = 0;

    Field p(Grid), r(Grid);

    // --- initialise: v_1 = src / ||src|| ---
    Field v(Grid);
    v = src;
    RealD nrm = std::sqrt(norm2(v));
    v = (1.0 / nrm) * v;
    V.push_back(v);

    for (int k = 0; k < Nmax; ++k) {

      // p = A v_k
      Linop.Op(V[k], p);

      // p = p - beta_k * u_{k-1}   (remove previous left vector)
      if (k > 0) {
        p = p - beta[k-1] * U[k-1];
      }

      // alpha_k = ||p||
      RealD ak = std::sqrt(norm2(p));
      if (ak < 1.0e-14) {
        std::cout << GridLogMessage
                  << "LanczosBidiagonalization: lucky breakdown at step "
                  << k << " (alpha = " << ak << ")" << std::endl;
        break;
      }
      alpha.push_back(ak);

      // u_k = p / alpha_k
      Field u(Grid);
      u = (1.0 / ak) * p;

      // full reorthogonalisation of u against previous U
      if (reorth) {
        for (int j = 0; j < (int)U.size(); ++j) {
          ComplexD ip = innerProduct(U[j], u);
          u = u - ip * U[j];
        }
        RealD unrm = std::sqrt(norm2(u));
        if (unrm > 1.0e-14) u = (1.0 / unrm) * u;
      }
      U.push_back(u);

      // r = A^dag u_k - alpha_k * v_k
      Linop.AdjOp(U[k], r);
      r = r - ak * V[k];

      // full reorthogonalisation of r against previous V
      if (reorth) {
        for (int j = 0; j < (int)V.size(); ++j) {
          ComplexD ip = innerProduct(V[j], r);
          r = r - ip * V[j];
        }
      }

      // beta_{k+1} = ||r||
      RealD bk = std::sqrt(norm2(r));
      beta.push_back(bk);

      Nm = k + 1;

      std::cout << GridLogMessage
                << "LanczosBidiagonalization step " << k
                << "  alpha = " << ak
                << "  beta  = " << bk << std::endl;

      // convergence: residual beta / alpha small enough
      if (bk / ak < Tolerance) {
        std::cout << GridLogMessage
                  << "LanczosBidiagonalization converged at step " << k
                  << "  (beta/alpha = " << bk / ak << ")" << std::endl;
        break;
      }

      if (k == Nmax - 1) break;   // no v_{k+2} needed after last step

      // v_{k+1} = r / beta_{k+1}
      Field vnext(Grid);
      vnext = (1.0 / bk) * r;
      V.push_back(vnext);
    }
  }

  /**
   * Compute the SVD of the bidiagonal matrix B using Eigen.
   * Singular values are stored in descending order.
   */
  void computeSVD()
  {
    int m = Nm;
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m, m);

    for (int k = 0; k < m; ++k) {
      B(k, k) = alpha[k];
      if (k + 1 < m && k < (int)beta.size())
        B(k, k+1) = beta[k];
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(B,
        Eigen::ComputeThinU | Eigen::ComputeThinV);

    singularValues = svd.singularValues();   // already sorted descending
    leftSVecs      = svd.matrixU();
    rightSVecs     = svd.matrixV();

    std::cout << GridLogMessage
              << "LanczosBidiagonalization: singular values of B_" << m
              << std::endl;
    for (int k = 0; k < m; ++k)
      std::cout << GridLogMessage << "  sigma[" << k << "] = "
                << singularValues(k) << std::endl;
  }

  /**
   * Return the k-th approximate left singular vector of A in the full
   * lattice space.  computeSVD() must have been called first.
   */
  Field leftSingularVector(int k)
  {
    assert(k < (int)leftSVecs.cols());
    Field svec(Grid);
    svec = Zero();
    for (int j = 0; j < Nm; ++j)
      svec = svec + leftSVecs(j, k) * U[j];
    return svec;
  }

  /**
   * Return the k-th approximate right singular vector of A in the full
   * lattice space.  computeSVD() must have been called first.
   */
  Field rightSingularVector(int k)
  {
    assert(k < (int)rightSVecs.cols());
    Field svec(Grid);
    svec = Zero();
    for (int j = 0; j < Nm; ++j)
      svec = svec + rightSVecs(j, k) * V[j];
    return svec;
  }

  /**
   * Verify the bidiagonalization: returns max residual
   *   max_k || A v_k - alpha_k u_k - beta_k u_{k-1} ||
   */
  RealD verify()
  {
    Field tmp(Grid);
    RealD maxres = 0.0;
    for (int k = 0; k < Nm; ++k) {
      Linop.Op(V[k], tmp);
      tmp = tmp - alpha[k] * U[k];
      if (k > 0 && k-1 < (int)beta.size())
        tmp = tmp - beta[k-1] * U[k-1];
      RealD res = std::sqrt(norm2(tmp));
      if (res > maxres) maxres = res;
      std::cout << GridLogMessage
                << "LanczosBidiagonalization verify step " << k
                << "  ||A v_k - alpha_k u_k - beta_{k-1} u_{k-1}|| = "
                << res << std::endl;
    }
    return maxres;
  }

  /* Getters */
  int                       getNm()           const { return Nm; }
  const std::vector<Field>& getV()            const { return V; }
  const std::vector<Field>& getU()            const { return U; }
  const std::vector<RealD>& getAlpha()        const { return alpha; }
  const std::vector<RealD>& getBeta()         const { return beta; }
  Eigen::VectorXd           getSingularValues() const { return singularValues; }
};

NAMESPACE_END(Grid);
#endif
