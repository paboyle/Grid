/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: lib/algorithms/iterative/GeneralisedMinimalResidual.h

Copyright (C) 2015
Copyright (C) 2016


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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef GRID_GENERALISED_MINIMAL_RESIDUAL_H
#define GRID_GENERALISED_MINIMAL_RESIDUAL_H

// from Y. Saad - Iterative Methods for Sparse Linear Systems, PP 172
// Compute r0 = b − Ax0 , β := ||r0||2 , and v1 := r0 /β
// For j = 1, 2, ..., m Do:
//   Compute wj := Avj
//   For i = 1, ..., j Do:
//     hij := (wj , vi)
//     wj := wj − hij vi
//   EndDo
//   hj+1,j = ||wj||2 . If hj+1,j = 0 set m := j and go to HERE
//   vj+1 = wj /hj+1,j
// EndDo
// Define the (m + 1) × m Hessenberg matrix H̄m = {hij}1≤i≤m+1,1≤j≤m. [HERE]
// Compute ym the minimizer of ||βe1 − H̄m y||2 and xm = x0 + Vm ym.
///////////////////////////////////////////////////////////////////////////////////////////////////////

// want to solve Ax = b -> A = LinOp, psi = x, b = src

namespace Grid {

template<class Field>
class GeneralisedMinimalResidual : public OperatorFunction<Field> {
 public:
  bool ErrorOnNoConverge; // Throw an assert when GMRES fails to converge,
                          // defaults to True.
  RealD   Tolerance;
  Integer MaxIterations;
  Integer RestartLength;
  Integer IterationsToComplete; // Number of iterations the GMRES took to
                                // finish. Filled in upon completion

  GeneralisedMinimalResidual(RealD   tol,
                             Integer maxit,
                             Integer restart_length,
                             bool    err_on_no_conv = true)
    : Tolerance(tol), MaxIterations(maxit), RestartLength(restart_length), ErrorOnNoConverge(err_on_no_conv){};

  // want to solve Ax = b -> A = LinOp, psi = x, b = src

  /* void */
  /* operator()(LinearOperatorBase<Field> &LinOp, const Field &src, Field &psi)
   * { */
  /*   typedef typename Eigen::MatrixXcd MyMatrix; */
  /*   typedef typename Eigen::VectorXcd MyVector; */

  /*   Field r(src); */
  /*   Field w(src); */
  /*   Field mmv(src); */

  /*   std::vector<Field>                V(MaxIterations + 1, src); */
  /*   std::vector<std::complex<double>> y(MaxIterations + 1, 0.); */
  /*   std::vector<std::complex<double>> gamma(MaxIterations + 1, 0.); */
  /*   std::vector<std::complex<double>> c(MaxIterations + 1, 0.); */
  /*   std::vector<std::complex<double>> s(MaxIterations + 1, 0.); */

  /*   int m = MaxIterations; */

  /*   RealD gamma0{}; */

  /*   MyMatrix H = Eigen::MatrixXcd::Zero(MaxIterations + 1, MaxIterations); */

  /*   RealD normPsiSq   = norm2(psi); */
  /*   RealD normSrcSq   = norm2(src); */
  /*   RealD TargetResSq = Tolerance * Tolerance * normSrcSq; */

  /*   LinOp.Op(psi, mmv); */

  /*   r        = src - mmv; */
  /*   gamma[0] = norm2(r); */
  /*   std::cout << gamma[0] << std::endl; */
  /*   gamma0 = std::real(gamma[0]); */
  /*   V[0]   = (1. / gamma[0]) * r; */

  /*   std::cout << GridLogMessage << std::setprecision(4) */
  /*             << "GeneralisedMinimalResidual:    psi " << normPsiSq */
  /*             << std::endl; */
  /*   std::cout << GridLogMessage << std::setprecision(4) */
  /*             << "GeneralisedMinimalResidual:    src " << normSrcSq */
  /*             << std::endl; */
  /*   std::cout << GridLogMessage << std::setprecision(4) */
  /*             << "GeneralisedMinimalResidual: target " << TargetResSq */
  /*             << std::endl; */
  /*   std::cout << GridLogMessage << std::setprecision(4) */
  /*             << "GeneralisedMinimalResidual:      r " << gamma0 <<
   * std::endl; */

  /*   std::cout */
  /*     << GridLogIterative << std::setprecision(4) */
  /*     << "GeneralisedMinimalResidual: before starting to iterate residual "
   */
  /*     << gamma0 << " target " << TargetResSq << std::endl; */

  /*   for(auto j = 0; j < m; ++j) { */
  /*     LinOp.Op(V[j], w); */

  /*     for(auto i = 0; i <= j; ++i) { */
  /*       H(i, j) = innerProduct(V[i], w); */
  /*       w = w - H(i, j) * V[i]; */
  /*     } */

  /*     H(j + 1, j) = norm2(w); */
  /*     V[j + 1] = (1. / H(j + 1, j)) * w; */

  /*     if(std::abs(H(j + 1, j)) > 1e-15) { */
  /*       qrUpdate(gamma, c, s, H, j); */
  /*     } */

  /*     /\* std::cout << GridLogMessage << "GeneralisedMinimalResidual: H( "
   * *\/ */
  /*     /\*           << j + 1 << "," << j << " ) = " << H( j + 1, j ) *\/ */
  /*     /\*           << std::endl; *\/ */

  /*     std::cout << GridLogIterative << "GeneralisedMinimalResidual: Iteration
   * " */
  /*               << j << " residual " << std::abs(gamma[j + 1]) << " target "
   */
  /*               << TargetResSq << std::endl; */
  /*     if(std::abs(gamma[j + 1]) / gamma0 < Tolerance) { */
  /*       IterationsToComplete = j; */
  /*       break; */
  /*     } */
  /*   } */
  /*   computeSolution(y, gamma, H, V, psi, IterationsToComplete); */
  /*   std::cout << GridLogMessage */
  /*             << "GeneralisedMinimalResidual: End of operator() after " */
  /*             << IterationsToComplete << " iterations" << std::endl; */

  /*   RealD normSrc       = sqrt(normSrcSq); */
  /*   RealD resnorm       = sqrt(norm2(mmv)); */
  /*   RealD true_residual = resnorm / srcnorm; */
  /*   Field result        = mmv; */
  /*   Field Dx(src); */
  /*   Field tmp(src); */

  /*   // Test the correctness */
  /*   LinOp.Op(result, Dx); */

  /*   tmp = Dx - src; */

  /*   std::cout << norm2(tmp) << " " << norm2(tmp) / gamma0 << std::endl; */
  /* } */

  void operator()(LinearOperatorBase<Field> &LinOp, const Field &src, Field &psi) {

    std::cout << GridLogIterative << "GMRES: Start of operator()" << std::endl;

    psi.checkerboard = src.checkerboard;
    conformable(psi, src);

    int m = MaxIterations;

    Field r(src);
    Field w(src);
    Field Dpsi(src);
    Field Dv(src);

    std::vector<Field> v(m + 1, src);

    Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(m + 1, m);

    std::vector<std::complex<double>> y(m + 1, 0.);
    std::vector<std::complex<double>> gamma(m + 1, 0.);
    std::vector<std::complex<double>> c(m + 1, 0.);
    std::vector<std::complex<double>> s(m + 1, 0.);

    // Initial residual computation & set up
    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    RealD ssq = norm2(src); // flopcount.addSiteFlops(4*Nc*Ns,s); // stands for "source squared"
    RealD rsd_sq = Tolerance * Tolerance * ssq; // flopcount.addSiteFlops(4*Nc*Ns,s); // stands for "residual squared"

    LinOp.Op(psi, Dpsi);
    r = src - Dpsi;

    RealD cp = norm2(r); // cp = beta in WMG nomenclature, in WMG there is no norm2 but a sqrt(norm2) here
    gamma[0]  = cp;

    std::cout << GridLogIterative << "cp " << cp << std::endl;

    v[0] = (1. / cp) * r;

    std::cout << GridLogIterative << std::setprecision(4) << "GeneralizedMinimalResidual: guess " << guess << std::endl;
    std::cout << GridLogIterative << std::setprecision(4) << "GeneralizedMinimalResidual:   src " << ssq << std::endl;
    // std::cout << GridLogIterative << std::setprecision(4) << "GeneralizedMinimalResidual:    mp " << d << std::endl;
    std::cout << GridLogIterative << std::setprecision(4) << "GeneralizedMinimalResidual:  cp,r " << cp << std::endl;

    if (cp <= rsd_sq) {
      return;
    }

    std::cout << GridLogIterative << std::setprecision(4)
              << "GeneralizedMinimalResidual: k=0 residual " << cp << " target " << rsd_sq << std::endl;

    GridStopWatch SolverTimer;
    GridStopWatch MatrixTimer;

    SolverTimer.Start();
    for(auto j = 0; j < m; ++j) {

      // std::cout << GridLogIterative << "GeneralizedMinimalResidual: Start of outer loop with index j = " << j << std::endl;

      MatrixTimer.Start();
      LinOp.Op(v[j], Dv);
      MatrixTimer.Stop();

      w = Dv;

      for(auto i = 0; i <= j; ++i) {
        H(i, j) = innerProduct(v[i], w);
        w = w - H(i, j) * v[i];
      }

      H(j + 1, j) = norm2(w);
      v[j + 1] = (1. / H(j + 1, j)) * w;

      // end of arnoldi process, begin of givens rotations
      // apply old Givens rotation
      for(auto i = 0; i < j ; ++i) {
        auto tmp = -s[i] * H(i, j) + c[i] * H(i + 1, j);
        H(i, j)     = std::conj(c[i]) * H(i, j) + std::conj(s[i]) * H(i + 1, j);
        H(i + 1, j) = tmp;
      }

      // compute new Givens Rotation
      ComplexD nu = sqrt(std::norm(H(j, j)) + std::norm(H(j + 1, j)));
      c[j]        = H(j, j) / nu;
      s[j]        = H(j + 1, j) / nu;
      std::cout << GridLogIterative << "GeneralizedMinimalResidual: nu" << nu << std::endl;
      std::cout << GridLogIterative << "GeneralizedMinimalResidual: H("<<j<<","<<j<<")" << H(j,j) << std::endl;
      std::cout << GridLogIterative << "GeneralizedMinimalResidual: H("<<j+1<<","<<j<<")" << H(j+1,j) << std::endl;

      // apply new Givens rotation
      H(j, j)     = nu;
      H(j + 1, j) = 0.;

      /* ORDERING??? */
      gamma[j + 1] = -s[j] * gamma[j];
      gamma[j]     = std::conj(c[j]) * gamma[j];

      /* for(auto k = 0; k <= j+1 ; ++k) */
      /*   std::cout << GridLogIterative << "k " << k << "nu " << nu << " c["<<k<<"]" << c[k]<< " s["<<k<<"]" << s[k] << " gamma["<<k<<"]" << gamma[k] << std::endl; */

      std::cout << GridLogIterative << "GeneralisedMinimalResidual: Iteration "
                << j << " residual " << std::abs(gamma[j + 1]) << std::endl; //" target "
                /* << TargetResSq << std::endl; */
      if(std::abs(gamma[j + 1]) / sqrt(cp) < Tolerance) {
        SolverTimer.Stop();

        std::cout << GridLogMessage << "GeneralizedMinimalResidual Converged on iteration " << j << std::endl;
        // std::cout << GridLogMessage << "\tComputed residual " << sqrt(cp / ssq) << std::endl;
        // std::cout << GridLogMessage << "\tTrue residual " << true_residual << std::endl;
        std::cout << GridLogMessage << "\tTarget " << Tolerance << std::endl;

        std::cout << GridLogMessage << "Time breakdown " << std::endl;
        std::cout << GridLogMessage << "\tElapsed    " << SolverTimer.Elapsed() << std::endl;
        std::cout << GridLogMessage << "\tMatrix     " << MatrixTimer.Elapsed() << std::endl;
        // std::cout << GridLogMessage << "\tLinalg     " << LinalgTimer.Elapsed() << std::endl;

        IterationsToComplete = j;

        break;
      }
    }

    // backward substitution
    computeSolution(y, gamma, H, v, psi, IterationsToComplete);

    std::cout << GridLogIterative << "GeneralizedMinimalResidual: End of operator()" << std::endl;
  }

  void alternativeOperatorImplementation()(LinearOperatorBase<Field> &LinOp, const Field &src, Field &psi) {

    psi.checkerboard = src.checkerboard;
    psi(conformable, src);

    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    RealD cp;
    RealD ssq    = norm2(src);
    RealD rsd_sq = Tolerance * Tolerance * ssq;

    Field r(src._grid);
    Field Dpsi(src._grid);

    PrecTimer.Reset();
    MatTimer.Reset();
    LinalgTimer.Reset();

    GridStopWatch SolverTimer;
    SolverTimer.Start();

    int iterations = 0;
    for (int k=0; k<MaxIterations; k++) {

      cp = outerLoopBody();

      // Stopping condition
      if (cp <= rsd_sq) {

        SolverTimer.Stop();

        Linop.Op(psi, Dpsi); // maybe can improve these two lines
        r = src - Dpsi;      // by technique used in VPGCR

        RealD srcnorm       = sqrt(ssq);
        RealD resnorm       = sqrt(norm2(r));
        RealD true_residual = resnorm / srcnorm;

        std::cout << GridLogMessage << "GeneralizedMinimalResidual: Converged on iteration " << k              << std::endl;
        std::cout << GridLogMessage << "\tComputed residual "                                << sqrt(cp / ssq) << std::endl;
        std::cout << GridLogMessage << "\tTrue residual "                                    << true_residual  << std::endl;
        std::cout << GridLogMessage << "\tTarget "                                           << Tolerance      << std::endl;

        std::cout << GridLogMessage << "GeneralizedMinimalResidual Time breakdown" << std::endl;
        std::cout << GridLogMessage << "\tElapsed " << SolverTimer.Elapsed()       << std::endl;
        std::cout << GridLogMessage << "\tPrecon "  << PrecTimer.Elapsed()         << std::endl;
        std::cout << GridLogMessage << "\tMatrix "  << MatTimer.Elapsed()          << std::endl;
        std::cout << GridLogMessage << "\tLinalg "  << LinalgTimer.Elapsed()       << std::endl;
        return;
      }
    }

    std::cout << GridLogMessage << "GeneralizedMinimalResidual did NOT converge" << std::endl;

    if (ErrorOnNoConverge)
      assert(0);
  }

  RealD outerLoopBody() {
  }

  void Step() {

    int m = MaxIterations;

    Field r(src);
    Field w(src);
    Field Dpsi(src);
    Field Dv(src);
    std::vector<Field> v(m + 1, src);

    Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(m + 1, m);

    std::vector<std::complex<double>> y(m + 1, 0.);
    std::vector<std::complex<double>> gamma(m + 1, 0.);
    std::vector<std::complex<double>> c(m + 1, 0.);
    std::vector<std::complex<double>> s(m + 1, 0.);

    // Initial residual computation & set up
    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    RealD ssq = norm2(src); // flopcount.addSiteFlops(4*Nc*Ns,s); // stands for "source squared"
    RealD rsd_sq = Tolerance * Tolerance * ssq; // flopcount.addSiteFlops(4*Nc*Ns,s); // stands for "residual squared"

    LinOp.Op(psi, Dpsi);
    r = src - Dpsi;

    RealD cp = norm2(r); // cp = beta in WMG nomenclature, in WMG there is no norm2 but a sqrt(norm2) here
    gamma[0]  = cp;

    std::cout << GridLogIterative << "cp " << cp << std::endl;

    v[0] = (1. / cp) * r;

    std::cout << GridLogIterative << std::setprecision(4) << "GeneralizedMinimalResidual: guess " << guess << std::endl;
    std::cout << GridLogIterative << std::setprecision(4) << "GeneralizedMinimalResidual:   src " << ssq << std::endl;
    // std::cout << GridLogIterative << std::setprecision(4) << "GeneralizedMinimalResidual:    mp " << d << std::endl;
    std::cout << GridLogIterative << std::setprecision(4) << "GeneralizedMinimalResidual:  cp,r " << cp << std::endl;

    if (cp <= rsd_sq) {
      return;
    }

    std::cout << GridLogIterative << std::setprecision(4)
              << "GeneralizedMinimalResidual: k=0 residual " << cp << " target " << rsd_sq << std::endl;

    GridStopWatch SolverTimer;
    GridStopWatch MatrixTimer;

    SolverTimer.Start();
    for(auto j = 0; j < m; ++j) {

      // std::cout << GridLogIterative << "GeneralizedMinimalResidual: Start of outer loop with index j = " << j << std::endl;

      MatrixTimer.Start();
      LinOp.Op(v[j], Dv);
      MatrixTimer.Stop();

      w = Dv;

      for(auto i = 0; i <= j; ++i) {
        H(i, j) = innerProduct(v[i], w);
        w = w - H(i, j) * v[i];
      }

      H(j + 1, j) = norm2(w);
      v[j + 1] = (1. / H(j + 1, j)) * w;

      // end of arnoldi process, begin of givens rotations
      // apply old Givens rotation
      for(auto i = 0; i < j ; ++i) {
        auto tmp = -s[i] * H(i, j) + c[i] * H(i + 1, j);
        H(i, j)     = std::conj(c[i]) * H(i, j) + std::conj(s[i]) * H(i + 1, j);
        H(i + 1, j) = tmp;
      }

      // compute new Givens Rotation
      ComplexD nu = sqrt(std::norm(H(j, j)) + std::norm(H(j + 1, j)));
      c[j]        = H(j, j) / nu;
      s[j]        = H(j + 1, j) / nu;
      std::cout << GridLogIterative << "GeneralizedMinimalResidual: nu" << nu << std::endl;
      std::cout << GridLogIterative << "GeneralizedMinimalResidual: H("<<j<<","<<j<<")" << H(j,j) << std::endl;
      std::cout << GridLogIterative << "GeneralizedMinimalResidual: H("<<j+1<<","<<j<<")" << H(j+1,j) << std::endl;

      // apply new Givens rotation
      H(j, j)     = nu;
      H(j + 1, j) = 0.;

      /* ORDERING??? */
      gamma[j + 1] = -s[j] * gamma[j];
      gamma[j]     = std::conj(c[j]) * gamma[j];

      /* for(auto k = 0; k <= j+1 ; ++k) */
      /*   std::cout << GridLogIterative << "k " << k << "nu " << nu << " c["<<k<<"]" << c[k]<< " s["<<k<<"]" << s[k] << " gamma["<<k<<"]" << gamma[k] << std::endl; */

      std::cout << GridLogIterative << "GeneralisedMinimalResidual: Iteration "
                << j << " residual " << std::abs(gamma[j + 1]) << std::endl; //" target "
                /* << TargetResSq << std::endl; */
      if(std::abs(gamma[j + 1]) / sqrt(cp) < Tolerance) {
        SolverTimer.Stop();

        std::cout << GridLogMessage << "GeneralizedMinimalResidual Converged on iteration " << j << std::endl;
        // std::cout << GridLogMessage << "\tComputed residual " << sqrt(cp / ssq) << std::endl;
        // std::cout << GridLogMessage << "\tTrue residual " << true_residual << std::endl;
        std::cout << GridLogMessage << "\tTarget " << Tolerance << std::endl;

        std::cout << GridLogMessage << "Time breakdown " << std::endl;
        std::cout << GridLogMessage << "\tElapsed    " << SolverTimer.Elapsed() << std::endl;
        std::cout << GridLogMessage << "\tMatrix     " << MatrixTimer.Elapsed() << std::endl;
        // std::cout << GridLogMessage << "\tLinalg     " << LinalgTimer.Elapsed() << std::endl;

        IterationsToComplete = j;

        break;
      }
    }

    // backward substitution
    computeSolution(y, gamma, H, v, psi, IterationsToComplete);

    std::cout << GridLogIterative << "GeneralizedMinimalResidual: End of operator()" << std::endl;
  }

  private:
  /*  void qrUpdate(std::vector<std::complex<double>> &gamma, */
  /*                std::vector<std::complex<double>> &c, */
  /*                std::vector<std::complex<double>> &s, */
  /*                Eigen::MatrixXcd &                 H, */
  /*                int                                j) { */
  /*    ComplexD cp{}; */
  /*    // update QR factorization */
  /*    // apply previous Givens rotation */
  /*    for(auto i = 0; i < j; i++) { */
  /*      cp = -s[i] * H(i, j) + c[i] * H(i + 1, j); */
  /*      H(i, j)     = std::conj(c[i]) * H(i, j) + std::conj(s[i]) * H(i + 1,
   * j); */
  /*      H(i + 1, j) = cp; */
  /*    } */

  /*    // compute current Givens rotation */
  /*    cp = sqrt(std::norm(H(j, j)) + std::norm(H(j + 1, j))); */
  /*    s[j] = H(j + 1, j) / cp; */
  /*    c[j] = H(j, j) / cp; */
  /*    /\* std::cout << GridLogIterative << "cp= " << cp << std::endl; *\/ */
  /*    /\* std::cout << GridLogIterative << "s[j]= " << s[ j ] << std::endl; *\/ */
  /*    /\* std::cout << GridLogIterative << "c[j]= " << c[ j ] << std::endl; *\/ */

  /*    /\* std::cout << GridLogIterative << "gamma[j+1]= " << gamma[ j + 1 ] << std::endl; *\/ */
  /*    /\* std::cout << GridLogIterative << "gamma[j]= " << gamma[ j ] << std::endl; *\/ */
  /*    // update right column */
  /*    gamma[j + 1] = -s[j] * gamma[j]; */
  /*    gamma[j]     = std::conj(c[j]) * gamma[j]; */
  /*    /\* std::cout << GridLogIterative << "gamma[j+1]= " << gamma[ j + 1 ] << std::endl; *\/ */
  /*    /\* std::cout << GridLogIterative << "gamma[j]= " << gamma[ j ] << std::endl; *\/ */

  /*    // apply current Givens rotation */
  /*    H(j, j)     = cp; */
  /*    H(j + 1, j) = 0.; */
  /*    /\* std::cout << GridLogIterative << "H(j,j)= " << H( j, j ) << std::endl; *\/ */
  /*    /\* std::cout << GridLogIterative << "H(j+1,j)= " << H( j + 1, j ) << std::endl; *\/ */
  /*  } */

  void computeSolution(std::vector<std::complex<double>> &      y,
                       std::vector<std::complex<double>> const &gamma,
                       Eigen::MatrixXcd const &                 H,
                       std::vector<Field> const &               v,
                       Field &                                  x,
                       int                                      j) {
    for(auto i = j; i >= 0; i--) {
      y[i] = gamma[i];
      for(auto k = i + 1; k <= j; k++)
        y[i] -= H(i, k) * y[k];
      y[i] /= H(i, i);
    }

    /* if(true) // TODO ??? */
    /* { */
    /*   for(auto i = 0; i <= j; i++) */
    /*     x = x + v[i] * y[i]; */
    /* } else { */
      x = y[0] * v[0];
      for(auto i = 1; i <= j; i++)
        x = x + v[i] * y[i];
    /* } */
  }
};
}
#endif
