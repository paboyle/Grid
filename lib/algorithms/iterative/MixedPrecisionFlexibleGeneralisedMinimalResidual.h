/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/MixedPrecisionFlexibleGeneralisedMinimalResidual.h

Copyright (C) 2015

Author: Daniel Richtmann <daniel.richtmann@ur.de>

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
#ifndef GRID_MIXED_PRECISION_FLEXIBLE_GENERALISED_MINIMAL_RESIDUAL_H
#define GRID_MIXED_PRECISION_FLEXIBLE_GENERALISED_MINIMAL_RESIDUAL_H

namespace Grid {

template<class FieldD, class FieldF, typename std::enable_if<getPrecision<FieldD>::value == 2, int>::type = 0, typename std::enable_if< getPrecision<FieldF>::value == 1, int>::type = 0>
class MixedPrecisionFlexibleGeneralisedMinimalResidual : public OperatorFunction<FieldD> {
 public:
  bool ErrorOnNoConverge; // Throw an assert when MPFGMRES fails to converge,
                          // defaults to true

  RealD   Tolerance;

  Integer MaxIterations;
  Integer RestartLength;
  Integer MaxNumberOfRestarts;
  Integer IterationCount; // Number of iterations the MPFGMRES took to finish,
                          // filled in upon completion

  GridStopWatch MatrixTimer;
  GridStopWatch PrecTimer;
  GridStopWatch LinalgTimer;
  GridStopWatch QrTimer;
  GridStopWatch CompSolutionTimer;
  GridStopWatch ChangePrecTimer;

  Eigen::MatrixXcd H;

  std::vector<std::complex<double>> y;
  std::vector<std::complex<double>> gamma;
  std::vector<std::complex<double>> c;
  std::vector<std::complex<double>> s;

  GridBase* SinglePrecGrid;

  LinearFunction<FieldF> &Preconditioner;

  MixedPrecisionFlexibleGeneralisedMinimalResidual(RealD   tol,
                                                   Integer maxit,
                                                   GridBase * sp_grid,
                                                   LinearFunction<FieldF> &Prec,
                                                   Integer restart_length,
                                                   bool    err_on_no_conv = true)
      : Tolerance(tol)
      , MaxIterations(maxit)
      , RestartLength(restart_length)
      , MaxNumberOfRestarts(MaxIterations/RestartLength + ((MaxIterations%RestartLength == 0) ? 0 : 1))
      , ErrorOnNoConverge(err_on_no_conv)
      , H(Eigen::MatrixXcd::Zero(RestartLength, RestartLength + 1)) // sizes taken from DD-αAMG code base
      , y(RestartLength + 1, 0.)
      , gamma(RestartLength + 1, 0.)
      , c(RestartLength + 1, 0.)
      , s(RestartLength + 1, 0.)
      , SinglePrecGrid(sp_grid)
      , Preconditioner(Prec) {};

  void operator()(LinearOperatorBase<FieldD> &LinOp, const FieldD &src, FieldD &psi) {

    psi.checkerboard = src.checkerboard;
    conformable(psi, src);

    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    RealD cp;
    RealD ssq = norm2(src);
    RealD rsq = Tolerance * Tolerance * ssq;

    FieldD r(src._grid);

    std::cout << std::setprecision(4) << std::scientific;
    std::cout << GridLogIterative << "MPFGMRES: guess " << guess << std::endl;
    std::cout << GridLogIterative << "MPFGMRES:   src " << ssq   << std::endl;

    PrecTimer.Reset();
    MatrixTimer.Reset();
    LinalgTimer.Reset();
    QrTimer.Reset();
    CompSolutionTimer.Reset();
    ChangePrecTimer.Reset();

    GridStopWatch SolverTimer;
    SolverTimer.Start();

    IterationCount = 0;

    for (int k=0; k<MaxNumberOfRestarts; k++) {

      cp = outerLoopBody(LinOp, src, psi, rsq);

      // Stopping condition
      if (cp <= rsq) {

        SolverTimer.Stop();

        LinOp.Op(psi,r);
        axpy(r,-1.0,src,r);

        RealD srcnorm       = sqrt(ssq);
        RealD resnorm       = sqrt(norm2(r));
        RealD true_residual = resnorm / srcnorm;

        std::cout << GridLogMessage        << "MPFGMRES: Converged on iteration " << IterationCount
                  << " computed residual " << sqrt(cp / ssq)
                  << " true residual "     << true_residual
                  << " target "            << Tolerance << std::endl;

        std::cout << GridLogMessage << "MPFGMRES Time elapsed: Total      " <<       SolverTimer.Elapsed() << std::endl;
        std::cout << GridLogMessage << "MPFGMRES Time elapsed: Precon     " <<         PrecTimer.Elapsed() << std::endl;
        std::cout << GridLogMessage << "MPFGMRES Time elapsed: Matrix     " <<       MatrixTimer.Elapsed() << std::endl;
        std::cout << GridLogMessage << "MPFGMRES Time elapsed: Linalg     " <<       LinalgTimer.Elapsed() << std::endl;
        std::cout << GridLogMessage << "MPFGMRES Time elapsed: QR         " <<           QrTimer.Elapsed() << std::endl;
        std::cout << GridLogMessage << "MPFGMRES Time elapsed: CompSol    " << CompSolutionTimer.Elapsed() << std::endl;
        std::cout << GridLogMessage << "MPFGMRES Time elapsed: PrecChange " <<   ChangePrecTimer.Elapsed() << std::endl;
        return;
      }
    }

    std::cout << GridLogMessage << "MPFGMRES did NOT converge" << std::endl;

    if (ErrorOnNoConverge)
      assert(0);
  }

  RealD outerLoopBody(LinearOperatorBase<FieldD> &LinOp, const FieldD &src, FieldD &psi, RealD rsq) {

    RealD cp = 0;

    FieldD w(src._grid);
    FieldD r(src._grid);

    std::vector<FieldD> v(RestartLength + 1, src._grid); // these should probably be made class members
    std::vector<FieldD> z(RestartLength + 1, src._grid); // so that they are only allocated once, not in every restart

    MatrixTimer.Start();
    LinOp.Op(psi, w);
    MatrixTimer.Stop();

    LinalgTimer.Start();
    r = src - w;

    gamma[0] = sqrt(norm2(r));

    v[0] = (1. / gamma[0]) * r;
    LinalgTimer.Stop();

    for (int i=0; i<RestartLength; i++) {

      IterationCount++;

      arnoldiStep(LinOp, v, z, w, i);

      qrUpdate(i);

      cp = std::norm(gamma[i+1]);

      std::cout << GridLogIterative << "MPFGMRES: Iteration " << IterationCount
                << " residual " << cp << " target " << rsq << std::endl;

      if ((i == RestartLength - 1) || (IterationCount == MaxIterations) || (cp <= rsq)) {

        computeSolution(z, psi, i);

        return cp;
      }
    }

    assert(0); // Never reached
    return cp;
  }

  void arnoldiStep(LinearOperatorBase<FieldD> &LinOp, std::vector<FieldD> &v, std::vector<FieldD> &z, FieldD &w, int iter) {

    FieldF v_f(SinglePrecGrid);
    FieldF z_f(SinglePrecGrid);

    ChangePrecTimer.Start();
    precisionChange(v_f, v[iter]);
    precisionChange(z_f, z[iter]);
    ChangePrecTimer.Stop();

    PrecTimer.Start();
    Preconditioner(v_f, z_f);
    PrecTimer.Stop();

    ChangePrecTimer.Start();
    precisionChange(z[iter], z_f);
    ChangePrecTimer.Stop();

    MatrixTimer.Start();
    LinOp.Op(z[iter], w);
    MatrixTimer.Stop();

    LinalgTimer.Start();
    for (int i = 0; i <= iter; ++i) {
      H(iter, i) = innerProduct(v[i], w);
      w = w - H(iter, i) * v[i];
    }

    H(iter, iter + 1) = sqrt(norm2(w));
    v[iter + 1] = (1. / H(iter, iter + 1)) * w;
    LinalgTimer.Stop();
  }

  void qrUpdate(int iter) {

    QrTimer.Start();
    for (int i = 0; i < iter ; ++i) {
      auto tmp       = -s[i] * H(iter, i) + c[i] * H(iter, i + 1);
      H(iter, i)     = std::conj(c[i]) * H(iter, i) + std::conj(s[i]) * H(iter, i + 1);
      H(iter, i + 1) = tmp;
    }

    // Compute new Givens Rotation
    ComplexD nu = sqrt(std::norm(H(iter, iter)) + std::norm(H(iter, iter + 1)));
    c[iter]     = H(iter, iter) / nu;
    s[iter]     = H(iter, iter + 1) / nu;

    // Apply new Givens rotation
    H(iter, iter)     = nu;
    H(iter, iter + 1) = 0.;

    gamma[iter + 1] = -s[iter] * gamma[iter];
    gamma[iter]     = std::conj(c[iter]) * gamma[iter];
    QrTimer.Stop();
  }

  void computeSolution(std::vector<FieldD> const &z, FieldD &psi, int iter) {

    CompSolutionTimer.Start();
    for (int i = iter; i >= 0; i--) {
      y[i] = gamma[i];
      for (int k = i + 1; k <= iter; k++)
        y[i] = y[i] - H(k, i) * y[k];
      y[i] = y[i] / H(i, i);
    }

    for (int i = 0; i <= iter; i++)
      psi = psi + z[i] * y[i];
    CompSolutionTimer.Stop();
  }
};
}
#endif
