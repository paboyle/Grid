/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/FlexibleGeneralisedMinimalResidual.h

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
#ifndef GRID_FLEXIBLE_GENERALISED_MINIMAL_RESIDUAL_H
#define GRID_FLEXIBLE_GENERALISED_MINIMAL_RESIDUAL_H

namespace Grid {

template<class Field>
class FlexibleGeneralisedMinimalResidual : public OperatorFunction<Field> {
 public:
  bool ErrorOnNoConverge; // Throw an assert when FGMRES fails to converge,
                          // defaults to true

  RealD   Tolerance;

  Integer MaxIterations;
  Integer RestartLength;
  Integer IterationCount; // Number of iterations the FGMRES took to finish,
                          // filled in upon completion

  GridStopWatch MatrixTimer;
  GridStopWatch PrecTimer;
  GridStopWatch LinalgTimer;
  GridStopWatch QrTimer;
  GridStopWatch CompSolutionTimer;

  Eigen::MatrixXcd H;

  std::vector<std::complex<double>> y;
  std::vector<std::complex<double>> gamma;
  std::vector<std::complex<double>> c;
  std::vector<std::complex<double>> s;

  LinearFunction<Field> &Preconditioner;

  FlexibleGeneralisedMinimalResidual(RealD   tol,
                                     Integer maxit,
                                     LinearFunction<Field> &Prec,
                                     Integer restart_length,
                                     bool    err_on_no_conv = true)
      : Tolerance(tol)
      , MaxIterations(maxit)
      , RestartLength(restart_length)
      , ErrorOnNoConverge(err_on_no_conv)
      , H(Eigen::MatrixXcd::Zero(RestartLength, RestartLength + 1)) // sizes taken from DD-αAMG code base
      , y(RestartLength + 1, 0.)
      , gamma(RestartLength + 1, 0.)
      , c(RestartLength + 1, 0.)
      , s(RestartLength + 1, 0.)
      , Preconditioner(Prec) {};

  void operator()(LinearOperatorBase<Field> &LinOp, const Field &src, Field &psi) {

    psi.checkerboard = src.checkerboard;
    conformable(psi, src);

    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    RealD cp;
    RealD ssq = norm2(src);
    RealD rsq = Tolerance * Tolerance * ssq;

    Field r(src._grid);

    std::cout << std::setprecision(4) << std::scientific << std::endl;
    std::cout << GridLogIterative << "FlexibleGeneralisedMinimalResidual: guess " << guess << std::endl;
    std::cout << GridLogIterative << "FlexibleGeneralisedMinimalResidual:   src " << ssq   << std::endl;


    PrecTimer.Reset();
    MatrixTimer.Reset();
    LinalgTimer.Reset();
    QrTimer.Reset();
    CompSolutionTimer.Reset();

    GridStopWatch SolverTimer;
    SolverTimer.Start();

    IterationCount = 0;
    for (int k=0; k<MaxIterations; k++) {

      cp = outerLoopBody(LinOp, src, psi, rsq);

      // Stopping condition
      if (cp <= rsq) {

        SolverTimer.Stop();

        LinOp.Op(psi,r);
        axpy(r,-1.0,src,r);

        RealD srcnorm       = sqrt(ssq);
        RealD resnorm       = sqrt(norm2(r));
        RealD true_residual = resnorm / srcnorm;

        std::cout << GridLogMessage << "FlexibleGeneralisedMinimalResidual: Converged on iteration " << IterationCount << std::endl;
        std::cout << GridLogMessage << "\tComputed residual "                                        << sqrt(cp / ssq) << std::endl;
        std::cout << GridLogMessage << "\tTrue residual "                                            << true_residual  << std::endl;
        std::cout << GridLogMessage << "\tTarget "                                                   << Tolerance      << std::endl;

        std::cout << GridLogMessage << "FlexibleGeneralisedMinimalResidual Time breakdown" << std::endl;
        std::cout << GridLogMessage << "\tElapsed " << SolverTimer.Elapsed()               << std::endl;
        std::cout << GridLogMessage << "\tPrecon "  << PrecTimer.Elapsed()                 << std::endl;
        std::cout << GridLogMessage << "\tMatrix "  << MatrixTimer.Elapsed()               << std::endl;
        std::cout << GridLogMessage << "\tLinalg "  << LinalgTimer.Elapsed()               << std::endl;
        std::cout << GridLogMessage << "\tQR "      << QrTimer.Elapsed()                   << std::endl;
        std::cout << GridLogMessage << "\tCompSol " << CompSolutionTimer.Elapsed()         << std::endl;
        return;
      }
    }

    std::cout << GridLogMessage << "FlexibleGeneralisedMinimalResidual did NOT converge" << std::endl;

    if (ErrorOnNoConverge)
      assert(0);
  }

  RealD outerLoopBody(LinearOperatorBase<Field> &LinOp, const Field &src, Field &psi, RealD rsq) {

    RealD cp = 0;

    Field w(src._grid);
    Field r(src._grid);
    Field z(src._grid);

    std::vector<Field> v(RestartLength + 1, src._grid);

    MatrixTimer.Start();
    LinOp.Op(psi, z);
    MatrixTimer.Stop();

    PrecTimer.Start();
    Preconditioner(z, w);
    PrecTimer.Stop();

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

      std::cout << GridLogIterative << "FlexibleGeneralisedMinimalResidual: Iteration " << IterationCount
                << " residual " << cp << " target " << rsq << std::endl;

      if ((i == RestartLength - 1) || (cp <= rsq)) {

        computeSolution(v, psi, i);

        return cp;
      }
    }

    assert(0); // Never reached
    return cp;
  }

  void arnoldiStep(LinearOperatorBase<Field> &LinOp, std::vector<Field> &v, Field &z, Field &w, int iter) {

    MatrixTimer.Start();
    LinOp.Op(v[iter], z);
    MatrixTimer.Stop();

    PrecTimer.Start();
    Preconditioner(z, w);
    PrecTimer.Stop();

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

  void computeSolution(std::vector<Field> const &v, Field &psi, int iter) {

    CompSolutionTimer.Start();
    for (int i = iter; i >= 0; i--) {
      y[i] = gamma[i];
      for (int k = i + 1; k <= iter; k++)
        y[i] = y[i] - H(k, i) * y[k];
      y[i] = y[i] / H(i, i);
    }

    if (true) {
      for (int i = 0; i <= iter; i++)
        psi = psi + v[i] * y[i];
    }
    else {
      psi = y[0] * v[0];
      for (int i = 1; i <= iter; i++)
        psi = psi + v[i] * y[i];
    }
    CompSolutionTimer.Stop();
  }
};
}
#endif
