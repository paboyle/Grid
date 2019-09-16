/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/FlexibleCommunicationAvoidingGeneralisedMinimalResidual.h

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
#ifndef GRID_FLEXIBLE_COMMUNICATION_AVOIDING_GENERALISED_MINIMAL_RESIDUAL_H
#define GRID_FLEXIBLE_COMMUNICATION_AVOIDING_GENERALISED_MINIMAL_RESIDUAL_H

namespace Grid {

template<class Field>
class FlexibleCommunicationAvoidingGeneralisedMinimalResidual : public OperatorFunction<Field> {
 public:
  using OperatorFunction<Field>::operator();

  bool ErrorOnNoConverge; // Throw an assert when FCAGMRES fails to converge,
                          // defaults to true

  RealD   Tolerance;

  Integer MaxIterations;
  Integer RestartLength;
  Integer MaxNumberOfRestarts;
  Integer IterationCount; // Number of iterations the FCAGMRES took to finish,
                          // filled in upon completion

  GridStopWatch MatrixTimer;
  GridStopWatch PrecTimer;
  GridStopWatch LinalgTimer;
  GridStopWatch QrTimer;
  GridStopWatch CompSolutionTimer;

  Eigen::MatrixXcd H;

  std::vector<ComplexD> y;
  std::vector<ComplexD> gamma;
  std::vector<ComplexD> c;
  std::vector<ComplexD> s;

  LinearFunction<Field> &Preconditioner;

  FlexibleCommunicationAvoidingGeneralisedMinimalResidual(RealD   tol,
                                                          Integer maxit,
                                                          LinearFunction<Field> &Prec,
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
      , Preconditioner(Prec) {};

  void operator()(LinearOperatorBase<Field> &LinOp, const Field &src, Field &psi) {

    std::cout << GridLogWarning << "This algorithm currently doesn't differ from regular FGMRES" << std::endl;

    psi.Checkerboard() = src.Checkerboard();
    conformable(psi, src);

    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    RealD cp;
    RealD ssq = norm2(src);
    RealD rsq = Tolerance * Tolerance * ssq;

    Field r(src.Grid());

    std::cout << std::setprecision(4) << std::scientific;
    std::cout << GridLogIterative << "FlexibleCommunicationAvoidingGeneralisedMinimalResidual: guess " << guess << std::endl;
    std::cout << GridLogIterative << "FlexibleCommunicationAvoidingGeneralisedMinimalResidual:   src " << ssq   << std::endl;

    PrecTimer.Reset();
    MatrixTimer.Reset();
    LinalgTimer.Reset();
    QrTimer.Reset();
    CompSolutionTimer.Reset();

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

        std::cout << GridLogMessage        << "FlexibleCommunicationAvoidingGeneralisedMinimalResidual: Converged on iteration " << IterationCount
                  << " computed residual " << sqrt(cp / ssq)
                  << " true residual "     << true_residual
                  << " target "            << Tolerance << std::endl;

        std::cout << GridLogMessage << "FCAGMRES Time elapsed: Total   " <<       SolverTimer.Elapsed() << std::endl;
        std::cout << GridLogMessage << "FCAGMRES Time elapsed: Precon  " <<         PrecTimer.Elapsed() << std::endl;
        std::cout << GridLogMessage << "FCAGMRES Time elapsed: Matrix  " <<       MatrixTimer.Elapsed() << std::endl;
        std::cout << GridLogMessage << "FCAGMRES Time elapsed: Linalg  " <<       LinalgTimer.Elapsed() << std::endl;
        std::cout << GridLogMessage << "FCAGMRES Time elapsed: QR      " <<           QrTimer.Elapsed() << std::endl;
        std::cout << GridLogMessage << "FCAGMRES Time elapsed: CompSol " << CompSolutionTimer.Elapsed() << std::endl;
        return;
      }
    }

    std::cout << GridLogMessage << "FlexibleCommunicationAvoidingGeneralisedMinimalResidual did NOT converge" << std::endl;

    if (ErrorOnNoConverge)
      assert(0);
  }

  RealD outerLoopBody(LinearOperatorBase<Field> &LinOp, const Field &src, Field &psi, RealD rsq) {

    RealD cp = 0;

    Field w(src.Grid());
    Field r(src.Grid());

    // these should probably be made class members so that they are only allocated once, not in every restart
    std::vector<Field> v(RestartLength + 1, src.Grid()); for (auto &elem : v) elem = Zero();
    std::vector<Field> z(RestartLength + 1, src.Grid()); for (auto &elem : z) elem = Zero();

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

      cp = norm(gamma[i+1]);

      std::cout << GridLogIterative << "FlexibleCommunicationAvoidingGeneralisedMinimalResidual: Iteration " << IterationCount
                << " residual " << cp << " target " << rsq << std::endl;

      if ((i == RestartLength - 1) || (IterationCount == MaxIterations) || (cp <= rsq)) {

        computeSolution(z, psi, i);

        return cp;
      }
    }

    assert(0); // Never reached
    return cp;
  }

  void arnoldiStep(LinearOperatorBase<Field> &LinOp, std::vector<Field> &v, std::vector<Field> &z, Field &w, int iter) {

    PrecTimer.Start();
    Preconditioner(v[iter], z[iter]);
    PrecTimer.Stop();

    MatrixTimer.Start();
    LinOp.Op(z[iter], w);
    MatrixTimer.Stop();

    LinalgTimer.Start();
    for (int i = 0; i <= iter; ++i) {
      H(iter, i) = innerProduct(v[i], w);
      w = w - ComplexD(H(iter, i)) * v[i];
    }

    H(iter, iter + 1) = sqrt(norm2(w));
    v[iter + 1] = ComplexD(1. / H(iter, iter + 1)) * w;
    LinalgTimer.Stop();
  }

  void qrUpdate(int iter) {

    QrTimer.Start();
    for (int i = 0; i < iter ; ++i) {
      auto tmp       = -s[i] * ComplexD(H(iter, i)) + c[i] * ComplexD(H(iter, i + 1));
      H(iter, i)     = conjugate(c[i]) * ComplexD(H(iter, i)) + conjugate(s[i]) * ComplexD(H(iter, i + 1));
      H(iter, i + 1) = tmp;
    }

    // Compute new Givens Rotation
    auto nu = sqrt(std::norm(H(iter, iter)) + std::norm(H(iter, iter + 1)));
    c[iter]     = H(iter, iter) / nu;
    s[iter]     = H(iter, iter + 1) / nu;

    // Apply new Givens rotation
    H(iter, iter)     = nu;
    H(iter, iter + 1) = 0.;

    gamma[iter + 1] = -s[iter] * gamma[iter];
    gamma[iter]     = conjugate(c[iter]) * gamma[iter];
    QrTimer.Stop();
  }

  void computeSolution(std::vector<Field> const &z, Field &psi, int iter) {

    CompSolutionTimer.Start();
    for (int i = iter; i >= 0; i--) {
      y[i] = gamma[i];
      for (int k = i + 1; k <= iter; k++)
        y[i] = y[i] - ComplexD(H(k, i)) * y[k];
      y[i] = y[i] / ComplexD(H(i, i));
    }

    for (int i = 0; i <= iter; i++)
      psi = psi + z[i] * y[i];
    CompSolutionTimer.Stop();
  }
};
}
#endif
