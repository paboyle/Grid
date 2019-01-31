/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/MinimalResidual.h

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
#ifndef GRID_MINIMAL_RESIDUAL_H
#define GRID_MINIMAL_RESIDUAL_H

namespace Grid {

template<class Field> class MinimalResidual : public OperatorFunction<Field> {
 public:
  bool ErrorOnNoConverge; // throw an assert when the MR fails to converge.
                          // Defaults true.
  RealD   Tolerance;
  Integer MaxIterations;
  RealD   overRelaxParam;
  Integer IterationsToComplete; // Number of iterations the MR took to finish.
                                // Filled in upon completion

  MinimalResidual(RealD tol, Integer maxit, Real ovrelparam = 1.0, bool err_on_no_conv = true)
    : Tolerance(tol), MaxIterations(maxit), overRelaxParam(ovrelparam), ErrorOnNoConverge(err_on_no_conv){};

  void operator()(LinearOperatorBase<Field> &Linop, const Field &src, Field &psi) {

    psi.checkerboard = src.checkerboard;
    conformable(psi, src);

    Complex a, c;
    Real    d;

    Field Mr(src);
    Field r(src);

    // Initial residual computation & set up
    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    RealD ssq = norm2(src);
    RealD rsq = Tolerance * Tolerance * ssq;

    Linop.Op(psi, Mr);

    r = src - Mr;

    RealD cp = norm2(r);

    std::cout << std::setprecision(4) << std::scientific;
    std::cout << GridLogIterative << "MinimalResidual: guess " << guess << std::endl;
    std::cout << GridLogIterative << "MinimalResidual:   src " << ssq << std::endl;
    std::cout << GridLogIterative << "MinimalResidual:    mp " << d << std::endl;
    std::cout << GridLogIterative << "MinimalResidual:  cp,r " << cp << std::endl;

    if (cp <= rsq) {
      return;
    }

    std::cout << GridLogIterative << "MinimalResidual: k=0 residual " << cp << " target " << rsq << std::endl;

    GridStopWatch LinalgTimer;
    GridStopWatch MatrixTimer;
    GridStopWatch SolverTimer;

    SolverTimer.Start();
    int k;
    for (k = 1; k <= MaxIterations; k++) {

      MatrixTimer.Start();
      Linop.Op(r, Mr);
      MatrixTimer.Stop();

      LinalgTimer.Start();

      c = innerProduct(Mr, r);

      d = norm2(Mr);

      a = c / d;

      a = a * overRelaxParam;

      psi = psi + r * a;

      r = r - Mr * a;

      cp = norm2(r);

      LinalgTimer.Stop();

      std::cout << GridLogIterative << "MinimalResidual: Iteration " << k
                << " residual " << cp << " target " << rsq << std::endl;
      std::cout << GridLogDebug << "a = " << a << " c = " << c << " d = " << d << std::endl;

      // Stopping condition
      if (cp <= rsq) {
        SolverTimer.Stop();

        Linop.Op(psi, Mr);
        r = src - Mr;

        RealD srcnorm       = sqrt(ssq);
        RealD resnorm       = sqrt(norm2(r));
        RealD true_residual = resnorm / srcnorm;

        std::cout << GridLogMessage        << "MinimalResidual Converged on iteration " << k
                  << " computed residual " << sqrt(cp / ssq)
                  << " true residual "     << true_residual
                  << " target "            << Tolerance << std::endl;

        std::cout << GridLogMessage << "MR Time elapsed: Total   " << SolverTimer.Elapsed() << std::endl;
        std::cout << GridLogMessage << "MR Time elapsed: Matrix  " << MatrixTimer.Elapsed() << std::endl;
        std::cout << GridLogMessage << "MR Time elapsed: Linalg  " << LinalgTimer.Elapsed() << std::endl;

        if (ErrorOnNoConverge)
          assert(true_residual / Tolerance < 10000.0);

        IterationsToComplete = k;

        return;
      }
    }

    std::cout << GridLogMessage << "MinimalResidual did NOT converge"
              << std::endl;

    if (ErrorOnNoConverge)
      assert(0);

    IterationsToComplete = k;
  }
};
} // namespace Grid
#endif
