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

/////////////////////////////////////////////////////////////
// Base classes for iterative processes based on operators
// single input vec, single output vec.
/////////////////////////////////////////////////////////////

template<class Field> class MinimalResidual : public OperatorFunction<Field> {
 public:
  bool ErrorOnNoConverge; // throw an assert when the MR fails to converge.
                          // Defaults true.
  RealD   Tolerance;
  Integer MaxIterations;
  Integer IterationsToComplete; // Number of iterations the MR took to finish.
                                // Filled in upon completion

  MinimalResidual(RealD tol, Integer maxit, bool err_on_no_conv = true)
    : Tolerance(tol), MaxIterations(maxit), ErrorOnNoConverge(err_on_no_conv){};

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

    RealD ssq = norm2(src); // flopcount.addSiteFlops(4*Nc*Ns,s);
    RealD rsq = Tolerance * Tolerance * ssq; // flopcount.addSiteFlops(4*Nc*Ns,s);

    Linop.Op(psi, Mr); // flopcount.addFlops(M.nFlops());

    r = src - Mr; // flopcount.addSiteFlops(2*Nc*Ns,s);

    RealD cp = norm2(r); //  Cp = |r[0]|^2 // 2 Nc Ns  flops // flopcount.addSiteFlops(4*Nc*Ns, s);

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
    for (k = 1; k <= MaxIterations; k++) { //  a[k-1] := < M.r[k-1], r[k-1] >/ < M.r[k-1], M.r[k-1] >

      MatrixTimer.Start();
      Linop.Op(r, Mr); //  Mr = M * r // flopcount.addFlops(M.nFlops());
      MatrixTimer.Stop();

      LinalgTimer.Start();

      c = innerProduct(Mr, r); //  c = < M.r, r > // // flopcount.addSiteFlops(4*Nc*Ns,s);

      d = norm2(Mr); //  d = | M.r | ** 2  // // flopcount.addSiteFlops(4*Nc*Ns,s);

      a = c / d;

      // a = a * MRovpar; //  a[k-1] *= MRovpar // from chroma code. TODO: check what to do with this

      psi = psi + r * a; //  Psi[k] += a[k-1] r[k-1] ; // flopcount.addSiteFlops(4*Nc*Ns,s);

      r = r - Mr * a; //  r[k] -= a[k-1] M . r[k-1] ; // flopcount.addSiteFlops(4*Nc*Ns,s);

      cp = norm2(r); //  cp  =  | r[k] |**2 // flopcount.addSiteFlops(4*Nc*Ns,s);

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
