/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/MinimalResidual.h

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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

template <class Field>
class MinimalResidual : public OperatorFunction<Field> {
 public:
  bool ErrorOnNoConverge;  // throw an assert when the MR fails to converge.
                           // Defaults true.
  RealD Tolerance;
  Integer MaxIterations;
  Integer IterationsToComplete; //Number of iterations the MR took to finish. Filled in upon completion
  
  MinimalResidual(RealD tol, Integer maxit, bool err_on_no_conv = true)
      : Tolerance(tol),
        MaxIterations(maxit),
        ErrorOnNoConverge(err_on_no_conv){};

  void operator()(LinearOperatorBase<Field> &Linop, const Field &src,
                  Field &psi) {
    psi.checkerboard = src.checkerboard; // Check
    conformable(psi, src);

    Field p {src};
    Field matrixTimesPsi {src};
    Field r {src};

    RealD alpha {};

    // Initial residual computation & set up
    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    Linop.HermOp(psi, matrixTimesPsi);

    r = src - matrixTimesPsi;

    Linop.HermOp(r, p);

    alpha = innerProduct(p,r) / innerProduct(p,p);
    psi = psi + alpha * r;
    r   = r   - alpha * p;

    Linop.HermOp(r, p);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

    // RealD cp, c, a, d, b, ssq, qq, b_pred;

    Field p(src);
    Field matrixTimesPsi(src);
    // Field r(src);

    // Initial residual computation & set up
    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    
    Linop.HermOpAndNorm(psi, matrixTimesPsi, d, b);
    

    r = src - matrixTimesPsi;
    p = matrixTimesPsi;

    a = norm2(p);
    cp = a;
    ssq = norm2(src);

    std::cout << GridLogIterative << std::setprecision(4)
              << "MinimalResidual: guess " << guess << std::endl;
    std::cout << GridLogIterative << std::setprecision(4)
              << "MinimalResidual:   src " << ssq << std::endl;
    std::cout << GridLogIterative << std::setprecision(4)
              << "MinimalResidual:    mp " << d << std::endl;
    std::cout << GridLogIterative << std::setprecision(4)
              << "MinimalResidual:   matrixTimesPsi " << b << std::endl;
    std::cout << GridLogIterative << std::setprecision(4)
              << "MinimalResidual:  cp,r " << cp << std::endl;
    std::cout << GridLogIterative << std::setprecision(4)
              << "MinimalResidual:     p " << a << std::endl;

    RealD rsq = Tolerance * Tolerance * ssq;

    // Check if guess is really REALLY good :)
    if (cp <= rsq) {
      return;
    }

    std::cout << GridLogIterative << std::setprecision(4)
              << "MinimalResidual: k=0 residual " << cp << " target " << rsq
              << std::endl;

    GridStopWatch LinalgTimer;
    GridStopWatch MatrixTimer;
    GridStopWatch SolverTimer;

    SolverTimer.Start();
    int k;
    for (k = 1; k <= MaxIterations; k++) {
      c = cp;

      MatrixTimer.Start();
      Linop.HermOpAndNorm(p, matrixTimesPsi, d, qq);
      MatrixTimer.Stop();

      LinalgTimer.Start();
      //  RealD    qqck = norm2(matrixTimesPsi);
      //  ComplexD dck  = innerProduct(p,matrixTimesPsi);

      a = c / d;
      b_pred = a * (a * qq - d) / c;

      cp = axpy_norm(r, -a, matrixTimesPsi, r);
      b = cp / c;

      // Fuse these loops ; should be really easy
      psi = a * p + psi;
      p = p * b + r;

      LinalgTimer.Stop();
      std::cout << GridLogIterative << "MinimalResidual: Iteration " << k
                << " residual " << cp << " target " << rsq << std::endl;

      // Stopping condition
      if (cp <= rsq) {
        SolverTimer.Stop();
        Linop.HermOpAndNorm(psi, matrixTimesPsi, d, qq);
        p = matrixTimesPsi - src;

        RealD matrixTimesPsiNorm = sqrt(norm2(matrixTimesPsi));
        RealD psinorm = sqrt(norm2(psi));
        RealD srcnorm = sqrt(norm2(src));
        RealD resnorm = sqrt(norm2(p));
        RealD true_residual = resnorm / srcnorm;

        std::cout << GridLogMessage
                  << "MinimalResidual: Converged on iteration " << k << std::endl;
        std::cout << GridLogMessage << "Computed residual " << sqrt(cp / ssq)
                  << " true residual " << true_residual << " target "
                  << Tolerance << std::endl;
        std::cout << GridLogMessage << "Time elapsed: Iterations "
                  << SolverTimer.Elapsed() << " Matrix  "
                  << MatrixTimer.Elapsed() << " Linalg "
                  << LinalgTimer.Elapsed();
        std::cout << std::endl;

        if (ErrorOnNoConverge) assert(true_residual / Tolerance < 10000.0);
	IterationsToComplete = k;	
        return;
      }
    }
    std::cout << GridLogMessage << "MinimalResidual did NOT converge"
              << std::endl;
    if (ErrorOnNoConverge) assert(0);
    IterationsToComplete = k;
  }
};
}
#endif
