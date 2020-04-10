/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/ConjugateGradient.h

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
#ifndef GRID_CONJUGATE_GRADIENT_H
#define GRID_CONJUGATE_GRADIENT_H

NAMESPACE_BEGIN(Grid);

/////////////////////////////////////////////////////////////
// Base classes for iterative processes based on operators
// single input vec, single output vec.
/////////////////////////////////////////////////////////////

template <class Field>
class ConjugateGradient : public OperatorFunction<Field> {
public:

  using OperatorFunction<Field>::operator();

  bool ErrorOnNoConverge;  // throw an assert when the CG fails to converge.
                           // Defaults true.
  RealD Tolerance;
  Integer MaxIterations;
  Integer IterationsToComplete; //Number of iterations the CG took to finish. Filled in upon completion
  
  ConjugateGradient(RealD tol, Integer maxit, bool err_on_no_conv = true)
    : Tolerance(tol),
      MaxIterations(maxit),
      ErrorOnNoConverge(err_on_no_conv){};

  void operator()(LinearOperatorBase<Field> &Linop, const Field &src, Field &psi) {

    psi.Checkerboard() = src.Checkerboard();

    conformable(psi, src);

    RealD cp, c, a, d, b, ssq, qq;
    //RealD b_pred;

    Field p(src);
    Field mmp(src);
    Field r(src);

    // Initial residual computation & set up
    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);
    
    Linop.HermOpAndNorm(psi, mmp, d, b);
    
    r = src - mmp;
    p = r;

    a = norm2(p);
    cp = a;
    ssq = norm2(src);

    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient: guess " << guess << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient:   src " << ssq << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient:    mp " << d << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient:   mmp " << b << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient:  cp,r " << cp << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient:     p " << a << std::endl;

    RealD rsq = Tolerance * Tolerance * ssq;

    // Check if guess is really REALLY good :)
    if (cp <= rsq) {
      std::cout << GridLogMessage << "ConjugateGradient guess is converged already " << std::endl;
      IterationsToComplete = 0;	
      return;
    }

    std::cout << GridLogIterative << std::setprecision(8)
              << "ConjugateGradient: k=0 residual " << cp << " target " << rsq << std::endl;

    GridStopWatch LinalgTimer;
    GridStopWatch InnerTimer;
    GridStopWatch AxpyNormTimer;
    GridStopWatch LinearCombTimer;
    GridStopWatch MatrixTimer;
    GridStopWatch SolverTimer;

    SolverTimer.Start();
    int k;
    for (k = 1; k <= MaxIterations; k++) {
      c = cp;

      MatrixTimer.Start();
      Linop.HermOp(p, mmp);
      MatrixTimer.Stop();

      LinalgTimer.Start();

      InnerTimer.Start();
      ComplexD dc  = innerProduct(p,mmp);
      InnerTimer.Stop();
      d = dc.real();
      a = c / d;

      AxpyNormTimer.Start();
      cp = axpy_norm(r, -a, mmp, r);
      AxpyNormTimer.Stop();
      b = cp / c;

      LinearCombTimer.Start();
      auto psi_v = psi.View();
      auto p_v   = p.View();
      auto r_v   = r.View();
      accelerator_for(ss,p_v.size(), Field::vector_object::Nsimd(),{
	  coalescedWrite(psi_v[ss], a      *  p_v(ss) + psi_v(ss));
	  coalescedWrite(p_v[ss]  , b      *  p_v(ss) + r_v  (ss));
      });
      LinearCombTimer.Stop();
      LinalgTimer.Stop();

      std::cout << GridLogIterative << "ConjugateGradient: Iteration " << k
                << " residual^2 " << sqrt(cp/ssq) << " target " << Tolerance << std::endl;

      // Stopping condition
      if (cp <= rsq) {
        SolverTimer.Stop();
        Linop.HermOpAndNorm(psi, mmp, d, qq);
        p = mmp - src;

        RealD srcnorm = std::sqrt(norm2(src));
        RealD resnorm = std::sqrt(norm2(p));
        RealD true_residual = resnorm / srcnorm;

        std::cout << GridLogMessage << "ConjugateGradient Converged on iteration " << k 
		  << "\tComputed residual " << std::sqrt(cp / ssq)
		  << "\tTrue residual " << true_residual
		  << "\tTarget " << Tolerance << std::endl;

        std::cout << GridLogIterative << "Time breakdown "<<std::endl;
	std::cout << GridLogIterative << "\tElapsed    " << SolverTimer.Elapsed() <<std::endl;
	std::cout << GridLogIterative << "\tMatrix     " << MatrixTimer.Elapsed() <<std::endl;
	std::cout << GridLogIterative << "\tLinalg     " << LinalgTimer.Elapsed() <<std::endl;
	std::cout << GridLogIterative << "\tInner      " << InnerTimer.Elapsed() <<std::endl;
	std::cout << GridLogIterative << "\tAxpyNorm   " << AxpyNormTimer.Elapsed() <<std::endl;
	std::cout << GridLogIterative << "\tLinearComb " << LinearCombTimer.Elapsed() <<std::endl;

        if (ErrorOnNoConverge) assert(true_residual / Tolerance < 10000.0);

	IterationsToComplete = k;	

        return;
      }
    }
    std::cout << GridLogMessage << "ConjugateGradient did NOT converge "<<k<<" / "<< MaxIterations<< std::endl;

    if (ErrorOnNoConverge) assert(0);
    IterationsToComplete = k;

  }
};
NAMESPACE_END(Grid);
#endif
