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
  RealD TrueResidual;
  
  ConjugateGradient(RealD tol, Integer maxit, bool err_on_no_conv = true)
    : Tolerance(tol),
      MaxIterations(maxit),
      ErrorOnNoConverge(err_on_no_conv)
  {};

  void operator()(LinearOperatorBase<Field> &Linop, const Field &src, Field &psi) {

    GRID_TRACE("ConjugateGradient");
    GridStopWatch PreambleTimer;
    PreambleTimer.Start();
    psi.Checkerboard() = src.Checkerboard();

    conformable(psi, src);

    RealD cp, c, a, d, b, ssq, qq;
    //RealD b_pred;

    // Was doing copies
    Field p(src.Grid());
    Field mmp(src.Grid());
    Field r(src.Grid());

    // Initial residual computation & set up
    ssq = norm2(src);
    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);
    if ( guess == 0.0 ) {
      r = src;
      p = r;
      a = ssq;
    } else { 
      Linop.HermOpAndNorm(psi, mmp, d, b);
      r = src - mmp;
      p = r;
      a = norm2(p);
    }
    cp = a;

    // Handle trivial case of zero src
    if (ssq == 0.){
      psi = Zero();
      IterationsToComplete = 1;
      TrueResidual = 0.;
      return;
    }

    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient: guess " << guess << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient:   src " << ssq << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient:    mp " << d << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient:   mmp " << b << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient:  cp,r " << cp << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "ConjugateGradient:     p " << a << std::endl;

    RealD rsq = Tolerance * Tolerance * ssq;

    // Check if guess is really REALLY good :)
    if (cp <= rsq) {
      TrueResidual = std::sqrt(a/ssq);
      std::cout << GridLogMessage << "ConjugateGradient guess is converged already " << std::endl;
      IterationsToComplete = 0;	
      return;
    }

    std::cout << GridLogIterative << std::setprecision(8)
              << "ConjugateGradient: k=0 residual " << cp << " target " << rsq << std::endl;

    PreambleTimer.Stop();
    GridStopWatch LinalgTimer;
    GridStopWatch InnerTimer;
    GridStopWatch AxpyNormTimer;
    GridStopWatch LinearCombTimer;
    GridStopWatch MatrixTimer;
    GridStopWatch SolverTimer;

    RealD usecs = -usecond();
    SolverTimer.Start();
    int k;
    for (k = 1; k <= MaxIterations; k++) {

      GridStopWatch IterationTimer;
      IterationTimer.Start();
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
      {
	autoView( psi_v , psi, AcceleratorWrite);
	autoView( p_v   , p,   AcceleratorWrite);
	autoView( r_v   , r,   AcceleratorWrite);
	accelerator_for(ss,p_v.size(), Field::vector_object::Nsimd(),{
	    coalescedWrite(psi_v[ss], a      *  p_v(ss) + psi_v(ss));
	    coalescedWrite(p_v[ss]  , b      *  p_v(ss) + r_v  (ss));
	});
      }
      LinearCombTimer.Stop();
      LinalgTimer.Stop();

      IterationTimer.Stop();
      if ( (k % 500) == 0 ) {
	std::cout << GridLogMessage << "ConjugateGradient: Iteration " << k
                << " residual " << sqrt(cp/ssq) << " target " << Tolerance << std::endl;
      } else { 
	std::cout << GridLogIterative << "ConjugateGradient: Iteration " << k
		  << " residual " << sqrt(cp/ssq) << " target " << Tolerance << " took " << IterationTimer.Elapsed() << std::endl;
      }

      // Stopping condition
      if (cp <= rsq) {
	usecs +=usecond();
        SolverTimer.Stop();
        Linop.HermOpAndNorm(psi, mmp, d, qq);
        p = mmp - src;
	GridBase *grid = src.Grid();
	RealD DwfFlops = (1452. )*grid->gSites()*4*k
   	               + (8+4+8+4+4)*12*grid->gSites()*k; // CG linear algebra
        RealD srcnorm = std::sqrt(norm2(src));
        RealD resnorm = std::sqrt(norm2(p));
        RealD true_residual = resnorm / srcnorm;
        std::cout << GridLogMessage << "ConjugateGradient Converged on iteration " << k 
		  << "\tComputed residual " << std::sqrt(cp / ssq)
		  << "\tTrue residual " << true_residual
		  << "\tTarget " << Tolerance << std::endl;

	//	std::cout << GridLogMessage << "\tPreamble   " << PreambleTimer.Elapsed() <<std::endl;
	std::cout << GridLogMessage << "\tSolver Elapsed    " << SolverTimer.Elapsed() <<std::endl;
        std::cout << GridLogPerformance << "Time breakdown "<<std::endl;
	std::cout << GridLogPerformance << "\tMatrix     " << MatrixTimer.Elapsed() <<std::endl;
	std::cout << GridLogPerformance << "\tLinalg     " << LinalgTimer.Elapsed() <<std::endl;
	std::cout << GridLogPerformance << "\t\tInner      " << InnerTimer.Elapsed() <<std::endl;
	std::cout << GridLogPerformance << "\t\tAxpyNorm   " << AxpyNormTimer.Elapsed() <<std::endl;
	std::cout << GridLogPerformance << "\t\tLinearComb " << LinearCombTimer.Elapsed() <<std::endl;

	std::cout << GridLogDebug << "\tMobius flop rate " << DwfFlops/ usecs<< " Gflops " <<std::endl;

        if (ErrorOnNoConverge) assert(true_residual / Tolerance < 10000.0);

	IterationsToComplete = k;	
	TrueResidual = true_residual;

        return;
      }
    }
    // Failed. Calculate true residual before giving up                                                         
    // Linop.HermOpAndNorm(psi, mmp, d, qq);
    //    p = mmp - src;
    //TrueResidual = sqrt(norm2(p)/ssq);
    //    TrueResidual = 1;

    std::cout << GridLogMessage << "ConjugateGradient did NOT converge "<<k<<" / "<< MaxIterations
    	      <<" residual "<< std::sqrt(cp / ssq)<< std::endl;
    SolverTimer.Stop();
    std::cout << GridLogMessage << "\tPreamble   " << PreambleTimer.Elapsed() <<std::endl;
    std::cout << GridLogMessage << "\tSolver     " << SolverTimer.Elapsed() <<std::endl;
    std::cout << GridLogMessage << "Solver breakdown "<<std::endl;
    std::cout << GridLogMessage << "\tMatrix     " << MatrixTimer.Elapsed() <<std::endl;
    std::cout << GridLogMessage<< "\tLinalg     " << LinalgTimer.Elapsed() <<std::endl;
    std::cout << GridLogPerformance << "\t\tInner      " << InnerTimer.Elapsed() <<std::endl;
    std::cout << GridLogPerformance << "\t\tAxpyNorm   " << AxpyNormTimer.Elapsed() <<std::endl;
    std::cout << GridLogPerformance << "\t\tLinearComb " << LinearCombTimer.Elapsed() <<std::endl;

    if (ErrorOnNoConverge) assert(0);
    IterationsToComplete = k;

  }
};
NAMESPACE_END(Grid);
#endif
