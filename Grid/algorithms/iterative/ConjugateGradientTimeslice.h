/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/ConjugateGradientTimeslice.h

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
#ifndef GRID_CONJUGATE_GRADIENT_TIMESLICE_H
#define GRID_CONJUGATE_GRADIENT_TIMESLICE_H

NAMESPACE_BEGIN(Grid);

/////////////////////////////////////////////////////////////
// Base classes for iterative processes based on operators
// single input vec, single output vec.
/////////////////////////////////////////////////////////////

/**
 * Simple modification of conjugate gradient that outputs the residual as a function 
 * of time, in order to study the large wavelength behavior of the solver. 
 */


template <class Field>
class ConjugateGradientTimeslice : public OperatorFunction<Field> {
public:

  using OperatorFunction<Field>::operator();
  
  bool ErrorOnNoConverge;  // throw an assert when the CG fails to converge.
                           // Defaults true.
  RealD Tolerance;
  Integer MaxIterations;
  Integer IterationsToComplete; //Number of iterations the CG took to finish. Filled in upon completion
  RealD TrueResidual;
  
  ConjugateGradientTimeslice(RealD tol, Integer maxit, bool err_on_no_conv = true)
    : Tolerance(tol),
      MaxIterations(maxit),
      ErrorOnNoConverge(err_on_no_conv)
  {};

  virtual void LogIteration(int k,RealD a,RealD b){
    //    std::cout << "ConjugageGradient::LogIteration() "<<std::endl;
  };
  virtual void LogBegin(void){
    std::cout << "ConjugageGradient::LogBegin() "<<std::endl;
  };

    void operator()(LinearOperatorBase<Field> &Linop, const Field &src, Field &psi) {

      this->LogBegin();

      GRID_TRACE("ConjugateGradientTimeslice");
    GridStopWatch PreambleTimer;
    GridStopWatch ConstructTimer;
    GridStopWatch NormTimer;
    GridStopWatch AssignTimer;
    PreambleTimer.Start();
    psi.Checkerboard() = src.Checkerboard();

    conformable(psi, src);

    RealD cp, c, a, d, b, ssq, qq;
    //RealD b_pred;

    // Was doing copies
    ConstructTimer.Start();
    Field p  (src.Grid());
    Field mmp(src.Grid());
    Field r  (src.Grid());
    ConstructTimer.Stop();

    // Initial residual computation & set up
    NormTimer.Start();
    ssq = norm2(src);                 // Norm of source vector ||b||^2

    ssqtx = localNorm2(src);          // Norm |b(x, t)|^2 as a field
    std::vector<RealD> ssqt;          // Norm of source not summed over time slices, ssq(t) = \sum_x |b(x, t)|^2
    sliceSum(ssqtx, ssqt, Tdir);      // TODO make sure Tdir is globally defined

    RealD guess = norm2(psi);         // Norm of initial guess ||psi||^2
    NormTimer.Stop();
    assert(std::isnan(guess) == 0);
    AssignTimer.Start();
    if ( guess == 0.0 ) {
      r = src;
      p = r;
      a = ssq;
    } else { 
      Linop.HermOpAndNorm(psi, mmp, d, b);        // 
      r = src - mmp;      // Initial residual r0 = b - A guess
      p = r;              // initial conj vector p0 = r0
      a = norm2(p);
    }
    cp = a;
    AssignTimer.Stop();

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
      Linop.HermOp(p, mmp);         // Computes mmp = Ap
      MatrixTimer.Stop();

      LinalgTimer.Start();

      InnerTimer.Start();
      ComplexD dc  = innerProduct(p,mmp);         // p^\dagger A p
      InnerTimer.Stop();
      d = dc.real();
      a = c / d;

      // What is axpy? Some accelerator or something? Check Lattice_arith.h
      AxpyNormTimer.Start();

      // axpy_norm computes ax+by for vectors x and y compatible with a GPU. Here b is set to 1 (see the function in Lattice_reduction.h). 
      // The first argument passes r by reference, so it stores r --> -a * Ap + 1 * r, i.e. it performs an update on 
      // r_k --> r_{k+1} = r_k - \alpha_k A p_k. The function returns the norm squared of the first variable, i.e. ||r_{k+1}||^2.
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
      LogIteration(k,a,b);

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

        // GridLogMessage logs the message to the terminal output; GridLogPerformance probably writes to a log file?
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
    std::cout << GridLogMessage << "\tConstruct  " << ConstructTimer.Elapsed() <<std::endl;
    std::cout << GridLogMessage << "\tNorm       " << NormTimer.Elapsed() <<std::endl;
    std::cout << GridLogMessage << "\tAssign     " << AssignTimer.Elapsed() <<std::endl;
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
