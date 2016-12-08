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

namespace Grid {

struct CG_state {
  bool do_repro;
  std::vector<RealD> residuals;

  CG_state() {reset();}

  void reset(){
    do_repro = false;
    residuals.clear();
  }
};

/////////////////////////////////////////////////////////////
// Base classes for iterative processes based on operators
// single input vec, single output vec.
/////////////////////////////////////////////////////////////

template <class Field>
class ConjugateGradient : public OperatorFunction<Field> {
 public:
  bool ErrorOnNoConverge;  // throw an assert when the CG fails to converge.
                           // Defaults true.
  RealD Tolerance;
  Integer MaxIterations;

  // Reproducibility controls
  bool ReproTest;
  CG_state CGState; //to check reproducibility by repeating the CG
  ReproducibilityState<typename Field::vector_object> ReprTest; // for the inner proucts

  ConjugateGradient(RealD tol, Integer maxit, bool err_on_no_conv = true,
        bool ReproducibilityTest = false)
      : Tolerance(tol),
        MaxIterations(maxit),
        ErrorOnNoConverge(err_on_no_conv),
        ReproTest(ReproducibilityTest){
        	if(ReproducibilityTest == true && err_on_no_conv == true){
        		std::cout << GridLogMessage << "CG: Reproducibility test ON "<<
        					"and error on convergence ON are incompatible options" << std::endl;
        	exit(1);
        	}
        	
        };


  void operator()(LinearOperatorBase<Field> &Linop, const Field &src,
                  Field &psi) {
    psi.checkerboard = src.checkerboard;
    conformable(psi, src);

    RealD cp, c, a, d, b, ssq, qq, b_pred;

    Field p(src);
    Field mmp(src);
    Field r(src);
    Field psi_start(psi);// save for the repro test

    if (CGState.do_repro)
        std::cout << GridLogMessage << "Starting reproducibility test" << std::endl;

    // Initial residual computation & set up
    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    
    Linop.HermOpAndNorm(psi, mmp, d, b);
    
    if(!ReprTest.do_check)
        ReprTest.reset();
    ReprTest.enable_reprocheck=ReproTest;


    r = src - mmp;
    p = r;

    a = norm2(p);
    cp = a;
    ssq = norm2(src);

    std::cout << GridLogIterative << std::setprecision(4)
              << "ConjugateGradient: guess " << guess << std::endl;
    std::cout << GridLogIterative << std::setprecision(4)
              << "ConjugateGradient:   src " << ssq << std::endl;
    std::cout << GridLogIterative << std::setprecision(4)
              << "ConjugateGradient:    mp " << d << std::endl;
    std::cout << GridLogIterative << std::setprecision(4)
              << "ConjugateGradient:   mmp " << b << std::endl;
    std::cout << GridLogIterative << std::setprecision(4)
              << "ConjugateGradient:  cp,r " << cp << std::endl;
    std::cout << GridLogIterative << std::setprecision(4)
              << "ConjugateGradient:     p " << a << std::endl;

    RealD rsq = Tolerance * Tolerance * ssq;

    // Check if guess is really REALLY good :)
    if (cp <= rsq) {
      return;
    }

    std::cout << GridLogIterative << std::setprecision(4)
              << "ConjugateGradient: k=0 residual " << cp << " target " << rsq
              << std::endl;

    GridStopWatch LinalgTimer;
    GridStopWatch MatrixTimer;
    GridStopWatch SolverTimer;

    SolverTimer.Start();
    int k;
    for (k = 1; k <= MaxIterations; k++) {
      c = cp;// old residual

      MatrixTimer.Start();
      Linop.HermOpAndNorm(p, mmp, d, qq);// mmp = Ap, d=pAp
      MatrixTimer.Stop();

      LinalgTimer.Start();
      //  RealD    qqck = norm2(mmp);
      //  ComplexD dck  = innerProduct(p,mmp);

      a = c / d;
      b_pred = a * (a * qq - d) / c;// a check


      axpy(r, -a, mmp, r);// new residual r = r_old - a * Ap
      cp = norm2(r, ReprTest);// 
      if (ReproTest && !CGState.do_repro) {
        CGState.residuals.push_back(cp);  // save residuals state
                std::cout << GridLogIterative << "ReproTest: Saving state" << std::endl;
        }
      if (ReproTest && CGState.do_repro){
        // check that the residual agrees with the previous run
        std::cout << GridLogIterative << "ReproTest: Checking state k=" << k << std::endl;
        if (cp != CGState.residuals[k-1]){
                std::cout << GridLogMessage << "Failing reproducibility test";
                std::cout << GridLogMessage << " at k=" << k << std::endl;
                std::cout << GridLogMessage << "saved residual = " << CGState.residuals[k-1] 
                        << " cp = " << cp << std::endl;
                exit(-1);
        }
      }
      b = cp / c;

      // Fuse these loops ; should be really easy
      psi = a * p + psi; // update solution
      p = p * b + r;  // update search direction

      LinalgTimer.Stop();
      std::cout << GridLogIterative << "ConjugateGradient: Iteration " << k
                << " residual " << cp << " target " << rsq << std::endl;

      // Stopping condition
      if (cp <= rsq) {
        SolverTimer.Stop();
        Linop.HermOpAndNorm(psi, mmp, d, qq);
        p = mmp - src;

        RealD mmpnorm = sqrt(norm2(mmp));
        RealD psinorm = sqrt(norm2(psi));
        RealD srcnorm = sqrt(norm2(src));
        RealD resnorm = sqrt(norm2(p));
        RealD true_residual = resnorm / srcnorm;

        std::cout << GridLogMessage
                  << "ConjugateGradient: Converged on iteration " << k << std::endl;
        std::cout << GridLogMessage << "Computed residual " << sqrt(cp / ssq)
                  << " true residual " << true_residual << " target "
                  << Tolerance << std::endl;
        std::cout << GridLogMessage << "Time elapsed: Iterations "
                  << SolverTimer.Elapsed() << " Matrix  "
                  << MatrixTimer.Elapsed() << " Linalg "
                  << LinalgTimer.Elapsed();
        std::cout << std::endl;

        if (ErrorOnNoConverge) assert(true_residual / Tolerance < 10000.0);

        if (!CGState.do_repro && ReproTest){
                CGState.do_repro = true;
                ReprTest.do_check = true;
                ReprTest.reset_counter();
                this->operator()(Linop, src, psi_start);// run the repro test
                if (ReprTest.success)
                	std::cout << GridLogMessage << "Reproducibility test passed" << std::endl;
                else{
                	std::cout << GridLogMessage << "Reproducibility test failed" << std::endl;
                	exit(1);
                }
        }

        // Clear state
        CGState.reset();
        ReprTest.reset();
        return;
      }
    }
    std::cout << GridLogMessage << "ConjugateGradient did NOT converge"
              << std::endl;
    if (ErrorOnNoConverge) assert(0);
  }
};
}
#endif
