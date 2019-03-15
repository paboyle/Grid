    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/ConjugateGradientReliableUpdate.h

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@phys.columbia.edu>

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_CONJUGATE_GRADIENT_RELIABLE_UPDATE_H
#define GRID_CONJUGATE_GRADIENT_RELIABLE_UPDATE_H

namespace Grid {

  template<class FieldD,class FieldF, typename std::enable_if< getPrecision<FieldD>::value == 2, int>::type = 0,typename std::enable_if< getPrecision<FieldF>::value == 1, int>::type = 0> 
  class ConjugateGradientReliableUpdate : public LinearFunction<FieldD> {
  public:
    bool ErrorOnNoConverge;  // throw an assert when the CG fails to converge.
    // Defaults true.
    RealD Tolerance;
    Integer MaxIterations;
    Integer IterationsToComplete; //Number of iterations the CG took to finish. Filled in upon completion
    Integer ReliableUpdatesPerformed;

    bool DoFinalCleanup; //Final DP cleanup, defaults to true
    Integer IterationsToCleanup; //Final DP cleanup step iterations
    
    LinearOperatorBase<FieldF> &Linop_f;
    LinearOperatorBase<FieldD> &Linop_d;
    GridBase* SinglePrecGrid;
    RealD Delta; //reliable update parameter

    //Optional ability to switch to a different linear operator once the tolerance reaches a certain point. Useful for single/half -> single/single
    LinearOperatorBase<FieldF> *Linop_fallback;
    RealD fallback_transition_tol;

    
    ConjugateGradientReliableUpdate(RealD tol, Integer maxit, RealD _delta, GridBase* _sp_grid, LinearOperatorBase<FieldF> &_Linop_f, LinearOperatorBase<FieldD> &_Linop_d, bool err_on_no_conv = true)
      : Tolerance(tol),
        MaxIterations(maxit),
	Delta(_delta),
	Linop_f(_Linop_f),
	Linop_d(_Linop_d),
	SinglePrecGrid(_sp_grid),
        ErrorOnNoConverge(err_on_no_conv),
	DoFinalCleanup(true),
	Linop_fallback(NULL)
    {};

    void setFallbackLinop(LinearOperatorBase<FieldF> &_Linop_fallback, const RealD _fallback_transition_tol){
      Linop_fallback = &_Linop_fallback;
      fallback_transition_tol = _fallback_transition_tol;      
    }
    
    void operator()(const FieldD &src, FieldD &psi) {
      LinearOperatorBase<FieldF> *Linop_f_use = &Linop_f;
      bool using_fallback = false;
      
      psi.checkerboard = src.checkerboard;
      conformable(psi, src);

      RealD cp, c, a, d, b, ssq, qq, b_pred;

      FieldD p(src);
      FieldD mmp(src);
      FieldD r(src);

      // Initial residual computation & set up
      RealD guess = norm2(psi);
      assert(std::isnan(guess) == 0);
    
      Linop_d.HermOpAndNorm(psi, mmp, d, b);
    
      r = src - mmp;
      p = r;

      a = norm2(p);
      cp = a;
      ssq = norm2(src);

      std::cout << GridLogIterative << std::setprecision(4) << "ConjugateGradientReliableUpdate: guess " << guess << std::endl;
      std::cout << GridLogIterative << std::setprecision(4) << "ConjugateGradientReliableUpdate:   src " << ssq << std::endl;
      std::cout << GridLogIterative << std::setprecision(4) << "ConjugateGradientReliableUpdate:    mp " << d << std::endl;
      std::cout << GridLogIterative << std::setprecision(4) << "ConjugateGradientReliableUpdate:   mmp " << b << std::endl;
      std::cout << GridLogIterative << std::setprecision(4) << "ConjugateGradientReliableUpdate:  cp,r " << cp << std::endl;
      std::cout << GridLogIterative << std::setprecision(4) << "ConjugateGradientReliableUpdate:     p " << a << std::endl;

      RealD rsq = Tolerance * Tolerance * ssq;

      // Check if guess is really REALLY good :)
      if (cp <= rsq) {
	std::cout << GridLogMessage << "ConjugateGradientReliableUpdate guess was REALLY good\n";
	std::cout << GridLogMessage << "\tComputed residual " << sqrt(cp / ssq)<<std::endl;
	return;
      }

      //Single prec initialization
      FieldF r_f(SinglePrecGrid);
      r_f.checkerboard = r.checkerboard;
      precisionChange(r_f, r);

      FieldF psi_f(r_f);
      psi_f = zero;

      FieldF p_f(r_f);
      FieldF mmp_f(r_f);

      RealD MaxResidSinceLastRelUp = cp; //initial residual    
    
      std::cout << GridLogIterative << std::setprecision(4)
		<< "ConjugateGradient: k=0 residual " << cp << " target " << rsq << std::endl;

      GridStopWatch LinalgTimer;
      GridStopWatch MatrixTimer;
      GridStopWatch SolverTimer;

      SolverTimer.Start();
      int k = 0;
      int l = 0;
    
      for (k = 1; k <= MaxIterations; k++) {
	c = cp;

	MatrixTimer.Start();
	Linop_f_use->HermOpAndNorm(p_f, mmp_f, d, qq);
	MatrixTimer.Stop();

	LinalgTimer.Start();

	a = c / d;
	b_pred = a * (a * qq - d) / c;

	cp = axpy_norm(r_f, -a, mmp_f, r_f);
	b = cp / c;

	// Fuse these loops ; should be really easy
	psi_f = a * p_f + psi_f;
	//p_f = p_f * b + r_f;

	LinalgTimer.Stop();

	std::cout << GridLogIterative << "ConjugateGradientReliableUpdate: Iteration " << k
		  << " residual " << cp << " target " << rsq << std::endl;
	std::cout << GridLogDebug << "a = "<< a << " b_pred = "<< b_pred << "  b = "<< b << std::endl;
	std::cout << GridLogDebug << "qq = "<< qq << " d = "<< d << "  c = "<< c << std::endl;

	if(cp > MaxResidSinceLastRelUp){
	  std::cout << GridLogIterative << "ConjugateGradientReliableUpdate: updating MaxResidSinceLastRelUp : " << MaxResidSinceLastRelUp << " -> " << cp << std::endl;
	  MaxResidSinceLastRelUp = cp;
	}
	  
	// Stopping condition
	if (cp <= rsq) {
	  //Although not written in the paper, I assume that I have to add on the final solution
	  precisionChange(mmp, psi_f);
	  psi = psi + mmp;
	
	
	  SolverTimer.Stop();
	  Linop_d.HermOpAndNorm(psi, mmp, d, qq);
	  p = mmp - src;

	  RealD srcnorm = sqrt(norm2(src));
	  RealD resnorm = sqrt(norm2(p));
	  RealD true_residual = resnorm / srcnorm;

	  std::cout << GridLogMessage << "ConjugateGradientReliableUpdate Converged on iteration " << k << " after " << l << " reliable updates" << std::endl;
	  std::cout << GridLogMessage << "\tComputed residual " << sqrt(cp / ssq)<<std::endl;
	  std::cout << GridLogMessage << "\tTrue residual " << true_residual<<std::endl;
	  std::cout << GridLogMessage << "\tTarget " << Tolerance << std::endl;

	  std::cout << GridLogMessage << "Time breakdown "<<std::endl;
	  std::cout << GridLogMessage << "\tElapsed    " << SolverTimer.Elapsed() <<std::endl;
	  std::cout << GridLogMessage << "\tMatrix     " << MatrixTimer.Elapsed() <<std::endl;
	  std::cout << GridLogMessage << "\tLinalg     " << LinalgTimer.Elapsed() <<std::endl;

	  IterationsToComplete = k;	
	  ReliableUpdatesPerformed = l;
	  
	  if(DoFinalCleanup){
	    //Do a final CG to cleanup
	    std::cout << GridLogMessage << "ConjugateGradientReliableUpdate performing final cleanup.\n";
	    ConjugateGradient<FieldD> CG(Tolerance,MaxIterations);
	    CG.ErrorOnNoConverge = ErrorOnNoConverge;
	    CG(Linop_d,src,psi);
	    IterationsToCleanup = CG.IterationsToComplete;
	  }
	  else if (ErrorOnNoConverge) assert(true_residual / Tolerance < 10000.0);

	  std::cout << GridLogMessage << "ConjugateGradientReliableUpdate complete.\n";
	  return;
	}
	else if(cp < Delta * MaxResidSinceLastRelUp) { //reliable update
	  std::cout << GridLogMessage << "ConjugateGradientReliableUpdate "
		    << cp << "(residual) < " << Delta << "(Delta) * " << MaxResidSinceLastRelUp << "(MaxResidSinceLastRelUp) on iteration " << k << " : performing reliable update\n";
	  precisionChange(mmp, psi_f);
	  psi = psi + mmp;

	  Linop_d.HermOpAndNorm(psi, mmp, d, qq);
	  r = src - mmp;

	  psi_f = zero;
	  precisionChange(r_f, r);
	  cp = norm2(r);
	  MaxResidSinceLastRelUp = cp;

	  b = cp/c;
	  
	  std::cout << GridLogMessage << "ConjugateGradientReliableUpdate new residual " << cp << std::endl;
	  
	  l = l+1;
	}

	p_f = p_f * b + r_f; //update search vector after reliable update appears to help convergence

	if(!using_fallback && Linop_fallback != NULL && cp < fallback_transition_tol){
	  std::cout << GridLogMessage << "ConjugateGradientReliableUpdate switching to fallback linear operator on iteration " << k << " at residual " << cp << std::endl;
	  Linop_f_use = Linop_fallback;
	  using_fallback = true;
	}

	
      }
      std::cout << GridLogMessage << "ConjugateGradientReliableUpdate did NOT converge"
		<< std::endl;
      
      if (ErrorOnNoConverge) assert(0);
      IterationsToComplete = k;
      ReliableUpdatesPerformed = l;      
    }    
  };


};



#endif
