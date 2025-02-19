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

  virtual void LogIteration(int k,RealD a,RealD b){
    //    std::cout << "ConjugageGradient::LogIteration() "<<std::endl;
  };
  virtual void LogBegin(void){
    std::cout << "ConjugageGradient::LogBegin() "<<std::endl;
  };

    void operator()(LinearOperatorBase<Field> &Linop, const Field &src, Field &psi) {

      this->LogBegin();

      GRID_TRACE("ConjugateGradient");
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
    ssq = norm2(src);
    RealD guess = norm2(psi);
    NormTimer.Stop();
    assert(std::isnan(guess) == 0);
    AssignTimer.Start();
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


template <class Field>
class ConjugateGradientPolynomial : public ConjugateGradient<Field> {
public:
  // Optionally record the CG polynomial
  std::vector<double> ak;
  std::vector<double> bk;
  std::vector<double> poly_p;
  std::vector<double> poly_r;
  std::vector<double> poly_Ap;
  std::vector<double> polynomial;

public:
  ConjugateGradientPolynomial(RealD tol, Integer maxit, bool err_on_no_conv = true)
    : ConjugateGradient<Field>(tol,maxit,err_on_no_conv)
  { };
  void PolyHermOp(LinearOperatorBase<Field> &Linop, const Field &src, Field &psi)
  {
    Field tmp(src.Grid());
    Field AtoN(src.Grid());
    AtoN = src;
    psi=AtoN*polynomial[0];
    for(int n=1;n<polynomial.size();n++){
      tmp = AtoN;
      Linop.HermOp(tmp,AtoN);
      psi = psi + polynomial[n]*AtoN;
    }
  }
  void CGsequenceHermOp(LinearOperatorBase<Field> &Linop, const Field &src, Field &x)
  {
    Field Ap(src.Grid());
    Field r(src.Grid());
    Field p(src.Grid());
    p=src;
    r=src;
    x=Zero();
    x.Checkerboard()=src.Checkerboard();
    for(int k=0;k<ak.size();k++){
      x = x + ak[k]*p;
      Linop.HermOp(p,Ap);
      r = r - ak[k] * Ap;
      p = r + bk[k] * p;
    }
  }
  void Solve(LinearOperatorBase<Field> &Linop, const Field &src, Field &psi)
  {
    psi=Zero();
    this->operator ()(Linop,src,psi);
  }
  virtual void LogBegin(void)
  {
    std::cout << "ConjugageGradientPolynomial::LogBegin() "<<std::endl;
    ak.resize(0);
    bk.resize(0);
    polynomial.resize(0);
    poly_Ap.resize(0);
    poly_Ap.resize(0);
    poly_p.resize(1);
    poly_r.resize(1);
    poly_p[0]=1.0;
    poly_r[0]=1.0;
  };
  virtual void LogIteration(int k,RealD a,RealD b)
  {
    // With zero guess,
    // p = r = src
    //
    // iterate:
    //   x =  x + a p
    //   r =  r - a A p
    //   p =  r + b p
    //
    // [0]
    // r = x
    // p = x
    // Ap=0
    //
    // [1]
    // Ap = A x + 0  ==> shift poly P right by 1 and add 0.
    // x  = x + a p  ==> add polynomials term by term 
    // r  = r - a A p  ==> add polynomials term by term
    // p  = r + b p  ==> add polynomials term by term
    //
    std::cout << "ConjugageGradientPolynomial::LogIteration() "<<k<<std::endl;
    ak.push_back(a);
    bk.push_back(b);
    //  Ap= right_shift(p)
    poly_Ap.resize(k+1);
    poly_Ap[0]=0.0;
    for(int i=0;i<k;i++){
      poly_Ap[i+1]=poly_p[i];
    }

    //  x = x + a p
    polynomial.resize(k);
    polynomial[k-1]=0.0;
    for(int i=0;i<k;i++){
      polynomial[i] = polynomial[i] + a * poly_p[i];
    }
    
    //  r = r - a Ap
    //  p = r + b p
    poly_r.resize(k+1);
    poly_p.resize(k+1);
    poly_r[k] = poly_p[k] = 0.0;
    for(int i=0;i<k+1;i++){
      poly_r[i] = poly_r[i] - a * poly_Ap[i];
      poly_p[i] = poly_r[i] + b * poly_p[i];
    }
  }
};

NAMESPACE_END(Grid);
#endif
