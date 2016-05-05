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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_CONJUGATE_GRADIENT_H
#define GRID_CONJUGATE_GRADIENT_H

namespace Grid {

    /////////////////////////////////////////////////////////////
    // Base classes for iterative processes based on operators
    // single input vec, single output vec.
    /////////////////////////////////////////////////////////////

  template<class Field> 
    class ConjugateGradient : public OperatorFunction<Field> {
public:                                                
    RealD   Tolerance;
    Integer MaxIterations;
    ConjugateGradient(RealD tol,Integer maxit) : Tolerance(tol), MaxIterations(maxit) { 
    };


    void operator() (LinearOperatorBase<Field> &Linop,const Field &src, Field &psi){

      psi.checkerboard = src.checkerboard;
      conformable(psi,src);

      RealD cp,c,a,d,b,ssq,qq,b_pred;
      
      Field   p(src);
      Field mmp(src);
      Field   r(src);
      
      //Initial residual computation & set up
      RealD guess = norm2(psi);
      assert(std::isnan(guess)==0);

      Linop.HermOpAndNorm(psi,mmp,d,b);
      
      r= src-mmp;
      p= r;
      
      a  =norm2(p);
      cp =a;
      ssq=norm2(src);

      std::cout<<GridLogIterative <<std::setprecision(4)<< "ConjugateGradient: guess "<<guess<<std::endl;
      std::cout<<GridLogIterative <<std::setprecision(4)<< "ConjugateGradient:   src "<<ssq  <<std::endl;
      std::cout<<GridLogIterative <<std::setprecision(4)<< "ConjugateGradient:    mp "<<d    <<std::endl;
      std::cout<<GridLogIterative <<std::setprecision(4)<< "ConjugateGradient:   mmp "<<b    <<std::endl;
      std::cout<<GridLogIterative <<std::setprecision(4)<< "ConjugateGradient:  cp,r "<<cp   <<std::endl;
      std::cout<<GridLogIterative <<std::setprecision(4)<< "ConjugateGradient:     p "<<a    <<std::endl;

      RealD rsq =  Tolerance* Tolerance*ssq;
      
      //Check if guess is really REALLY good :)
      if ( cp <= rsq ) {
	return;
      }
      
      std::cout<<GridLogIterative << std::setprecision(4)<< "ConjugateGradient: k=0 residual "<<cp<<" target "<<rsq<<std::endl;

      GridStopWatch LinalgTimer;
      GridStopWatch MatrixTimer;
      GridStopWatch SolverTimer;

      SolverTimer.Start();
      int k;
      for (k=1;k<=MaxIterations;k++){
	
	c=cp;

	MatrixTimer.Start();
	Linop.HermOpAndNorm(p,mmp,d,qq);
	MatrixTimer.Stop();

	LinalgTimer.Start();
	//	RealD    qqck = norm2(mmp);
	//	ComplexD dck  = innerProduct(p,mmp);
      
	a      = c/d;
	b_pred = a*(a*qq-d)/c;

	cp = axpy_norm(r,-a,mmp,r);
	b = cp/c;
	
	// Fuse these loops ; should be really easy
	psi= a*p+psi;
	p  = p*b+r;
	  
	LinalgTimer.Stop();
	std::cout<<GridLogIterative<<"ConjugateGradient: Iteration " <<k<<" residual "<<cp<< " target "<< rsq<<std::endl;
	
	// Stopping condition
	if ( cp <= rsq ) { 
	  
	  SolverTimer.Stop();
	  Linop.HermOpAndNorm(psi,mmp,d,qq);
	  p=mmp-src;
	  
	  RealD mmpnorm = sqrt(norm2(mmp));
	  RealD psinorm = sqrt(norm2(psi));
	  RealD srcnorm = sqrt(norm2(src));
	  RealD resnorm = sqrt(norm2(p));
	  RealD true_residual = resnorm/srcnorm;

	  std::cout<<GridLogMessage<<"ConjugateGradient: Converged on iteration " <<k
		   <<" computed residual "<<sqrt(cp/ssq)
		   <<" true residual "    <<true_residual
		   <<" target "<<Tolerance<<std::endl;
	  std::cout<<GridLogMessage<<"Time elapsed: Total "<< SolverTimer.Elapsed() << " Matrix  "<<MatrixTimer.Elapsed() << " Linalg "<<LinalgTimer.Elapsed();
	  std::cout<<std::endl;
	  
	  assert(true_residual/Tolerance < 1000.0);

	  return;
	}
      }
      std::cout<<GridLogMessage<<"ConjugateGradient did NOT converge"<<std::endl;
      assert(0);
    }
  };
}
#endif
