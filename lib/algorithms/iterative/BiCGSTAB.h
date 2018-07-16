/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/BiCGSTAB.h

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: juettner <juettner@soton.ac.uk>

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
#ifndef GRID_BICGSTAB_H
#define GRID_BICGSTAB_H

namespace Grid {

/////////////////////////////////////////////////////////////
// Base classes for iterative processes based on operators
// single input vec, single output vec.
// 
// Implementation of BiCGSTAB method following Van Der Vorts's
// article SIAM J. Sci. and Stat. Comput., 13(2), 631–644
// Read More: https://epubs.siam.org/doi/10.1137/0913035"
// "Bi-CGSTAB: A Fast and Smoothly Converging Variant of 
// Bi-CG for the Solution of Nonsymmetric Linear Systems"
//
/////////////////////////////////////////////////////////////

template <class Field>
class BiCGSTAB : public OperatorFunction<Field> {
 public:
  bool ErrorOnNoConverge;  // throw an assert when the BiCGSTAB fails to converge.
                           // Defaults true.
  RealD Tolerance;
  Integer MaxIterations;
  Integer IterationsToComplete; //Number of iterations the BiCGSTAB took to finish. Filled in upon completion
  
  BiCGSTAB(RealD tol, Integer maxit, bool err_on_no_conv = true)
      : Tolerance(tol),
        MaxIterations(maxit),
        ErrorOnNoConverge(err_on_no_conv){};

  void operator()(LinearOperatorBase<Field> &Linop, const Field &src, Field &psi) {


    psi.checkerboard = src.checkerboard;
    conformable(psi, src);

    
    RealD cp,rho,oldrho,alpha,beta,omega;
    RealD a,bo;
    RealD d,b,ssq,qq;
    rho=1;
    alpha=1;
    omega=1;

    Field p(src);
    Field r(src);
    Field rhat(src);
    Field v(src);
    Field s(src);
    Field t(src);
    
    v=0;
    p=0;

    // Initial residual computation & set up
    RealD guess = norm2(psi);
    assert(std::isnan(guess) == 0);

    
    Linop.HermOpAndNorm(psi, v, d, b);
    
    r = src - v; // r0=b-Ax0
    rhat = r;
    a = norm2(r);
    ssq = norm2(src);

    std::cout << GridLogIterative << std::setprecision(8) << "BiCGSTAB: guess " << guess << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "BiCGSTAB:   src " << ssq << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "BiCGSTAB:    mp " << d << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "BiCGSTAB:   mmp " << b << std::endl;
    std::cout << GridLogIterative << std::setprecision(8) << "BiCGSTAB:     r " << a << std::endl;

    RealD rsq = Tolerance * Tolerance * ssq;

    // Check if guess is really REALLY good :)
    if (cp <= rsq) {
      return;
    }

    std::cout << GridLogIterative << std::setprecision(8)
              << "BiCGSTAB: k=0 residual " << a << " target " << rsq << std::endl;

    GridStopWatch LinalgTimer;
    GridStopWatch InnerTimer;
    GridStopWatch AxpyNormTimer;
    GridStopWatch LinearCombTimer;
    GridStopWatch MatrixTimer;
    GridStopWatch SolverTimer;

    SolverTimer.Start();
    int k;
    for (k = 1; k <= MaxIterations; k++) {

      oldrho=rho;

      LinalgTimer.Start();
      InnerTimer.Start();
      ComplexD Crho = innerProduct(rhat,r);
      InnerTimer.Stop();
      rho=Crho.real();

      beta = (rho/oldrho)*(alpha/omega);

      LinearCombTimer.Start();
      bo = +beta*omega;
      parallel_for(int ss=0;ss<src._grid->oSites();ss++){
	      vstream(p[ss], +beta * p[ss] + r[ss]);
        vstream(p[ss], -bo   * v[ss] + p[ss]);
      }
      LinearCombTimer.Stop();
      LinalgTimer.Stop();

      MatrixTimer.Start();
      Linop.HermOp(p, v);
      MatrixTimer.Stop();

      LinalgTimer.Start();
      InnerTimer.Start();
      ComplexD Calpha = innerProduct(rhat,v);
      InnerTimer.Stop();
      alpha = Calpha.real();

      alpha = rho/alpha;

      LinearCombTimer.Start();
      parallel_for(int ss=0;ss<src._grid->oSites();ss++){
	      vstream(s[ss], -alpha * v[ss] + r[ss]);
      }
      LinearCombTimer.Stop();
      LinalgTimer.Stop();

      MatrixTimer.Start();
      Linop.HermOp(s, t);
      MatrixTimer.Stop();


      LinalgTimer.Start();
      InnerTimer.Start();
      ComplexD Comega = innerProduct(t,s);
      InnerTimer.Stop();
      omega = Comega.real()/norm2(t);

      LinearCombTimer.Start();
      parallel_for(int ss=0;ss<src._grid->oSites();ss++){
	      vstream(psi[ss]   , + alpha * p[ss] + psi[ss]);
        vstream(psi[ss]   , + omega * s[ss] + psi[ss]);
      }
      LinearCombTimer.Stop();

      LinearCombTimer.Start();
      parallel_for(int ss=0;ss<src._grid->oSites();ss++){
	      vstream(r   [ss]   , -omega *  t[ss] + s   [ss]);        
      }
      LinearCombTimer.Stop();

      cp = norm2(r);
      LinalgTimer.Stop();

      std::cout << GridLogIterative << "BiCGSTAB: Iteration " << k
                << " residual " << cp << " target " << rsq << std::endl;

      // Stopping condition
      if (cp <= rsq) {
        SolverTimer.Stop();
        Linop.HermOpAndNorm(psi, v, d, qq);
        p = v - src;

        RealD srcnorm = sqrt(norm2(src));
        RealD resnorm = sqrt(norm2(p));
        RealD true_residual = resnorm / srcnorm;

        std::cout << GridLogMessage << "BiCGSTAB Converged on iteration " << k << std::endl;
        std::cout << GridLogMessage << "\tComputed residual " << sqrt(cp / ssq)<<std::endl;
	std::cout << GridLogMessage << "\tTrue residual " << true_residual<<std::endl;
	std::cout << GridLogMessage << "\tTarget " << Tolerance << std::endl;

        std::cout << GridLogMessage << "Time breakdown "<<std::endl;
	std::cout << GridLogMessage << "\tElapsed    " << SolverTimer.Elapsed() <<std::endl;
	std::cout << GridLogMessage << "\tMatrix     " << MatrixTimer.Elapsed() <<std::endl;
	std::cout << GridLogMessage << "\tLinalg     " << LinalgTimer.Elapsed() <<std::endl;
	std::cout << GridLogMessage << "\tInner      " << InnerTimer.Elapsed() <<std::endl;
	std::cout << GridLogMessage << "\tAxpyNorm   " << AxpyNormTimer.Elapsed() <<std::endl;
	std::cout << GridLogMessage << "\tLinearComb " << LinearCombTimer.Elapsed() <<std::endl;

        if (ErrorOnNoConverge) assert(true_residual / Tolerance < 10000.0);

	IterationsToComplete = k;	

        return;
      }
    }
    std::cout << GridLogMessage << "BiCGSTAB did NOT converge"
              << std::endl;

    if (ErrorOnNoConverge) assert(0);
    IterationsToComplete = k;

  }
};
}
#endif
