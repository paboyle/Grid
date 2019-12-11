/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/BiCGSTAB.h

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: juettner <juettner@soton.ac.uk>
Author: David Murphy <djmurphy@mit.edu>

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

NAMESPACE_BEGIN(Grid);

/////////////////////////////////////////////////////////////
// Base classes for iterative processes based on operators
// single input vec, single output vec.
/////////////////////////////////////////////////////////////

template <class Field>
class BiCGSTAB : public OperatorFunction<Field> 
{
  public:
    using OperatorFunction<Field>::operator();
    
    bool ErrorOnNoConverge;  // throw an assert when the CG fails to converge.
                             // Defaults true.
    RealD Tolerance;
    Integer MaxIterations;
    Integer IterationsToComplete; //Number of iterations the CG took to finish. Filled in upon completion
  
    BiCGSTAB(RealD tol, Integer maxit, bool err_on_no_conv = true) : 
      Tolerance(tol), MaxIterations(maxit), ErrorOnNoConverge(err_on_no_conv){};

    void operator()(LinearOperatorBase<Field>& Linop, const Field& src, Field& psi) 
    {
      psi.Checkerboard() = src.Checkerboard();
      conformable(psi, src);

      RealD cp(0), rho(1), rho_prev(0), alpha(1), beta(0), omega(1);
      RealD a(0), bo(0), b(0), ssq(0);

      Field p(src);
      Field r(src);
      Field rhat(src);
      Field v(src);
      Field s(src);
      Field t(src);
      Field h(src);

      v = Zero();
      p = Zero();

      // Initial residual computation & set up
      RealD guess = norm2(psi);
      assert(std::isnan(guess) == 0);
    
      Linop.Op(psi, v);
      b = norm2(v);

      r = src - v;
      rhat = r;
      a = norm2(r);
      ssq = norm2(src);

      std::cout << GridLogIterative << std::setprecision(8) << "BiCGSTAB: guess " << guess << std::endl;
      std::cout << GridLogIterative << std::setprecision(8) << "BiCGSTAB:   src " << ssq << std::endl;
      std::cout << GridLogIterative << std::setprecision(8) << "BiCGSTAB:    mp " << b << std::endl;
      std::cout << GridLogIterative << std::setprecision(8) << "BiCGSTAB:     r " << a << std::endl;

      RealD rsq = Tolerance * Tolerance * ssq;

      // Check if guess is really REALLY good :)
      if(a <= rsq){ return; }

      std::cout << GridLogIterative << std::setprecision(8) << "BiCGSTAB: k=0 residual " << a << " target " << rsq << std::endl;

      GridStopWatch LinalgTimer;
      GridStopWatch InnerTimer;
      GridStopWatch AxpyNormTimer;
      GridStopWatch LinearCombTimer;
      GridStopWatch MatrixTimer;
      GridStopWatch SolverTimer;

      SolverTimer.Start();
      int k;
      for (k = 1; k <= MaxIterations; k++) 
      {
        rho_prev = rho;

        LinalgTimer.Start();
        InnerTimer.Start();
        ComplexD Crho  = innerProduct(rhat,r);
        InnerTimer.Stop();
        rho = Crho.real();

        beta = (rho / rho_prev) * (alpha / omega);

        LinearCombTimer.Start();
        bo = beta * omega;
        auto p_v = p.View();
        auto r_v = r.View();
        auto v_v = v.View();
        accelerator_for(ss, p_v.size(), Field::vector_object::Nsimd(),{
          coalescedWrite(p_v[ss], beta*p_v(ss) - bo*v_v(ss) + r_v(ss));
        });
        LinearCombTimer.Stop();
        LinalgTimer.Stop();

        MatrixTimer.Start();
        Linop.Op(p,v);
        MatrixTimer.Stop();

        LinalgTimer.Start();
        InnerTimer.Start();
        ComplexD Calpha = innerProduct(rhat,v);
        InnerTimer.Stop();
        alpha = rho / Calpha.real();

        LinearCombTimer.Start();
        auto h_v = h.View();
        auto psi_v = psi.View();
        accelerator_for(ss, h_v.size(), Field::vector_object::Nsimd(),{
          coalescedWrite(h_v[ss], alpha*p_v(ss) + psi_v(ss));
        });
        
        auto s_v = s.View();
        accelerator_for(ss, s_v.size(), Field::vector_object::Nsimd(),{
          coalescedWrite(s_v[ss], -alpha*v_v(ss) + r_v(ss));
        });
        LinearCombTimer.Stop();
        LinalgTimer.Stop();

        MatrixTimer.Start();
        Linop.Op(s,t);
        MatrixTimer.Stop();

        LinalgTimer.Start();
        InnerTimer.Start();
        ComplexD Comega = innerProduct(t,s);
        InnerTimer.Stop();
        omega = Comega.real() / norm2(t);

        LinearCombTimer.Start();
        auto t_v = t.View();
        accelerator_for(ss, psi_v.size(), Field::vector_object::Nsimd(),{
          coalescedWrite(psi_v[ss], h_v(ss) + omega * s_v(ss));
          coalescedWrite(r_v[ss], -omega * t_v(ss) + s_v(ss));
        });
        LinearCombTimer.Stop();

        cp = norm2(r);
        LinalgTimer.Stop();

        std::cout << GridLogIterative << "BiCGSTAB: Iteration " << k << " residual " << sqrt(cp/ssq) << " target " << Tolerance << std::endl;

        // Stopping condition
        if(cp <= rsq) 
        {
          SolverTimer.Stop();
          Linop.Op(psi, v);
          p = v - src;

          RealD srcnorm = sqrt(norm2(src));
          RealD resnorm = sqrt(norm2(p));
          RealD true_residual = resnorm / srcnorm;

          std::cout << GridLogMessage << "BiCGSTAB Converged on iteration " << k << std::endl;
          std::cout << GridLogMessage << "\tComputed residual " << sqrt(cp/ssq) << std::endl;
          std::cout << GridLogMessage << "\tTrue residual " << true_residual << std::endl;
          std::cout << GridLogMessage << "\tTarget " << Tolerance << std::endl;

          std::cout << GridLogMessage << "Time breakdown " << std::endl;
          std::cout << GridLogMessage << "\tElapsed    " << SolverTimer.Elapsed() << std::endl;
          std::cout << GridLogMessage << "\tMatrix     " << MatrixTimer.Elapsed() << std::endl;
          std::cout << GridLogMessage << "\tLinalg     " << LinalgTimer.Elapsed() << std::endl;
          std::cout << GridLogMessage << "\tInner      " << InnerTimer.Elapsed() << std::endl;
          std::cout << GridLogMessage << "\tAxpyNorm   " << AxpyNormTimer.Elapsed() << std::endl;
          std::cout << GridLogMessage << "\tLinearComb " << LinearCombTimer.Elapsed() << std::endl;

          if(ErrorOnNoConverge){ assert(true_residual / Tolerance < 10000.0); }

          IterationsToComplete = k;	

          return;
        }
      }
      
      std::cout << GridLogMessage << "BiCGSTAB did NOT converge" << std::endl;

      if(ErrorOnNoConverge){ assert(0); }
      IterationsToComplete = k;
    }
};

NAMESPACE_END(Grid);

#endif
