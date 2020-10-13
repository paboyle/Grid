/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithmsf/iterative/QuasiMinimalResidual.h

Copyright (C) 2019

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#pragma once

NAMESPACE_BEGIN(Grid);

template<class Field> 
RealD innerG5ProductReal(Field &l, Field &r)
{
  Gamma G5(Gamma::Algebra::Gamma5);
  Field tmp(l.Grid());
  //  tmp = G5*r;
  G5R5(tmp,r);
  ComplexD ip =innerProduct(l,tmp);
  std::cout << "innerProductRealG5R5 "<<ip<<std::endl;
  return ip.real();
}

template<class Field>
class QuasiMinimalResidual : public OperatorFunction<Field> {
 public:
  using OperatorFunction<Field>::operator();

  bool ErrorOnNoConverge; 
  RealD   Tolerance;
  Integer MaxIterations;
  Integer IterationCount;

  QuasiMinimalResidual(RealD   tol,
		       Integer maxit,
		       bool    err_on_no_conv = true)
      : Tolerance(tol)
      , MaxIterations(maxit)
      , ErrorOnNoConverge(err_on_no_conv) 
  {};

#if 1
  void operator()(LinearOperatorBase<Field> &LinOp, const Field &b, Field &x) 
  {
    RealD resid;
    IterationCount=0;

    RealD  rho, rho_1, xi, gamma, gamma_1, theta, theta_1;
    RealD  eta, delta, ep, beta; 

    GridBase *Grid = b.Grid();
    Field r(Grid), d(Grid), s(Grid);
    Field v(Grid), w(Grid), y(Grid),  z(Grid);
    Field v_tld(Grid), w_tld(Grid), y_tld(Grid), z_tld(Grid);
    Field p(Grid), q(Grid), p_tld(Grid);

    Real normb = norm2(b);

    LinOp.Op(x,r); r = b - r;

    assert(normb> 0.0);

    resid = norm2(r)/normb;
    if (resid <= Tolerance) {
      return;
    }

    v_tld = r;
    y = v_tld;
    rho = norm2(y);

    // Take Gamma5 conjugate
    //    Gamma G5(Gamma::Algebra::Gamma5);
    //    G5R5(w_tld,r);
    //    w_tld = G5* v_tld;
    w_tld=v_tld;
    z = w_tld;
    xi = norm2(z);

    gamma = 1.0;
    eta   = -1.0;
    theta = 0.0;

    for (int i = 1; i <= MaxIterations; i++) {

      // Breakdown tests
      assert( rho != 0.0);
      assert( xi  != 0.0);

      v = (1. / rho) * v_tld;
      y = (1. / rho) * y;

      w = (1. / xi) * w_tld;
      z = (1. / xi) * z;

      ComplexD Zdelta = innerProduct(z, y); // Complex?
      std::cout << "Zdelta "<<Zdelta<<std::endl;
      delta = Zdelta.real();

      y_tld = y; 
      z_tld = z;

      if (i > 1) {
	p = y_tld - (xi  * delta / ep) * p;
	q = z_tld - (rho * delta / ep) * q;
      } else {
	p = y_tld;
	q = z_tld;
      }

      LinOp.Op(p,p_tld);      //     p_tld = A * p;
      ComplexD Zep = innerProduct(q, p_tld);
      ep=Zep.real();
      std::cout << "Zep "<<Zep <<std::endl;
      // Complex Audit
      assert(abs(ep)>0);

      beta = ep / delta;
      assert(abs(beta)>0);

      v_tld = p_tld - beta * v;
      y = v_tld;

      rho_1 = rho;
      rho   = norm2(y);
      LinOp.AdjOp(q,w_tld);
      w_tld = w_tld - beta * w;
      z = w_tld;

      xi = norm2(z);

      gamma_1 = gamma;
      theta_1 = theta;

      theta   = rho / (gamma_1 * beta);
      gamma   = 1.0 / sqrt(1.0 + theta * theta);
      std::cout << "theta "<<theta<<std::endl;
      std::cout << "gamma "<<gamma<<std::endl;

      assert(abs(gamma)> 0.0);

      eta = -eta * rho_1 * gamma* gamma / (beta * gamma_1 * gamma_1);

      if (i > 1) {
	d = eta * p + (theta_1 * theta_1 * gamma * gamma) * d;
	s = eta * p_tld + (theta_1 * theta_1 * gamma * gamma) * s;
      } else {
	d = eta * p;
	s = eta * p_tld;
      }

      x =x+d;                            // update approximation vector
      r =r-s;                            // compute residual

      if ((resid = norm2(r) / normb) <= Tolerance) {
	return;
      }
      std::cout << "Iteration "<<i<<" resid " << resid<<std::endl;
    }
    assert(0);
    return;                            // no convergence
  }
#else
  // QMRg5 SMP thesis
  void operator()(LinearOperatorBase<Field> &LinOp, const Field &b, Field &x) 
  {
    // Real scalars
    GridBase *grid = b.Grid();

    Field    r(grid);
    Field    p_m(grid), p_m_minus_1(grid), p_m_minus_2(grid);
    Field    v_m(grid), v_m_minus_1(grid), v_m_plus_1(grid);
    Field    tmp(grid);

    RealD    w;
    RealD    z1, z2;
    RealD    delta_m, delta_m_minus_1;
    RealD    c_m_plus_1, c_m, c_m_minus_1;
    RealD    s_m_plus_1, s_m, s_m_minus_1;
    RealD    alpha, beta, gamma, epsilon;
    RealD    mu, nu, rho, theta, xi, chi;
    RealD    mod2r, mod2b;
    RealD    tau2, target2;

    mod2b=norm2(b);

    /////////////////////////
    // Initial residual
    /////////////////////////
    LinOp.Op(x,tmp);
    r = b - tmp;

    /////////////////////////
    // \mu = \rho = |r_0|
    /////////////////////////
    mod2r = norm2(r);
    rho = sqrt( mod2r);
    mu=rho;
    
    std::cout << "QuasiMinimalResidual rho "<< rho<<std::endl;
    /////////////////////////
    // Zero negative history
    /////////////////////////
    v_m_plus_1  = Zero();
    v_m_minus_1 = Zero();
    p_m_minus_1 = Zero();
    p_m_minus_2 = Zero();

    // v0
    v_m = (1.0/rho)*r;

    /////////////////////////
    // Initial coeffs
    /////////////////////////
    delta_m_minus_1 = 1.0;
    c_m_minus_1     = 1.0;
    c_m             = 1.0;
    s_m_minus_1     = 0.0;
    s_m             = 0.0;

    /////////////////////////
    // Set up convergence check
    /////////////////////////
    tau2    = mod2r;
    target2 = mod2b * Tolerance*Tolerance;
 
    for(int iter = 0 ; iter < MaxIterations; iter++){

      /////////////////////////
      // \delta_m = (v_m, \gamma_5 v_m) 
      /////////////////////////
      delta_m = innerG5ProductReal(v_m,v_m);
      std::cout << "QuasiMinimalResidual delta_m "<< delta_m<<std::endl;

      /////////////////////////
      // tmp = A v_m
      /////////////////////////
      LinOp.Op(v_m,tmp);

      /////////////////////////
      // \alpha = (v_m, \gamma_5 temp) / \delta_m 
      /////////////////////////
      alpha = innerG5ProductReal(v_m,tmp);
      alpha = alpha/delta_m ;
      std::cout << "QuasiMinimalResidual alpha "<< alpha<<std::endl;

      /////////////////////////
      // \beta = \rho \delta_m / \delta_{m-1}
      /////////////////////////
      beta = rho * delta_m / delta_m_minus_1;
      std::cout << "QuasiMinimalResidual beta "<< beta<<std::endl;

      /////////////////////////
      // \tilde{v}_{m+1} = temp - \alpha v_m - \beta v_{m-1}
      /////////////////////////
      v_m_plus_1 = tmp - alpha*v_m - beta*v_m_minus_1;

      ///////////////////////////////
      // \rho = || \tilde{v}_{m+1} ||
      ///////////////////////////////
      rho = sqrt( norm2(v_m_plus_1) );
      std::cout << "QuasiMinimalResidual rho "<< rho<<std::endl;

      ///////////////////////////////
      //      v_{m+1} = \tilde{v}_{m+1}
      ///////////////////////////////
      v_m_plus_1 = (1.0 / rho) * v_m_plus_1;

      ////////////////////////////////
      // QMR recurrence coefficients.
      ////////////////////////////////
      theta      = s_m_minus_1 * beta;
      gamma      = c_m_minus_1 * beta;
      epsilon    =  c_m * gamma + s_m * alpha;
      xi         = -s_m * gamma + c_m * alpha;
      nu         = sqrt( xi*xi + rho*rho );
      c_m_plus_1 = fabs(xi) / nu;
      if ( xi == 0.0 ) {
	s_m_plus_1 = 1.0;
      } else {
	s_m_plus_1 = c_m_plus_1 * rho / xi;
      }
      chi = c_m_plus_1 * xi + s_m_plus_1 * rho;

      std::cout << "QuasiMinimalResidual coeffs "<< theta <<" "<<gamma<<" "<< epsilon<<" "<< xi<<" "<< nu<<std::endl;
      std::cout << "QuasiMinimalResidual coeffs "<< chi   <<std::endl;

      ////////////////////////////////
      //p_m=(v_m - \epsilon p_{m-1} - \theta p_{m-2}) / \chi
      ////////////////////////////////
      p_m = (1.0/chi) * v_m - (epsilon/chi) * p_m_minus_1 - (theta/chi) * p_m_minus_2;

      ////////////////////////////////////////////////////////////////
      //      \psi = \psi + c_{m+1} \mu p_m	
      ////////////////////////////////////////////////////////////////
      x = x + ( c_m_plus_1 * mu ) * p_m;

      ////////////////////////////////////////
      //
      ////////////////////////////////////////
      mu              = -s_m_plus_1 * mu;
      delta_m_minus_1 = delta_m;
      c_m_minus_1     = c_m;
      c_m             = c_m_plus_1;
      s_m_minus_1     = s_m;
      s_m             = s_m_plus_1;

      ////////////////////////////////////
      // Could use pointer swizzle games.
      ////////////////////////////////////
      v_m_minus_1 = v_m;
      v_m         = v_m_plus_1;
      p_m_minus_2 = p_m_minus_1;
      p_m_minus_1 = p_m;


      /////////////////////////////////////
      // Convergence checks
      /////////////////////////////////////
      z1 = RealD(iter+1.0);
      z2 = z1 + 1.0;
      tau2 = tau2 *( z2 / z1 ) * s_m * s_m;
      std::cout << " QuasiMinimumResidual iteration "<< iter<<std::endl;
      std::cout << " QuasiMinimumResidual tau bound "<< tau2<<std::endl;

      // Compute true residual
      mod2r = tau2;
      if ( 1 || (tau2 < (100.0 * target2)) ) {
	LinOp.Op(x,tmp);
	r = b - tmp;
	mod2r = norm2(r);
	std::cout << " QuasiMinimumResidual true residual is "<< mod2r<<std::endl;
      }


      if ( mod2r < target2 ) { 

	std::cout << " QuasiMinimumResidual has converged"<<std::endl;
	return;

      }

    }


  }
#endif
};

NAMESPACE_END(Grid);
