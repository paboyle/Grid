    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/PrecConjugateResidual.h

    Copyright (C) 2015

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_PREC_CONJUGATE_RESIDUAL_H
#define GRID_PREC_CONJUGATE_RESIDUAL_H

namespace Grid {

    /////////////////////////////////////////////////////////////
    // Base classes for iterative processes based on operators
    // single input vec, single output vec.
    /////////////////////////////////////////////////////////////

  template<class Field> 
    class PrecConjugateResidual : public OperatorFunction<Field> {
  public:                                                
    RealD   Tolerance;
    Integer MaxIterations;
    int verbose;
    LinearFunction<Field> &Preconditioner;

    PrecConjugateResidual(RealD tol,Integer maxit,LinearFunction<Field> &Prec) : Tolerance(tol), MaxIterations(maxit),      Preconditioner(Prec)
    { 
      verbose=1;
    };

    void operator() (LinearOperatorBase<Field> &Linop,const Field &src, Field &psi){

      RealD a, b, c, d;
      RealD cp, ssq,rsq;
      
      RealD rAr, rAAr, rArp;
      RealD pAp, pAAp;

      GridBase *grid = src._grid;
      Field r(grid),  p(grid), Ap(grid), Ar(grid), z(grid);
      
      psi=zero;
      r  = src;
      Preconditioner(r,p);

      

      Linop.HermOpAndNorm(p,Ap,pAp,pAAp);
      Ar=Ap;
      rAr=pAp;
      rAAr=pAAp;

      cp =norm2(r);
      ssq=norm2(src);
      rsq=Tolerance*Tolerance*ssq;

      if (verbose) std::cout<<GridLogMessage<<"PrecConjugateResidual: iteration " <<0<<" residual "<<cp<< " target"<< rsq<<std::endl;

      for(int k=0;k<MaxIterations;k++){


	Preconditioner(Ap,z);
	RealD rq= real(innerProduct(Ap,z)); 

	a = rAr/rq;

   	axpy(psi,a,p,psi);
   cp = axpy_norm(r,-a,z,r);

	rArp=rAr;

	Linop.HermOpAndNorm(r,Ar,rAr,rAAr);

	b   =rAr/rArp;
 
	axpy(p,b,p,r);
	pAAp=axpy_norm(Ap,b,Ap,Ar);
	
	if(verbose) std::cout<<GridLogMessage<<"PrecConjugateResidual: iteration " <<k<<" residual "<<cp<< " target"<< rsq<<std::endl;

	if(cp<rsq) {
	  Linop.HermOp(psi,Ap);
	  axpy(r,-1.0,src,Ap);
	  RealD true_resid = norm2(r)/ssq;
	  std::cout<<GridLogMessage<<"PrecConjugateResidual: Converged on iteration " <<k
		   << " computed residual "<<sqrt(cp/ssq)
	           << " true residual "<<sqrt(true_resid)
	           << " target "       <<Tolerance <<std::endl;
	  return;
	}

      }

      std::cout<<GridLogMessage<<"PrecConjugateResidual did NOT converge"<<std::endl;
      assert(0);
    }
  };
}
#endif
