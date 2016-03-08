    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/ConjugateGradientMultiShift.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
#ifndef GRID_CONJUGATE_MULTI_SHIFT_GRADIENT_H
#define GRID_CONJUGATE_MULTI_SHIFT_GRADIENT_H

namespace Grid {

    /////////////////////////////////////////////////////////////
    // Base classes for iterative processes based on operators
    // single input vec, single output vec.
    /////////////////////////////////////////////////////////////

  template<class Field> 
    class ConjugateGradientMultiShift : public OperatorMultiFunction<Field>,
                                        public OperatorFunction<Field>
    {
public:                                                
    RealD   Tolerance;
    Integer MaxIterations;
    int verbose;
    MultiShiftFunction shifts;

    ConjugateGradientMultiShift(Integer maxit,MultiShiftFunction &_shifts) : 
	MaxIterations(maxit),
	shifts(_shifts)
    { 
      verbose=1;
    }

void operator() (LinearOperatorBase<Field> &Linop, const Field &src, Field &psi)
{
  GridBase *grid = src._grid;
  int nshift = shifts.order;
  std::vector<Field> results(nshift,grid);
  (*this)(Linop,src,results,psi);
}
void operator() (LinearOperatorBase<Field> &Linop, const Field &src, std::vector<Field> &results, Field &psi)
{
  int nshift = shifts.order;

  (*this)(Linop,src,results);
  
  psi = shifts.norm*src;
  for(int i=0;i<nshift;i++){
    psi = psi + shifts.residues[i]*results[i];
  }

  return;
}

void operator() (LinearOperatorBase<Field> &Linop, const Field &src, std::vector<Field> &psi)
{
  
  GridBase *grid = src._grid;
  
  ////////////////////////////////////////////////////////////////////////
  // Convenience references to the info stored in "MultiShiftFunction"
  ////////////////////////////////////////////////////////////////////////
  int nshift = shifts.order;

  std::vector<RealD> &mass(shifts.poles); // Make references to array in "shifts"
  std::vector<RealD> &mresidual(shifts.tolerances);
  std::vector<RealD> alpha(nshift,1.0);
  std::vector<Field>   ps(nshift,grid);// Search directions

  assert(psi.size()==nshift);
  assert(mass.size()==nshift);
  assert(mresidual.size()==nshift);
  
  // dynamic sized arrays on stack; 2d is a pain with vector
  RealD  bs[nshift];
  RealD  rsq[nshift];
  RealD  z[nshift][2];
  int     converged[nshift];
  
  const int       primary =0;
  
  //Primary shift fields CG iteration
  RealD a,b,c,d;
  RealD cp,bp,qq; //prev
  
  // Matrix mult fields
  Field r(grid);
  Field p(grid);
  Field tmp(grid);
  Field mmp(grid);
  
  // Check lightest mass
  for(int s=0;s<nshift;s++){
    assert( mass[s]>= mass[primary] );
    converged[s]=0;
  }
  
  // Wire guess to zero
  // Residuals "r" are src
  // First search direction "p" is also src
  cp = norm2(src);
  for(int s=0;s<nshift;s++){
    rsq[s] = cp * mresidual[s] * mresidual[s];
    std::cout<<GridLogMessage<<"ConjugateGradientMultiShift: shift "<<s
	     <<" target resid "<<rsq[s]<<std::endl;
    ps[s] = src;
  }
  // r and p for primary
  r=src;
  p=src;
  
  //MdagM+m[0]
  Linop.HermOpAndNorm(p,mmp,d,qq);
  axpy(mmp,mass[0],p,mmp);
  RealD rn = norm2(p);
  d += rn*mass[0];
  
  // have verified that inner product of 
  // p and mmp is equal to d after this since
  // the d computation is tricky
  //  qq = real(innerProduct(p,mmp));
  //  std::cout<<GridLogMessage << "debug equal ?  qq "<<qq<<" d "<< d<<std::endl;
  
  b = -cp /d;
  
  // Set up the various shift variables
  int       iz=0;
  z[0][1-iz] = 1.0;
  z[0][iz]   = 1.0;
  bs[0]      = b;
  for(int s=1;s<nshift;s++){
    z[s][1-iz] = 1.0;
    z[s][iz]   = 1.0/( 1.0 - b*(mass[s]-mass[0]));
    bs[s]      = b*z[s][iz]; 
  }
  
  // r += b[0] A.p[0]
  // c= norm(r)
  c=axpy_norm(r,b,mmp,r);
  
  for(int s=0;s<nshift;s++) {
    axpby(psi[s],0.,-bs[s]*alpha[s],src,src);
  }
  
  
  // Iteration loop
  int k;
  
  for (k=1;k<=MaxIterations;k++){
    
    a = c /cp;
    axpy(p,a,p,r);
    
    // Note to self - direction ps is iterated seperately
    // for each shift. Does not appear to have any scope
    // for avoiding linear algebra in "single" case.
    // 
    // However SAME r is used. Could load "r" and update
    // ALL ps[s]. 2/3 Bandwidth saving
    // New Kernel: Load r, vector of coeffs, vector of pointers ps
    for(int s=0;s<nshift;s++){
      if ( ! converged[s] ) { 
	if (s==0){
	  axpy(ps[s],a,ps[s],r);
	} else{
	  RealD as =a *z[s][iz]*bs[s] /(z[s][1-iz]*b);
	  axpby(ps[s],z[s][iz],as,r,ps[s]);
	}
      }
    }
    
    cp=c;
    
    Linop.HermOpAndNorm(p,mmp,d,qq);
    axpy(mmp,mass[0],p,mmp);
    RealD rn = norm2(p);
    d += rn*mass[0];
    
    bp=b;
    b=-cp/d;
    
    c=axpy_norm(r,b,mmp,r);

    // Toggle the recurrence history
    bs[0] = b;
    iz = 1-iz;
    for(int s=1;s<nshift;s++){
      if((!converged[s])){
	RealD z0 = z[s][1-iz];
	RealD z1 = z[s][iz];
	z[s][iz] = z0*z1*bp
	  / (b*a*(z1-z0) + z1*bp*(1- (mass[s]-mass[0])*b)); 
	bs[s] = b*z[s][iz]/z0; // NB sign  rel to Mike
      }
    }
    
    for(int s=0;s<nshift;s++){
      int ss = s;
      // Scope for optimisation here in case of "single".
      // Could load psi[0] and pull all ps[s] in.
      //      if ( single ) ss=primary;
      // Bandwith saving in single case is Ls * 3 -> 2+Ls, so ~ 3x saving
      // Pipelined CG gain:
      //
      // New Kernel: Load r, vector of coeffs, vector of pointers ps
      // New Kernel: Load psi[0], vector of coeffs, vector of pointers ps
      // If can predict the coefficient bs then we can fuse these and avoid write reread cyce
      //  on ps[s].
      // Before:  3 x npole  + 3 x npole
      // After :  2 x npole (ps[s])        => 3x speed up of multishift CG.
      
      if( (!converged[s]) ) { 
	axpy(psi[ss],-bs[s]*alpha[s],ps[s],psi[ss]);
      }
    }
    
    // Convergence checks
    int all_converged = 1;
    for(int s=0;s<nshift;s++){
      
      if ( (!converged[s]) ){
	
	RealD css  = c * z[s][iz]* z[s][iz];
	
	if(css<rsq[s]){
	  if ( ! converged[s] )
	    std::cout<<GridLogMessage<<"ConjugateGradientMultiShift k="<<k<<" Shift "<<s<<" has converged"<<std::endl;
	      converged[s]=1;
	} else {
	  all_converged=0;
	}

      }
    }
    
    if ( all_converged ){

      std::cout<<GridLogMessage<< "CGMultiShift: All shifts have converged iteration "<<k<<std::endl;
      std::cout<<GridLogMessage<< "CGMultiShift: Checking solutions"<<std::endl;
      
      // Check answers 
      for(int s=0; s < nshift; s++) { 
	Linop.HermOpAndNorm(psi[s],mmp,d,qq);
	axpy(tmp,mass[s],psi[s],mmp);
	axpy(r,-alpha[s],src,tmp);
	RealD rn = norm2(r);
	RealD cn = norm2(src);
	std::cout<<GridLogMessage<<"CGMultiShift: shift["<<s<<"] true residual "<<std::sqrt(rn/cn)<<std::endl;
      }
      return;
    }
  }
  // ugly hack
  std::cout<<GridLogMessage<<"CG multi shift did not converge"<<std::endl;
//  assert(0);
}

  };
}
#endif
