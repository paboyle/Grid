/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/ConjugateGradientMultiShift.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Christopher Kelly <ckelly@bnl.gov>

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
#pragma once

NAMESPACE_BEGIN(Grid);

//CK 2020: A variant of the multi-shift conjugate gradient with the matrix multiplication in single precision. 
//The residual is stored in single precision, but the search directions and solution are stored in double precision. 
//Every update_freq iterations the residual is corrected in double precision. 
//For safety the a final regular CG is applied to clean up if necessary

//PB Pure single, then double fixup

template<class FieldD, class FieldF,
	 typename std::enable_if< getPrecision<FieldD>::value == 2, int>::type = 0,
	 typename std::enable_if< getPrecision<FieldF>::value == 1, int>::type = 0> 
class ConjugateGradientMultiShiftMixedPrecCleanup : public OperatorMultiFunction<FieldD>,
					     public OperatorFunction<FieldD>
{
public:                                                

  using OperatorFunction<FieldD>::operator();

  RealD   Tolerance;
  Integer MaxIterationsMshift;
  Integer MaxIterations;
  Integer IterationsToComplete; //Number of iterations the CG took to finish. Filled in upon completion
  std::vector<int> IterationsToCompleteShift;  // Iterations for this shift
  int verbose;
  MultiShiftFunction shifts;
  std::vector<RealD> TrueResidualShift;

  int ReliableUpdateFreq; //number of iterations between reliable updates

  GridBase* SinglePrecGrid; //Grid for single-precision fields
  LinearOperatorBase<FieldF> &Linop_f; //single precision

  ConjugateGradientMultiShiftMixedPrecCleanup(Integer maxit, const MultiShiftFunction &_shifts,
				       GridBase* _SinglePrecGrid, LinearOperatorBase<FieldF> &_Linop_f,
				       int _ReliableUpdateFreq) : 
    MaxIterationsMshift(maxit),  shifts(_shifts), SinglePrecGrid(_SinglePrecGrid), Linop_f(_Linop_f), ReliableUpdateFreq(_ReliableUpdateFreq),
    MaxIterations(20000)
  { 
    verbose=1;
    IterationsToCompleteShift.resize(_shifts.order);
    TrueResidualShift.resize(_shifts.order);
  }

  void operator() (LinearOperatorBase<FieldD> &Linop, const FieldD &src, FieldD &psi)
  {
    GridBase *grid = src.Grid();
    int nshift = shifts.order;
    std::vector<FieldD> results(nshift,grid);
    (*this)(Linop,src,results,psi);
  }
  void operator() (LinearOperatorBase<FieldD> &Linop, const FieldD &src, std::vector<FieldD> &results, FieldD &psi)
  {
    int nshift = shifts.order;

    (*this)(Linop,src,results);
  
    psi = shifts.norm*src;
    for(int i=0;i<nshift;i++){
      psi = psi + shifts.residues[i]*results[i];
    }

    return;
  }

  void operator() (LinearOperatorBase<FieldD> &Linop_d, const FieldD &src_d, std::vector<FieldD> &psi_d)
  { 
    GRID_TRACE("ConjugateGradientMultiShiftMixedPrecCleanup");
    GridBase *DoublePrecGrid = src_d.Grid();

    ////////////////////////////////////////////////////////////////////////
    // Convenience references to the info stored in "MultiShiftFunction"
    ////////////////////////////////////////////////////////////////////////
    int nshift = shifts.order;

    std::vector<RealD> &mass(shifts.poles); // Make references to array in "shifts"
    std::vector<RealD> &mresidual(shifts.tolerances);
    std::vector<RealD> alpha(nshift,1.0);

    //Double precision search directions
    FieldD p_d(DoublePrecGrid);
    std::vector<FieldF> ps_f (nshift, SinglePrecGrid);// Search directions (single precision)
    std::vector<FieldF> psi_f(nshift, SinglePrecGrid);// solutions (single precision)

    FieldD tmp_d(DoublePrecGrid);
    FieldD r_d(DoublePrecGrid);
    FieldF r_f(SinglePrecGrid);
    FieldD mmp_d(DoublePrecGrid);

    assert(psi_d.size()==nshift);
    assert(mass.size()==nshift);
    assert(mresidual.size()==nshift);
  
    // dynamic sized arrays on stack; 2d is a pain with vector
    std::vector<RealD>  bs(nshift);
    std::vector<RealD>  rsq(nshift);
    std::vector<RealD>  rsqf(nshift);
    std::vector<std::array<RealD,2> >  z(nshift);
    std::vector<int>     converged(nshift);
  
    const int       primary =0;
  
    //Primary shift fields CG iteration
    RealD a,b,c,d;
    RealD cp,bp,qq; //prev
  
    // Matrix mult fields
    FieldF p_f(SinglePrecGrid);
    FieldF mmp_f(SinglePrecGrid);

    // Check lightest mass
    for(int s=0;s<nshift;s++){
      assert( mass[s]>= mass[primary] );
      converged[s]=0;
    }
  
    // Wire guess to zero
    // Residuals "r" are src
    // First search direction "p" is also src
    cp = norm2(src_d);

    // Handle trivial case of zero src.
    if( cp == 0. ){
      for(int s=0;s<nshift;s++){
	psi_d[s] = Zero();
	psi_f[s] = Zero();
	IterationsToCompleteShift[s] = 1;
	TrueResidualShift[s] = 0.;
      }
      return;
    }

    for(int s=0;s<nshift;s++){
      rsq[s] = cp * mresidual[s] * mresidual[s];
      rsqf[s] =rsq[s];
      std::cout<<GridLogMessage<<"ConjugateGradientMultiShiftMixedPrecCleanup: shift "<< s <<" target resid "<<rsq[s]<<std::endl;
      //      ps_d[s] = src_d;
      precisionChange(ps_f[s],src_d);
    }
    // r and p for primary
    p_d = src_d; //primary copy --- make this a reference to ps_d to save axpys
    r_d = p_d;
    
    //MdagM+m[0]
    precisionChange(p_f,p_d);
    Linop_f.HermOpAndNorm(p_f,mmp_f,d,qq); // mmp = MdagM p        d=real(dot(p, mmp)),  qq=norm2(mmp)
    precisionChange(tmp_d,mmp_f);
    Linop_d.HermOpAndNorm(p_d,mmp_d,d,qq); // mmp = MdagM p        d=real(dot(p, mmp)),  qq=norm2(mmp)
    tmp_d = tmp_d - mmp_d;
    std::cout << " Testing operators match "<<norm2(mmp_d)<<" f "<<norm2(mmp_f)<<" diff "<< norm2(tmp_d)<<std::endl;
    //    assert(norm2(tmp_d)< 1.0e-4);

    axpy(mmp_d,mass[0],p_d,mmp_d);
    RealD rn = norm2(p_d);
    d += rn*mass[0];

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
    c=axpy_norm(r_d,b,mmp_d,r_d);
  
    for(int s=0;s<nshift;s++) {
      axpby(psi_d[s],0.,-bs[s]*alpha[s],src_d,src_d);
      precisionChange(psi_f[s],psi_d[s]);
    }
  
    ///////////////////////////////////////
    // Timers
    ///////////////////////////////////////
    GridStopWatch AXPYTimer, ShiftTimer, QRTimer, MatrixTimer, SolverTimer, PrecChangeTimer, CleanupTimer;

    SolverTimer.Start();
  
    // Iteration loop
    int k;
  
    for (k=1;k<=MaxIterationsMshift;k++){    

      a = c /cp;
      AXPYTimer.Start();
      axpy(p_d,a,p_d,r_d); 
      AXPYTimer.Stop();

      PrecChangeTimer.Start();
      precisionChange(r_f, r_d);
      PrecChangeTimer.Stop();

      AXPYTimer.Start();
      for(int s=0;s<nshift;s++){
	if ( ! converged[s] ) { 
	  if (s==0){
	    axpy(ps_f[s],a,ps_f[s],r_f);
	  } else{
	    RealD as =a *z[s][iz]*bs[s] /(z[s][1-iz]*b);
	    axpby(ps_f[s],z[s][iz],as,r_f,ps_f[s]);
	  }
	}
      }
      AXPYTimer.Stop();

      cp=c;
      PrecChangeTimer.Start();
      precisionChange(p_f, p_d); //get back single prec search direction for linop
      PrecChangeTimer.Stop();
      MatrixTimer.Start();  
      Linop_f.HermOp(p_f,mmp_f);
      MatrixTimer.Stop();  
      PrecChangeTimer.Start();
      precisionChange(mmp_d, mmp_f); // From Float to Double
      PrecChangeTimer.Stop();

      d=real(innerProduct(p_d,mmp_d));    
      axpy(mmp_d,mass[0],p_d,mmp_d);
      RealD rn = norm2(p_d);
      d += rn*mass[0];
    
      bp=b;
      b=-cp/d;

      // Toggle the recurrence history
      bs[0] = b;
      iz = 1-iz;
      ShiftTimer.Start();
      for(int s=1;s<nshift;s++){
	if((!converged[s])){
	  RealD z0 = z[s][1-iz];
	  RealD z1 = z[s][iz];
	  z[s][iz] = z0*z1*bp
	    / (b*a*(z1-z0) + z1*bp*(1- (mass[s]-mass[0])*b)); 
	  bs[s] = b*z[s][iz]/z0; // NB sign  rel to Mike
	}
      }
      ShiftTimer.Stop();

      //Update single precision solutions
      AXPYTimer.Start();
      for(int s=0;s<nshift;s++){
	int ss = s;
	if( (!converged[s]) ) { 
	  axpy(psi_f[ss],-bs[s]*alpha[s],ps_f[s],psi_f[ss]);
	}
      }
      c = axpy_norm(r_d,b,mmp_d,r_d);
      AXPYTimer.Stop();
    
      // Convergence checks
      int all_converged = 1;
      for(int s=0;s<nshift;s++){
      
	if ( (!converged[s]) ){
	  IterationsToCompleteShift[s] = k;
	
	  RealD css  = c * z[s][iz]* z[s][iz];
	
	  if(css<rsqf[s]){
	    if ( ! converged[s] )
	      std::cout<<GridLogMessage<<"ConjugateGradientMultiShiftMixedPrecCleanup k="<<k<<" Shift "<<s<<" has converged"<<std::endl;
	    converged[s]=1;
	  } else {
	    all_converged=0;
	  }

	}
      }

      if ( all_converged || k == MaxIterationsMshift-1){

	SolverTimer.Stop();

	for(int s=0;s<nshift;s++){
	  precisionChange(psi_d[s],psi_f[s]);
	}

	
	if ( all_converged ){
	  std::cout<<GridLogMessage<< "ConjugateGradientMultiShiftMixedPrecCleanup: All shifts have converged iteration "<<k<<std::endl;
	  std::cout<<GridLogMessage<< "ConjugateGradientMultiShiftMixedPrecCleanup: Checking solutions"<<std::endl;
	} else {
	  std::cout<<GridLogMessage<< "ConjugateGradientMultiShiftMixedPrecCleanup: Not all shifts have converged iteration "<<k<<std::endl;
	}
	
	// Check answers 
	for(int s=0; s < nshift; s++) { 
	  Linop_d.HermOpAndNorm(psi_d[s],mmp_d,d,qq);
	  axpy(tmp_d,mass[s],psi_d[s],mmp_d);
	  axpy(r_d,-alpha[s],src_d,tmp_d);
	  RealD rn = norm2(r_d);
	  RealD cn = norm2(src_d);
	  TrueResidualShift[s] = std::sqrt(rn/cn);
	  std::cout<<GridLogMessage<<"ConjugateGradientMultiShiftMixedPrecCleanup: shift["<<s<<"] true residual "<< TrueResidualShift[s] << " target " << mresidual[s] << std::endl;

	  //If we have not reached the desired tolerance, do a (mixed precision) CG cleanup
	  if(rn >= rsq[s]){
	    CleanupTimer.Start();
	    std::cout<<GridLogMessage<<"ConjugateGradientMultiShiftMixedPrecCleanup: performing cleanup step for shift " << s << std::endl;

	    //Setup linear operators for final cleanup
	    ConjugateGradientMultiShiftMixedPrecSupport::ShiftedLinop<FieldD> Linop_shift_d(Linop_d, mass[s]);
	    ConjugateGradientMultiShiftMixedPrecSupport::ShiftedLinop<FieldF> Linop_shift_f(Linop_f, mass[s]);
					       
	    MixedPrecisionConjugateGradient<FieldD,FieldF> cg(mresidual[s], MaxIterations, MaxIterations, SinglePrecGrid, Linop_shift_f, Linop_shift_d); 
	    cg(src_d, psi_d[s]);
	    
	    TrueResidualShift[s] = cg.TrueResidual;
	    CleanupTimer.Stop();
	  }
	}

	std::cout << GridLogMessage << "ConjugateGradientMultiShiftMixedPrecCleanup: Time Breakdown for body"<<std::endl;
	std::cout << GridLogMessage << "\tSolver    " << SolverTimer.Elapsed()     <<std::endl;
	std::cout << GridLogMessage << "\t\tAXPY    " << AXPYTimer.Elapsed()     <<std::endl;
	std::cout << GridLogMessage << "\t\tMatrix    " << MatrixTimer.Elapsed()     <<std::endl;
	std::cout << GridLogMessage << "\t\tShift    " << ShiftTimer.Elapsed()     <<std::endl;
	std::cout << GridLogMessage << "\t\tPrecision Change " << PrecChangeTimer.Elapsed()     <<std::endl;
	std::cout << GridLogMessage << "\tFinal Cleanup " << CleanupTimer.Elapsed()     <<std::endl;
	std::cout << GridLogMessage << "\tSolver+Cleanup " << SolverTimer.Elapsed() + CleanupTimer.Elapsed() << std::endl;

	IterationsToComplete = k;	

	return;
      }
   
    }
    std::cout<<GridLogMessage<<"CG multi shift did not converge"<<std::endl;
    assert(0);
  }

};
NAMESPACE_END(Grid);

