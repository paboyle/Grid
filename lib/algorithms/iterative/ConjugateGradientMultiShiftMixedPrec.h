    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/ConjugateGradientMultiShiftMixedPrec.h

    Copyright (C) 2015

Author: Chulwoo Jung <chulwoo@quark.phy.bnl.gov>

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
    /*  END/ LEGAL */
#ifndef GRID_CONJUGATE_GRADIENT_MULTI_MIXED_PREC_H
#define GRID_CONJUGATE_GRADIENT_MULTI_MIXED_PREC_H

namespace Grid {

  //Mixed precision restarted defect correction CG
  template<class FieldD,class FieldF
//, typename std::enable_if< getPrecision<FieldD>::value == 2, int>::type = 0
//, typename std::enable_if< getPrecision<FieldF>::value == 1, int>::type = 0
> 
  class MixedPrecisionConjugateGradientMultiShift : public LinearFunction<FieldD> {
  public:                                                
//    RealD   Tolerance;
    Integer MaxInnerIterations;
    Integer MaxOuterIterations;
    GridBase* SinglePrecGrid; //Grid for single-precision fields
    RealD OuterLoopNormMult; //Stop the outer loop and move to a final double prec solve when the residual is OuterLoopNormMult * Tolerance
    LinearOperatorBase<FieldF> &Linop_f;
    LinearOperatorBase<FieldD> &Linop_d;
    MultiShiftFunction shifts;
    Integer iter;

    //Option to speed up *inner single precision* solves using a LinearFunction that produces a guess
//    LinearFunction<FieldF> *guesser;
    
    MixedPrecisionConjugateGradientMultiShift(GridBase* _sp_grid, LinearOperatorBase<FieldF> &_Linop_f, LinearOperatorBase<FieldD> &_Linop_d, 
Integer maxinnerit,	MultiShiftFunction &_shifts ) :
      Linop_f(_Linop_f), Linop_d(_Linop_d),
      MaxInnerIterations(maxinnerit), SinglePrecGrid(_sp_grid),
      OuterLoopNormMult(100.), shifts(_shifts) {};

  
    void operator() (const FieldD &src_d_in, FieldD &sol_d){
	assert(0); // not yet implemented
    }
    void operator() (const FieldD &src_d_in, std::vector<FieldD> &sol_d){
      GridStopWatch TotalTimer;
      TotalTimer.Start();
    
      int cb = src_d_in.checkerboard;

      int nshift = shifts.order;
      assert(nshift == sol_d.size());
      for(int i=0;i<nshift;i++) sol_d[i].checkerboard = cb;
    
      RealD src_norm = norm2(src_d_in);
//      RealD stop = src_norm * Tolerance*Tolerance;

      GridBase* DoublePrecGrid = src_d_in._grid;
      FieldD tmp_d(DoublePrecGrid); tmp_d.checkerboard = cb;
    
      FieldD tmp2_d(DoublePrecGrid); tmp2_d.checkerboard = cb;
    
      FieldD src_d(DoublePrecGrid);
      src_d = src_d_in; //source for next inner iteration, computed from residual during operation
    
//      RealD inner_tol = Tolerance;
  	FieldD psi_d(DoublePrecGrid);psi_d.checkerboard = cb;
    
      FieldF src_f(SinglePrecGrid);
      src_f.checkerboard = cb;
    
      std::vector<FieldF> sol_f(nshift,SinglePrecGrid);
      for(int i=0;i<nshift;i++) sol_f[i].checkerboard = cb;
    
//      ConjugateGradientShifted<FieldF> CG_f(inner_tol, MaxInnerIterations);
      ConjugateGradientMultiShift<FieldF> MSCG(MaxInnerIterations,shifts);
//      CG_f.ErrorOnNoConverge = false;

      GridStopWatch InnerCGtimer;

      GridStopWatch PrecChangeTimer;
    
{
//	std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradient: Outer iteration " <<outer_iter<<" residual "<< norm<< " target "<< stop<<std::endl;

//	if(norm < OuterLoopNormMult * stop){
//	  std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradient: Outer iteration converged on iteration " <<outer_iter <<std::endl;
//	  break;
//	}
//	while(norm * inner_tol * inner_tol < stop) inner_tol *= 2;  // inner_tol = sqrt(stop/norm) ??

	PrecChangeTimer.Start();
	precisionChange(src_f, src_d);
	PrecChangeTimer.Stop();
      
//	zeroit(sol_f);


	//Inner CG
	InnerCGtimer.Start();
  int if_relup = 0;
#if 0
        MSCG(Linop_f,src_f,sol_f);
#else
{
  
  GridBase *grid = SinglePrecGrid;
  
  ////////////////////////////////////////////////////////////////////////
  // Convenience references to the info stored in "MultiShiftFunction"
  ////////////////////////////////////////////////////////////////////////
  int nshift = shifts.order;


  std::vector<RealD> &mass(shifts.poles); // Make references to array in "shifts"
  std::vector<RealD> &mresidual(shifts.tolerances);
  std::vector<RealD> alpha(nshift,1.);
  std::vector<FieldF>   ps(nshift,grid);// Search directions

  assert(sol_f.size()==nshift);
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
  
  int cb=src_f.checkerboard;
  // Matrix mult fields
  FieldF r(grid); r.checkerboard = src_f.checkerboard;
  FieldF p(grid); p.checkerboard = src_f.checkerboard;
  FieldF tmp(grid); tmp.checkerboard = src_f.checkerboard;
  FieldF mmp(grid);mmp.checkerboard = src_f.checkerboard;
  FieldF psi(grid);psi.checkerboard = src_f.checkerboard;
    std::cout.precision(12);
    std::cout<<GridLogMessage<<"norm2(psi_d)= "<<norm2(psi_d)<<std::endl;
    std::cout<<GridLogMessage<<"norm2(psi)= "<<norm2(psi)<<std::endl;
  
  
  // Check lightest mass
  for(int s=0;s<nshift;s++){
    assert( mass[s]>= mass[primary] );
    converged[s]=0;
  }
  
  // Wire guess to zero
  // Residuals "r" are src
  // First search direction "p" is also src
  cp = norm2(src_f);
  Real c_relup = cp;
  for(int s=0;s<nshift;s++){
    rsq[s] = cp * mresidual[s] * mresidual[s];
    std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradientMultiShift: shift "<<s
	     <<" target resid "<<rsq[s]<<std::endl;
    ps[s] = src_f;
  }
  // r and p for primary
  r=src_f;
  p=src_f;
  
  //MdagM+m[0]
  std::cout << "p.checkerboard " << p.checkerboard
  << "mmp.checkerboard " << mmp.checkerboard << std::endl;

  Linop_f.HermOpAndNorm(p,mmp,d,qq);
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
  
 axpby(psi,0.,-bs[0],src_f,src_f);
  for(int s=0;s<nshift;s++) {
    axpby(sol_f[s],0.,-bs[s]*alpha[s],src_f,src_f);
  }
  
  
  // Iteration loop
  int k;
 // inefficient zeroing, please replace!
//  RealD sol_norm = axpy_norm(sol_d[0],-1.,sol_d[0],sol_d[0]);
  zeroit(sol_d[0]);
  std::cout<<GridLogMessage<<"norm(sol_d[0])= "<<norm2(sol_d[0])<<std::endl;
  

  int all_converged = 1;
	RealD tmp1,tmp2;
  for (k=1;k<=MaxOuterIterations;k++){
    
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
    
    Linop_f.HermOpAndNorm(p,mmp,d,qq);
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
    
    axpy(psi,-bs[0],ps[0],psi);
    for(int s=0;s<nshift;s++){
      int ss = s;
      // Scope for optimisation here in case of "single".
      // Could load sol_f[0] and pull all ps[s] in.
      //      if ( single ) ss=primary;
      // Bandwith saving in single case is Ls * 3 -> 2+Ls, so ~ 3x saving
      // Pipelined CG gain:
      //
      // New Kernel: Load r, vector of coeffs, vector of pointers ps
      // New Kernel: Load sol_f[0], vector of coeffs, vector of pointers ps
      // If can predict the coefficient bs then we can fuse these and avoid write reread cyce
      //  on ps[s].
      // Before:  3 x npole  + 3 x npole
      // After :  2 x npole (ps[s])        => 3x speed up of multishift CG.
      
      if( (!converged[s]) ) { 
	axpy(sol_f[ss],-bs[s]*alpha[s],ps[s],sol_f[ss]);
      }
    }


    if (k%MaxInnerIterations==0){
//    if (c < 1e-4*c_relup){
       RealD c_f=c;
       precisionChange(tmp_d,psi);
       RealD sol_norm =axpy_norm (psi_d,1.,tmp_d,psi_d);
       tmp1 = norm2(psi);
       zeroit(psi);
       tmp2 = norm2(psi);
       std::cout<<GridLogMessage<<"k= "<<k<<" norm2(sol)= "<<sol_norm<<" "<<tmp1<<" "<<tmp2<<std::endl;
//       precisionChange(sol_d[0],sol_f[0]);
       Linop_d.HermOpAndNorm(psi_d,tmp_d,tmp1,tmp2);
       axpy(tmp2_d,mass[0],psi_d,tmp_d);
       axpy(tmp_d,-1.,tmp2_d,src_d);
       precisionChange(r,tmp_d);
	c_relup = norm2(r);
       std::cout<<GridLogMessage<<"k= "<<k<<" norm2(r)= "<<c<<" "<<c_relup<<" "<<c_f<<std::endl;
	if_relup=1;
    }
    
    // Convergence checks
  all_converged=1;
    for(int s=0;s<nshift;s++){
      
      if ( (!converged[s]) ){
	
	RealD css  = c * z[s][iz]* z[s][iz];
	
	if(css<rsq[s]){
	  if ( ! converged[s] )
	    std::cout<<GridLogMessage<<"ConjugateGradientMultiShift k="<<k<<" Shift "<<s<<" has converged"<<std::endl;
	      converged[s]=1;
	} else {
		if (k%MaxInnerIterations==0)
	    std::cout<<GridLogMessage<<"ConjugateGradientMultiShift k="<<k<<" Shift "<<s<<" has not converged "<<css<<"<"<<rsq[s]<<std::endl;
	  all_converged=0;
	}

      }
    }
    
#if 0
    if ( all_converged ){
      std::cout<<GridLogMessage<< "CGMultiShift: All shifts have converged iteration "<<k<<std::endl;
#else
    if ( converged[0] ){
      std::cout<<GridLogMessage<< "CGMultiShift: Shift 0 have converged iteration, terminating  "<<k<<std::endl;
#endif
      
#if 1
      for(int s=1; s < nshift; s++) { 
	Linop_f.HermOpAndNorm(sol_f[s],mmp,d,qq);
	axpy(tmp,mass[s],sol_f[s],mmp);
	axpy(r,-alpha[s],src_f,tmp);
	RealD rn = norm2(r);
	RealD cn = norm2(src_f);
	std::cout<<GridLogMessage<<"CGMultiShift: shift["<<s<<"] true residual "<<std::sqrt(rn/cn)<<std::endl;
      }
#endif
     iter = k;
      break;
    }
  }
  // ugly hack
  if ( !all_converged )
  std::cout<<GridLogMessage<<"CG multi shift did not converge"<<std::endl;
//  assert(0);
}
	
#endif
	InnerCGtimer.Stop();
      
	//Convert sol back to double and add to double prec solution
	PrecChangeTimer.Start();
	sol_d[0]=psi_d;
	for(int i=1;i<nshift;i++)precisionChange(sol_d[i], sol_f[i]);
      std::cout<<GridLogMessage<< "CGMultiShift: Checking solutions"<<std::endl;
      // Check answers 
      for(int s=0; s < nshift; s++) { 
	RealD tmp1,tmp2;
       Linop_d.HermOpAndNorm(sol_d[s],tmp_d,tmp1,tmp2);
       axpy(tmp2_d,shifts.poles[s],sol_d[s],tmp_d);
       axpy(tmp_d,-1.,src_d,tmp2_d);
	std::cout<<GridLogMessage<<"CGMultiShift: shift["<<s<<"] true residual "<<std::sqrt(norm2(tmp_d)/norm2(src_d))<<std::endl;
      }
	PrecChangeTimer.Stop();
      
}
    
      //Final trial CG
 //     std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradient: Starting final patch-up double-precision solve"<<std::endl;
    
      TotalTimer.Stop();
      std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradient: Total " << TotalTimer.Elapsed() << " Precision change " << PrecChangeTimer.Elapsed() << " Inner CG total " << InnerCGtimer.Elapsed() << std::endl;
    }
  };

}

#endif
