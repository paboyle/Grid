    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/ConjugateGradientMixedPrec.h

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
#ifndef GRID_CONJUGATE_GRADIENT_MIXED_PREC_H
#define GRID_CONJUGATE_GRADIENT_MIXED_PREC_H

namespace Grid {

  //Mixed precision restarted defect correction CG
  template<class FieldD,class FieldF, typename std::enable_if< getPrecision<FieldD>::value == 2, int>::type = 0,typename std::enable_if< getPrecision<FieldF>::value == 1, int>::type = 0> 
  class MixedPrecisionConjugateGradient : public LinearFunction<FieldD> {
  public:                                                
    RealD   Tolerance;
    RealD   InnerTolerance; //Initial tolerance for inner CG. Defaults to Tolerance but can be changed
    Integer MaxInnerIterations;
    Integer MaxOuterIterations;
    GridBase* SinglePrecGrid; //Grid for single-precision fields
    RealD OuterLoopNormMult; //Stop the outer loop and move to a final double prec solve when the residual is OuterLoopNormMult * Tolerance
    LinearOperatorBase<FieldF> &Linop_f;
    LinearOperatorBase<FieldD> &Linop_d;

    Integer TotalInnerIterations; //Number of inner CG iterations
    Integer TotalOuterIterations; //Number of restarts
    Integer TotalFinalStepIterations; //Number of CG iterations in final patch-up step

    //Option to speed up *inner single precision* solves using a LinearFunction that produces a guess
    LinearFunction<FieldF> *guesser;
    
    MixedPrecisionConjugateGradient(RealD tol, Integer maxinnerit, Integer maxouterit, GridBase* _sp_grid, LinearOperatorBase<FieldF> &_Linop_f, LinearOperatorBase<FieldD> &_Linop_d) :
      Linop_f(_Linop_f), Linop_d(_Linop_d),
      Tolerance(tol), InnerTolerance(tol), MaxInnerIterations(maxinnerit), MaxOuterIterations(maxouterit), SinglePrecGrid(_sp_grid),
      OuterLoopNormMult(100.), guesser(NULL){ };

    void useGuesser(LinearFunction<FieldF> &g){
      guesser = &g;
    }
  
    void operator() (const FieldD &src_d_in, FieldD &sol_d){

      TotalInnerIterations = 0;
	
      GridStopWatch TotalTimer;
      TotalTimer.Start();
    
      int cb = src_d_in.checkerboard;
      sol_d.checkerboard = cb;
    
      RealD src_norm = norm2(src_d_in);
      RealD stop = src_norm * Tolerance*Tolerance;

      GridBase* DoublePrecGrid = src_d_in._grid;
      FieldD tmp_d(DoublePrecGrid);
      tmp_d.checkerboard = cb;
    
      FieldD tmp2_d(DoublePrecGrid);
      tmp2_d.checkerboard = cb;
    
      FieldD src_d(DoublePrecGrid);
      src_d = src_d_in; //source for next inner iteration, computed from residual during operation
    
      RealD inner_tol = InnerTolerance;
    
      FieldF src_f(SinglePrecGrid);
      src_f.checkerboard = cb;
    
      FieldF sol_f(SinglePrecGrid);
      sol_f.checkerboard = cb;
    
      ConjugateGradient<FieldF> CG_f(inner_tol, MaxInnerIterations);
      CG_f.ErrorOnNoConverge = false;

      GridStopWatch InnerCGtimer;

      GridStopWatch PrecChangeTimer;
    
      Integer &outer_iter = TotalOuterIterations; //so it will be equal to the final iteration count
      
      for(outer_iter = 0; outer_iter < MaxOuterIterations; outer_iter++){
	//Compute double precision rsd and also new RHS vector.
	Linop_d.HermOp(sol_d, tmp_d);
	RealD norm = axpy_norm(src_d, -1., tmp_d, src_d_in); //src_d is residual vector
      
	std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradient: Outer iteration " <<outer_iter<<" residual "<< norm<< " target "<< stop<<std::endl;

	if(norm < OuterLoopNormMult * stop){
	  std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradient: Outer iteration converged on iteration " <<outer_iter <<std::endl;
	  break;
	}
	while(norm * inner_tol * inner_tol < stop) inner_tol *= 2;  // inner_tol = sqrt(stop/norm) ??

	PrecChangeTimer.Start();
	precisionChange(src_f, src_d);
	PrecChangeTimer.Stop();
      
	zeroit(sol_f);

	//Optionally improve inner solver guess (eg using known eigenvectors)
	if(guesser != NULL)
	  (*guesser)(src_f, sol_f);

	//Inner CG
	CG_f.Tolerance = inner_tol;
	InnerCGtimer.Start();
	CG_f(Linop_f, src_f, sol_f);
	InnerCGtimer.Stop();
	TotalInnerIterations += CG_f.IterationsToComplete;
      
	//Convert sol back to double and add to double prec solution
	PrecChangeTimer.Start();
	precisionChange(tmp_d, sol_f);
	PrecChangeTimer.Stop();
      
	axpy(sol_d, 1.0, tmp_d, sol_d);
      }
    
      //Final trial CG
      std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradient: Starting final patch-up double-precision solve"<<std::endl;
    
      ConjugateGradient<FieldD> CG_d(Tolerance, MaxInnerIterations);
      CG_d(Linop_d, src_d_in, sol_d);
      TotalFinalStepIterations = CG_d.IterationsToComplete;

      TotalTimer.Stop();
      std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradient: Inner CG iterations " << TotalInnerIterations << " Restarts " << TotalOuterIterations << " Final CG iterations " << TotalFinalStepIterations << std::endl;
      std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradient: Total time " << TotalTimer.Elapsed() << " Precision change " << PrecChangeTimer.Elapsed() << " Inner CG total " << InnerCGtimer.Elapsed() << std::endl;
    }
  };

}

#endif
