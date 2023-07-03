/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/ConjugateGradientMixedPrecBatched.h

    Copyright (C) 2015

    Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk>

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
#ifndef GRID_CONJUGATE_GRADIENT_MIXED_PREC_BATCHED_H
#define GRID_CONJUGATE_GRADIENT_MIXED_PREC_BATCHED_H

NAMESPACE_BEGIN(Grid);

//Mixed precision restarted defect correction CG
template<class FieldD,class FieldF, 
  typename std::enable_if< getPrecision<FieldD>::value == 2, int>::type = 0,
  typename std::enable_if< getPrecision<FieldF>::value == 1, int>::type = 0> 
class MixedPrecisionConjugateGradientBatched : public LinearFunction<FieldD> {
public:
  using LinearFunction<FieldD>::operator();
  RealD   Tolerance;
  RealD   InnerTolerance; //Initial tolerance for inner CG. Defaults to Tolerance but can be changed
  Integer MaxInnerIterations;
  Integer MaxOuterIterations;
  Integer MaxPatchupIterations;
  GridBase* SinglePrecGrid; //Grid for single-precision fields
  RealD OuterLoopNormMult; //Stop the outer loop and move to a final double prec solve when the residual is OuterLoopNormMult * Tolerance
  LinearOperatorBase<FieldF> &Linop_f;
  LinearOperatorBase<FieldD> &Linop_d;

  //Option to speed up *inner single precision* solves using a LinearFunction that produces a guess
  LinearFunction<FieldF> *guesser;
  bool updateResidual;
  
  MixedPrecisionConjugateGradientBatched(RealD tol, 
          Integer maxinnerit, 
          Integer maxouterit, 
          Integer maxpatchit,
          GridBase* _sp_grid, 
          LinearOperatorBase<FieldF> &_Linop_f, 
          LinearOperatorBase<FieldD> &_Linop_d,
          bool _updateResidual=true) :
    Linop_f(_Linop_f), Linop_d(_Linop_d),
    Tolerance(tol), InnerTolerance(tol), MaxInnerIterations(maxinnerit), MaxOuterIterations(maxouterit), MaxPatchupIterations(maxpatchit), SinglePrecGrid(_sp_grid),
    OuterLoopNormMult(100.), guesser(NULL), updateResidual(_updateResidual) { };

  void useGuesser(LinearFunction<FieldF> &g){
    guesser = &g;
  }
  
  void operator() (const FieldD &src_d_in, FieldD &sol_d){
    std::vector<FieldD> srcs_d_in{src_d_in};
    std::vector<FieldD> sols_d{sol_d};

    (*this)(srcs_d_in,sols_d);

    sol_d = sols_d[0];
  }

  void operator() (const std::vector<FieldD> &src_d_in, std::vector<FieldD> &sol_d){
    assert(src_d_in.size() == sol_d.size());
    int NBatch = src_d_in.size();

    std::cout << GridLogMessage << "NBatch = " << NBatch << std::endl;

    Integer TotalOuterIterations = 0; //Number of restarts
    std::vector<Integer> TotalInnerIterations(NBatch,0);     //Number of inner CG iterations
    std::vector<Integer> TotalFinalStepIterations(NBatch,0); //Number of CG iterations in final patch-up step
  
    GridStopWatch TotalTimer;
    TotalTimer.Start();

    GridStopWatch InnerCGtimer;
    GridStopWatch PrecChangeTimer;
    
    int cb = src_d_in[0].Checkerboard();
    
    std::vector<RealD> src_norm;
    std::vector<RealD> norm;
    std::vector<RealD> stop;
    
    GridBase* DoublePrecGrid = src_d_in[0].Grid();
    FieldD tmp_d(DoublePrecGrid);
    tmp_d.Checkerboard() = cb;
    
    FieldD tmp2_d(DoublePrecGrid);
    tmp2_d.Checkerboard() = cb;

    std::vector<FieldD> src_d;
    std::vector<FieldF> src_f;
    std::vector<FieldF> sol_f;

    for (int i=0; i<NBatch; i++) {
      sol_d[i].Checkerboard() = cb;

      src_norm.push_back(norm2(src_d_in[i]));
      norm.push_back(0.);
      stop.push_back(src_norm[i] * Tolerance*Tolerance);

      src_d.push_back(src_d_in[i]); //source for next inner iteration, computed from residual during operation

      src_f.push_back(SinglePrecGrid);
      src_f[i].Checkerboard() = cb;

      sol_f.push_back(SinglePrecGrid);
      sol_f[i].Checkerboard() = cb;
    }
    
    RealD inner_tol = InnerTolerance;
    
    ConjugateGradient<FieldF> CG_f(inner_tol, MaxInnerIterations);
    CG_f.ErrorOnNoConverge = false;
    
    Integer &outer_iter = TotalOuterIterations; //so it will be equal to the final iteration count
      
    for(outer_iter = 0; outer_iter < MaxOuterIterations; outer_iter++){
      std::cout << GridLogMessage << std::endl;
      std::cout << GridLogMessage << "Outer iteration " << outer_iter << std::endl;
      
      bool allConverged = true;
      
      for (int i=0; i<NBatch; i++) {
        //Compute double precision rsd and also new RHS vector.
        Linop_d.HermOp(sol_d[i], tmp_d);
        norm[i] = axpy_norm(src_d[i], -1., tmp_d, src_d_in[i]); //src_d is residual vector
        
        std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradientBatched: Outer iteration " << outer_iter <<" solve " << i << " residual "<< norm[i] << " target "<< stop[i] <<std::endl;

        PrecChangeTimer.Start();
        precisionChange(src_f[i], src_d[i]);
        PrecChangeTimer.Stop();
        
        sol_f[i] = Zero();
      
        if(norm[i] > OuterLoopNormMult * stop[i]) {
          allConverged = false;
        }
      }
      if (allConverged) break;

      if (updateResidual) {
        RealD normMax = *std::max_element(std::begin(norm), std::end(norm));
        RealD stopMax = *std::max_element(std::begin(stop), std::end(stop));
        while( normMax * inner_tol * inner_tol < stopMax) inner_tol *= 2;  // inner_tol = sqrt(stop/norm) ??
        CG_f.Tolerance = inner_tol;
      }

      //Optionally improve inner solver guess (eg using known eigenvectors)
      if(guesser != NULL) {
        (*guesser)(src_f, sol_f);
      }

      for (int i=0; i<NBatch; i++) {
        //Inner CG
        InnerCGtimer.Start();
        CG_f(Linop_f, src_f[i], sol_f[i]);
        InnerCGtimer.Stop();
        TotalInnerIterations[i] += CG_f.IterationsToComplete;
        
        //Convert sol back to double and add to double prec solution
        PrecChangeTimer.Start();
        precisionChange(tmp_d, sol_f[i]);
        PrecChangeTimer.Stop();
        
        axpy(sol_d[i], 1.0, tmp_d, sol_d[i]);
      }

    }
    
    //Final trial CG
    std::cout << GridLogMessage << std::endl;
    std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradientBatched: Starting final patch-up double-precision solve"<<std::endl;
    
    for (int i=0; i<NBatch; i++) {
      ConjugateGradient<FieldD> CG_d(Tolerance, MaxPatchupIterations);
      CG_d(Linop_d, src_d_in[i], sol_d[i]);
      TotalFinalStepIterations[i] += CG_d.IterationsToComplete;
    }

    TotalTimer.Stop();

    std::cout << GridLogMessage << std::endl;
    for (int i=0; i<NBatch; i++) {
      std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradientBatched: solve " << i << " Inner CG iterations " << TotalInnerIterations[i] << " Restarts " << TotalOuterIterations << " Final CG iterations " << TotalFinalStepIterations[i] << std::endl;
    }
    std::cout << GridLogMessage << std::endl;
    std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradientBatched: Total time " << TotalTimer.Elapsed() << " Precision change " << PrecChangeTimer.Elapsed() << " Inner CG total " << InnerCGtimer.Elapsed() << std::endl;
    
  }
};

NAMESPACE_END(Grid);

#endif
