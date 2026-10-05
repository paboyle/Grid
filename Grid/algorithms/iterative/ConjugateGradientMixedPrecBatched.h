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
  
  // Inner solves on independent partitions of the communicator (--batched-solver-split).
  // BatchedSplit is the partition MPI layout, empty for no split; BatchedSplitNode selects
  // one partition per node. Fixed at construction, default from the command line.
  const Coordinate BatchedSplit;
  const bool       BatchedSplitNode;

  // Copy of Linop_f on the partitions, made once at construction and owned; nullptr when
  // there is no split or the operator cannot be split. It holds the gauge field as it was at
  // construction: a solver must not outlive a change of the operator's gauge field. The
  // outer defect correction uses Linop_d, so a stale copy would cost convergence, not
  // correctness.
  SplitOperator<FieldF> *split = nullptr;
  int                    Partitions = 1;

  MixedPrecisionConjugateGradientBatched(RealD tol, 
          Integer maxinnerit, 
          Integer maxouterit, 
          Integer maxpatchit,
          GridBase* _sp_grid, 
          LinearOperatorBase<FieldF> &_Linop_f, 
          LinearOperatorBase<FieldD> &_Linop_d,
          bool _updateResidual=true,
          const Coordinate &_split=GridDefaultBatchedSolverSplit(),
          bool _split_node=GridDefaultBatchedSolverSplitNode()) :
    Linop_f(_Linop_f), Linop_d(_Linop_d),
    Tolerance(tol), InnerTolerance(tol), MaxInnerIterations(maxinnerit), MaxOuterIterations(maxouterit), MaxPatchupIterations(maxpatchit), SinglePrecGrid(_sp_grid),
    OuterLoopNormMult(100.), guesser(NULL), updateResidual(_updateResidual),
    BatchedSplit(_split), BatchedSplitNode(_split_node)
  {
    Coordinate layout = BatchedSolverSplitLayout(SinglePrecGrid,BatchedSplit,BatchedSplitNode,Partitions);
    if ( Partitions > 1 ) {
      double t0 = usecond();
      split = Linop_f.SplitClone(layout);
      double t1 = usecond();
      if ( split == nullptr ) {
        std::cout << GridLogMessage << "MixedPrecisionConjugateGradientBatched: operator cannot be split; serial inner solves" << std::endl;
      } else {
        std::cout << GridLogMessage << "MixedPrecisionConjugateGradientBatched: split clone " << (t1-t0)/1.0e6 << " s" << std::endl;
        HostMemoryReport(SinglePrecGrid,GridLogMessage,"MixedPrecisionConjugateGradientBatched: after split clone");
      }
    }
  };

  // Owns split
  MixedPrecisionConjugateGradientBatched(const MixedPrecisionConjugateGradientBatched &) = delete;

  MixedPrecisionConjugateGradientBatched &operator=(const MixedPrecisionConjugateGradientBatched &) = delete;

  ~MixedPrecisionConjugateGradientBatched(void)
  {
    delete split;
  }

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
    GRID_ASSERT(src_d_in.size() == sol_d.size());
    int NBatch = src_d_in.size();

    std::cout << GridLogMessage << "NBatch = " << NBatch << std::endl;

    Integer TotalOuterIterations = 0; //Number of restarts
    std::vector<Integer> TotalInnerIterations(NBatch,0);     //Number of inner CG iterations
    std::vector<Integer> TotalFinalStepIterations(NBatch,0); //Number of CG iterations in final patch-up step
  
    GridStopWatch TotalTimer;
    TotalTimer.Start();

    GridStopWatch InnerCGtimer;
    GridStopWatch PrecChangeTimer;
    GridStopWatch OuterResidualTimer;
    GridStopWatch PatchupTimer;
    
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
    
    SplitTimers splitTimers;
    if ( split != nullptr && split->Partitions > NBatch ) {
      std::cout << GridLogWarning << "MixedPrecisionConjugateGradientBatched: " << split->Partitions << " partitions for "
                << NBatch << " right-hand sides; the extra partitions only solve zero padding" << std::endl;
    }
    
    Integer &outer_iter = TotalOuterIterations; //so it will be equal to the final iteration count
      
    for(outer_iter = 0; outer_iter < MaxOuterIterations; outer_iter++){
      std::cout << GridLogMessage << std::endl;
      std::cout << GridLogMessage << "Outer iteration " << outer_iter << std::endl;
      HostMemoryReport(DoublePrecGrid,GridLogMessage,"MixedPrecisionConjugateGradientBatched: outer iteration "+std::to_string(outer_iter));
      
      bool allConverged = true;
      
      for (int i=0; i<NBatch; i++) {
        //Compute double precision rsd and also new RHS vector.
        OuterResidualTimer.Start();
        Linop_d.HermOp(sol_d[i], tmp_d);
        norm[i] = axpy_norm(src_d[i], -1., tmp_d, src_d_in[i]); //src_d is residual vector
        OuterResidualTimer.Stop();
        
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

      if ( split != nullptr ) {
        InnerSplitSolves(*split, CG_f, src_f, sol_f, TotalInnerIterations, InnerCGtimer, splitTimers);
      }

      for (int i=0; i<NBatch; i++) {
        //Inner CG
        if ( split == nullptr ) {
          InnerCGtimer.Start();
          CG_f(Linop_f, src_f[i], sol_f[i]);
          InnerCGtimer.Stop();
          TotalInnerIterations[i] += CG_f.IterationsToComplete;
        }
        
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
    
    PatchupTimer.Start();
    for (int i=0; i<NBatch; i++) {
      ConjugateGradient<FieldD> CG_d(Tolerance, MaxPatchupIterations);
      CG_d(Linop_d, src_d_in[i], sol_d[i]);
      TotalFinalStepIterations[i] += CG_d.IterationsToComplete;
    }
    PatchupTimer.Stop();

    TotalTimer.Stop();

    std::cout << GridLogMessage << std::endl;
    for (int i=0; i<NBatch; i++) {
      std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradientBatched: solve " << i << " Inner CG iterations " << TotalInnerIterations[i] << " Restarts " << TotalOuterIterations << " Final CG iterations " << TotalFinalStepIterations[i] << std::endl;
    }
    std::cout << GridLogMessage << std::endl;
    std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradientBatched: Total time " << TotalTimer.Elapsed() << " Precision change " << PrecChangeTimer.Elapsed() << " Inner CG total " << InnerCGtimer.Elapsed() << std::endl;
    std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradientBatched: Outer residual " << OuterResidualTimer.Elapsed() << " Patch-up " << PatchupTimer.Elapsed() << std::endl;
    if ( split != nullptr ) {
      std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradientBatched: Grid_split " << splitTimers.Split.Elapsed() << " (" << splitTimers.SplitCalls << " calls)"
               << " Grid_unsplit " << splitTimers.Unsplit.Elapsed() << " (" << splitTimers.UnsplitCalls << " calls)"
               << " Staging copies " << splitTimers.Staging.Elapsed()
               << " Host to device " << splitTimers.H2D.Elapsed()
               << " Iteration count sum " << splitTimers.Reduce.Elapsed() << std::endl;
    }
    double accounted = PrecChangeTimer.useconds() + InnerCGtimer.useconds() + OuterResidualTimer.useconds()
                     + PatchupTimer.useconds() + splitTimers.Split.useconds()
                     + splitTimers.Unsplit.useconds() + splitTimers.Staging.useconds() + splitTimers.H2D.useconds()
                     + splitTimers.Reduce.useconds();
    std::cout<<GridLogMessage<<"MixedPrecisionConjugateGradientBatched: Unaccounted " << (TotalTimer.useconds()-accounted)/1.0e6 << " s" << std::endl;
    
  }

private:

  // Coarse breakdown of the split path outside the inner CG
  struct SplitTimers {
    GridStopWatch Split;
    GridStopWatch Unsplit;
    GridStopWatch Staging;
    GridStopWatch H2D;
    GridStopWatch Reduce;
    int SplitCalls   = 0;
    int UnsplitCalls = 0;
  };

  ////////////////////////////////////////////////////////////////////////////////////////
  // One inner solve per partition, in groups of Partitions right-hand sides. The last
  // group is zero-padded; a zero source returns at once from CG. Collective.
  ////////////////////////////////////////////////////////////////////////////////////////
  void InnerSplitSolves(SplitOperator<FieldF> &split,
                        ConjugateGradient<FieldF> &CG_f,
                        std::vector<FieldF> &src_f,
                        std::vector<FieldF> &sol_f,
                        std::vector<Integer> &TotalInnerIterations,
                        GridStopWatch &InnerCGtimer,
                        SplitTimers &timers)
  {
    int NBatch = src_f.size();
    int P      = split.Partitions;
    int cb     = src_f[0].Checkerboard();

    FieldF s_src(split.FieldGrid);
    FieldF s_sol(split.FieldGrid);
    std::vector<FieldF> group_src(P,SinglePrecGrid);
    std::vector<FieldF> group_sol(P,SinglePrecGrid);
    std::vector<uint64_t> iters(P);

    for(int g=0;g<NBatch;g+=P){

      // Gather the group; zero-pad past the end of the batch
      timers.Staging.Start();
      for(int p=0;p<P;p++){
        group_src[p].Checkerboard() = cb;
        group_sol[p].Checkerboard() = cb;
        if ( g+p < NBatch ) {
          group_src[p] = src_f[g+p];
          group_sol[p] = sol_f[g+p];
        } else {
          group_src[p] = Zero();
          group_sol[p] = Zero();
        }
      }
      timers.Staging.Stop();

      // The initial guess (e.g. from the guesser) travels with the source
      timers.Split.Start();
      Grid_split(group_src,s_src);
      Grid_split(group_sol,s_sol);
      timers.Split.Stop();
      timers.SplitCalls += 2;

      // Grid_split leaves its output on the host; move it now so the copy is not in the CG time
      timers.H2D.Start();
      {
        autoView(s_src_v, s_src, AcceleratorRead);
        autoView(s_sol_v, s_sol, AcceleratorRead);
      }
      timers.H2D.Stop();

      InnerCGtimer.Start();
      CG_f(*split.Linop,s_src,s_sol);
      InnerCGtimer.Stop();

      timers.Unsplit.Start();
      Grid_unsplit(group_sol,s_sol);
      timers.Unsplit.Stop();
      timers.UnsplitCalls += 1;

      // One iteration count per partition, contributed by the partition's rank 0 only
      timers.Reduce.Start();
      for(int p=0;p<P;p++){
        iters[p] = 0;
      }
      if ( split.FieldGrid->ThisRank() == 0 ) {
        iters[split.Partition] = CG_f.IterationsToComplete;
      }
      SinglePrecGrid->GlobalSumVector(&iters[0],P);
      timers.Reduce.Stop();

      // Includes the host to device copy of group_sol left by Grid_unsplit
      timers.Staging.Start();
      for(int p=0;p<P && g+p<NBatch;p++){
        sol_f[g+p] = group_sol[p];
        TotalInnerIterations[g+p] += iters[p];
      }
      timers.Staging.Stop();
    }
  }
};

NAMESPACE_END(Grid);

#endif
