/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/multigrid/MrhsPreconditioner.h

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See the full license in the file "LICENSE" in the top level distribution
    directory
*************************************************************************************/
/*  END LEGAL */
#pragma once

NAMESPACE_BEGIN(Grid);

//////////////////////////////////////////////////////////////////////
// The mrhs interfaces.  This header is in the Algorithms.h include
// chain (AdefMrhs.h needs it); the concrete mrhs classes are in
// MrhsMultiGrid.h, reached through PVdagMMultiGrid.h / HDCGMultiGrid.h.
//
// MrhsLinearFunction: vector-of-fields in, vector out.  The outer level
// carries mrhs as std::vector<Field>; below it mrhs is PACKED into a
// single D+1 field (rhs = dim 0) and the coarse classes are plain
// LinearFunctions on that.
//////////////////////////////////////////////////////////////////////
template<class Field>
class MrhsLinearFunction {
public:
  virtual void operator()(std::vector<Field> &in, std::vector<Field> &out) = 0;
  virtual ~MrhsLinearFunction(){};
};

//////////////////////////////////////////////////////////////////////
// An mrhs preconditioner for the two-level CG family (AdefMrhs.h):
// operator() is M1; Vstart chooses x_0 from the source (Tang's ADEF-2
// start x_0 = Q b; the default is the zero start).  The timers are the
// M1 breakdown the outer solver reports.
//////////////////////////////////////////////////////////////////////
template<class Field>
class MrhsPreconditioner : public MrhsLinearFunction<Field> {
public:
  GridStopWatch ProjectTimer, PromoteTimer, DeflateTimer, CoarseTimer, FineTimer, SmoothTimer;
  virtual void Vstart(std::vector<Field> &x, std::vector<Field> &src)
  {
    for(int r=0;r<(int)x.size();r++) x[r] = Zero();
  }
  virtual void ResetTimers(void)
  {
    ProjectTimer.Reset(); PromoteTimer.Reset(); DeflateTimer.Reset();
    CoarseTimer.Reset();  FineTimer.Reset();    SmoothTimer.Reset();
  }
  virtual void ReportTimers(const std::string &prefix)
  {
    std::cout<<GridLogMessage<<"**** M1 breakdown:"<<std::endl;
    std::cout<<GridLogMessage<<prefix<<"Project "<<ProjectTimer.Elapsed()<<std::endl;
    std::cout<<GridLogMessage<<prefix<<"Promote "<<PromoteTimer.Elapsed()<<std::endl;
    std::cout<<GridLogMessage<<prefix<<"Deflate "<<DeflateTimer.Elapsed()<<std::endl;
    std::cout<<GridLogMessage<<prefix<<"Coarse  "<<CoarseTimer.Elapsed()<<std::endl;
    std::cout<<GridLogMessage<<prefix<<"Fine    "<<FineTimer.Elapsed()<<std::endl;
    std::cout<<GridLogMessage<<prefix<<"Smooth  "<<SmoothTimer.Elapsed()<<std::endl;
  }
  virtual ~MrhsPreconditioner(){};
};

NAMESPACE_END(Grid);
