/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/iterative/GCRCoefficients.h

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

//////////////////////////////////////////////////////////////////////////////
// Recorded GCR coefficients: per-step means over calls of the step length
// a_k and the orthogonalisation coefficients b_kj (already scaled and
// signed as applied: p_{k+1} = r + sum_j b_kj p_{k-j}).
//////////////////////////////////////////////////////////////////////////////
struct GCRCoefficients {
  int mmax = 0;
  // Every recorded call is kept: calls[c] = list of (a_k, [b_kj]) per step.
  // A(k)/B(k,j) return the coefficients of the SELECTED call: by default the
  // last complete one.  Selection "mean" averages coefficients over calls --
  // kept for comparison only: the mean of the coefficients of a nonlinear
  // recurrence is not the mean of the polynomials, and in practice (Frontier
  // M3, 2026-08-26) it was worse than every individual call.
  enum Select { Last=0, First=1, Index=2, Mean=3 };
  Select select = Last;
  int    index  = 0;
  typedef std::vector<std::pair<ComplexD,std::vector<ComplexD> > > Call;
  std::vector<Call> calls;
  Call current;
  void RecordA(int k, ComplexD a){
    if ( k==0 && current.size() ) { calls.push_back(current); current.clear(); }
    if ( (int)current.size() <= k ) current.resize(k+1);
    current[k].first = a;
  }
  void RecordB(int k, const std::vector<ComplexD> &b){
    if ( (int)current.size() <= k ) current.resize(k+1);
    current[k].second = b;
  }
  void Flush(void){ if ( current.size() ) { calls.push_back(current); current.clear(); } }
  int  Calls(void) const { return calls.size() + (current.size() ? 1 : 0); }
  const Call & Chosen(void) const {
    GRID_ASSERT( calls.size() || current.size() );
    if ( calls.empty() ) return current;
    if ( select==First ) return calls.front();
    if ( select==Index ) { GRID_ASSERT(index>=0 && index<(int)calls.size()); return calls[index]; }
    return calls.back();
  }
  int  Steps(void) const { return select==Mean ? MeanSteps() : Chosen().size(); }
  int  NB(int k) const   { return select==Mean ? MeanNB(k)   : Chosen()[k].second.size(); }
  ComplexD A(int k) const    { return select==Mean ? MeanA(k)   : Chosen()[k].first; }
  ComplexD B(int k,int j) const { return select==Mean ? MeanB(k,j) : Chosen()[k].second[j]; }
  // mean over calls (comparison only)
  int MeanSteps(void) const { int m=0; for(auto &c:calls) m = std::max(m,(int)c.size()); return m; }
  int MeanNB(int k) const { int m=0; for(auto &c:calls) if(k<(int)c.size()) m = std::max(m,(int)c[k].second.size()); return m; }
  ComplexD MeanA(int k) const { ComplexD s(0.0); int n=0; for(auto &c:calls) if(k<(int)c.size()){ s+=c[k].first; n++; } return s/(double)n; }
  ComplexD MeanB(int k,int j) const { ComplexD s(0.0); int n=0; for(auto &c:calls) if(k<(int)c.size() && j<(int)c[k].second.size()){ s+=c[k].second[j]; n++; } return s/(double)n; }
  void Report(const std::string &name) const {
    const char *sel[4]={"last","first","index","mean"};
    std::cout << GridLogMessage << "GCRCoefficients " << name << ": " << Calls() << " calls, " << Steps() << " steps, mmax " << mmax
              << ", selection " << sel[select] << std::endl;
    for(int k=0;k<Steps();k++){
      std::cout << GridLogMessage << "  step " << k << " a=(" << real(A(k)) << "," << imag(A(k)) << ")";
      for(int j=0;j<NB(k);j++) std::cout << " b[" << j << "]=(" << real(B(k,j)) << "," << imag(B(k,j)) << ")";
      // spread of a_k across calls: how different the individual polynomials are
      if ( calls.size()>1 ) {
        RealD lo=1e300, hi=0; for(auto &c:calls) if(k<(int)c.size()){ RealD x=real(c[k].first), y=imag(c[k].first); RealD m=std::sqrt(x*x+y*y); lo=std::min(lo,m); hi=std::max(hi,m); }
        std::cout << "   |a| over calls [" << lo << "," << hi << "]";
      }
      std::cout << std::endl;
    }
  }
};

NAMESPACE_END(Grid);
