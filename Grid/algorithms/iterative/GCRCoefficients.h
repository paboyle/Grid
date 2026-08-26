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
  std::vector<ComplexD>               a_sum;   // [k]
  std::vector<int>                    a_n;
  std::vector<std::vector<ComplexD> > b_sum;   // [k][j]
  std::vector<std::vector<int> >      b_n;
  void RecordA(int k, ComplexD a){
    if ( (int)a_sum.size() <= k ) { a_sum.resize(k+1,ComplexD(0.0)); a_n.resize(k+1,0); }
    a_sum[k] += a; a_n[k]++;
  }
  void RecordB(int k, const std::vector<ComplexD> &b){
    if ( (int)b_sum.size() <= k ) { b_sum.resize(k+1); b_n.resize(k+1); }
    if ( b_sum[k].size() < b.size() ) { b_sum[k].resize(b.size(),ComplexD(0.0)); b_n[k].resize(b.size(),0); }
    for(int j=0;j<(int)b.size();j++){ b_sum[k][j] += b[j]; b_n[k][j]++; }
  }
  int  Steps(void) const { return a_sum.size(); }
  int  Calls(void) const { return a_n.size() ? a_n[0] : 0; }
  ComplexD A(int k) const { return a_sum[k]/(double)a_n[k]; }
  int  NB(int k) const { return (k<(int)b_sum.size()) ? b_sum[k].size() : 0; }
  ComplexD B(int k,int j) const { return b_sum[k][j]/(double)b_n[k][j]; }
  void Report(const std::string &name) const {
    std::cout << GridLogMessage << "GCRCoefficients " << name << ": " << Calls() << " calls, " << Steps() << " steps, mmax " << mmax << std::endl;
    for(int k=0;k<Steps();k++){
      std::cout << GridLogMessage << "  step " << k << " a=(" << real(A(k)) << "," << imag(A(k)) << ")";
      for(int j=0;j<NB(k);j++) std::cout << " b[" << j << "]=(" << real(B(k,j)) << "," << imag(B(k,j)) << ")";
      std::cout << std::endl;
    }
  }
};

NAMESPACE_END(Grid);
