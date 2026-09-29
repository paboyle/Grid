/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/algorithms/iterative/PrecGeneralisedConjugateResidual.h

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
#ifndef GRID_PREC_GCR_H
#define GRID_PREC_GCR_H
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>

NAMESPACE_BEGIN(Grid);

///////////////////////////////////////////////////////////////////////////////////////////////////////
// Hermitian spelling of PrecGeneralisedConjugateResidualNonHermitian.
//
// The algorithm is the same; only the operator entry point differs.  The solver
// drives Linop.Op(), while a caller of this class means "invert HermOp", so the
// operator is wrapped in HermOpAdaptor.  The wrap is what makes the two classes
// one: for an operator whose Op() already IS its HermOp() (HermitianLinearOperator,
// Gamma5R5HermitianLinearOperator, HermOpAdaptor) it changes nothing, and for one
// where they differ (MdagMLinearOperator: Op()=M, HermOp()=MdagM) it preserves the
// Hermitian meaning the caller asked for.
//
// A non-zero initial guess is respected.  Callers that guarantee x0 = 0 can call
// SetZeroGuess(1) to skip the first r0 = src - A x0 apply.
///////////////////////////////////////////////////////////////////////////////////////////////////////
template<class Field>
class PrecGeneralisedConjugateResidual : public LinearFunction<Field> {
  // Declaration order is initialisation order: the adaptor is complete before
  // the solver binds its reference to it.
  HermOpAdaptor<Field>                                HermLinop;
  PrecGeneralisedConjugateResidualNonHermitian<Field> GCR;
public:
  using LinearFunction<Field>::operator();

  PrecGeneralisedConjugateResidual(RealD tol,Integer maxit,LinearOperatorBase<Field> &Linop,
                                   LinearFunction<Field> &Prec,int mmax,int nstep)
    : HermLinop(Linop),
      GCR(tol,maxit,HermLinop,Prec,mmax,nstep) {};

  void operator() (const Field &src, Field &psi) { GCR(src,psi); };

  void Name(std::string n)                       { GCR.Name(n); };
  void Level(int n)                              { GCR.Level(n); };
  void SetZeroGuess(int z)                       { GCR.SetZeroGuess(z); };
  int  Steps(void) const                         { return GCR.Steps(); };
  void LogCoefficients(int l)                    { GCR.LogCoefficients(l); };
  void SetCoefficientRecorder(GCRCoefficients *r){ GCR.SetCoefficientRecorder(r); };
  void ReleaseHistory(void)                      { GCR.ReleaseHistory(); };
};

NAMESPACE_END(Grid);
#endif
