/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/SplitGridBlockKrylovSchur.h

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Chulwoo Jung <chulwoo@bnl.gov>

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
#ifndef GRID_SPLIT_BLOCK_KRYLOV_SCHUR_H
#define GRID_SPLIT_BLOCK_KRYLOV_SCHUR_H

NAMESPACE_BEGIN(Grid);

/**
 * Block Krylov-Schur eigensolver with Split-Grid batched operator application.
 *
 * Derives from BlockKrylovSchur<Field> and overrides only the operator
 * application: instead of calling Linop.Op() once per vector, mrhs vectors
 * are packed onto a smaller split grid with Grid_split, the polynomial-
 * filtered operator poly(SLinop) is applied once to the combined field, and
 * the results are unpacked with Grid_unsplit.  All other algorithmic logic
 * (Arnoldi orthogonalisation, Schur restart, convergence test, etc.) is
 * inherited unchanged from the base class.
 *
 * Constructor extras (beyond BlockKrylovSchur)
 * -------------------------------------------
 * SLinop  : split-grid linear operator used inside poly
 * poly    : polynomial filter OperatorFunction applied on the split grid
 * SGrid   : split-grid GridBase (mrhs fields packed here)
 * mrhs    : RHS batched per split-grid call; Nblock must be divisible by mrhs
 *
 * Notes
 * -----
 * - The Rayleigh quotient H is built from poly(A), so its Ritz values are
 *   eigenvalues of poly(A).  Convergence is assessed against the full-grid
 *   Linop (inherited from base), giving true A-eigenvalues.
 * - Grid_ (inherited) is the full grid; SGrid is the split grid.
 */
template<class Field>
class SplitGridBlockKrylovSchur : public BlockKrylovSchur<Field> {

  using Base = BlockKrylovSchur<Field>;

  // Bring protected base members into scope
  using Base::Nblock;
  using Base::Grid_;
  using Base::basis;

  // Split-grid extras
  LinearOperatorBase<Field>& SLinop;
  OperatorFunction<Field>&   poly;
  GridBase*                  SGrid;
  int                        mrhs;

public:

  SplitGridBlockKrylovSchur(LinearOperatorBase<Field>& _Linop,
                             LinearOperatorBase<Field>& _SLinop,
                             OperatorFunction<Field>&   _poly,
                             GridBase*                  _FGrid,
                             GridBase*                  _SGrid,
                             int                        _mrhs,
                             RealD                      _Tolerance,
                             RitzFilter                 _rf = EvalReSmall)
    : Base(_Linop, _FGrid, _Tolerance, _rf),
      SLinop(_SLinop), poly(_poly), SGrid(_SGrid), mrhs(_mrhs)
  {}

protected:

  // Validate mrhs divisibility before the Arnoldi loop starts
  void preRun() override
  {
    assert(Nblock % mrhs == 0 && "Nblock must be divisible by mrhs");
  }

  // Apply poly(A) to basis[kBase .. kBase+Nblock-1] via Grid_split / Grid_unsplit,
  // batching mrhs vectors per split-grid call.
  void applyBlock(std::vector<Field>& W, int kBase) override
  {
    std::vector<Field> in(mrhs, Field(Grid_));
    Field s_in(SGrid);
    Field s_out(SGrid);

    int k_start = 0;
    while (k_start < Nblock) {
      for (int u = 0; u < mrhs; u++)
        in[u] = basis[kBase + k_start + u];

      Grid_split(in, s_in);
      poly(SLinop, s_in, s_out);
      Grid_unsplit(in, s_out);

      for (int u = 0; u < mrhs; u++)
        W[k_start + u] = in[u];

      k_start += mrhs;
    }
  }
};

NAMESPACE_END(Grid);

#endif  // GRID_SPLIT_BLOCK_KRYLOV_SCHUR_H
