/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/algorithms/multigrid/MrhsPromotedOperator.h

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

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

///////////////////////////////////////////////////////////////////////////////
// Present a D dimensional operator as a D+1 dimensional operator with Nrhs in
// dimension 0.  Field type is unchanged; only the Grid differs, so this is a
// LinearOperatorBase and callers need no template parameter: a native mrhs
// operator derives from the same base and substitutes without a call site
// change.
//
// Slices in and out around the wrapped operator.  No arithmetic beyond the
// wrapped call, but one ExtractSliceFast/InsertSliceFast pair per right hand
// side: data motion, not work.
//
// AdjOp is carried so that A^dag may be coarsened as a separate coarse
// operator when needed, rather than doubling coarse storage.
///////////////////////////////////////////////////////////////////////////////
template<class Field>
class MrhsPromotedOperator : public LinearOperatorBase<Field>
{
private:

  LinearOperatorBase<Field> &_LinOp;
  GridBase                  *_LowGrid;
  int                        _Nrhs;

public:

  MrhsPromotedOperator(LinearOperatorBase<Field> &LinOp,GridBase *LowGrid,int Nrhs)
    : _LinOp(LinOp), _LowGrid(LowGrid), _Nrhs(Nrhs)
  {
    GRID_ASSERT(_Nrhs>=1);
  }

  GridBase *LowGrid(void) { return _LowGrid; }
  int       Nrhs(void)    { return _Nrhs;    }

  // Reset on each call; retrieve and accumulate in the caller
  RealD tslice;
  RealD top;

  void OpDiag (const Field &in, Field &out)
  {
    SliceLoop(in,out,[&](Field &i,Field &o){ _LinOp.OpDiag(i,o); });
  }

  void Op     (const Field &in, Field &out)
  {
    SliceLoop(in,out,[&](Field &i,Field &o){ _LinOp.Op(i,o); });
  }

  void AdjOp  (const Field &in, Field &out)
  {
    SliceLoop(in,out,[&](Field &i,Field &o){ _LinOp.AdjOp(i,o); });
  }

  void HermOp (const Field &in, Field &out)
  {
    SliceLoop(in,out,[&](Field &i,Field &o){ _LinOp.HermOp(i,o); });
  }

  void OpDir  (const Field &in, Field &out,int dir,int disp)
  {
    SliceLoop(in,out,[&](Field &i,Field &o){ _LinOp.OpDir(i,o,dir,disp); });
  }

  //////////////////////////////////////////////////////////////////
  // Norms of the D+1 field are the sums over slices
  //////////////////////////////////////////////////////////////////
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2)
  {
    Conformable(in,out);
    Field lo_in (_LowGrid);
    Field lo_out(_LowGrid);
    n1=0.0;
    n2=0.0;
    for(int r=0;r<_Nrhs;r++){
      RealD r1,r2;
      ExtractSliceFast(lo_in,in,r,0);
      _LinOp.HermOpAndNorm(lo_in,lo_out,r1,r2);
      InsertSliceFast(lo_out,out,r,0);
      n1=n1+r1;
      n2=n2+r2;
    }
  }

  void OpDirAll(const Field &in, std::vector<Field> &out)
  {
    int npoint = out.size();
    Field lo_in(_LowGrid);
    std::vector<Field> lo_out(npoint,_LowGrid);
    for(int r=0;r<_Nrhs;r++){
      ExtractSliceFast(lo_in,in,r,0);
      _LinOp.OpDirAll(lo_in,lo_out);
      for(int p=0;p<npoint;p++){
        InsertSliceFast(lo_out[p],out[p],r,0);
      }
    }
  }

private:

  void Conformable(const Field &in,const Field &out)
  {
    conformable(in.Grid(),out.Grid());
    GRID_ASSERT(in.Grid()->_ndimension == _LowGrid->_ndimension+1);
    GRID_ASSERT(in.Grid()->_fdimensions[0] == _Nrhs);
  }

  template<class Kernel>
  void SliceLoop(const Field &in,Field &out,Kernel K)
  {
    Conformable(in,out);
    Field lo_in (_LowGrid);
    Field lo_out(_LowGrid);
    tslice=0.0;
    top=0.0;
    for(int r=0;r<_Nrhs;r++){
      tslice-=usecond();
      ExtractSliceFast(lo_in,in,r,0);
      tslice+=usecond();
      top-=usecond();
      K(lo_in,lo_out);
      top+=usecond();
      tslice-=usecond();
      InsertSliceFast(lo_out,out,r,0);
      tslice+=usecond();
    }
  }

};

NAMESPACE_END(Grid);
