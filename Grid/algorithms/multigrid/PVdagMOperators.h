/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/multigrid/PVdagMOperators.h

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
// A = PV^dag M (non-Hermitian).  The operator the PVdagM multigrid
// coarsens and solves.
//
// SloppyComms delegates to BOTH underlying fermion operators, so one
// call controls every fine halo in the chain -- including any shifted
// wrapper built on the same two matrices.  Halo-precision POLICY:
// reduced-precision (fp32 wire) halos belong in the PRECONDITIONER --
// the smoother, the V-cycle's own residuals, the coarsening -- and
// NEVER in the outer Krylov, whose applications define what
// "converged" means.  The operators default to EXACT; the V-cycle
// turns sloppiness on for its own scope only (MrhsTwoLevelMG).
//////////////////////////////////////////////////////////////////////
template<class Matrix,class Field>
class PVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
public:
  PVdagMLinearOperator(Matrix &Mat,Matrix &PV):
    _Mat(Mat),_PV(PV)
  {};
  void SloppyComms(int sloppy)
  {
    _Mat.SloppyComms(sloppy);
    _PV.SloppyComms(sloppy);
  }
  void OpDiag (const Field &in, Field &out)
  {
    GRID_ASSERT(0);
  }
  void OpDir  (const Field &in, Field &out,int dir,int disp)
  {
    GRID_ASSERT(0);
  }
  void OpDirAll  (const Field &in, std::vector<Field> &out)
  {
    GRID_ASSERT(0);
  };
  void Op     (const Field &in, Field &out)
  {
    Field tmp(in.Grid());
    _Mat.M(in,tmp);
    _PV.Mdag(tmp,out);
  }
  void AdjOp  (const Field &in, Field &out)
  {
    Field tmp(in.Grid());
    _PV.M(in,tmp);
    _Mat.Mdag(tmp,out);
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2)
  {
    HermOp(in,out);
    ComplexD d=innerProduct(in,out);
    n1=real(d);
    n2=norm2(out);
  }
  void HermOp(const Field &in, Field &out)
  {
    Field tmp(in.Grid());
    Op(in,tmp);
    AdjOp(tmp,out);
  }
};

//////////////////////////////////////////////////////////////////////
// A + shift, built directly on the two matrices (fine level: avoids
// an extra fine-field pass through a generic shifted wrapper).
//////////////////////////////////////////////////////////////////////
template<class Matrix,class Field>
class ShiftedPVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
public:
  RealD shift;
  ShiftedPVdagMLinearOperator(RealD _shift,Matrix &Mat,Matrix &PV): shift(_shift),_Mat(Mat),_PV(PV){};
  void OpDiag (const Field &in, Field &out)
  {
    GRID_ASSERT(0);
  }
  void OpDir  (const Field &in, Field &out,int dir,int disp)
  {
    GRID_ASSERT(0);
  }
  void OpDirAll  (const Field &in, std::vector<Field> &out)
  {
    GRID_ASSERT(0);
  };
  void Op     (const Field &in, Field &out)
  {
    Field tmp(in.Grid());
    _Mat.M(in,tmp);
    _PV.Mdag(tmp,out);
    out = out + shift*in;
  }
  void AdjOp  (const Field &in, Field &out)
  {
    Field tmp(in.Grid());
    _PV.M(in,tmp);
    _Mat.Mdag(tmp,out);
    out = out + shift*in;
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2)
  {
    GRID_ASSERT(0);
  }
  void HermOp(const Field &in, Field &out){
    Field tmp(in.Grid());
    Op(in,tmp);
    AdjOp(tmp,out);
  }
};

//////////////////////////////////////////////////////////////////////
// Op + shift for any operator (coarse levels).
//////////////////////////////////////////////////////////////////////
template<class Field>
class ShiftedLinearOperator : public LinearOperatorBase<Field> {
  LinearOperatorBase<Field> &_Op;
  RealD shift;
public:
  ShiftedLinearOperator(RealD _shift, LinearOperatorBase<Field> &Op) : _Op(Op), shift(_shift) {}
  void OpDiag  (const Field &in, Field &out)
  {
    GRID_ASSERT(0);
  }
  void OpDir   (const Field &in, Field &out,int dir,int disp)
  {
    GRID_ASSERT(0);
  }
  void OpDirAll (const Field &in, std::vector<Field> &out)
  {
    GRID_ASSERT(0);
  }
  void Op      (const Field &in, Field &out)
  {
    _Op.Op(in,out);
    out = out + shift*in;
  }
  void AdjOp   (const Field &in, Field &out) {
    _Op.AdjOp(in,out);
    out = out + shift*in;
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2)
  {
    GRID_ASSERT(0);
  }
  void HermOp  (const Field &in, Field &out)
  {
    Field tmp(in.Grid());
    Op(in,tmp);
    AdjOp(tmp,out);
  }
};

NAMESPACE_END(Grid);
