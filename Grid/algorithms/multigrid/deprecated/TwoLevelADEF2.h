    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./Grid/algorithms/multigrid/deprecated/TwoLevelADEF2.h

    Copyright (C) 2015

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
#pragma once

// DEPRECATED with the V1 coarse operators: the single-RHS ADEF-2 on an
// Aggregation (ProjectToSubspace/PromoteFromSubspace) and a D-dimensional
// coarse operator.  The mrhs solver in AdefMrhs.h covers one right-hand
// side through its LinearFunction interface, on the V2 coarse operator.

NAMESPACE_BEGIN(Grid);

template<class Field, class CoarseField, class Aggregation>
class TwoLevelADEF2 : public TwoLevelCG<Field>
{
 public:
  ///////////////////////////////////////////////////////////////////////////////////
  // Need something that knows how to get from Coarse to fine and back again
  //  void ProjectToSubspace(CoarseVector &CoarseVec,const FineField &FineVec){
  //  void PromoteFromSubspace(const CoarseVector &CoarseVec,FineField &FineVec){
  ///////////////////////////////////////////////////////////////////////////////////
  GridBase *coarsegrid;
  Aggregation &_Aggregates;                    
  LinearFunction<CoarseField> &_CoarseSolver;
  LinearFunction<CoarseField> &_CoarseSolverPrecise;
  ///////////////////////////////////////////////////////////////////////////////////
  
  // more most opertor functions
  TwoLevelADEF2(RealD tol,
		Integer maxit,
		LinearOperatorBase<Field>    &FineLinop,
		LinearFunction<Field>        &Smoother,
		LinearFunction<CoarseField>  &CoarseSolver,
		LinearFunction<CoarseField>  &CoarseSolverPrecise,
		Aggregation &Aggregates
		) :
      TwoLevelCG<Field>(tol,maxit,FineLinop,Smoother,Aggregates.FineGrid),
      _CoarseSolver(CoarseSolver),
      _CoarseSolverPrecise(CoarseSolverPrecise),
      _Aggregates(Aggregates)
  {
    coarsegrid = Aggregates.CoarseGrid;
  };

  virtual void PcgM1(Field & in, Field & out)
  {
    GRID_TRACE("MultiGridPreconditioner ");
    // [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]

    Field tmp(this->grid);
    Field Min(this->grid);
    CoarseField PleftProj(this->coarsegrid);
    CoarseField PleftMss_proj(this->coarsegrid);

    this->SmoothTimer.Start();
    this->_Smoother(in,Min);
    this->SmoothTimer.Stop();
    this->SmoothCalls++;

    this->MatrixTimer.Start();
    this->_FineLinop.HermOp(Min,out);
    this->MatrixTimer.Stop();
    this->MatrixCalls++;
    axpy(tmp,-1.0,out,in);          // tmp  = in - A Min

    this->ProjectTimer.Start();
    this->_Aggregates.ProjectToSubspace(PleftProj,tmp);
    this->ProjectTimer.Stop();
    this->ProjectCalls++;
    this->CoarseTimer.Start();
    this->_CoarseSolver(PleftProj,PleftMss_proj); // Ass^{-1} [in - A Min]_s
    this->CoarseTimer.Stop();
    this->CoarseCalls++;
    this->PromoteTimer.Start();
    this->_Aggregates.PromoteFromSubspace(PleftMss_proj,tmp);// tmp = Q[in - A Min]
    this->PromoteTimer.Stop();
    this->PromoteCalls++;

    axpy(out,1.0,Min,tmp); // Min+tmp
  }

  virtual void Vstart(Field & x,const Field & src)
  {
    std::cout << GridLogMessage<<"HDCG: fPcg Vstart "<<std::endl;
    ///////////////////////////////////
    // Choose x_0 such that 
    // x_0 = guess +  (A_ss^inv) r_s = guess + Ass_inv [src -Aguess]
    //                               = [1 - Ass_inv A] Guess + Assinv src
    //                               = P^T guess + Assinv src 
    //                               = Vstart  [Tang notation]
    // This gives:
    // W^T (src - A x_0) = src_s - A guess_s - r_s
    //                   = src_s - (A guess)_s - src_s  + (A guess)_s 
    //                   = 0 
    ///////////////////////////////////
    Field r(this->grid);
    Field mmp(this->grid);
    CoarseField PleftProj(this->coarsegrid);
    CoarseField PleftMss_proj(this->coarsegrid);

    std::cout << GridLogMessage<<"HDCG: fPcg Vstart projecting "<<std::endl;
    this->_Aggregates.ProjectToSubspace(PleftProj,src);     
    std::cout << GridLogMessage<<"HDCG: fPcg Vstart coarse solve "<<std::endl;
    this->_CoarseSolverPrecise(PleftProj,PleftMss_proj); // Ass^{-1} r_s
    std::cout << GridLogMessage<<"HDCG: fPcg Vstart promote "<<std::endl;
    this->_Aggregates.PromoteFromSubspace(PleftMss_proj,x);  

  }

};

NAMESPACE_END(Grid);
