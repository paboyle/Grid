/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/SchurFactoredFermionOperator.h

    Copyright (C) 2021

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

NAMESPACE_BEGIN(Grid);

  ////////////////////////////////////////////////////////
  // Some explanation of class structure for domain decomposition:
  //
  // Need a dirichlet operator for two flavour determinant - acts on both Omega and OmegaBar.
  //
  // Possible gain if the global sums and CG are run independently?? Could measure this.
  //
  // Types of operations
  //
  // 1) assemble local det dOmega det dOmegaBar pseudofermion
  //
  // - DirichletFermionOperator - can either do a global solve, or independent/per cell coefficients.
  //
  // 2) assemble dOmegaInverse and dOmegaBarInverse in R
  //
  // - DirichletFermionOperator - can also be used to 
  //                                       - need two or more cells per node. Options
  //                                       - a) solve one cell at a time, no new code, CopyRegion and reduced /split Grids
  //                                       - b) solve multiple cells in parallel. predicated dslash implementation
  //
  //                                       - b) has more parallelism, experience with block solver suggest might not be aalgorithmically inefficient
  //                                         a) has more cache friendly and easier code.
  //                                         b) is easy to implement in a "trial" or inefficient code with projection.
  //
  // 3)  Additional functionality for domain operations
  //
  // - SchurFactoredFermionOperator  - Need a DDHMC utility - whether used in two flavour or one flavour 
  //
  // - dBoundary - needs non-dirichlet operator
  // - Contains one Dirichlet Op, and one non-Dirichlet op. Implements dBoundary etc...
  // - The Dirichlet ops can be passed to dOmega(Bar) solvers etc...
  //
  ////////////////////////////////////////////////////////

template<class Impl>
class SchurFactoredFermionOperator : public Impl
{
  INHERIT_IMPL_TYPES(Impl);

public:
  
  FermionOperator<Impl> & DirichletFermOp;
  FermionOperator<Impl> & FermOp; 
  OperatorFunction<FermionField> &OmegaSolver;
  OperatorFunction<FermionField> &OmegaDagSolver;
  OperatorFunction<FermionField> &DSolver;
  OperatorFunction<FermionField> &DdagSolver;
  Coordinate Block;

  SchurFactoredFermionOperator(FermionOperator<Impl> & _FermOp,
			       FermionOperator<Impl> & _DirichletFermOp,
			       OperatorFunction<FermionField> &_OmegaSolver,
			       OperatorFunction<FermionField> &_OmegaDagSolver,
			       OperatorFunction<FermionField> &_DSolver,
			       OperatorFunction<FermionField> &_DdagSolver,
			       Coordinate &_Block)
    : Block(_Block),
      FermOp(_FermOp),
      DirichletFermOp(_DirichletFermOp),
      OmegaSolver(_OmegaSolver),
      OmegaDagSolver(_OmegaDagSolver),
      DSolver(_DSolver),
      DdagSolver(_DdagSolver)
  {
    // Pass in Dirichlet FermOp because we really need two dirac operators
    // as double stored gauge fields differ and they will otherwise overwrite
    assert(_FermOp.FermionGrid() == _DirichletFermOp.FermionGrid()); // May not be true in future if change communicator scheme
  };

  enum Domain { Omega=0, OmegaBar=1 };

  void ImportGauge(const GaugeField &Umu)
  {
    FermOp.ImportGauge(Umu);
    DirichletFermOp.ImportGauge(Umu);
  }
  void ProjectBoundaryBothDomains (FermionField &f,int sgn)
  {
    assert((sgn==1)||(sgn==-1));
    Real rsgn = sgn;

    Gamma::Algebra Gmu [] = {
      Gamma::Algebra::GammaX,
      Gamma::Algebra::GammaY,
      Gamma::Algebra::GammaZ,
      Gamma::Algebra::GammaT
    };

    GridBase *grid = f.Grid();
    LatticeInteger  coor(grid);
    LatticeInteger  face(grid);
    LatticeInteger  one(grid); one = 1;
    LatticeInteger  zero(grid); zero = 0;
    LatticeInteger nface(grid); nface=Zero();
    
    ComplexField zz(grid); zz=Zero();
    FermionField projected(grid); projected=Zero();
    FermionField sp_proj  (grid);

    int dims = grid->Nd();
    int isDWF= (dims==Nd+1);
    assert((dims==Nd)||(dims==Nd+1));
    Coordinate Global=grid->GlobalDimensions();

    for(int mu=0;mu<Nd;mu++){

      if ( Block[mu] <= Global[mu+isDWF] ) {
	// need to worry about DWF 5th dim first
	LatticeCoordinate(coor,mu+isDWF); 
      
	face = where(mod(coor,Block[mu]) == Integer(0),one,zero );
	nface = nface + face;

	Gamma G(Gmu[mu]);
	// Lower face receives (1-gamma)/2 in normal forward hopping term
	sp_proj  = 0.5*(f-G*f*rsgn);
	projected= where(face,sp_proj,projected);
	//projected= where(face,f,projected);
      
	face = where(mod(coor,Block[mu]) == Integer(Block[mu]-1) ,one,zero );
	nface = nface + face;

	// Upper face receives (1+gamma)/2 in normal backward hopping term
	sp_proj = 0.5*(f+G*f*rsgn);
	projected= where(face,sp_proj,projected);
	//projected= where(face,f,projected);
      }
      
    }
    // Initial Zero() where nface==0.
    // Keep the spin projected faces where nface==1
    // Full spinor where nface>=2
    projected = where(nface>Integer(1),f,projected);
    f=projected;
  }
  void ProjectDomain(FermionField &f,int domain)
  {    
    GridBase *grid = f.Grid();
    int dims = grid->Nd();
    int isDWF= (dims==Nd+1);
    assert((dims==Nd)||(dims==Nd+1));

    FermionField zz(grid); zz=Zero();
    LatticeInteger coor(grid);
    LatticeInteger domaincb(grid); domaincb=Zero();
    for(int d=0;d<Nd;d++){
      LatticeCoordinate(coor,d+isDWF);
      domaincb = domaincb + div(coor,Block[d]);
    }
    f = where(mod(domaincb,2)==Integer(domain),f,zz);
  };
  void ProjectOmegaBar   (FermionField &f) {ProjectDomain(f,OmegaBar);}
  void ProjectOmega      (FermionField &f) {ProjectDomain(f,Omega);}
  // See my notes(!).
  // Notation: Following Luscher, we introduce projectors $\hPdb$ with both spinor and space structure
  // projecting all spinor elements in $\Omega$ connected by $\Ddb$ to $\bar{\Omega}$,
  void ProjectBoundaryBar(FermionField &f)
  {
    ProjectBoundaryBothDomains(f,1);
    ProjectOmega(f);
  }
  // and $\hPd$ projecting all spinor elements in $\bar{\Omega}$ connected by $\Dd$ to $\Omega$.
  void ProjectBoundary   (FermionField &f)
  {
    ProjectBoundaryBothDomains(f,1);
    ProjectOmegaBar(f);
    //    DumpSliceNorm("ProjectBoundary",f,f.Grid()->Nd()-1);
  };

  void dBoundary    (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmegaBar(tmp);
    FermOp.M(tmp,out);
    ProjectOmega(out);
  };
  void dBoundaryDag (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmega(tmp);
    FermOp.Mdag(tmp,out);
    ProjectOmegaBar(out);
  };
  void dBoundaryBar (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmega(tmp);
    FermOp.M(tmp,out);
    ProjectOmegaBar(out);
  };
  void dBoundaryBarDag (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmegaBar(tmp);
    FermOp.Mdag(tmp,out);
    ProjectOmega(out);
  };
  void dOmega       (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmega(tmp);
    DirichletFermOp.M(tmp,out);
    ProjectOmega(out);
  };
  void dOmegaBar    (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmegaBar(tmp);
    DirichletFermOp.M(tmp,out);
    ProjectOmegaBar(out);
  };
  void dOmegaDag       (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmega(tmp);
    DirichletFermOp.Mdag(tmp,out);
    ProjectOmega(out);
  };
  void dOmegaBarDag    (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmegaBar(tmp);
    DirichletFermOp.Mdag(tmp,out);
    ProjectOmegaBar(out);
  };
  void dOmegaInv   (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmega(tmp);
    dOmegaInvAndOmegaBarInv(tmp,out); // Inefficient warning
    ProjectOmega(out);
  };
  void dOmegaBarInv(FermionField &in,FermionField &out)
  {    
    FermionField tmp(in);
    ProjectOmegaBar(tmp);
    dOmegaInvAndOmegaBarInv(tmp,out);
    ProjectOmegaBar(out);
  };
  void dOmegaDagInv   (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmega(tmp);
    dOmegaDagInvAndOmegaBarDagInv(tmp,out);
    ProjectOmega(out);
  };
  void dOmegaBarDagInv(FermionField &in,FermionField &out)
  {    
    FermionField tmp(in);
    ProjectOmegaBar(tmp);
    dOmegaDagInvAndOmegaBarDagInv(tmp,out);
    ProjectOmegaBar(out);
  };
  void dOmegaInvAndOmegaBarInv(FermionField &in,FermionField &out)
  {
    SchurRedBlackDiagMooeeSolve<FermionField> PrecSolve(OmegaSolver);
    PrecSolve(DirichletFermOp,in,out);
  };
  void dOmegaDagInvAndOmegaBarDagInv(FermionField &in,FermionField &out)
  {
    SchurRedBlackDiagMooeeDagSolve<FermionField> PrecSolve(OmegaDagSolver);
    PrecSolve(DirichletFermOp,in,out);
  };

  // Rdag = Pdbar - DdbarDag DomegabarDagInv  DdDag DomegaDagInv Pdbar 
  void RDag(FermionField &in,FermionField &out)
  {
    FermionField tmp1(FermOp.FermionGrid());
    FermionField tmp2(FermOp.FermionGrid());
    out = in;
    ProjectBoundaryBar(out);
    dOmegaDagInv(out,tmp1);   
    dBoundaryDag(tmp1,tmp2);   
    dOmegaBarDagInv(tmp2,tmp1);
    dBoundaryBarDag(tmp1,tmp2); 
    out = out - tmp2;
  };

  // R = Pdbar - Pdbar DomegaInv Dd DomegabarInv Ddbar
  void R(FermionField &in,FermionField &out)
  {
    FermionField tmp1(FermOp.FermionGrid());
    FermionField tmp2(FermOp.FermionGrid());
    out = in;
    ProjectBoundaryBar(out);
    dBoundaryBar(out,tmp1); 
    dOmegaBarInv(tmp1,tmp2);
    dBoundary(tmp2,tmp1);   
    dOmegaInv(tmp1,tmp2);   
    out = in - tmp2 ;       
    ProjectBoundaryBar(out);
    //    DumpSliceNorm("R",out,out.Grid()->Nd()-1);
  };
  
  // R = Pdbar - Pdbar Dinv Ddbar 
  void RInv(FermionField &in,FermionField &out)
  {
    FermionField tmp1(FermOp.FermionGrid());
    dBoundaryBar(in,out);
    Dinverse(out,tmp1);  
    out =in -tmp1; 
    ProjectBoundaryBar(out);
  };
  // R = Pdbar - DdbarDag DinvDag Pdbar 
  void RDagInv(FermionField &in,FermionField &out)
  {
    FermionField tmp(FermOp.FermionGrid());
    FermionField Pin(FermOp.FermionGrid());
    Pin = in; ProjectBoundaryBar(Pin);
    DinverseDag(Pin,out);  
    dBoundaryBarDag(out,tmp);
    out =Pin -tmp; 
  };
  // Non-dirichlet inverter using red-black preconditioning
  void Dinverse(FermionField &in,FermionField &out)
  {
    SchurRedBlackDiagMooeeSolve<FermionField> Solve(DSolver);
    Solve(FermOp,in,out);
  }
  void DinverseDag(FermionField &in,FermionField &out)
  {
    SchurRedBlackDiagMooeeDagSolve<FermionField> Solve(DdagSolver);
    Solve(FermOp,in,out);
  }
};

NAMESPACE_END(Grid);

