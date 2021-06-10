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

#include <Grid/qcd/utils/MixedPrecisionOperatorFunction.h>
#include <Grid/qcd/action/domains/Domains.h>

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


template<class ImplD,class ImplF>
class SchurFactoredFermionOperator : public ImplD
{
  INHERIT_IMPL_TYPES(ImplD);
  
  typedef typename ImplF::FermionField FermionFieldF;
  typedef typename ImplD::FermionField FermionFieldD;

  typedef SchurDiagMooeeOperator<FermionOperator<ImplD>,FermionFieldD> LinearOperatorD;
  typedef SchurDiagMooeeOperator<FermionOperator<ImplF>,FermionFieldF> LinearOperatorF;
  typedef SchurDiagMooeeDagOperator<FermionOperator<ImplD>,FermionFieldD> LinearOperatorDagD;
  typedef SchurDiagMooeeDagOperator<FermionOperator<ImplF>,FermionFieldF> LinearOperatorDagF;

  typedef MixedPrecisionConjugateGradientOperatorFunction<FermionOperator<ImplD>,
							  FermionOperator<ImplF>,
							  LinearOperatorD,
							  LinearOperatorF> MxPCG;

  typedef MixedPrecisionConjugateGradientOperatorFunction<FermionOperator<ImplD>,
							  FermionOperator<ImplF>,
							  LinearOperatorDagD,
							  LinearOperatorDagF> MxDagPCG;
public:

  GridBase *FermionGrid(void) { return PeriodicFermOpD.FermionGrid(); };
  GridBase *GaugeGrid(void)   { return PeriodicFermOpD.GaugeGrid(); };
  
  FermionOperator<ImplD> & DirichletFermOpD;
  FermionOperator<ImplF> & DirichletFermOpF;
  FermionOperator<ImplD> & PeriodicFermOpD; 
  FermionOperator<ImplF> & PeriodicFermOpF; 

  LinearOperatorD DirichletLinOpD;
  LinearOperatorF DirichletLinOpF;
  LinearOperatorD PeriodicLinOpD;
  LinearOperatorF PeriodicLinOpF;

  LinearOperatorDagD DirichletLinOpDagD;
  LinearOperatorDagF DirichletLinOpDagF;
  LinearOperatorDagD PeriodicLinOpDagD;
  LinearOperatorDagF PeriodicLinOpDagF;

  // Can tinker with these in the pseudofermion for force vs. action solves
  Integer maxinnerit;
  Integer maxouterit;
  RealD tol;
  RealD tolinner;
  
  Coordinate Block;

  DomainDecomposition Domains;

  SchurFactoredFermionOperator(FermionOperator<ImplD>  & _PeriodicFermOpD,
			       FermionOperator<ImplF>  & _PeriodicFermOpF,
			       FermionOperator<ImplD>  & _DirichletFermOpD,
			       FermionOperator<ImplF>  & _DirichletFermOpF,
			       Coordinate &_Block)
   : Domains(Block),

      PeriodicFermOpD(_PeriodicFermOpD),
      PeriodicFermOpF(_PeriodicFermOpF),
      DirichletFermOpD(_DirichletFermOpD),
      DirichletFermOpF(_DirichletFermOpF),
      DirichletLinOpD(DirichletFermOpD),
      DirichletLinOpF(DirichletFermOpF),
      PeriodicLinOpD(PeriodicFermOpD),
      PeriodicLinOpF(PeriodicFermOpF),
      DirichletLinOpDagD(DirichletFermOpD),
      DirichletLinOpDagF(DirichletFermOpF),
      PeriodicLinOpDagD(PeriodicFermOpD),
      PeriodicLinOpDagF(PeriodicFermOpF)
  {
    tol=1.0e-10;
    tolinner=1.0e-6;
    maxinnerit=1000;
    maxouterit=10;
    assert(PeriodicFermOpD.FermionGrid() == DirichletFermOpD.FermionGrid());
    assert(PeriodicFermOpF.FermionGrid() == DirichletFermOpF.FermionGrid());
  };

  enum Domain { Omega=0, OmegaBar=1 };

  void ImportGauge(const GaugeField &Umu)
  {
    // Single precision will update in the mixed prec CG
    PeriodicFermOpD.ImportGauge(Umu);
    GaugeField dUmu(Umu.Grid());
    dUmu=Umu;
    //    DirchletBCs(dUmu);
    DirichletFilter<GaugeField> Filter(Block);
    Filter.applyFilter(dUmu);
    DirichletFermOpD.ImportGauge(dUmu);
  }

/*
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
*/
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
    LatticeInteger  one(grid);   one = 1;
    LatticeInteger  zero(grid); zero = 0;
    LatticeInteger  omega(grid);
    LatticeInteger  omegabar(grid);
    LatticeInteger  tmp(grid);

    omega=one;    Domains.ProjectDomain(omega,0);
    omegabar=one; Domains.ProjectDomain(omegabar,1);
    
    LatticeInteger nface(grid); nface=Zero();
    
    FermionField projected(grid); projected=Zero();
    FermionField sp_proj  (grid);

    int dims = grid->Nd();
    int isDWF= (dims==Nd+1);
    assert((dims==Nd)||(dims==Nd+1));
    Coordinate Global=grid->GlobalDimensions();

    for(int mmu=0;mmu<Nd;mmu++){
      Gamma G(Gmu[mmu]);

      // need to worry about DWF 5th dim first
      int mu = mmu+isDWF;
      if ( Block[mmu] && (Block[mmu] <= Global[mu]) ) {

	// Lower face receives (1-gamma)/2 in normal forward hopping term
 	tmp = Cshift(omegabar,mu,-1);
	tmp = tmp + omega;
	face = where(tmp == Integer(2),one,zero );

 	tmp = Cshift(omega,mu,-1);
	tmp = tmp + omegabar;
	face = where(tmp == Integer(2),one,face );

	nface = nface + face;

	sp_proj  = 0.5*(f-G*f*rsgn);
	projected= where(face,sp_proj,projected);

	// Upper face receives (1+gamma)/2 in normal backward hopping term
 	tmp = Cshift(omegabar,mu,1);
	tmp = tmp + omega;
	face = where(tmp == Integer(2),one,zero );

 	tmp = Cshift(omega,mu,1);
	tmp = tmp + omegabar;
	face = where(tmp == Integer(2),one,face );

	nface = nface + face;

	sp_proj = 0.5*(f+G*f*rsgn);
	projected= where(face,sp_proj,projected);
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
/*
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
*/
    Domains.ProjectDomain(f,domain);

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
    PeriodicFermOpD.M(tmp,out);
    ProjectOmega(out);
  };
  void dBoundaryDag (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmega(tmp);
    PeriodicFermOpD.Mdag(tmp,out);
    ProjectOmegaBar(out);
  };
  void dBoundaryBar (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmega(tmp);
    PeriodicFermOpD.M(tmp,out);
    ProjectOmegaBar(out);
  };
  void dBoundaryBarDag (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmegaBar(tmp);
    PeriodicFermOpD.Mdag(tmp,out);
    ProjectOmega(out);
  };
  void dOmega       (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmega(tmp);
    DirichletFermOpD.M(tmp,out);
    ProjectOmega(out);
  };
  void dOmegaBar    (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmegaBar(tmp);
    DirichletFermOpD.M(tmp,out);
    ProjectOmegaBar(out);
  };
  void dOmegaDag       (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmega(tmp);
    DirichletFermOpD.Mdag(tmp,out);
    ProjectOmega(out);
  };
  void dOmegaBarDag    (FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmegaBar(tmp);
    DirichletFermOpD.Mdag(tmp,out);
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
    MxPCG OmegaSolver(tol,
		      tolinner,
		      maxinnerit,
		      maxouterit,
		      DirichletFermOpF.FermionRedBlackGrid(),
		      DirichletFermOpF,
		      DirichletFermOpD,
		      DirichletLinOpF,
		      DirichletLinOpD);
    SchurRedBlackDiagMooeeSolve<FermionField> PrecSolve(OmegaSolver);
    PrecSolve(DirichletFermOpD,in,out);
  };
  void dOmegaDagInvAndOmegaBarDagInv(FermionField &in,FermionField &out)
  {
    MxDagPCG OmegaDagSolver(tol,
			    tolinner,
			    maxinnerit,
			    maxouterit,
			    DirichletFermOpF.FermionRedBlackGrid(),
			    DirichletFermOpF,
			    DirichletFermOpD,
			    DirichletLinOpDagF,
			    DirichletLinOpDagD);
    SchurRedBlackDiagMooeeDagSolve<FermionField> PrecSolve(OmegaDagSolver);
    PrecSolve(DirichletFermOpD,in,out);
  };

  // Rdag = Pdbar - DdbarDag DomegabarDagInv  DdDag DomegaDagInv Pdbar 
  void RDag(FermionField &in,FermionField &out)
  {
    FermionField tmp1(PeriodicFermOpD.FermionGrid());
    FermionField tmp2(PeriodicFermOpD.FermionGrid());
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
    FermionField tmp1(PeriodicFermOpD.FermionGrid());
    FermionField tmp2(PeriodicFermOpD.FermionGrid());
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
    FermionField tmp1(PeriodicFermOpD.FermionGrid());
    dBoundaryBar(in,out);
    Dinverse(out,tmp1);  
    out =in -tmp1; 
    ProjectBoundaryBar(out);
  };
  // R = Pdbar - DdbarDag DinvDag Pdbar 
  void RDagInv(FermionField &in,FermionField &out)
  {
    FermionField tmp(PeriodicFermOpD.FermionGrid());
    FermionField Pin(PeriodicFermOpD.FermionGrid());
    Pin = in; ProjectBoundaryBar(Pin);
    DinverseDag(Pin,out);  
    dBoundaryBarDag(out,tmp);
    out =Pin -tmp; 
  };
  // Non-dirichlet inverter using red-black preconditioning
  void Dinverse(FermionField &in,FermionField &out)
  {
    MxPCG DSolver(tol,
		  tolinner,
		  maxinnerit,
		  maxouterit,
		  PeriodicFermOpF.FermionRedBlackGrid(),
		  PeriodicFermOpF,
		  PeriodicFermOpD,
		  PeriodicLinOpF,
		  PeriodicLinOpD);
    SchurRedBlackDiagMooeeSolve<FermionField> Solve(DSolver);
    Solve(PeriodicFermOpD,in,out);
  }
  void DinverseDag(FermionField &in,FermionField &out)
  {
    MxDagPCG DdagSolver(tol,
			tolinner,
			maxinnerit,
			maxouterit,
			PeriodicFermOpF.FermionRedBlackGrid(),
			PeriodicFermOpF,
			PeriodicFermOpD,
			PeriodicLinOpDagF,
			PeriodicLinOpDagD);
    SchurRedBlackDiagMooeeDagSolve<FermionField> Solve(DdagSolver);
    Solve(PeriodicFermOpD,in,out);
  }
};

NAMESPACE_END(Grid);

