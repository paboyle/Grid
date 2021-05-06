/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/pseudofermion/TwoFlavourRatio.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>

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

///////////////////////////////////////
// Two flavour ratio
///////////////////////////////////////
template<class Impl>
class DomainDecomposedBoundary {
public:
  INHERIT_IMPL_TYPES(Impl);

  typedef typename GaugeField::vector_type vector_type; //SIMD-vectorized complex type
  typedef typename GaugeField::scalar_type scalar_type; //scalar complex type

  typedef iVector<iScalar<iScalar<vector_type> >, Nd > LorentzScalarType; //complex phase for each site/direction
  typedef iScalar<iScalar<iScalar<vector_type> > >            ScalarType; //complex phase for each site
  typedef Lattice<LorentzScalarType> LatticeLorentzScalarType;
  typedef Lattice<ScalarType> LatticeScalarType;

  Coordinate Block;
  DDHMCFilter Filter;
  const int Omega=0;
  const int OmegaBar=1;

  void ProjectBoundaryBothDomains (FermionField &f,int sgn)
  {
    assert((sgn==1)||(sgn==-1));

    Gamma::Algebra Gmu [] = {
			     Gamma::Algebra::GammaX,
			     Gamma::Algebra::GammaY,
			     Gamma::Algebra::GammaZ,
			     Gamma::Algebra::GammaT
    };

    GridBase *grid = f.Grid();
    LatticeInteger  coor(grid);
    LatticeInteger  face(grid);
    LatticeInteger nface(grid); nface=Zero();
    
    ComplexField zz(grid); zz=Zero();

    FermionField projected(grid); projected=Zero();
    FermionField sp_proj  (grid);

    int dims = grid->Nd();
    int isDWF= (dims==Nd+1);
    assert((dims==Nd)||(dims==Nd+1));

    for(int mu=0;mu<Nd;mu++){ 
      // need to worry about DWF 5th dim first
      // Could extend to domain decompose in FIFTH dimension.
      // With chiral projectors here
      LatticeCoordinate(coor,mu+isDWF); 
      
      face = (mod(coor,Block[mu]) == 0 );
      nface = nface + face;
      
      // Lower face receives (1-gamma)/2 in normal forward hopping term
      sp_proj  = 0.5*(f-sgn*Gamma(Gmu[mu])*f)
      projected= where(face==cb,f,projected);
       
      face = (mod(coor,Block[mu]) == Block[mu]-1 );
      nface = nface + face;

      // Upper face receives (1+gamma)/2 in normal backward hopping term
      sp_proj = 0.5*(f+sgn*Gamma(Gmu[mu])*f)
      projected= where(face==cb,f,projected);

    }
    // Keep the spin projected faces where nface==1 and initial Zero() where nface==0.
    projected = where(nface>1,f,projected);
  }
  void ProjectDomain(FermionField &f,int cb)
  {    
    GridBase *grid = f.Grid();
    ComplexField zz(grid); zz=Zero();
    LatticeInteger coor(grid);
    LatticeInteger domaincb(grid); domaincb=Zero();
    for(int d=0;d<grid->Nd();d++){
      LatticeCoordinate(coor,mu);
      domaincb = domaincb + div(coor,Block[d]);
    }
    f = where(mod(domaincb,2)==cb,f,zz);
  };

  void ProjectOmegaBar   (FermionField &f) {ProjectDomain(f,OmegaBar);}
  void ProjectOmega      (FermionField &f) {ProjectDomain(f,Omega);}

  // See my notes(!).
  // Notation: Following Luscher, we introduce projectors $\hPdb$ with both spinor and space structure
  // projecting all spinor elements in $\Omega$ connected by $\Ddb$ to $\bar{\Omega}$,
  void ProjectBoundaryBar(FermionField &f)
  {
    ProjectBoundaryBothDomains(f);
    ProjectOmega(f);
  }
  // and $\hPd$ projecting all spinor elements in $\bar{\Omega}$ connected by $\Dd$ to $\Omega$.
  void ProjectBoundary   (FermionField &f)
  {
    ProjectBoundaryBothDomains(f);
    ProjectOmegaBar(f);
  };

  void dBoundary    (FermionOperator<Impl>  &Op,FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmegaBar(tmp);
    Op.M(tmp,out);
    ProjectOmega(out);
  };
  void dBoundaryBar (FermionOperator<Impl>  &Op,FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmega(tmp);
    Op.M(tmp,out);
    ProjectOmegaBar(out);
  };
  void dOmega       (FermionOperator<Impl>  &Op,FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmega(tmp);
    Op.M(tmp,out);
    ProjectOmega(out);
  };
  void dOmegaBar    (FermionOperator<Impl>  &Op,FermionField &in,FermionField &out)
  {
    FermionField tmp(in);
    ProjectOmegaBar(tmp);
    Op.M(tmp,out);
    ProjectOmegaBar(out);
  };

  void SolveOmega   (FermionOperator<Impl>  &Op,FermionField &in,FermionField &out){    assert(0);  };
  void SolveOmegaBar(FermionOperator<Impl>  &Op,FermionField &in,FermionField &out){    assert(0);  };
  void SolveOmegaAndOmegaBar(FermionOperator<Impl>  &Op,FermionField &in,FermionField &out){    assert(0);  };
  void dInverse     (FermionOperator<Impl>  &Op,FermionField &in,FermionField &out){    assert(0);  };

  // R = Pdbar - Pdbar DomegaInv Dd DomegabarInv Ddbar
  void R(FermionOperator<Impl>  &Op,FermionOperator<Impl>  &OpDirichlet,FermionField &in,FermionField &out)
  {
    FermionField tmp1(Op.FermionGrid());
    FermionField tmp2(Op.FermionGrid());
    dBoundaryBar(Op,in,tmp1);
    SolveOmegaBar(OpDirichlet,tmp1,tmp2); // 1/2 cost
    dBoundary(Op,tmp2,tmp1);
    SolveOmega(OpDirichlet,tmp1,tmp2); // 1/2 cost
    out = in - tmp2 ;
    ProjectBoundaryBar(out);
  };
  
  // R = Pdbar - Pdbar Dinv Ddbar 
  void Rinverse(FermionField &in,FermionField &out)
  {
    FermionField tmp1(NumOp.FermionGrid());
    out = in;
    ProjectBoundaryBar(out);
    dInverse(out,tmp1);
    ProjectBoundaryBar(tmp1);
    out = out -tmp1;
  };
  
}
  
template<class Impl>
class DomainDecomposedBoundaryPseudoFermionAction : public Action<typename Impl::GaugeField> {
public:
  INHERIT_IMPL_TYPES(Impl);

private:
  FermionOperator<Impl> & NumOp;// the basic operator
  FermionOperator<Impl> & DenOp;// the basic operator
  FermionOperator<Impl> & NumOpDirichlet;// the basic operator
  FermionOperator<Impl> & DenOpDirichlet;// the basic operator

  OperatorFunction<FermionField> &DerivativeSolver;
  OperatorFunction<FermionField> &ActionSolver;

  FermionField Phi; // the pseudo fermion field for this trajectory
  

public:
  DomainBoundaryPseudoFermionAction(FermionOperator<Impl>  &_NumOp, 
				    FermionOperator<Impl>  &_DenOp,
				    FermionOperator<Impl>  &_NumOpDirichlet, 
				    FermionOperator<Impl>  &_DenOpDirichlet, 
				    OperatorFunction<FermionField> & DS,
				    OperatorFunction<FermionField> & AS,
				    Coordinate &_Block
				    ) : NumOp(_NumOp),
					DenOp(_DenOp),
					DerivativeSolver(DS),
					ActionSolver(AS),
					Phi(_NumOp.FermionGrid()),
					Block(_Block)
					//					LinkFilter(Block)
  {};
      
  virtual std::string action_name(){return "DomainBoundaryPseudoFermionRatioAction";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "["<<action_name()<<"] Block "<<_Block << std::endl;
    return sstream.str();
  }  

  
  virtual void refresh(const GaugeField &U, GridParallelRNG& pRNG)
  {
    // P(phi) = e^{- phi^dag V (MdagM)^-1 Vdag phi}
    //
    // NumOp == V
    // DenOp == M
    //
    // Take phi = Vdag^{-1} Mdag eta  ; eta = Mdag^{-1} Vdag Phi
    //
    // P(eta) = e^{- eta^dag eta}
    //
    // e^{x^2/2 sig^2} => sig^2 = 0.5.
    // 
    // So eta should be of width sig = 1/sqrt(2) and must multiply by 0.707....
    //
    RealD scale = std::sqrt(0.5);

    FermionField eta(NumOp.FermionGrid());
    FermionField tmp(NumOp.FermionGrid());

    gaussian(pRNG,eta);

    ProjectBoundary(eta);
    
    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    // Note: this hard codes normal equations type solvers; alternate implementation needed for 
    // non-herm style solvers.
    MdagMLinearOperator<FermionOperator<Impl> ,FermionField> MdagMOp(NumOp);

    DenOp.Mdag(eta,Phi);            // Mdag eta
    tmp = Zero();
    ActionSolver(MdagMOp,Phi,tmp);  // (VdagV)^-1 Mdag eta = V^-1 Vdag^-1 Mdag eta
    NumOp.M(tmp,Phi);               // Vdag^-1 Mdag eta

    Phi=Phi*scale;
	
  };

  //////////////////////////////////////////////////////
  // S = phi^dag V (Mdag M)^-1 Vdag phi
  //////////////////////////////////////////////////////
  virtual RealD S(const GaugeField &U) {

    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    FermionField X(NumOp.FermionGrid());
    FermionField Y(NumOp.FermionGrid());
	
    MdagMLinearOperator<FermionOperator<Impl> ,FermionField> MdagMOp(DenOp);

    NumOp.Mdag(Phi,Y);              // Y= Vdag phi
    X=Zero();
    ActionSolver(MdagMOp,Y,X);      // X= (MdagM)^-1 Vdag phi
    DenOp.M(X,Y);                  // Y=  Mdag^-1 Vdag phi

    RealD action = norm2(Y);

    return action;
  };

  //////////////////////////////////////////////////////
  // dS/du = phi^dag dV (Mdag M)^-1 V^dag  phi
  //       - phi^dag V (Mdag M)^-1 [ Mdag dM + dMdag M ]  (Mdag M)^-1 V^dag  phi
  //       + phi^dag V (Mdag M)^-1 dV^dag  phi
  //////////////////////////////////////////////////////
  virtual void deriv(const GaugeField &U,GaugeField & dSdU) {

    NumOp.ImportGauge(U);
    DenOp.ImportGauge(U);

    MdagMLinearOperator<FermionOperator<Impl> ,FermionField> MdagMOp(DenOp);

    FermionField  X(NumOp.FermionGrid());
    FermionField  Y(NumOp.FermionGrid());

    GaugeField   force(NumOp.GaugeGrid());	


    //Y=Vdag phi
    //X = (Mdag M)^-1 V^dag phi
    //Y = (Mdag)^-1 V^dag  phi
    NumOp.Mdag(Phi,Y);              // Y= Vdag phi
    X=Zero();
    DerivativeSolver(MdagMOp,Y,X);      // X= (MdagM)^-1 Vdag phi
    DenOp.M(X,Y);                  // Y=  Mdag^-1 Vdag phi

    // phi^dag V (Mdag M)^-1 dV^dag  phi
    NumOp.MDeriv(force , X, Phi, DaggerYes );  dSdU=force;
  
    // phi^dag dV (Mdag M)^-1 V^dag  phi
    NumOp.MDeriv(force , Phi, X ,DaggerNo  );  dSdU=dSdU+force;

    //    -    phi^dag V (Mdag M)^-1 Mdag dM   (Mdag M)^-1 V^dag  phi
    //    -    phi^dag V (Mdag M)^-1 dMdag M   (Mdag M)^-1 V^dag  phi
    DenOp.MDeriv(force,Y,X,DaggerNo);   dSdU=dSdU-force;
    DenOp.MDeriv(force,X,Y,DaggerYes);  dSdU=dSdU-force;

    dSdU *= -1.0;
    //dSdU = - Ta(dSdU);

  };
};

NAMESPACE_END(Grid);


