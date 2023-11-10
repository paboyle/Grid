/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/integrators/Integrator.h

Copyright (C) 2015

Author: Guido Cossu <guido.cossu@ed.ac.uk>
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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
			   /*  END LEGAL */
			   //--------------------------------------------------------------------
#pragma once

NAMESPACE_BEGIN(Grid);

template <typename Field> 
class Metric{
protected:
  int triv;
public:
  Metric(){this->triv=1;}
  int Trivial(){ return triv;}
//printf("Metric::Trivial=%d\n",triv); ;
  virtual void ImportGauge(const Field&)   = 0;
  virtual void M(const Field&, Field&)     = 0;
  virtual void Minv(const Field&, Field&)  = 0;
  virtual void MSquareRoot(Field&) = 0;
  virtual void MInvSquareRoot(Field&) = 0;
  virtual void MDeriv(const Field&, Field&) = 0;
  virtual void MDeriv(const Field&, const Field&, Field&) = 0;
  virtual void MinvDeriv(const Field&, Field&) = 0;
//  virtual void MinvDeriv(const Field&, const Field&, Field&) = 0;
};


// Need a trivial operator
template <typename Field>
class TrivialMetric : public Metric<Field>{
public:
//  TrivialMetric(){this->triv=1;printf("TrivialMetric::triv=%d\n",this->Trivial());}
  virtual void ImportGauge(const Field&){};
  virtual void M(const Field& in, Field& out){
//    printf("M:norm=%0.15e\n",norm2(in));
    std::cout << GridLogIntegrator << " M:norm(in)= " << std::sqrt(norm2(in)) << std::endl;
    out = in;
  }
  virtual void Minv(const Field& in, Field& out){
    std::cout << GridLogIntegrator << " Minv:norm(in)= " << std::sqrt(norm2(in)) << std::endl;
    out = in;
  }
  virtual void MSquareRoot(Field& P){
    std::cout << GridLogIntegrator << " MSquareRoot:norm(P)= " << std::sqrt(norm2(P)) << std::endl;
    // do nothing
  }
  virtual void MInvSquareRoot(Field& P){
    std::cout << GridLogIntegrator << " MInvSquareRoot:norm(P)= " << std::sqrt(norm2(P)) << std::endl;
    // do nothing
  }
  virtual void MDeriv(const Field& in, Field& out){
    std::cout << GridLogIntegrator << " MDeriv:norm(in)= " << std::sqrt(norm2(in)) << std::endl;
    out = Zero();
  }
  virtual void MinvDeriv(const Field& in, Field& out){
    std::cout << GridLogIntegrator << " MinvDeriv:norm(in)= " << std::sqrt(norm2(in)) << std::endl;
    out = Zero();
  }
  virtual void MDeriv(const Field& left, const Field& right, Field& out){
    std::cout << GridLogIntegrator << " MDeriv:norm(left)= " << std::sqrt(norm2(left)) << std::endl;
    std::cout << GridLogIntegrator << " MDeriv:norm(right)= " << std::sqrt(norm2(right)) << std::endl;
    out = Zero();
  }

};

///////////////////////////////
// Generalised momenta
///////////////////////////////

template <typename Implementation>
class GeneralisedMomenta{
public:
  typedef typename Implementation::Field MomentaField;  //for readability
  typedef typename Implementation::GaugeLinkField MomentaLinkField;  //for readability
  Metric<MomentaField>& M;
  MomentaField Mom;

  // Auxiliary fields
  // not hard coded but inherit the type from the metric
  // created Nd new fields
  // hide these in the metric?
  //typedef Lattice<iVector<iScalar<iMatrix<vComplex, Nc> >, Nd/2 > > AuxiliaryMomentaType;
  MomentaField AuxMom;
  MomentaField AuxField;

  GeneralisedMomenta(GridBase* grid, Metric<MomentaField>& M): M(M), Mom(grid), AuxMom(grid), AuxField(grid){}

  // Correct
  void MomentaDistribution(GridSerialRNG & sRNG, GridParallelRNG& pRNG){
    // Generate a distribution for
    // P^dag G P
    // where G = M^-1

    // Generate gaussian momenta
    Implementation::generate_momenta(Mom, sRNG, pRNG);
    // Modify the distribution with the metric
//    if(M.Trivial()) return;
    M.MSquareRoot(Mom);

    if (1) {
      // Auxiliary momenta
      // do nothing if trivial, so hide in the metric
      MomentaField AuxMomTemp(Mom.Grid());
      Implementation::generate_momenta(AuxMom, sRNG,pRNG);
      Implementation::generate_momenta(AuxField, sRNG,pRNG);
      // Modify the distribution with the metric
      // Aux^dag M Aux
      M.MInvSquareRoot(AuxMom);  // AuxMom = M^{-1/2} AuxMomTemp
    }
  }

  // Correct
  RealD MomentaAction(){
    static RealD Saux=0.,Smom=0.;
    MomentaField inv(Mom.Grid());
    inv = Zero();
    M.Minv(Mom, inv);
    LatticeComplex Hloc(Mom.Grid()); Hloc = Zero();
    LatticeComplex Hloc2(Mom.Grid()); Hloc2 = Zero();
    for (int mu = 0; mu < Nd; mu++) {
      // This is not very general
      // hide in the metric
      auto Mom_mu = PeekIndex<LorentzIndex>(Mom, mu);
      auto inv_mu = PeekIndex<LorentzIndex>(inv, mu);
      Hloc += trace(Mom_mu * inv_mu);
    }
    auto Htmp1 = TensorRemove(sum(Hloc));
    std::cout << GridLogMessage << "S:dSmom = " << Htmp1.real()-Smom << "\n";
    Smom=Htmp1.real()/HMC_MOMENTUM_DENOMINATOR;
    

    

//    if(!M.Trivial()) 
    {
      // Auxiliary Fields
      // hide in the metric
      M.M(AuxMom, inv);
      for (int mu = 0; mu < Nd; mu++) {
        // This is not very general
        // hide in the operators
        auto inv_mu = PeekIndex<LorentzIndex>(inv, mu);
        auto am_mu = PeekIndex<LorentzIndex>(AuxMom, mu);
        auto af_mu = PeekIndex<LorentzIndex>(AuxField, mu);
        Hloc += trace(am_mu * inv_mu);
        Hloc2 += trace(af_mu * af_mu);
      }
    }
    auto Htmp2 = TensorRemove(sum(Hloc))-Htmp1;
    std::cout << GridLogMessage << "S:dSaux = " << Htmp2.real()-Saux << "\n";
    Saux=Htmp2.real();

    auto Hsum = TensorRemove(sum(Hloc))/HMC_MOMENTUM_DENOMINATOR;
    auto Hsum2 = TensorRemove(sum(Hloc2));
    std::cout << GridLogIntegrator << "MomentaAction: " <<  Hsum.real()+Hsum2.real() << std::endl;
    return Hsum.real()+Hsum2.real();
  }

  // Correct
  void DerivativeU(MomentaField& in, MomentaField& der){

    // Compute the derivative of the kinetic term
    // with respect to the gauge field
    MomentaField MDer(in.Grid());
    MomentaField X(in.Grid());
    X = Zero();
    M.MinvDeriv(in, MDer);  // MDer = U * dS/dU
    der = -1.0* Implementation::projectForce(MDer);  // Ta if gauge fields
//    std::cout << GridLogIntegrator << " DerivativeU: norm(in)= " << std::sqrt(norm2(in)) << std::endl;
//    std::cout << GridLogIntegrator << " DerivativeU: norm(der)= " << std::sqrt(norm2(der)) << std::endl;
    
  }

  void AuxiliaryFieldsDerivative(MomentaField& der){
    der = Zero();
//    if(!M.Trivial()) 
    {
      // Auxiliary fields
      MomentaField der_temp(der.Grid());
      MomentaField X(der.Grid());
      X=Zero();
      //M.M(AuxMom, X); // X = M Aux
      // Two derivative terms
      // the Mderiv need separation of left and right terms
    std::cout << GridLogIntegrator << " AuxiliaryFieldsDerivative:norm(AuxMom)= " << std::sqrt(norm2(AuxMom)) << std::endl;
      M.MDeriv(AuxMom, der); 


      // this one should not be necessary (identical to the previous one)
      //M.MDeriv(X, AuxMom, der_temp); der += der_temp;

      der = -1.0*Implementation::projectForce(der);
      std::cout << GridLogIntegrator << " AuxiliaryFieldsDerivative:norm(der)= " << std::sqrt(norm2(der)) << std::endl;
    }
  }

  void DerivativeP(MomentaField& der){
    der = Zero();
    M.Minv(Mom, der);
    // is the projection necessary here?
    // no for fields in the algebra
    der = Implementation::projectForce(der); 
    std::cout << GridLogIntegrator << " DerivativeP:norm(der)= " << std::sqrt(norm2(der)) << std::endl;
  }

  void update_auxiliary_momenta(RealD ep){
//    if(!M.Trivial()) 
    {
      AuxMom -= ep * AuxField * HMC_MOMENTUM_DENOMINATOR;
      std::cout << GridLogIntegrator << "AuxMom update_auxiliary_fields: " << std::sqrt(norm2(AuxMom)) << std::endl;
    }
  }

  void update_auxiliary_fields(RealD ep){
//    if(!M.Trivial()) 
    {
      MomentaField tmp(AuxMom.Grid());
      MomentaField tmp2(AuxMom.Grid());
      M.M(AuxMom, tmp);
      // M.M(tmp, tmp2);
      AuxField += ep * tmp;  // M^2 AuxMom
      // factor of 2?
      std::cout << GridLogIntegrator << "AuxField update_auxiliary_fields: " << std::sqrt(norm2(AuxField)) << std::endl;
    }
  }

};

NAMESPACE_END(Grid);

