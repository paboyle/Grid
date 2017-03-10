/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/integrators/Integrator.h

Copyright (C) 2015

Author: Guido Cossu <guido.cossu@ed.ac.uk>

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
#ifndef METRIC_H
#define METRIC_H

namespace Grid{
namespace QCD{

template <typename Field> 
class Metric{
public:
  virtual void ImportGauge(const Field&)   = 0;
  virtual void M(const Field&, Field&)     = 0;
  virtual void Minv(const Field&, Field&)  = 0;
  virtual void MSquareRoot(Field&) = 0;
  virtual void MInvSquareRoot(Field&) = 0;
  virtual void MDeriv(const Field&, Field&) = 0;
  virtual void MDeriv(const Field&, const Field&, Field&) = 0;
};


// Need a trivial operator
template <typename Field>
class TrivialMetric : public Metric<Field>{
public:
  virtual void ImportGauge(const Field&){};
  virtual void M(const Field& in, Field& out){
    out = in;
  }
  virtual void Minv(const Field& in, Field& out){
    out = in;
  }
  virtual void MSquareRoot(Field& P){
    // do nothing
  }
  virtual void MInvSquareRoot(Field& P){
    // do nothing
  }
  virtual void MDeriv(const Field& in, Field& out){
    out = zero;
  }
  virtual void MDeriv(const Field& left, const Field& right, Field& out){
    out = zero;
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
  void MomentaDistribution(GridParallelRNG& pRNG){
    // Generate a distribution for
    // P^dag G P
    // where G = M^-1

    // Generate gaussian momenta
    Implementation::generate_momenta(Mom, pRNG);
    // Modify the distribution with the metric
    M.MSquareRoot(Mom);

    if (1) {
      // Auxiliary momenta
      // do nothing if trivial, so hide in the metric
      MomentaField AuxMomTemp(Mom._grid);
      Implementation::generate_momenta(AuxMom, pRNG);
      Implementation::generate_momenta(AuxField, pRNG);
      // Modify the distribution with the metric
      // Aux^dag M Aux
      M.MInvSquareRoot(AuxMom);  // AuxMom = M^{-1/2} AuxMomTemp
    }
  }

  // Correct
  RealD MomentaAction(){
    MomentaField inv(Mom._grid);
    inv = zero;
    M.Minv(Mom, inv);
    LatticeComplex Hloc(Mom._grid);
    Hloc = zero;
    for (int mu = 0; mu < Nd; mu++) {
      // This is not very general
      // hide in the metric
      auto Mom_mu = PeekIndex<LorentzIndex>(Mom, mu);
      auto inv_mu = PeekIndex<LorentzIndex>(inv, mu);
      Hloc += trace(Mom_mu * inv_mu);
    }

    if (1) {
      // Auxiliary Fields
      // hide in the metric
      M.M(AuxMom, inv);
      for (int mu = 0; mu < Nd; mu++) {
        // This is not very general
        // hide in the operators
        auto inv_mu = PeekIndex<LorentzIndex>(inv, mu);
        auto am_mu = PeekIndex<LorentzIndex>(AuxMom, mu);
        auto af_mu = PeekIndex<LorentzIndex>(AuxField, mu);
        Hloc += trace(am_mu * inv_mu);// p M p
        Hloc += trace(af_mu * af_mu);
      }
    }

    Complex Hsum = sum(Hloc);
    return Hsum.real();
  }

  // Correct
  void DerivativeU(MomentaField& in, MomentaField& der){

    // Compute the derivative of the kinetic term
    // with respect to the gauge field
    MomentaField MDer(in._grid);
    MomentaField X(in._grid);
    X = zero;
    M.Minv(in, X);  // X = G in
    M.MDeriv(X, MDer);  // MDer = U * dS/dU
    der = Implementation::projectForce(MDer);  // Ta if gauge fields
    
  }

  void AuxiliaryFieldsDerivative(MomentaField& der){
    der = zero;
    if (1){
    // Auxiliary fields
    MomentaField der_temp(der._grid);
    MomentaField X(der._grid);
    X=zero;
    //M.M(AuxMom, X); // X = M Aux
    // Two derivative terms
    // the Mderiv need separation of left and right terms
    M.MDeriv(AuxMom, der); 


    // this one should not be necessary (identical to the previous one)
    //M.MDeriv(X, AuxMom, der_temp); der += der_temp;

    der = -1.0*Implementation::projectForce(der);
    }
  }

  void DerivativeP(MomentaField& der){
    der = zero;
    M.Minv(Mom, der);
    // is the projection necessary here?
    // no for fields in the algebra
    der = Implementation::projectForce(der); 
  }

  void update_auxiliary_momenta(RealD ep){
    if(1){
      AuxMom -= ep * AuxField;
    }
  }

  void update_auxiliary_fields(RealD ep){
    if (1) {
      MomentaField tmp(AuxMom._grid);
      MomentaField tmp2(AuxMom._grid);
      M.M(AuxMom, tmp);
      // M.M(tmp, tmp2);
      AuxField += ep * tmp;  // M^2 AuxMom
      // factor of 2?
    }
  }

};








}
}


#endif //METRIC_H