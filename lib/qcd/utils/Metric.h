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
  virtual void MDeriv(const Field&, Field&) = 0;
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
  virtual void MDeriv(const Field& in, Field& out){
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

  GeneralisedMomenta(GridBase* grid, Metric<MomentaField>& M): M(M), Mom(grid){}

  // Correct
  void MomentaDistribution(GridParallelRNG& pRNG){
    // Generate a distribution for
    // 1/2 P^dag G P
    // where G = M^-1

    // Generate gaussian momenta
    Implementation::generate_momenta(Mom, pRNG);
    // Modify the distribution with the metric
    M.MSquareRoot(Mom);
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
      // hide in the operators
      auto Mom_mu = PeekIndex<LorentzIndex>(Mom, mu);
      auto inv_mu = PeekIndex<LorentzIndex>(inv, mu);
      Hloc += trace(Mom_mu * inv_mu);
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


  //   
  void DerivativeP(MomentaField& der){
    der = zero;
    M.Minv(Mom, der);
    der = Implementation::projectForce(der); 
  }

};








}
}


#endif //METRIC_H