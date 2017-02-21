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
  virtual void MInvSquareRoot(Field&) = 0;
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
  virtual void MInvSquareRoot(Field& P){
    // do nothing
  }

};

///////////////////////////////
// Generalised momenta
///////////////////////////////

template <typename Implementation, typename Metric>
class GeneralisedMomenta{
public:
  typedef typename Implementation::Field MomentaField;  //for readability

  Metric M;
  MomentaField Mom;

  GeneralisedMomenta(GridBase* grid): Mom(grid){}


  void MomentaDistribution(GridParallelRNG& pRNG){
    // Generate a distribution for
    // 1/2 P^dag G P
    // where G = M^-1

    // Generate gaussian momenta
    Implementation::generate_momenta(Mom, pRNG);
    // Modify the distribution with the metric
    M.MInvSquareRoot(Mom);
  }

  void Derivative(MomentaField& in, MomentaField& der){
    // Compute the derivative of the kinetic term
    // with respect to the gauge field
    MomentaField MomDer(in._grid);
    MomentaField X(in._grid);

    M.Minv(in, X); // X = G in
    M.MDeriv(X, MomDer, DaggerNo); // MomDer = dM/dU X
    // MomDer is just the derivative
    MomDer = adj(X)* MomDer;
    // Traceless Antihermitian
    // assuming we are in the algebra
    der = Implementation::projectForce(MomDer);
  }

  void DerivativeP(MomentaField& der){
    M.Minv(Mom, der);
  }

};








}
}


#endif //METRIC_H