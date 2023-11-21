/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/gauge/WilsonGaugeAction.h

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: paboyle <paboyle@ph.ed.ac.uk>
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
#ifndef QCD_WILSON_GAUGE_ACTION_H
#define QCD_WILSON_GAUGE_ACTION_H

NAMESPACE_BEGIN(Grid);

////////////////////////////////////////////////////////////////////////
// Wilson Gauge Action .. should I template the Nc etc..
////////////////////////////////////////////////////////////////////////
template <class Gimpl>
class WilsonGaugeAction : public Action<typename Gimpl::GaugeField> {
public:  
  INHERIT_GIMPL_TYPES(Gimpl);
  typedef GaugeImplParams ImplParams;
  ImplParams Params;

  /////////////////////////// constructors
  explicit WilsonGaugeAction(RealD beta_,
		  const ImplParams &p = ImplParams()
		  ):beta(beta_),Params(p){};

  virtual std::string action_name() {return "WilsonGaugeAction";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "[WilsonGaugeAction] Beta: " << beta << std::endl;
    return sstream.str();
  }

  virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG &pRNG){};  // noop as no pseudoferms

// Umu<->U maximally confusing
  virtual void boundary(const GaugeField &Umu, GaugeField &Ub){
    typedef typename Simd::scalar_type scalar_type;
    assert(Params.boundary_phases.size() == Nd);
    GridBase *GaugeGrid=Umu.Grid();
    GaugeLinkField U(GaugeGrid);
    GaugeLinkField tmp(GaugeGrid);

    Lattice<iScalar<vInteger> > coor(GaugeGrid);
    for (int mu = 0; mu < Nd; mu++) {
	////////// boundary phase /////////////
      auto pha = Params.boundary_phases[mu];
      scalar_type phase( real(pha),imag(pha) );
      std::cout<< GridLogIterative << "[WilsonGaugeAction] boundary "<<mu<<" "<<phase<< std::endl; 

	int L   = GaugeGrid->GlobalDimensions()[mu];
        int Lmu = L - 1;

      LatticeCoordinate(coor, mu);

      U = PeekIndex<LorentzIndex>(Umu, mu);
      tmp = where(coor == Lmu, phase * U, U);
      PokeIndex<LorentzIndex>(Ub, tmp, mu);
//      PokeIndex<LorentzIndex>(Ub, U, mu);
//      PokeIndex<LorentzIndex>(Umu, tmp, mu);

    }
  };

  virtual RealD S(const GaugeField &U) {
    GaugeField Ub(U.Grid());
    this->boundary(U,Ub);
    static RealD lastG=0.;
    RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(Ub);
    RealD vol = Ub.Grid()->gSites();
    RealD action = beta * (1.0 - plaq) * (Nd * (Nd - 1.0)) * vol * 0.5;
    std::cout << GridLogMessage << "[WilsonGaugeAction] dH: " << action-lastG << std::endl;
    RealD plaq_o = WilsonLoops<Gimpl>::avgPlaquette(U);
    RealD action_o = beta * (1.0 - plaq_o) * (Nd * (Nd - 1.0)) * vol * 0.5;
    std::cout << GridLogMessage << "[WilsonGaugeAction] U: " << action_o <<" Ub: "<< action  << std::endl;
    lastG=action;
    return action;
  };

  virtual void deriv(const GaugeField &U, GaugeField &dSdU) {
    GaugeField Ub(U.Grid());
    this->boundary(U,Ub);
    // not optimal implementation FIXME
    // extend Ta to include Lorentz indexes

    RealD factor = 0.5 * beta / RealD(Nc);

    GaugeLinkField Umu(U.Grid());
    GaugeLinkField dSdU_mu(U.Grid());
    for (int mu = 0; mu < Nd; mu++) {

      Umu = PeekIndex<LorentzIndex>(Ub, mu);
      // Staple in direction mu
      WilsonLoops<Gimpl>::Staple(dSdU_mu, Ub, mu);
      dSdU_mu = Ta(Umu * dSdU_mu) * factor;
      
      PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);
    }
  }
private:
  RealD beta;  
 };

NAMESPACE_END(Grid);
#endif
