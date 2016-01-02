    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/gauge/WilsonGaugeAction.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
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
#ifndef QCD_WILSON_GAUGE_ACTION_H
#define QCD_WILSON_GAUGE_ACTION_H

namespace Grid{
  namespace QCD{
    
    ////////////////////////////////////////////////////////////////////////
    // Wilson Gauge Action .. should I template the Nc etc..
    ////////////////////////////////////////////////////////////////////////
    template<class Gimpl>
    class WilsonGaugeAction : public Action<typename Gimpl::GaugeField> {
    public:

      INHERIT_GIMPL_TYPES(Gimpl);

      //      typedef LorentzScalar<GaugeField> GaugeLinkField;

    private:
      RealD beta;
    public:
    WilsonGaugeAction(RealD b):beta(b){};
      
      virtual void refresh(const GaugeField &U, GridParallelRNG& pRNG) {}; // noop as no pseudoferms
      
      virtual RealD S(const GaugeField &U) {
	RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(U);
	RealD vol = U._grid->gSites();
	RealD action=beta*(1.0 -plaq)*(Nd*(Nd-1.0))*vol*0.5;
	return action;
      };

      virtual void deriv(const GaugeField &U,GaugeField & dSdU) {
	//not optimal implementation FIXME
	//extend Ta to include Lorentz indexes

	//RealD factor = 0.5*beta/RealD(Nc);
	RealD factor = 0.5*beta/RealD(Nc);

	GaugeLinkField Umu(U._grid);
	GaugeLinkField dSdU_mu(U._grid);
	for (int mu=0; mu < Nd; mu++){

	  Umu = PeekIndex<LorentzIndex>(U,mu);

	  // Staple in direction mu
	  WilsonLoops<Gimpl>::Staple(dSdU_mu,U,mu);
	  dSdU_mu = Ta(Umu*dSdU_mu)*factor;
	  PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);
	}
      };
    };
    
  }
}

#endif
