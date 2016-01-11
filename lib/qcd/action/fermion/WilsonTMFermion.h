    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/WilsonTMFermion.h

    Copyright (C) 2015

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
#ifndef  GRID_QCD_WILSON_TM_FERMION_H
#define  GRID_QCD_WILSON_TM_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    template<class Impl>
    class WilsonTMFermion : public WilsonFermion<Impl>
    {
    public:
      INHERIT_IMPL_TYPES(Impl);
    public:

      virtual void   Instantiatable(void) {};
      // Constructors
      WilsonTMFermion(GaugeField &_Umu,
		    GridCartesian         &Fgrid,
		    GridRedBlackCartesian &Hgrid, 
		    RealD _mass,
		    RealD _mu,
		    const ImplParams &p= ImplParams()
		      ) :
	WilsonFermion<Impl>(_Umu,
			    Fgrid,
			    Hgrid,
			    _mass,p)

      {
	mu = _mu;
      }


    // allow override for twisted mass and clover
    virtual void Mooee(const FermionField &in, FermionField &out) ;
    virtual void MooeeDag(const FermionField &in, FermionField &out) ;
    virtual void MooeeInv(const FermionField &in, FermionField &out) ;
    virtual void MooeeInvDag(const FermionField &in, FermionField &out) ;

  private:
     RealD mu; // TwistedMass parameter

  };

}}

#endif
