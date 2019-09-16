/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/OverlapWilsonCayleyTanhFermion.h

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
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
#ifndef OVERLAP_WILSON_CAYLEY_TANH_FERMION_H
#define OVERLAP_WILSON_CAYLEY_TANH_FERMION_H

#include <Grid/qcd/action/fermion/FermionCore.h>

NAMESPACE_BEGIN(Grid);

template<class Impl>
class OverlapWilsonCayleyTanhFermion : public MobiusFermion<Impl>
{
public:
  INHERIT_IMPL_TYPES(Impl);
public:

     void  MomentumSpacePropagator(FermionField &out,const FermionField &in,RealD _m,std::vector<double> twist) {
       this->MomentumSpacePropagatorHw(out,in,_m,twist);
  };

  // Constructors
  OverlapWilsonCayleyTanhFermion(GaugeField &_Umu,
				 GridCartesian         &FiveDimGrid,
				 GridRedBlackCartesian &FiveDimRedBlackGrid,
				 GridCartesian         &FourDimGrid,
				 GridRedBlackCartesian &FourDimRedBlackGrid,
				 RealD _mass,RealD _M5,
				 RealD scale,const ImplParams &p= ImplParams()) :
      
    // b+c=scale, b-c = 0 <=> b =c = scale/2
    MobiusFermion<Impl>(_Umu,
			FiveDimGrid,
			FiveDimRedBlackGrid,
			FourDimGrid,
			FourDimRedBlackGrid,_mass,_M5,0.5*scale,0.5*scale,p)
  {
  }
};

NAMESPACE_END(Grid);
#endif
