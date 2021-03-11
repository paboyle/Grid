/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/integrators/MomentumFilter.h

Copyright (C) 2015

Author: Christopher Kelly <ckelly@bnl.gov>
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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
//--------------------------------------------------------------------
#ifndef MOMENTUM_FILTER
#define MOMENTUM_FILTER

NAMESPACE_BEGIN(Grid);

//These filter objects allow the user to manipulate the conjugate momentum as part of the update / refresh

template<typename MomentaField>
struct MomentumFilterBase{
  virtual void applyFilter(MomentaField &P) const;
};

//Do nothing
template<typename MomentaField>
struct MomentumFilterNone: public MomentumFilterBase<MomentaField>{
  void applyFilter(MomentaField &P) const override{}
};

//Multiply each site/direction by a Lorentz vector complex number field
//Can be used to implement a mask, zeroing out sites
template<typename MomentaField>
struct MomentumFilterApplyPhase: public MomentumFilterBase<MomentaField>{
  typedef typename MomentaField::vector_type vector_type; //SIMD-vectorized complex type
  typedef typename MomentaField::scalar_type scalar_type; //scalar complex type
  typedef iVector<iScalar<iScalar<vector_type> >, Nd > LorentzScalarType; //complex phase for each site/direction
  typedef Lattice<LorentzScalarType> LatticeLorentzScalarType;
  
  LatticeLorentzScalarType phase;
 
  MomentumFilterApplyPhase(const LatticeLorentzScalarType _phase): phase(_phase){}

  //Default to uniform field of (1,0)
  MomentumFilterApplyPhase(GridBase* _grid): phase(_grid){
    LorentzScalarType one;
    for(int mu=0;mu<Nd;mu++)
      one(mu)()() = scalar_type(1.);
    
    phase = one;
  }

  void applyFilter(MomentaField &P) const override{
    conformable(P,phase);
    autoView( P_v , P, AcceleratorWrite);
    autoView( phase_v , phase, AcceleratorRead);

    accelerator_for(ss,P_v.size(),MomentaField::vector_type::Nsimd(),{
    	auto site_mom = P_v(ss);
    	auto site_phase = phase_v(ss);
	for(int mu=0;mu<Nd;mu++)
	  site_mom(mu) = site_mom(mu) * site_phase(mu);
    	coalescedWrite(P_v[ss], site_mom);
      });
    
  }


};




NAMESPACE_END(Grid);

#endif
