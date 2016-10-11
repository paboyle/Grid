    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/gauge/Photon.h

    Copyright (C) 2015

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
#ifndef QCD_PHOTON_ACTION_H
#define QCD_PHOTON_ACTION_H

namespace Grid{
  namespace QCD{

    template<class Gimpl>
    class Photon  {
      
    public:

      INHERIT_GIMPL_TYPES(Gimpl);

      enum PhotonType { FEYNMAN_L, FEYNMAN_TL };

      PhotonType GaugeBC;
      
      Photon(PhotonType _GaugeBC) : GaugeBC(_GaugeBC){};

      void  FreePropagator (const GaugeField &in,GaugeField &out) {
	FFT theFFT((GridCartesian *) in._grid);

	GaugeField in_k(in._grid);
	GaugeField prop_k(in._grid);

	theFFT.FFT_all_dim(in_k,in,FFT::forward);
        MomentumSpacePropagator(prop_k,in_k);
	theFFT.FFT_all_dim(out,prop_k,FFT::backward);
      }
      void  MomentumSpacePropagator(GaugeField &out,const GaugeField &in) {
	if ( GaugeBC == FEYNMAN_L ) { 
	  FeynmanGaugeMomentumSpacePropagator_L(out,in);
	} else if ( GaugeBC == FEYNMAN_TL ) { 
	  FeynmanGaugeMomentumSpacePropagator_TL(out,in);
	} else {  // No coulomb yet
	  assert(0);
	}
      } 
      void  FeynmanGaugeMomentumSpacePropagator_L(GaugeField &out, const GaugeField &in) { 

	FeynmanGaugeMomentumSpacePropagator_TL(out,in);

	GridBase *grid = out._grid;
	LatticeInteger     coor(grid);
	GaugeField zz(grid); zz=zero;

	// xyzt
	for(int d = 0; d < grid->_ndimension-1;d++){
	  LatticeCoordinate(coor,d);
	  out = where(coor==Integer(0),zz,out);
	}
      }

      void  FeynmanGaugeMomentumSpacePropagator_TL(GaugeField &out, const GaugeField &in) { 

	// what type LatticeComplex 
	GridBase *grid = out._grid;
	int nd = grid->_ndimension;

	typedef typename GaugeField::vector_type vector_type;
	typedef typename GaugeField::scalar_type ScalComplex;
	typedef Lattice<iSinglet<vector_type> > LatComplex;
    
	std::vector<int> latt_size   = grid->_fdimensions;

	LatComplex denom(grid); denom= zero;
	LatComplex   one(grid); one = ScalComplex(1.0,0.0);
	LatComplex   kmu(grid); 

	ScalComplex ci(0.0,1.0);
	// momphase = n * 2pi / L
	for(int mu=0;mu<Nd;mu++) {
	  
	  LatticeCoordinate(kmu,mu);
	  
	  RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];

	  kmu = TwoPiL * kmu ;
	  
	  denom = denom + 4.0*sin(kmu*0.5)*sin(kmu*0.5); // Wilson term
	}
	std::vector<int> zero_mode(nd,0);
	TComplexD Tone = ComplexD(1.0,0.0);
	TComplexD Tzero= ComplexD(0.0,0.0);

	pokeSite(Tone,denom,zero_mode); 

	denom= one/denom;
	
	pokeSite(Tzero,denom,zero_mode); 

	out = zero;
	out = in*denom;
      };

    };

}}
#endif
