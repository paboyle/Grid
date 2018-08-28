    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/DomainWallFermion.h

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Vera Guelpers <V.M.Guelpers@soton.ac.uk>

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
#ifndef  GRID_QCD_DOMAIN_WALL_FERMION_H
#define  GRID_QCD_DOMAIN_WALL_FERMION_H

#include <Grid/qcd/action/fermion/FermionCore.h>

namespace Grid {

  namespace QCD {

    template<class Impl>
    class DomainWallFermion : public CayleyFermion5D<Impl>
    {
    public:
     INHERIT_IMPL_TYPES(Impl);
    public:

      void FreePropagator(const FermionField &in,FermionField &out,RealD mass, std::vector<double> twist, bool fiveD) {
	FermionField in_k(in._grid);
	FermionField prop_k(in._grid);

	FFT theFFT((GridCartesian *) in._grid);

	//phase for boundary condition
	ComplexField coor(in._grid);
	ComplexField ph(in._grid);  ph = zero;
	FermionField in_buf(in._grid); in_buf = zero;
	Complex ci(0.0,1.0);
	assert(twist.size() == Nd);//check that twist is Nd
	int shift = 0;
	if(fiveD) shift = 1;
	for(unsigned int nu = 0; nu < Nd; nu++)
	{
	  // Shift coordinate lattice index by 1 to account for 5th dimension.
          LatticeCoordinate(coor, nu + shift);
	  ph = ph + twist[nu]*coor*((1./(in._grid->_fdimensions[nu+shift])));
	}
	in_buf = exp((Real)(2.0*M_PI)*ci*ph*(-1.0))*in;

	if(fiveD){//FFT only on temporal and spatial dimensions
          std::vector<int> mask(Nd+1,1); mask[0] = 0;
	  theFFT.FFT_dim_mask(in_k,in_buf,mask,FFT::forward);
          this->MomentumSpacePropagatorHt_5d(prop_k,in_k,mass,twist);
          theFFT.FFT_dim_mask(out,prop_k,mask,FFT::backward);
        }
	else{
	  theFFT.FFT_all_dim(in_k,in,FFT::forward);
          this->MomentumSpacePropagatorHt(prop_k,in_k,mass,twist);
	  theFFT.FFT_all_dim(out,prop_k,FFT::backward);
        }

	//phase for boundary condition
	out = out * exp((Real)(2.0*M_PI)*ci*ph);
      };

      virtual void FreePropagator(const FermionField &in,FermionField &out,RealD mass,std::vector<double> twist) {
        bool fiveD = true; //5d propagator by default
        FreePropagator(in,out,mass,twist,fiveD);
      };

      virtual void FreePropagator(const FermionField &in,FermionField &out,RealD mass, bool fiveD) {
	std::vector<double> twist(Nd,0.0); //default: periodic boundarys in all directions
        FreePropagator(in,out,mass,twist,fiveD);
      };

      virtual void FreePropagator(const FermionField &in,FermionField &out,RealD mass) {
        bool fiveD = true; //5d propagator by default
	std::vector<double> twist(Nd,0.0); //default: periodic boundarys in all directions
        FreePropagator(in,out,mass,twist,fiveD);
      };

      virtual void   Instantiatable(void) {};
      // Constructors
      DomainWallFermion(GaugeField &_Umu,
			GridCartesian         &FiveDimGrid,
			GridRedBlackCartesian &FiveDimRedBlackGrid,
			GridCartesian         &FourDimGrid,
			GridRedBlackCartesian &FourDimRedBlackGrid,
			RealD _mass,RealD _M5,const ImplParams &p= ImplParams()) : 


      CayleyFermion5D<Impl>(_Umu,
			    FiveDimGrid,
			    FiveDimRedBlackGrid,
			    FourDimGrid,
			    FourDimRedBlackGrid,_mass,_M5,p)

      {
	RealD eps = 1.0;

	Approx::zolotarev_data *zdata = Approx::higham(eps,this->Ls);// eps is ignored for higham
	assert(zdata->n==this->Ls);
	
	std::cout<<GridLogMessage << "DomainWallFermion with Ls="<<this->Ls<<std::endl;
	// Call base setter
	this->SetCoefficientsTanh(zdata,1.0,0.0);

	Approx::zolotarev_free(zdata);
      }

    };

  }
}

#endif
