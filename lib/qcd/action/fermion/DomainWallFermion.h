    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/DomainWallFermion.h

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
#ifndef  GRID_QCD_DOMAIN_WALL_FERMION_H
#define  GRID_QCD_DOMAIN_WALL_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    template<class Impl>
    class DomainWallFermion : public CayleyFermion5D<Impl>
    {
    public:
     INHERIT_IMPL_TYPES(Impl);
    public:

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
	
	//	std::cout<<GridLogMessage << "DomainWallFermion with Ls="<<this->Ls<<std::endl;
	// Call base setter
	this->SetCoefficientsTanh(zdata,1.0,0.0);

	Approx::zolotarev_free(zdata);
      }

    };

  }
}

#endif
