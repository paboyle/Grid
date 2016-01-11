    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/OverlapWilsonCayleyZolotarevFermion.h

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
#ifndef  OVERLAP_WILSON_CAYLEY_ZOLOTAREV_FERMION_H
#define  OVERLAP_WILSON_CAYLEY_ZOLOTAREV_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    template<class Impl>
    class OverlapWilsonCayleyZolotarevFermion : public MobiusZolotarevFermion<Impl>
    {
    public:
     INHERIT_IMPL_TYPES(Impl);
    public:

      // Constructors

    OverlapWilsonCayleyZolotarevFermion(GaugeField &_Umu,
					GridCartesian         &FiveDimGrid,
					GridRedBlackCartesian &FiveDimRedBlackGrid,
					GridCartesian         &FourDimGrid,
					GridRedBlackCartesian &FourDimRedBlackGrid,
					RealD _mass,RealD _M5,
					RealD lo, RealD hi,const ImplParams &p= ImplParams()) : 
      // b+c=1.0, b-c = 0 <=> b =c = 1/2
      MobiusZolotarevFermion<Impl>(_Umu,
				   FiveDimGrid,
				   FiveDimRedBlackGrid,
				   FourDimGrid,
				   FourDimRedBlackGrid,_mass,_M5,0.5,0.5,lo,hi,p)

      {}

    };

  }
}

#endif
