/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/ActionParams.h

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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

#ifndef GRID_QCD_ACTION_PARAMS_H
#define GRID_QCD_ACTION_PARAMS_H

namespace Grid {
namespace QCD {

  // These can move into a params header and be given MacroMagic serialisation
  struct GparityWilsonImplParams {
    bool overlapCommsCompute;
    std::vector<int> twists;
    GparityWilsonImplParams() : twists(Nd, 0), overlapCommsCompute(false){};
  };
  
  struct WilsonImplParams {
    bool overlapCommsCompute;
    std::vector<Real> twist_n_2pi_L;
    std::vector<Complex> boundary_phases;
    WilsonImplParams() : overlapCommsCompute(false) {
      boundary_phases.resize(Nd, 1.0);
      twist_n_2pi_L.resize(Nd, 0.0);
    };
    WilsonImplParams(const std::vector<Complex> phi) : boundary_phases(phi), overlapCommsCompute(false) {
      twist_n_2pi_L.resize(Nd, 0.0);
    }
  };

  struct StaggeredImplParams {
    StaggeredImplParams()  {};
  };
  
  struct OneFlavourRationalParams : Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(OneFlavourRationalParams, 
				    RealD, lo, 
				    RealD, hi, 
				    int,   MaxIter, 
				    RealD, tolerance, 
				    int,   degree, 
				    int,   precision);
    
    // MaxIter and tolerance, vectors??
    
    // constructor 
    OneFlavourRationalParams(	RealD _lo      = 0.0, 
				RealD _hi      = 1.0, 
				int _maxit     = 1000,
				RealD tol      = 1.0e-8, 
                           	int _degree    = 10,
				int _precision = 64)
      : lo(_lo),
	hi(_hi),
	MaxIter(_maxit),
	tolerance(tol),
	degree(_degree),
	precision(_precision){};
  };
  
  
}
}




#endif
