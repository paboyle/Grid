    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/ActionParams.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#ifndef GRID_QCD_ACTION_PARAMS_H
#define GRID_QCD_ACTION_PARAMS_H

namespace Grid {
namespace QCD {

    // These can move into a params header and be given MacroMagic serialisation
    struct GparityWilsonImplParams {
      bool overlapCommsCompute;
      std::vector<int> twists; 
      GparityWilsonImplParams () : twists(Nd,0), overlapCommsCompute(false) {};

    };

    struct WilsonImplParams {
      bool overlapCommsCompute;
      WilsonImplParams() : overlapCommsCompute(false) {};
    };

    struct OneFlavourRationalParams { 
      RealD  lo;
      RealD  hi;
      int MaxIter;   // Vector?
      RealD tolerance; // Vector? 
      int    degree=10;
      int precision=64;

      OneFlavourRationalParams (RealD _lo,RealD _hi,int _maxit,RealD tol=1.0e-8,int _degree = 10,int _precision=64) :
        lo(_lo), hi(_hi), MaxIter(_maxit), tolerance(tol), degree(_degree), precision(_precision)
      {};
    };

}}

#endif
