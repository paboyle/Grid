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

NAMESPACE_BEGIN(Grid);


struct GparityWilsonImplParams {
  Coordinate twists;
                     //mu=Nd-1 is assumed to be the time direction and a twist value of 1 indicates antiperiodic BCs
  Coordinate dirichlet; // Blocksize of dirichlet BCs
  int  partialDirichlet;
  GparityWilsonImplParams() : twists(Nd, 0) {
    dirichlet.resize(0);
    partialDirichlet=0;
  };
};
  
struct XconjWilsonImplParams {
  Coordinate dirichlet; // Blocksize of dirichlet BCs
  int  partialDirichlet;

  Coordinate twists; //Here the first Nd-1 directions are treated as "spatial", and a twist value of 1 indicates G-parity BCs in that direction. 
                     //mu=Nd-1 is assumed to be the time direction and a twist value of 1 indicates antiperiodic BCs
  ComplexD boundary_phase; //+1 for X-conjugate, -1 for Xbar-conjugate, or other
  XconjWilsonImplParams() : twists(Nd, 0), boundary_phase(1.0) { dirichlet.resize(0); partialDirichlet=0; };
};


  
struct WilsonImplParams {
  bool overlapCommsCompute;
  Coordinate dirichlet; // Blocksize of dirichlet BCs
  int  partialDirichlet;
  AcceleratorVector<Real,Nd> twist_n_2pi_L;
  AcceleratorVector<Complex,Nd> boundary_phases;
  WilsonImplParams()  {
    dirichlet.resize(0);
    partialDirichlet=0;
    boundary_phases.resize(Nd, 1.0);
      twist_n_2pi_L.resize(Nd, 0.0);
  };
  WilsonImplParams(const AcceleratorVector<Complex,Nd> phi) : boundary_phases(phi), overlapCommsCompute(false) {
    twist_n_2pi_L.resize(Nd, 0.0);
    partialDirichlet=0;
    dirichlet.resize(0);
  }
};

struct StaggeredImplParams {
  Coordinate dirichlet; // Blocksize of dirichlet BCs
  int  partialDirichlet;
  StaggeredImplParams()
  {
    partialDirichlet=0;
    dirichlet.resize(0);
  };
};
  
  struct OneFlavourRationalParams : Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(OneFlavourRationalParams, 
				    RealD, lo, 
				    RealD, hi, 
				    int,   MaxIter, 
				    RealD, tolerance, 
				    RealD, mdtolerance, 
				    int,   degree, 
				    int,   precision,
				    int,   BoundsCheckFreq,
				    RealD, BoundsCheckTol);
    
  // MaxIter and tolerance, vectors??
    
  // constructor 
  OneFlavourRationalParams(	RealD _lo      = 0.0, 
				RealD _hi      = 1.0, 
				int _maxit     = 1000,
				RealD tol      = 1.0e-8, 
                           	int _degree    = 10,
				int _precision = 64,
				int _BoundsCheckFreq=20,
				RealD mdtol    = 1.0e-6,
				double _BoundsCheckTol=1e-6)
      : lo(_lo),
	hi(_hi),
	MaxIter(_maxit),
	tolerance(tol),
        mdtolerance(mdtol),
	degree(_degree),
        precision(_precision),
        BoundsCheckFreq(_BoundsCheckFreq),
        BoundsCheckTol(_BoundsCheckTol){};
  };
  
  /*Action parameters for the generalized rational action
    The approximation is for (M^dag M)^{1/inv_pow}
    where inv_pow is the denominator of the fractional power.
    Default inv_pow=2 for square root, making this equivalent to 
    the OneFlavourRational action
  */
    struct RationalActionParams : Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(RationalActionParams, 
				    int, inv_pow, 
				    RealD, lo, //low eigenvalue bound of rational approx
				    RealD, hi, //high eigenvalue bound of rational approx
				    int,   MaxIter,  //maximum iterations in msCG
				    RealD, action_tolerance,  //msCG tolerance in action evaluation
				    int,   action_degree, //rational approx tolerance in action evaluation
				    RealD, md_tolerance,  //msCG tolerance in MD integration
				    int,   md_degree, //rational approx tolerance in MD integration
				    int,   precision, //precision of floating point arithmetic
				    int,   BoundsCheckFreq); //frequency the approximation is tested (with Metropolis degree/tolerance); 0 disables the check
  // constructor 
  RationalActionParams(int _inv_pow = 2,
		       RealD _lo      = 0.0, 
		       RealD _hi      = 1.0, 
		       int _maxit     = 1000,
		       RealD _action_tolerance      = 1.0e-8, 
		       int _action_degree    = 10,
		       RealD _md_tolerance      = 1.0e-8, 
		       int _md_degree    = 10,
		       int _precision = 64,
		       int _BoundsCheckFreq=20)
    : inv_pow(_inv_pow), 
      lo(_lo),
      hi(_hi),
      MaxIter(_maxit),
      action_tolerance(_action_tolerance),
      action_degree(_action_degree),
      md_tolerance(_md_tolerance),
      md_degree(_md_degree),
      precision(_precision),
      BoundsCheckFreq(_BoundsCheckFreq){};
  };


NAMESPACE_END(Grid);

#endif
