/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/scalar/CovariantLaplacian.h

Copyright (C) 2016

Author: Azusa Yamaguchi

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
#pragma once

namespace Grid {
namespace QCD {

template <class Gimpl> class CovariantSmearing : public Gimpl 
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);

  typedef typename Gimpl::GaugeLinkField GaugeMat;
  typedef typename Gimpl::GaugeField GaugeLorentz;

  template<typename T>
  static void GaussianSmear(const std::vector<LatticeColourMatrix>& U, 
			    T& chi, 
			    const Real& width, int Iterations, int orthog)
  {
    GridBase *grid = chi._grid;
    T psi(grid);

    ////////////////////////////////////////////////////////////////////////////////////
    // Follow Chroma conventions for width to keep compatibility with previous data
    // Free field iterates 
    //   chi = (1 - w^2/4N p^2)^N chi
    //
    //       ~ (e^(-w^2/4N p^2)^N chi
    //       ~ (e^(-w^2/4 p^2) chi
    //       ~ (e^(-w'^2/2 p^2) chi          [ w' = w/sqrt(2) ]
    //
    // Which in coordinate space is proportional to
    //
    //   e^(-x^2/w^2) = e^(-x^2/2w'^2) 
    //
    // The 4 is a bit unconventional from Gaussian width perspective, but... it's Chroma convention.
    // 2nd derivative approx d^2/dx^2  =  x+mu + x-mu - 2x
    //
    // d^2/dx^2 = - p^2
    //
    // chi = ( 1 + w^2/4N d^2/dx^2 )^N chi
    //
    ////////////////////////////////////////////////////////////////////////////////////
    Real coeff = (width*width) / Real(4*Iterations);
  
    int dims = Nd;
    if( orthog < Nd ) dims=Nd-1;

    for(int n = 0; n < Iterations; ++n) {
      psi = (-2.0*dims)*chi;
      for(int mu=0;mu<Nd;mu++) {
	if ( mu != orthog ) { 
	  psi = psi + Gimpl::CovShiftForward(U[mu],mu,chi);    
	  psi = psi + Gimpl::CovShiftBackward(U[mu],mu,chi);    
	}
      }
      chi = chi + coeff*psi;
    }
  }
};
}}
