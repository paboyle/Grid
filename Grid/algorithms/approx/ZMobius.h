/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/approx/ZMobius.h

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@phys.columbia.edu>

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
#ifndef GRID_ZMOBIUS_APPROX_H
#define GRID_ZMOBIUS_APPROX_H

#include <Grid/GridCore.h>

NAMESPACE_BEGIN(Grid);
NAMESPACE_BEGIN(Approx);

//Compute the Zmobius Omega parameters suitable for eigenvalue range   -lambda_bound <= lambda <= lambda_bound
//Note omega_i = 1/(b_i + c_i)   where b_i and c_i are the Mobius parameters
void computeZmobiusOmega(std::vector<ComplexD> &omega_out, const int Ls_out,
			 const std::vector<RealD> &omega_in, const int Ls_in,
			 const RealD lambda_bound);
  
//mobius_param = b+c   with b-c=1
void computeZmobiusOmega(std::vector<ComplexD> &omega_out, const int Ls_out, const RealD mobius_param, const int Ls_in, const RealD lambda_bound);

//ZMobius class takes  gamma_i = (b+c) omega_i as its input, where b, c are factored out
void computeZmobiusGamma(std::vector<ComplexD> &gamma_out, 
			 const RealD mobius_param_out, const int Ls_out, 
			 const RealD mobius_param_in, const int Ls_in,
			 const RealD lambda_bound);

//Assumes mobius_param_out == mobius_param_in
void computeZmobiusGamma(std::vector<ComplexD> &gamma_out, const int Ls_out, const RealD mobius_param, const int Ls_in, const RealD lambda_bound);

NAMESPACE_END(Approx);
NAMESPACE_END(Grid);

#endif
