/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/modules/polyakov_line.h

Copyright (C) 2017

Author: David Preti <david.preti@csic.es>

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

#pragma once

NAMESPACE_BEGIN(Grid);

// this is only defined for a gauge theory
template <class Impl>
class SpatialPolyakovLogger : public HmcObservable<typename Impl::Field> {
 public:
    // here forces the Impl to be of gauge fields
    // if not the compiler will complain
    INHERIT_GIMPL_TYPES(Impl);

     // necessary for HmcObservable compatibility
    typedef typename Impl::Field Field;

    void TrajectoryComplete(int traj,
                            Field &U,
                            GridSerialRNG &sRNG,
                            GridParallelRNG &pRNG) {

    // Save current numerical output precision
    int def_prec = std::cout.precision();

    // Assume that the dimensions are D=3+1
    int Ndim = 3;
    ComplexD polyakov;
   
    // Iterate over the spatial directions and print the average spatial polyakov loop
    // over them 
    for (int idx=0; idx<Ndim; idx++) {
        polyakov = WilsonLoops<Impl>::avgPolyakovLoop(U, idx);
    
        std::cout << GridLogMessage
            << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
            << "Polyakov Loop in the " << idx << " spatial direction : [ " << traj << " ] "<< polyakov << std::endl;

    }

    // Return to original output precision
    std::cout.precision(def_prec);

  }
};

NAMESPACE_END(Grid);
