/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/observables/hmc_observable.h

Copyright (C) 2017

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

#pragma once

NAMESPACE_BEGIN(Grid);

template <class Field>
class HmcObservable {
 public:
  virtual void TrajectoryComplete(int traj,
                                  Field &U,
                                  GridSerialRNG &sRNG,
                                  GridParallelRNG &pRNG) = 0;
};

NAMESPACE_END(Grid);

NAMESPACE_CHECK(HmcObs);
#include "plaquette.h"
NAMESPACE_CHECK(Plaq);
#include "topological_charge.h"
NAMESPACE_CHECK(Topo);
#include "polyakov_loop.h"
NAMESPACE_CHECK(Polyakov);


