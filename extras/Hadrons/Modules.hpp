/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules.hpp

Copyright (C) 2015
Copyright (C) 2016
Copyright (C) 2017

Author: Antonin Portelli <antonin.portelli@me.com>

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

#include <Grid/Hadrons/Modules/MAction/DWF.hpp>
#include <Grid/Hadrons/Modules/MAction/Wilson.hpp>
#include <Grid/Hadrons/Modules/MContraction/Baryon.hpp>
#include <Grid/Hadrons/Modules/MContraction/DiscLoop.hpp>
#include <Grid/Hadrons/Modules/MContraction/Gamma3pt.hpp>
#include <Grid/Hadrons/Modules/MContraction/Meson.hpp>
#include <Grid/Hadrons/Modules/MContraction/WardIdentity.hpp>
#include <Grid/Hadrons/Modules/MContraction/WeakHamiltonian.hpp>
#include <Grid/Hadrons/Modules/MContraction/WeakHamiltonianEye.hpp>
#include <Grid/Hadrons/Modules/MContraction/WeakHamiltonianNonEye.hpp>
#include <Grid/Hadrons/Modules/MContraction/WeakNeutral4ptDisc.hpp>
#include <Grid/Hadrons/Modules/MFermion/GaugeProp.hpp>
#include <Grid/Hadrons/Modules/MGauge/Load.hpp>
#include <Grid/Hadrons/Modules/MGauge/Random.hpp>
#include <Grid/Hadrons/Modules/MGauge/StochEm.hpp>
#include <Grid/Hadrons/Modules/MGauge/Unit.hpp>
#include <Grid/Hadrons/Modules/MLoop/NoiseLoop.hpp>
#include <Grid/Hadrons/Modules/MScalar/ChargedProp.hpp>
#include <Grid/Hadrons/Modules/MScalar/FreeProp.hpp>
#include <Grid/Hadrons/Modules/MScalar/Scalar.hpp>
#include <Grid/Hadrons/Modules/MSink/Point.hpp>
#include <Grid/Hadrons/Modules/MSink/Smear.hpp>
#include <Grid/Hadrons/Modules/MSolver/RBPrecCG.hpp>
#include <Grid/Hadrons/Modules/MSource/Point.hpp>
#include <Grid/Hadrons/Modules/MSource/SeqConserved.hpp>
#include <Grid/Hadrons/Modules/MSource/SeqGamma.hpp>
#include <Grid/Hadrons/Modules/MSource/Wall.hpp>
#include <Grid/Hadrons/Modules/MSource/Z2.hpp>
#include <Grid/Hadrons/Modules/MUtilities/TestSeqConserved.hpp>
#include <Grid/Hadrons/Modules/MUtilities/TestSeqGamma.hpp>
