/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/ActionCore.h

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>

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
#ifndef QCD_ACTION_CORE
#define QCD_ACTION_CORE

#include <Grid/qcd/action/ActionBase.h>
#include <Grid/qcd/action/ActionSet.h>
#include <Grid/qcd/action/ActionParams.h>

////////////////////////////////////////////
// Gauge Actions
////////////////////////////////////////////
#include <Grid/qcd/action/gauge/Gauge.h>

////////////////////////////////////////////
// Fermion prereqs
////////////////////////////////////////////
#include <Grid/qcd/action/fermion/FermionCore.h>

////////////////////////////////////////////
// Scalar Actions
////////////////////////////////////////////
#include <Grid/qcd/action/scalar/Scalar.h>

////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////
#include <Grid/qcd/utils/Metric.h>
#include <Grid/qcd/utils/CovariantLaplacian.h>




#endif
