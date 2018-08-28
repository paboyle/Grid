/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/HMC.h

Copyright (C) 2015

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
//--------------------------------------------------------------------
//--------------------------------------------------------------------
#ifndef HMC_AGGREGATE_INCLUDED
#define HMC_AGGREGATE_INCLUDED

#include <string>

#include <Grid/qcd/observables/hmc_observable.h>
#include <Grid/qcd/hmc/HMC.h>


// annoying location; should move this ?
#include <Grid/parallelIO/IldgIOtypes.h>
#include <Grid/parallelIO/IldgIO.h>
#include <Grid/parallelIO/NerscIO.h>

#include <Grid/qcd/hmc/checkpointers/CheckPointers.h>
#include <Grid/qcd/hmc/HMCModules.h>
#include <Grid/qcd/modules/mods.h>
#include <Grid/qcd/hmc/HMCResourceManager.h>
#include <Grid/qcd/hmc/GenericHMCrunner.h>
#include <Grid/qcd/hmc/HMCRunnerModule.h>

#endif
