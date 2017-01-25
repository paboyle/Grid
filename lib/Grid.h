    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/Grid.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: azusayamaguchi <ayamaguc@YAMAKAZE.local>
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
//
//  Grid.h
//  simd
//
//  Created by Peter Boyle on 09/05/2014.
//  Copyright (c) 2014 University of Edinburgh. All rights reserved.
//

#ifndef GRID_H
#define GRID_H

///////////////////
// Std C++ dependencies
///////////////////
#include <cassert>
#include <complex>
#include <vector>
#include <iostream>
#include <iomanip>
#include <random>
#include <functional>
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <ctime>
#include <sys/time.h>
#include <chrono>

///////////////////
// Grid headers
///////////////////
#include "Config.h"
#include <Grid/Timer.h>
#include <Grid/PerfCount.h>
#include <Grid/Log.h>
#include <Grid/AlignedAllocator.h>
#include <Grid/Simd.h>
#include <Grid/serialisation/Serialisation.h>
#include <Grid/Threads.h>
#include <Grid/Lexicographic.h>
#include <Grid/Init.h>
#include <Grid/Communicator.h> 
#include <Grid/Cartesian.h>    
#include <Grid/Tensors.h>      
#include <Grid/Lattice.h>      
#include <Grid/Cshift.h>       
#include <Grid/Stencil.h>      
#include <Grid/Algorithms.h>   

#include <Grid/FFT.h>

#include <Grid/qcd/QCD.h>
#include <Grid/parallelIO/IldgIOtypes.h>
#include <Grid/parallelIO/BinaryIO.h>
#include <Grid/parallelIO/IldgIO.h>
#include <Grid/parallelIO/NerscIO.h>

#include <Grid/qcd/hmc/checkpointers/CheckPointers.h>
#include <Grid/qcd/hmc/HMCModules.h>
#include <Grid/qcd/modules/mods.h>
#include <Grid/qcd/hmc/HMCResourceManager.h>
#include <Grid/qcd/hmc/HmcRunner.h>
#include <Grid/qcd/hmc/GenericHMCrunner.h>
#include <Grid/qcd/hmc/HMCRunnerModule.h>


#endif
