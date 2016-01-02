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
#include <serialisation/Serialisation.h>
#include <Config.h>
#include <Timer.h>
#include <Log.h>
#include <AlignedAllocator.h>
#include <Simd.h>
#include <Threads.h>
#include <Communicator.h> 
#include <Cartesian.h>    
#include <Tensors.h>      
#include <Lattice.h>      
#include <Cshift.h>       
#include <Stencil.h>      
#include <Algorithms.h>   
#include <parallelIO/BinaryIO.h>
#include <qcd/QCD.h>
#include <parallelIO/NerscIO.h>
#include <Init.h>

#include <qcd/hmc/NerscCheckpointer.h>
#include <qcd/hmc/HmcRunner.h>



#endif
