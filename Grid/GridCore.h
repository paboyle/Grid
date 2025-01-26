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

#ifndef GRID_BASE_H
#define GRID_BASE_H


#include <Grid/DisableWarnings.h>
#include <Grid/Namespace.h>
#include <Grid/GridStd.h>
#include <Grid/threads/Pragmas.h>
#include <Grid/perfmon/Timer.h>
//#include <Grid/perfmon/PerfCount.h>
#include <Grid/util/Util.h>
#include <Grid/log/Log.h>
#include <Grid/perfmon/Tracing.h>
#include <Grid/allocator/Allocator.h>
#include <Grid/simd/Simd.h>
#include <Grid/threads/ThreadReduction.h>
#include <Grid/serialisation/Serialisation.h>
#include <Grid/util/Sha.h>
#include <Grid/communicator/Communicator.h> 
#include <Grid/cartesian/Cartesian.h>    
#include <Grid/tensors/Tensors.h>      
#include <Grid/lattice/Lattice.h>      
#include <Grid/cshift/Cshift.h>       
#include <Grid/stencil/Stencil.h>      
#include <Grid/stencil/GeneralLocalStencil.h>      
#include <Grid/parallelIO/BinaryIO.h>
#include <Grid/algorithms/Algorithms.h>   
NAMESPACE_CHECK(GridCore)

#endif
