/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/Cshift.h

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

    See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef _GRID_CSHIFT_H_
#define _GRID_CSHIFT_H_

#include <Grid/cshift/Cshift_common.h>

#ifdef GRID_COMMS_NONE
#include <Grid/cshift/Cshift_none.h>
#endif

#ifdef GRID_COMMS_MPI
#include <Grid/cshift/Cshift_mpi.h>
#endif 

#ifdef GRID_COMMS_MPI3
#include <Grid/cshift/Cshift_mpi.h>
#endif 

#ifdef GRID_COMMS_MPIT
#include <Grid/cshift/Cshift_mpi.h>
#endif 

#ifdef GRID_COMMS_SHMEM
#include <Grid/cshift/Cshift_mpi.h> // uses same implementation of communicator
#endif 
#endif
