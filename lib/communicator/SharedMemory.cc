/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/communicator/SharedMemory.cc

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

#include <Grid/GridCore.h>

namespace Grid { 

// static data

uint64_t            GlobalSharedMemory::MAX_MPI_SHM_BYTES   = 1024LL*1024LL*1024LL; 
int                 GlobalSharedMemory::Hugepages = 0;
int                 GlobalSharedMemory::ShmSetup;

std::vector<void *> GlobalSharedMemory::WorldShmCommBufs;

Grid_MPI_Comm       GlobalSharedMemory::WorldShmComm;
int                 GlobalSharedMemory::WorldShmRank;
int                 GlobalSharedMemory::WorldShmSize;
std::vector<int>    GlobalSharedMemory::WorldShmRanks;

Grid_MPI_Comm       GlobalSharedMemory::WorldComm;
int                 GlobalSharedMemory::WorldSize;
int                 GlobalSharedMemory::WorldRank;

int                 GlobalSharedMemory::WorldNodes;
int                 GlobalSharedMemory::WorldNode;


}
