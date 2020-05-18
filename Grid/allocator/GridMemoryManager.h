/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/GridMemoryManager.h

    Copyright (C) 2020

Author: Christoph Lehner <christoph@lhnr.de>

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
#ifndef GRID_MEMORY_MANAGER_H
#define GRID_MEMORY_MANAGER_H

NAMESPACE_BEGIN(Grid);

void gridMemoryInit();
void gridMallocManaged(void** pp, size_t sz);
void gridMoveToHost(void** pp);
void gridAcceleratorPrefetch(void* p, size_t sz);
void gridMemGetInfo(size_t* pfree, size_t* ptotal);
void gridFree(void* p);

NAMESPACE_END(Grid);

#endif
