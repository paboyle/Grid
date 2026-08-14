/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/DeviceMemoryAllocator.h

    Copyright (C) 2025

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

NAMESPACE_BEGIN(Grid);

#ifdef GRID_DEVICE_MEMORY_ALLOCATOR
void *acceleratorAllocDevice(size_t bytes);
void acceleratorFreeDevice(void *ptr);
#endif

NAMESPACE_END(Grid);
