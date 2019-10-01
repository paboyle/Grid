/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/serialisation/Serialisation.h

    Copyright (C) 2015

Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef GRID_SERIALISATION_READER_H
#define GRID_SERIALISATION_READER_H

#include <stdint.h>

#include "MacroMagic.h"
#include "BaseIO.h"
#include "BinaryIO.h"
#include "TextIO.h"
#include "XmlIO.h"
#ifndef GRID_NVCC
#include "JSON_IO.h"
#endif

#ifdef HAVE_HDF5
#include "Hdf5IO.h"
#endif

//////////////////////////////////////////
// Todo:
//////////////////////////////////////////
//#include "YamlIO.h"

//////////////////////////////////////////
// Select the default serialiser use ifdef's
//////////////////////////////////////////
NAMESPACE_BEGIN(Grid);
typedef XmlReader DefaultReader;
typedef XmlWriter DefaultWriter;
NAMESPACE_END(Grid);
#endif
