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

#include <serialisation/MacroMagic.h>
#include <serialisation/BaseIO.h>
#include <stdint.h>

//////////////////////////////////////////
// Todo:
//////////////////////////////////////////
#include <serialisation/BinaryIO.h>
#include <serialisation/TextIO.h>
//#include <serialisation/JsonIO.h>
//#include <serialisation/YamlIO.h>
#include <serialisation/XmlIO.h>

//////////////////////////////////////////
// Select the default serialiser use ifdef's
//////////////////////////////////////////
namespace Grid {
  typedef XmlReader DefaultReader;
  typedef XmlWriter DefaultWriter;
}
#endif
