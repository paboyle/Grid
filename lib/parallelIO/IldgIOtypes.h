/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/parallelIO/IldgIO.h

Copyright (C) 2015

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
#ifndef GRID_ILDGTYPES_IO_H
#define GRID_ILDGTYPES_IO_H

#ifdef HAVE_LIME
extern "C" { // for linkage
#include "lime.h"
}

namespace Grid {

#define GRID_FORMAT      "grid-format"
#define ILDG_FORMAT      "ildg-format"
#define ILDG_BINARY_DATA "ildg-binary-data"
#define ILDG_DATA_LFN    "ildg-data-lfn"
#define USQCD_INFO       "usqcdInfo"
#define SCIDAC_CHECKSUM  "scidac-checksum"

/////////////////////////////////////////////////////////////////////////////////
// Data representation of records that enter ILDG and SciDac formats
/////////////////////////////////////////////////////////////////////////////////
struct ildgFormat : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(ildgFormat,
				  double, version,
				  std::string, field,
				  int, precision,
				  int, lx,
				  int, ly,
				  int, lz,
				  int, lt);
  ildgFormat() { 
    version=1.0; 
  };
};
struct usqcdInfo : Serializable { 
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(usqcdInfo,
				  double, version,
				  double, plaq,
				  double, linktr,
				  std::string, info);
  usqcdInfo() { 
    version=1.0; 
  };
};

struct usqcdPropFile : Serializable { 
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(usqcdPropFile,
				  double, version,
				  std::string, type,
				  std::string, info);
  usqcdPropFile() { 
    version=1.0; 
  };
};
struct usqcdSourceInfo : Serializable { 
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(usqcdSourceInfo,
				  double, version,
				  std::string, info);
  usqcdSourceInfo() { 
    version=1.0; 
  };
};
struct usqcdPropInfo : Serializable { 
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(usqcdPropInfo,
				  double, version,
				  int, spin,
				  int, color,
				  std::string, info);
  usqcdPropInfo() { 
    version=1.0; 
  };
};
struct scidacChecksum : Serializable { 
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(scidacChecksum,
				  double, version,
				  uint32_t, suma,
				  uint32_t, sumb);
  scidacChecksum() { 
    version=1.0; 
    suma=sumb=0;
  };
};
}
#endif
#endif
