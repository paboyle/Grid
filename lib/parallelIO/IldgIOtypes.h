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

extern "C" { // for linkage
#include "lime.h"
}

namespace Grid {

  struct ILDGtype{
    bool is_ILDG;
    LimeWriter* LW;
    LimeReader* LR;
    
    ILDGtype(bool is, LimeWriter* L):is_ILDG(is),LW(L),LR(NULL){}
    ILDGtype(bool is, LimeReader* L):is_ILDG(is),LW(NULL),LR(L){}
    ILDGtype():is_ILDG(false),LW(NULL),LR(NULL){}
  };
  



	class ILDGField {
	public:
    // header strings (not in order)
		std::vector<int> dimension;
		std::vector<std::string> boundary; 
		int data_start;
		std::string hdr_version;
		std::string storage_format;
    // Checks on data
		double link_trace;
		double plaquette;
		uint32_t checksum;
		unsigned int sequence_number;
		std::string data_type;
		std::string ensemble_id ;
		std::string ensemble_label ;
		std::string creator ;
		std::string creator_hardware ;
		std::string creation_date ;
		std::string archive_date ;
		std::string floating_point;
	};




}
#endif
