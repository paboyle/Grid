   /*************************************************************************************

     Grid physics library, www.github.com/paboyle/Grid 

     Source file: ./lib/Stencil.cc

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

void Gather_plane_table_compute (GridBase *grid,int dimension,int plane,int cbmask,
					int off,std::vector<std::pair<int,int> > & table)
{
  table.resize(0);

  if ( !grid->CheckerBoarded(dimension) ) {
    cbmask = 0x3;
  }
  int rd = grid->_rdimensions[dimension];
  int so= plane*grid->_ostride[dimension]; // base offset for start of plane 
  int e1=grid->_slice_nblock[dimension];
  int e2=grid->_slice_block[dimension];
  int stride=grid->_slice_stride[dimension];
  if ( cbmask == 0x3 ) { 
    table.resize(e1*e2);
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int o  = n*stride;
	int bo = n*e2;
	table[bo+b]=std::pair<int,int>(bo+b,o+b);
      }
    }
  } else { 
     int bo=0;
     table.resize(e1*e2/2);
     for(int n=0;n<e1;n++){
       for(int b=0;b<e2;b++){
	 int o  = n*stride;
	 int ocb=1<<grid->CheckerBoardFromOindexTable(o+b);
	 if ( ocb &cbmask ) {
	   table[bo]=std::pair<int,int>(bo,o+b); bo++;
	 }
       }
     }
  }
}

}
