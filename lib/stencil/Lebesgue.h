    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/stencil/Lebesgue.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#ifndef GRID_LEBESGUE_H
#define GRID_LEBESGUE_H

#include<vector>

// Lebesgue, Morton, Z-graph ordering assistance
namespace Grid {
  
  class LebesgueOrder { 
  public:

    typedef int32_t IndexInteger;
    
    static int UseLebesgueOrder;
    GridBase *grid;

  public:
    LebesgueOrder(GridBase *_grid);

    inline IndexInteger Reorder(IndexInteger ss) { 
      return _LebesgueReorder[ss] ;
    };

    ////////////////////////////
    // Space filling fractal for cache oblivious
    ////////////////////////////
    void ZGraph(void);
    IndexInteger alignup(IndexInteger n);

    /////////////////////////////////
    // Cartesian stencil blocking strategy
    /////////////////////////////////
    static std::vector<int> Block;
    void CartesianBlocking(void);
    void IterateO(int ND,int dim,
		  std::vector<IndexInteger> & xo,
		  std::vector<IndexInteger> & xi,
		  std::vector<IndexInteger> &dims);
    void IterateI(int ND,int dim,
		  std::vector<IndexInteger> & xo,
		  std::vector<IndexInteger> & xi,
		  std::vector<IndexInteger> &dims);

  private:
    std::vector<IndexInteger> _LebesgueReorder;

  };    
}
#endif
