    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/utils/SpaceTimeGrid.h

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
#ifndef GRID_QCD_SPACE_TIME_GRID_H
#define GRID_QCD_SPACE_TIME_GRID_H
namespace Grid {
namespace QCD {

class SpaceTimeGrid {
 public:

  static GridCartesian         *makeFourDimGrid(const std::vector<int> & latt,const std::vector<int> &simd,const std::vector<int> &mpi);
  static GridRedBlackCartesian *makeFourDimRedBlackGrid       (const GridCartesian *FourDimGrid);

  static GridCartesian         *makeFiveDimGrid        (int Ls,const GridCartesian *FourDimGrid);
  static GridRedBlackCartesian *makeFiveDimRedBlackGrid(int Ls,const GridCartesian *FourDimGrid);

  static GridCartesian         *makeFiveDimDWFGrid        (int Ls,const GridCartesian *FourDimGrid);
  static GridRedBlackCartesian *makeFiveDimDWFRedBlackGrid(int Ls,const GridCartesian *FourDimGrid);
  static GridCartesian         *makeFourDimDWFGrid        (const std::vector<int> & latt,const std::vector<int> &mpi);

};

}}

#endif
