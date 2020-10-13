/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/utils/SpaceTimeGrid.cc

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
#include <Grid/GridQCDcore.h>

NAMESPACE_BEGIN(Grid); 

/////////////////////////////////////////////////////////////////
// Public interface
/////////////////////////////////////////////////////////////////
GridCartesian *SpaceTimeGrid::makeFourDimGrid(const Coordinate & latt,const Coordinate &simd,const Coordinate &mpi)
{
  return new GridCartesian(latt,simd,mpi); 
}
GridRedBlackCartesian *SpaceTimeGrid::makeFourDimRedBlackGrid(const GridCartesian *FourDimGrid)
{
  return new GridRedBlackCartesian(FourDimGrid); 
}
GridCartesian *SpaceTimeGrid::makeFourDimDWFGrid(const Coordinate & latt,const Coordinate &mpi)
{
  Coordinate simd(4,1);
  return makeFourDimGrid(latt,simd,mpi);
}
GridCartesian         *SpaceTimeGrid::makeFiveDimGrid(int Ls,const GridCartesian *FourDimGrid)
{
  int N4=FourDimGrid->_ndimension;

  Coordinate latt5(1,Ls);
  Coordinate simd5(1,1);
  Coordinate  mpi5(1,1);
  
  for(int d=0;d<N4;d++){
    latt5.push_back(FourDimGrid->_fdimensions[d]);
    simd5.push_back(FourDimGrid->_simd_layout[d]);
    mpi5.push_back(FourDimGrid->_processors[d]);
  }
  return new GridCartesian(latt5,simd5,mpi5,*FourDimGrid); 
}


GridRedBlackCartesian *SpaceTimeGrid::makeFiveDimRedBlackGrid(int Ls,const GridCartesian *FourDimGrid)
{
  int N4=FourDimGrid->_ndimension;
  int cbd=1;
  Coordinate   cb5(1,0);
  for(int d=0;d<N4;d++){
    cb5.push_back(  1);
  }
  GridCartesian *tmp = makeFiveDimGrid(Ls,FourDimGrid);
  GridRedBlackCartesian *ret = new GridRedBlackCartesian(tmp,cb5,cbd); 
  delete tmp;
  return ret;
}


GridCartesian         *SpaceTimeGrid::makeFiveDimDWFGrid(int Ls,const GridCartesian *FourDimGrid)
{
  int N4    = FourDimGrid->_ndimension;
  int nsimd = FourDimGrid->Nsimd();

  Coordinate latt5(1,Ls);
  Coordinate simd5(1,nsimd);
  Coordinate  mpi5(1,1);
  
  for(int d=0;d<N4;d++){
    latt5.push_back(FourDimGrid->_fdimensions[d]);
    simd5.push_back(1);
    mpi5.push_back(FourDimGrid->_processors[d]);
  }
  return new GridCartesian(latt5,simd5,mpi5,*FourDimGrid); 
}
///////////////////////////////////////////////////
// Interface is inefficient and forces the deletion
// Pass in the non-redblack grid
///////////////////////////////////////////////////
GridRedBlackCartesian *SpaceTimeGrid::makeFiveDimDWFRedBlackGrid(int Ls,const GridCartesian *FourDimGrid)
{
  int N4=FourDimGrid->_ndimension;
  int cbd=1;
  Coordinate   cb5(1,0);
  for(int d=0;d<N4;d++){
    cb5.push_back(1);
  }
  GridCartesian *tmp         = makeFiveDimDWFGrid(Ls,FourDimGrid);
  GridRedBlackCartesian *ret = new GridRedBlackCartesian(tmp,cb5,cbd); 
  delete tmp;
  return ret;
}

NAMESPACE_END(Grid);
