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
#include <Grid.h>

namespace Grid { 
  namespace QCD {

/////////////////////////////////////////////////////////////////
// Public interface
/////////////////////////////////////////////////////////////////
GridCartesian *SpaceTimeGrid::makeFourDimGrid(const std::vector<int> & latt,const std::vector<int> &simd,const std::vector<int> &mpi)
{
  return new GridCartesian(latt,simd,mpi); 
}
GridRedBlackCartesian *SpaceTimeGrid::makeFourDimRedBlackGrid(const GridCartesian *FourDimGrid)
{
  return new GridRedBlackCartesian(FourDimGrid); 
}
GridCartesian *SpaceTimeGrid::makeFourDimDWFGrid(const std::vector<int> & latt,const std::vector<int> &mpi)
{
  std::vector<int> simd(4,1);
  return makeFourDimGrid(latt,simd,mpi);
}
GridCartesian         *SpaceTimeGrid::makeFiveDimGrid(int Ls,const GridCartesian *FourDimGrid)
{
  int N4=FourDimGrid->_ndimension;

  std::vector<int> latt5(1,Ls);
  std::vector<int> simd5(1,1);
  std::vector<int>  mpi5(1,1);
  
  for(int d=0;d<N4;d++){
    latt5.push_back(FourDimGrid->_fdimensions[d]);
    simd5.push_back(FourDimGrid->_simd_layout[d]);
     mpi5.push_back(FourDimGrid->_processors[d]);
  }
  return new GridCartesian(latt5,simd5,mpi5); 
}


GridRedBlackCartesian *SpaceTimeGrid::makeFiveDimRedBlackGrid(int Ls,const GridCartesian *FourDimGrid)
{
  int N4=FourDimGrid->_ndimension;
  int cbd=1;
  std::vector<int> latt5(1,Ls);
  std::vector<int> simd5(1,1);
  std::vector<int>  mpi5(1,1);
  std::vector<int>   cb5(1,0);
    
  for(int d=0;d<N4;d++){
    latt5.push_back(FourDimGrid->_fdimensions[d]);
    simd5.push_back(FourDimGrid->_simd_layout[d]);
     mpi5.push_back(FourDimGrid->_processors[d]);
      cb5.push_back(  1);
    }
  return new GridRedBlackCartesian(latt5,simd5,mpi5,cb5,cbd); 
}


GridCartesian         *SpaceTimeGrid::makeFiveDimDWFGrid(int Ls,const GridCartesian *FourDimGrid)
{
  int N4=FourDimGrid->_ndimension;
  int nsimd = FourDimGrid->Nsimd();

  std::vector<int> latt5(1,Ls);
  std::vector<int> simd5(1,nsimd);
  std::vector<int>  mpi5(1,1);
  
  for(int d=0;d<N4;d++){
    latt5.push_back(FourDimGrid->_fdimensions[d]);
    simd5.push_back(1);
     mpi5.push_back(FourDimGrid->_processors[d]);
  }
  return new GridCartesian(latt5,simd5,mpi5); 
}

GridRedBlackCartesian *SpaceTimeGrid::makeFiveDimDWFRedBlackGrid(int Ls,const GridCartesian *FourDimGrid)
{
  int N4=FourDimGrid->_ndimension;
  int nsimd = FourDimGrid->Nsimd();
  int cbd=0;
  std::vector<int> latt5(1,Ls);
  std::vector<int> simd5(1,nsimd);
  std::vector<int>  mpi5(1,1);
  std::vector<int>   cb5(1,1);
    
  for(int d=0;d<N4;d++){
    latt5.push_back(FourDimGrid->_fdimensions[d]);
    simd5.push_back(1);
     mpi5.push_back(FourDimGrid->_processors[d]);
      cb5.push_back(1);
    }
  return new GridRedBlackCartesian(latt5,simd5,mpi5,cb5,cbd); 
}


}}
