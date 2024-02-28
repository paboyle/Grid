/*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/PaddedCell.h

    Copyright (C) 2019

Author: Peter Boyle pboyle@bnl.gov

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
#pragma once

#include<Grid/cshift/Cshift.h>

NAMESPACE_BEGIN(Grid);

//Allow the user to specify how the C-shift is performed, e.g. to respect the appropriate boundary conditions
template<typename vobj>
struct CshiftImplBase{
  virtual Lattice<vobj> Cshift(const Lattice<vobj> &in, int dir, int shift) const = 0;
  virtual ~CshiftImplBase(){}
};
template<typename vobj>
struct CshiftImplDefault: public CshiftImplBase<vobj>{
  Lattice<vobj> Cshift(const Lattice<vobj> &in, int dir, int shift) const override{ return Grid::Cshift(in,dir,shift); }
};
template<typename Gimpl>
struct CshiftImplGauge: public CshiftImplBase<typename Gimpl::GaugeLinkField::vector_object>{
  typename Gimpl::GaugeLinkField Cshift(const typename Gimpl::GaugeLinkField &in, int dir, int shift) const override{ return Gimpl::CshiftLink(in,dir,shift); }
};  

class PaddedCell {
public:
  GridCartesian * unpadded_grid;
  int dims;
  int depth;
  std::vector<GridCartesian *> grids;

  ~PaddedCell()
  {
    DeleteGrids();
  }
  PaddedCell(int _depth,GridCartesian *_grid)
  {
    unpadded_grid = _grid;
    depth=_depth;
    dims=_grid->Nd();
    AllocateGrids();
    Coordinate local     =unpadded_grid->LocalDimensions();
    for(int d=0;d<dims;d++){
      assert(local[d]>=depth);
    }
  }
  void DeleteGrids(void)
  {
    for(int d=0;d<grids.size();d++){
      delete grids[d];
    }
    grids.resize(0);
  };
  void AllocateGrids(void)
  {
    Coordinate local     =unpadded_grid->LocalDimensions();
    Coordinate simd      =unpadded_grid->_simd_layout;
    Coordinate processors=unpadded_grid->_processors;
    Coordinate plocal    =unpadded_grid->LocalDimensions();
    Coordinate global(dims);

    // expand up one dim at a time
    for(int d=0;d<dims;d++){

      plocal[d] += 2*depth; 

      for(int d=0;d<dims;d++){
	global[d] = plocal[d]*processors[d];
      }

      grids.push_back(new GridCartesian(global,simd,processors));
    }
  };
  template<class vobj>
  inline Lattice<vobj> Extract(const Lattice<vobj> &in) const
  {
    Lattice<vobj> out(unpadded_grid);

    Coordinate local     =unpadded_grid->LocalDimensions();
    Coordinate fll(dims,depth); // depends on the MPI spread
    Coordinate tll(dims,0); // depends on the MPI spread
    localCopyRegion(in,out,fll,tll,local);
    return out;
  }
  template<class vobj>
  inline Lattice<vobj> Exchange(const Lattice<vobj> &in, const CshiftImplBase<vobj> &cshift = CshiftImplDefault<vobj>()) const
  {
    GridBase *old_grid = in.Grid();
    int dims = old_grid->Nd();
    Lattice<vobj> tmp = in;
    for(int d=0;d<dims;d++){
      tmp = Expand(d,tmp,cshift); // rvalue && assignment
    }
    return tmp;
  }
  // expand up one dim at a time
  template<class vobj>
  inline Lattice<vobj> Expand(int dim, const Lattice<vobj> &in, const CshiftImplBase<vobj> &cshift = CshiftImplDefault<vobj>()) const
  {
    GridBase *old_grid = in.Grid();
    GridCartesian *new_grid = grids[dim];//These are new grids
    Lattice<vobj>  padded(new_grid);
    Lattice<vobj> shifted(old_grid);    
    Coordinate local     =old_grid->LocalDimensions();
    Coordinate plocal    =new_grid->LocalDimensions();
    if(dim==0) conformable(old_grid,unpadded_grid);
    else       conformable(old_grid,grids[dim-1]);

    std::cout << " dim "<<dim<<" local "<<local << " padding to "<<plocal<<std::endl;

    double tins=0, tshift=0;
    
    // Middle bit
    double t = usecond();
    for(int x=0;x<local[dim];x++){
      InsertSliceLocal(in,padded,x,depth+x,dim);
    }
    tins += usecond() - t;
    
    // High bit
    t = usecond();
    shifted = cshift.Cshift(in,dim,depth);
    tshift += usecond() - t;

    t=usecond();
    for(int x=0;x<depth;x++){
      InsertSliceLocal(shifted,padded,local[dim]-depth+x,depth+local[dim]+x,dim);
    }
    tins += usecond() - t;
    
    // Low bit
    t = usecond();
    shifted = cshift.Cshift(in,dim,-depth);
    tshift += usecond() - t;
    
    t = usecond();
    for(int x=0;x<depth;x++){
      InsertSliceLocal(shifted,padded,x,x,dim);
    }
    tins += usecond() - t;

    std::cout << GridLogPerformance << "PaddedCell::Expand timings: cshift:" << tshift/1000 << "ms, insert-slice:" << tins/1000 << "ms" << std::endl;
    
    return padded;
  }

};
 

NAMESPACE_END(Grid);

