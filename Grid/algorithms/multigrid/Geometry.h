/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/GeneralCoarsenedMatrix.h

    Copyright (C) 2015

Author: Peter Boyle <pboyle@bnl.gov>

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

NAMESPACE_BEGIN(Grid);


/////////////////////////////////////////////////////////////////
// Geometry class in cartesian case
/////////////////////////////////////////////////////////////////

class Geometry {
public:
  int npoint;
  int base;
  std::vector<int> directions   ;
  std::vector<int> displacements;
  std::vector<int> points_dagger;

  Geometry(int _d)  {
    
    base = (_d==5) ? 1:0;

    // make coarse grid stencil for 4d , not 5d
    if ( _d==5 ) _d=4;

    npoint = 2*_d+1;
    directions.resize(npoint);
    displacements.resize(npoint);
    points_dagger.resize(npoint);
    for(int d=0;d<_d;d++){
      directions[d   ] = d+base;
      directions[d+_d] = d+base;
      displacements[d  ] = +1;
      displacements[d+_d]= -1;
      points_dagger[d   ] = d+_d;
      points_dagger[d+_d] = d;
    }
    directions   [2*_d]=0;
    displacements[2*_d]=0;
    points_dagger[2*_d]=2*_d;
  }

  int point(int dir, int disp) {
    assert(disp == -1 || disp == 0 || disp == 1);
    assert(base+0 <= dir && dir < base+4);

    // directions faster index = new indexing
    // 4d (base = 0):
    // point 0  1  2  3  4  5  6  7  8
    // dir   0  1  2  3  0  1  2  3  0
    // disp +1 +1 +1 +1 -1 -1 -1 -1  0
    // 5d (base = 1):
    // point 0  1  2  3  4  5  6  7  8
    // dir   1  2  3  4  1  2  3  4  0
    // disp +1 +1 +1 +1 -1 -1 -1 -1  0

    // displacements faster index = old indexing
    // 4d (base = 0):
    // point 0  1  2  3  4  5  6  7  8
    // dir   0  0  1  1  2  2  3  3  0
    // disp +1 -1 +1 -1 +1 -1 +1 -1  0
    // 5d (base = 1):
    // point 0  1  2  3  4  5  6  7  8
    // dir   1  1  2  2  3  3  4  4  0
    // disp +1 -1 +1 -1 +1 -1 +1 -1  0

    if(dir == 0 and disp == 0)
      return 8;
    else // New indexing
      return (1 - disp) / 2 * 4 + dir - base;
    // else // Old indexing
    //   return (4 * (dir - base) + 1 - disp) / 2;
  }
};

/////////////////////////////////////////////////////////////////
// Less local equivalent of Geometry class in cartesian case
/////////////////////////////////////////////////////////////////
class NonLocalStencilGeometry {
public:
  //  int depth;
  int skip;
  int hops;
  int npoint;
  std::vector<Coordinate> shifts;
  Coordinate stencil_size;
  Coordinate stencil_lo;
  Coordinate stencil_hi;
  GridCartesian *grid;
  GridCartesian *Grid() {return grid;};
  int Depth(void){return 1;};   // Ghost zone depth
  int Hops(void){return hops;}; // # of hops=> level of corner fill in in stencil
  int DimSkip(void){return skip;};

  virtual ~NonLocalStencilGeometry() {};

  int  Reverse(int point)
  {
    int Nd = Grid()->Nd();
    Coordinate shft = shifts[point];
    Coordinate rev(Nd);
    for(int mu=0;mu<Nd;mu++) rev[mu]= -shft[mu];
    for(int p=0;p<npoint;p++){
      if(rev==shifts[p]){
	return p;
      }
    }
    assert(0);
    return -1;
  }
  void BuildShifts(void)
  {
    this->shifts.resize(0);
    int Nd = this->grid->Nd();

    int dd = this->DimSkip();
    for(int s0=this->stencil_lo[dd+0];s0<=this->stencil_hi[dd+0];s0++){
    for(int s1=this->stencil_lo[dd+1];s1<=this->stencil_hi[dd+1];s1++){
    for(int s2=this->stencil_lo[dd+2];s2<=this->stencil_hi[dd+2];s2++){
    for(int s3=this->stencil_lo[dd+3];s3<=this->stencil_hi[dd+3];s3++){
      Coordinate sft(Nd,0);
      sft[dd+0] = s0;
      sft[dd+1] = s1;
      sft[dd+2] = s2;
      sft[dd+3] = s3;
      int nhops = abs(s0)+abs(s1)+abs(s2)+abs(s3);
      if(nhops<=this->hops) this->shifts.push_back(sft);
    }}}}
    this->npoint = this->shifts.size();
    std::cout << GridLogMessage << "NonLocalStencilGeometry has "<< this->npoint << " terms in stencil "<<std::endl;
  }
  
  NonLocalStencilGeometry(GridCartesian *_coarse_grid,int _hops,int _skip) : grid(_coarse_grid), hops(_hops), skip(_skip)
  {
    Coordinate latt = grid->GlobalDimensions();
    stencil_size.resize(grid->Nd());
    stencil_lo.resize(grid->Nd());
    stencil_hi.resize(grid->Nd());
    for(int d=0;d<grid->Nd();d++){
     if ( latt[d] == 1 ) {
      stencil_lo[d] = 0;
      stencil_hi[d] = 0;
      stencil_size[d]= 1;
     } else if ( latt[d] == 2 ) {
      stencil_lo[d] = -1;
      stencil_hi[d] = 0;
      stencil_size[d]= 2;
     } else if ( latt[d] > 2 ) {
       stencil_lo[d] = -1;
       stencil_hi[d] =  1;
       stencil_size[d]= 3;
     }
    }
    this->BuildShifts();
  };

};

// Need to worry about red-black now
class NonLocalStencilGeometry4D : public NonLocalStencilGeometry {
public:
  virtual int DerivedDimSkip(void) { return 0;};
  NonLocalStencilGeometry4D(GridCartesian *Coarse,int _hops) : NonLocalStencilGeometry(Coarse,_hops,0) { };
  virtual ~NonLocalStencilGeometry4D() {};
};
class NonLocalStencilGeometry5D : public NonLocalStencilGeometry {
public:
  virtual int DerivedDimSkip(void) { return 1; }; 
  NonLocalStencilGeometry5D(GridCartesian *Coarse,int _hops) : NonLocalStencilGeometry(Coarse,_hops,1)  { };
  virtual ~NonLocalStencilGeometry5D() {};
};
/*
 * Bunch of different options classes
 */
class NextToNextToNextToNearestStencilGeometry4D : public NonLocalStencilGeometry4D {
public:
  NextToNextToNextToNearestStencilGeometry4D(GridCartesian *Coarse) :  NonLocalStencilGeometry4D(Coarse,4)
  {
  };
};
class NextToNextToNextToNearestStencilGeometry5D : public  NonLocalStencilGeometry5D {
public:
  NextToNextToNextToNearestStencilGeometry5D(GridCartesian *Coarse) :  NonLocalStencilGeometry5D(Coarse,4)
  {
  };
};
class NextToNearestStencilGeometry4D : public  NonLocalStencilGeometry4D {
public:
  NextToNearestStencilGeometry4D(GridCartesian *Coarse) :  NonLocalStencilGeometry4D(Coarse,2)
  {
  };
};
class NextToNearestStencilGeometry5D : public  NonLocalStencilGeometry5D {
public:
  NextToNearestStencilGeometry5D(GridCartesian *Coarse) :  NonLocalStencilGeometry5D(Coarse,2)
  {
  };
};
class NearestStencilGeometry4D : public  NonLocalStencilGeometry4D {
public:
  NearestStencilGeometry4D(GridCartesian *Coarse) :  NonLocalStencilGeometry4D(Coarse,1)
  {
  };
};
class NearestStencilGeometry5D : public  NonLocalStencilGeometry5D {
public:
  NearestStencilGeometry5D(GridCartesian *Coarse) :  NonLocalStencilGeometry5D(Coarse,1)
  {
  };
};

NAMESPACE_END(Grid);
