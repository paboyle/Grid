    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/utils/LinalgUtils.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
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
#ifndef GRID_QCD_LINALG_UTILS_H
#define GRID_QCD_LINALG_UTILS_H

namespace Grid{
namespace QCD{
////////////////////////////////////////////////////////////////////////
//This file brings additional linear combination assist that is helpful
//to QCD such as chiral projectors and spin matrices applied to one of the inputs.
//These routines support five-D chiral fermions and contain s-subslice indexing 
//on the 5d (rb4d) checkerboarded lattices
////////////////////////////////////////////////////////////////////////

template<class vobj> 
void axpibg5x(Lattice<vobj> &z,const Lattice<vobj> &x,RealD a,RealD b)
{
  z.checkerboard = x.checkerboard;
  conformable(x,z);

  GridBase *grid=x._grid;

  Gamma G5(Gamma::Gamma5);
PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss++){
    vobj tmp;
    tmp = a*x._odata[ss];
    tmp = tmp + G5*(b*timesI(x._odata[ss]));
    vstream(z._odata[ss],tmp);
  }
}

template<class vobj> 
void axpby_ssp(Lattice<vobj> &z, RealD a,const Lattice<vobj> &x,RealD b,const Lattice<vobj> &y,int s,int sp)
{
  z.checkerboard = x.checkerboard;
  conformable(x,y);
  conformable(x,z);
  GridBase *grid=x._grid;
  int Ls = grid->_rdimensions[0];
PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss+=Ls){ // adds Ls
    vobj tmp = a*x._odata[ss+s]+b*y._odata[ss+sp];
    vstream(z._odata[ss+s],tmp);
  }
}

template<class vobj> 
void ag5xpby_ssp(Lattice<vobj> &z,RealD a,const Lattice<vobj> &x,RealD b,const Lattice<vobj> &y,int s,int sp)
{
  z.checkerboard = x.checkerboard;
  conformable(x,y);
  conformable(x,z);
  GridBase *grid=x._grid;
  int Ls = grid->_rdimensions[0];
  Gamma G5(Gamma::Gamma5);
PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss+=Ls){ // adds Ls
    vobj tmp;
    tmp = G5*x._odata[ss+s]*a;
    tmp = tmp + b*y._odata[ss+sp];
    vstream(z._odata[ss+s],tmp);
  }
}

template<class vobj> 
void axpbg5y_ssp(Lattice<vobj> &z,RealD a,const Lattice<vobj> &x,RealD b,const Lattice<vobj> &y,int s,int sp)
{
  z.checkerboard = x.checkerboard;
  conformable(x,y);
  conformable(x,z);
  GridBase *grid=x._grid;
  int Ls = grid->_rdimensions[0];
  Gamma G5(Gamma::Gamma5);
PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss+=Ls){ // adds Ls
    vobj tmp;
    tmp = G5*y._odata[ss+sp]*b;
    tmp = tmp + a*x._odata[ss+s];
    vstream(z._odata[ss+s],tmp);
  }
}

template<class vobj> 
void ag5xpbg5y_ssp(Lattice<vobj> &z,RealD a,const Lattice<vobj> &x,RealD b,const Lattice<vobj> &y,int s,int sp)
{
  z.checkerboard = x.checkerboard;
  conformable(x,y);
  conformable(x,z);
  GridBase *grid=x._grid;
  int Ls = grid->_rdimensions[0];
  Gamma G5(Gamma::Gamma5);
PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss+=Ls){ // adds Ls
    vobj tmp1;
    vobj tmp2;
    tmp1 = a*x._odata[ss+s]+b*y._odata[ss+sp];
    tmp2 = G5*tmp1;
    vstream(z._odata[ss+s],tmp2);
  }
}

template<class vobj> 
void axpby_ssp_pminus(Lattice<vobj> &z,RealD a,const Lattice<vobj> &x,RealD b,const Lattice<vobj> &y,int s,int sp)
{
  z.checkerboard = x.checkerboard;
  conformable(x,y);
  conformable(x,z);
  GridBase *grid=x._grid;
  int Ls = grid->_rdimensions[0];
PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss+=Ls){ // adds Ls
    vobj tmp;
    spProj5m(tmp,y._odata[ss+sp]);
    tmp = a*x._odata[ss+s]+b*tmp;
    vstream(z._odata[ss+s],tmp);
  }
}

template<class vobj> 
void axpby_ssp_pplus(Lattice<vobj> &z,RealD a,const Lattice<vobj> &x,RealD b,const Lattice<vobj> &y,int s,int sp)
{
  z.checkerboard = x.checkerboard;
  conformable(x,y);
  conformable(x,z);
  GridBase *grid=x._grid;
  int Ls = grid->_rdimensions[0];
PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss+=Ls){ // adds Ls
    vobj tmp;
    spProj5p(tmp,y._odata[ss+sp]);
    tmp = a*x._odata[ss+s]+b*tmp;
    vstream(z._odata[ss+s],tmp);
  }
}

template<class vobj> 
void G5R5(Lattice<vobj> &z,const Lattice<vobj> &x)
{
  GridBase *grid=x._grid;
  z.checkerboard = x.checkerboard;
  conformable(x,z);
  int Ls = grid->_rdimensions[0];
  Gamma G5(Gamma::Gamma5);
PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss+=Ls){ // adds Ls
    vobj tmp;
    for(int s=0;s<Ls;s++){
      int sp = Ls-1-s;
      tmp = G5*x._odata[ss+s];
      vstream(z._odata[ss+sp],tmp);
    }
  }
}

}}
#endif 
