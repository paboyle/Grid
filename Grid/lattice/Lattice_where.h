    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_where.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>

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
#ifndef GRID_LATTICE_WHERE_H
#define GRID_LATTICE_WHERE_H
namespace Grid {
// Must implement the predicate gating the 
// Must be able to reduce the predicate down to a single vInteger per site.
// Must be able to require the type be iScalar x iScalar x ....
//                              give a GetVtype method in iScalar
//                              and blow away the tensor structures.
//
template<class vobj,class iobj>
inline void whereWolf(Lattice<vobj> &ret,const Lattice<iobj> &predicate,Lattice<vobj> &iftrue,Lattice<vobj> &iffalse)
{
  conformable(iftrue,iffalse);
  conformable(iftrue,predicate);
  conformable(iftrue,ret);

  GridBase *grid=iftrue.Grid();

  typedef typename vobj::scalar_object scalar_object;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  typedef typename iobj::vector_type mask_type;

  const int Nsimd = grid->Nsimd();

  autoView(iftrue_v,iftrue,CpuRead);
  autoView(iffalse_v,iffalse,CpuRead);
  autoView(predicate_v,predicate,CpuRead);
  autoView(ret_v,ret,CpuWrite);
  Integer NN= grid->oSites();
  thread_for(ss,NN,{
    Integer mask;
    scalar_object trueval;
    scalar_object falseval;
    for(int l=0;l<Nsimd;l++){
      trueval =extractLane(l,iftrue_v[ss]);
      falseval=extractLane(l,iffalse_v[ss]);
      mask    =extractLane(l,predicate_v[ss]);
      if (mask) falseval=trueval;
      insertLane(l,ret_v[ss],falseval);
    }
  });
}

template<class vobj,class iobj>
inline Lattice<vobj> whereWolf(const Lattice<iobj> &predicate,Lattice<vobj> &iftrue,Lattice<vobj> &iffalse)
{
  conformable(iftrue,iffalse);
  conformable(iftrue,predicate);

  Lattice<vobj> ret(iftrue.Grid());

  whereWolf(ret,predicate,iftrue,iffalse);

  return ret;
}
}
#endif
