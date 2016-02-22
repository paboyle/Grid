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

  GridBase *grid=iftrue._grid;

  typedef typename vobj::scalar_object scalar_object;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  typedef typename iobj::vector_type mask_type;

  const int Nsimd = grid->Nsimd();

  std::vector<Integer> mask(Nsimd);
  std::vector<scalar_object> truevals (Nsimd);
  std::vector<scalar_object> falsevals(Nsimd);

PARALLEL_FOR_LOOP
  for(int ss=0;ss<iftrue._grid->oSites(); ss++){

    extract(iftrue._odata[ss]   ,truevals);
    extract(iffalse._odata[ss]  ,falsevals);
    extract<vInteger,Integer>(TensorRemove(predicate._odata[ss]),mask);

    for(int s=0;s<Nsimd;s++){
      if (mask[s]) falsevals[s]=truevals[s];
    }

    merge(ret._odata[ss],falsevals);
  }
}

template<class vobj,class iobj>
inline Lattice<vobj> whereWolf(const Lattice<iobj> &predicate,Lattice<vobj> &iftrue,Lattice<vobj> &iffalse)
{
  conformable(iftrue,iffalse);
  conformable(iftrue,predicate);

  Lattice<vobj> ret(iftrue._grid);

  where(ret,predicate,iftrue,iffalse);

  return ret;
}
}
#endif
