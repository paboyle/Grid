    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_overload.h

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
#ifndef GRID_LATTICE_OVERLOAD_H
#define GRID_LATTICE_OVERLOAD_H

namespace Grid {

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // unary negation
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class vobj>
  inline Lattice<vobj> operator -(const Lattice<vobj> &r)
  {
    Lattice<vobj> ret(r._grid);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<r._grid->oSites();ss++){
      vstream(ret._odata[ss], -r._odata[ss]);
    }
    return ret;
  } 
  /////////////////////////////////////////////////////////////////////////////////////
  // Lattice BinOp Lattice,
  //NB mult performs conformable check. Do not reapply here for performance.
  /////////////////////////////////////////////////////////////////////////////////////
  template<class left,class right>
    inline auto operator * (const Lattice<left> &lhs,const Lattice<right> &rhs)-> Lattice<decltype(lhs._odata[0]*rhs._odata[0])>
  {
    Lattice<decltype(lhs._odata[0]*rhs._odata[0])> ret(rhs._grid);
    mult(ret,lhs,rhs);
    return ret;
  }
  template<class left,class right>
    inline auto operator + (const Lattice<left> &lhs,const Lattice<right> &rhs)-> Lattice<decltype(lhs._odata[0]+rhs._odata[0])>
  {
    Lattice<decltype(lhs._odata[0]+rhs._odata[0])> ret(rhs._grid);
    add(ret,lhs,rhs);
    return ret;
  }
  template<class left,class right>
    inline auto operator - (const Lattice<left> &lhs,const Lattice<right> &rhs)-> Lattice<decltype(lhs._odata[0]-rhs._odata[0])>
  {
    Lattice<decltype(lhs._odata[0]-rhs._odata[0])> ret(rhs._grid);
    sub(ret,lhs,rhs);
    return ret;
  }
  
  // Scalar BinOp Lattice ;generate return type
  template<class left,class right>
  inline auto operator * (const left &lhs,const Lattice<right> &rhs) -> Lattice<decltype(lhs*rhs._odata[0])>
  {
    Lattice<decltype(lhs*rhs._odata[0])> ret(rhs._grid);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites(); ss++){
      decltype(lhs*rhs._odata[0]) tmp=lhs*rhs._odata[ss]; 
      vstream(ret._odata[ss],tmp);
	   //      ret._odata[ss]=lhs*rhs._odata[ss];
    }
    return ret;
  }
  template<class left,class right>
    inline auto operator + (const left &lhs,const Lattice<right> &rhs) -> Lattice<decltype(lhs+rhs._odata[0])>
    {
      Lattice<decltype(lhs+rhs._odata[0])> ret(rhs._grid);
PARALLEL_FOR_LOOP
      for(int ss=0;ss<rhs._grid->oSites(); ss++){
	decltype(lhs+rhs._odata[0]) tmp =lhs-rhs._odata[ss];  
	vstream(ret._odata[ss],tmp);
	//	ret._odata[ss]=lhs+rhs._odata[ss];
      }
        return ret;
    }
  template<class left,class right>
    inline auto operator - (const left &lhs,const Lattice<right> &rhs) -> Lattice<decltype(lhs-rhs._odata[0])>
  {
    Lattice<decltype(lhs-rhs._odata[0])> ret(rhs._grid);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites(); ss++){
      decltype(lhs-rhs._odata[0]) tmp=lhs-rhs._odata[ss];  
      vstream(ret._odata[ss],tmp);
      //      ret._odata[ss]=lhs-rhs._odata[ss];
    }
    return ret;
  }
    template<class left,class right>
      inline auto operator * (const Lattice<left> &lhs,const right &rhs) -> Lattice<decltype(lhs._odata[0]*rhs)>
    {
      Lattice<decltype(lhs._odata[0]*rhs)> ret(lhs._grid);
PARALLEL_FOR_LOOP
      for(int ss=0;ss<lhs._grid->oSites(); ss++){
	decltype(lhs._odata[0]*rhs) tmp =lhs._odata[ss]*rhs;
	vstream(ret._odata[ss],tmp);
	//            ret._odata[ss]=lhs._odata[ss]*rhs;
      }
      return ret;
    }
    template<class left,class right>
      inline auto operator + (const Lattice<left> &lhs,const right &rhs) -> Lattice<decltype(lhs._odata[0]+rhs)>
    {
        Lattice<decltype(lhs._odata[0]+rhs)> ret(lhs._grid);
PARALLEL_FOR_LOOP
        for(int ss=0;ss<rhs._grid->oSites(); ss++){
	  decltype(lhs._odata[0]+rhs) tmp=lhs._odata[ss]+rhs; 
	  vstream(ret._odata[ss],tmp);
	  //	  ret._odata[ss]=lhs._odata[ss]+rhs;
        }
        return ret;
    }
    template<class left,class right>
      inline auto operator - (const Lattice<left> &lhs,const right &rhs) -> Lattice<decltype(lhs._odata[0]-rhs)>
    {
      Lattice<decltype(lhs._odata[0]-rhs)> ret(lhs._grid);
PARALLEL_FOR_LOOP
      for(int ss=0;ss<rhs._grid->oSites(); ss++){
	  decltype(lhs._odata[0]-rhs) tmp=lhs._odata[ss]-rhs;
	  vstream(ret._odata[ss],tmp);
	  //	ret._odata[ss]=lhs._odata[ss]-rhs;
      }
      return ret;
    }


}
#endif
