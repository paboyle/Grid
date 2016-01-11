    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_unary.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
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
#ifndef GRID_LATTICE_UNARY_H
#define GRID_LATTICE_UNARY_H

namespace Grid {

  template<class obj> Lattice<obj> pow(const Lattice<obj> &rhs,RealD y){
    Lattice<obj> ret(rhs._grid);
    ret.checkerboard = rhs.checkerboard;
    conformable(ret,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites();ss++){
      ret._odata[ss]=pow(rhs._odata[ss],y);
    }
    return ret;
  }
  template<class obj> Lattice<obj> mod(const Lattice<obj> &rhs,Integer y){
    Lattice<obj> ret(rhs._grid);
    ret.checkerboard = rhs.checkerboard;
    conformable(ret,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites();ss++){
      ret._odata[ss]=mod(rhs._odata[ss],y);
    }
    return ret;
  }

  template<class obj> Lattice<obj> div(const Lattice<obj> &rhs,Integer y){
    Lattice<obj> ret(rhs._grid);
    ret.checkerboard = rhs.checkerboard;
    conformable(ret,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites();ss++){
      ret._odata[ss]=div(rhs._odata[ss],y);
    }
    return ret;
  }

  template<class obj> Lattice<obj> expMat(const Lattice<obj> &rhs, ComplexD alpha, Integer Nexp = DEFAULT_MAT_EXP){
    Lattice<obj> ret(rhs._grid);
    ret.checkerboard = rhs.checkerboard;
    conformable(ret,rhs);
PARALLEL_FOR_LOOP
    for(int ss=0;ss<rhs._grid->oSites();ss++){
      ret._odata[ss]=Exponentiate(rhs._odata[ss],alpha, Nexp);
    }
    return ret;
  }



}
#endif
