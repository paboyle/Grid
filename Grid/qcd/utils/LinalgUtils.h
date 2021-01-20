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
#pragma once

NAMESPACE_BEGIN(Grid);

////////////////////////////////////////////////////////////////////////
//This file brings additional linear combination assist that is helpful
//to QCD such as chiral projectors and spin matrices applied to one of the inputs.
//These routines support five-D chiral fermions and contain s-subslice indexing 
//on the 5d (rb4d) checkerboarded lattices
////////////////////////////////////////////////////////////////////////

template<class vobj,class Coeff>
void axpibg5x(Lattice<vobj> &z,const Lattice<vobj> &x,Coeff a,Coeff b)
{
  z.Checkerboard() = x.Checkerboard();
  conformable(x,z);

  GridBase *grid=x.Grid();

  Gamma G5(Gamma::Algebra::Gamma5);
  autoView(x_v, x, AcceleratorRead);
  autoView(z_v, z, AcceleratorWrite);
  accelerator_for( ss, x_v.size(),vobj::Nsimd(), {
    auto tmp = a*x_v(ss) + G5*(b*timesI(x_v(ss)));
    coalescedWrite(z_v[ss],tmp);
  });
}

template<class vobj,class Coeff> 
void axpby_ssp(Lattice<vobj> &z, Coeff a,const Lattice<vobj> &x,Coeff b,const Lattice<vobj> &y,int s,int sp)
{
  z.Checkerboard() = x.Checkerboard();
  conformable(x,y);
  conformable(x,z);
  GridBase *grid=x.Grid();
  int Ls = grid->_rdimensions[0];
  autoView( x_v, x, AcceleratorRead);
  autoView( y_v, y, AcceleratorRead);
  autoView( z_v, z, AcceleratorWrite);
  // FIXME -- need a new class of accelerator_loop to implement this
  //
  uint64_t nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,vobj::Nsimd(),{
    uint64_t ss = sss*Ls;
    auto tmp = a*x_v(ss+s)+b*y_v(ss+sp);
    coalescedWrite(z_v[ss+s],tmp);
  });
}

template<class vobj,class Coeff> 
void ag5xpby_ssp(Lattice<vobj> &z,Coeff a,const Lattice<vobj> &x,Coeff b,const Lattice<vobj> &y,int s,int sp)
{
  z.Checkerboard() = x.Checkerboard();
  conformable(x,y);
  conformable(x,z);
  GridBase *grid=x.Grid();
  int Ls = grid->_rdimensions[0];
  Gamma G5(Gamma::Algebra::Gamma5);
  autoView( x_v, x, AcceleratorRead);
  autoView( y_v, y, AcceleratorRead);
  autoView( z_v, z, AcceleratorWrite);
  uint64_t nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,vobj::Nsimd(),{
    uint64_t ss = sss*Ls;
    auto tmp = G5*x_v(ss+s)*a + b*y_v(ss+sp);
    coalescedWrite(z_v[ss+s],tmp);
  });
}

template<class vobj,class Coeff> 
void axpbg5y_ssp(Lattice<vobj> &z,Coeff a,const Lattice<vobj> &x,Coeff b,const Lattice<vobj> &y,int s,int sp)
{
  z.Checkerboard() = x.Checkerboard();
  conformable(x,y);
  conformable(x,z);
  GridBase *grid=x.Grid();
  int Ls = grid->_rdimensions[0];
  autoView( x_v, x, AcceleratorRead);
  autoView( y_v, y, AcceleratorRead);
  autoView( z_v, z, AcceleratorWrite);
  Gamma G5(Gamma::Algebra::Gamma5);
  uint64_t nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,vobj::Nsimd(),{
    uint64_t ss = sss*Ls;
    auto tmp = G5*y_v(ss+sp)*b + a*x_v(ss+s);
    coalescedWrite(z_v[ss+s],tmp);
  });
}

template<class vobj,class Coeff> 
void ag5xpbg5y_ssp(Lattice<vobj> &z,Coeff a,const Lattice<vobj> &x,Coeff b,const Lattice<vobj> &y,int s,int sp)
{
  z.Checkerboard() = x.Checkerboard();
  conformable(x,y);
  conformable(x,z);
  GridBase *grid=x.Grid();
  int Ls = grid->_rdimensions[0];

  autoView( x_v, x, AcceleratorRead);
  autoView( y_v, y, AcceleratorRead);
  autoView( z_v, z, AcceleratorWrite);
  Gamma G5(Gamma::Algebra::Gamma5);
  uint64_t nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,vobj::Nsimd(),{
    uint64_t ss = sss*Ls;
    auto tmp1 = a*x_v(ss+s)+b*y_v(ss+sp);
    auto tmp2 = G5*tmp1;
    coalescedWrite(z_v[ss+s],tmp2);
  });
}

template<class vobj,class Coeff> 
void axpby_ssp_pminus(Lattice<vobj> &z,Coeff a,const Lattice<vobj> &x,Coeff b,const Lattice<vobj> &y,int s,int sp)
{
  z.Checkerboard() = x.Checkerboard();
  conformable(x,y);
  conformable(x,z);
  GridBase *grid=x.Grid();
  int Ls = grid->_rdimensions[0];

  autoView( x_v, x, AcceleratorRead);
  autoView( y_v, y, AcceleratorRead);
  autoView( z_v, z, AcceleratorWrite);
  uint64_t nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,vobj::Nsimd(),{
    uint64_t ss = sss*Ls;
    decltype(coalescedRead(y_v[ss+sp])) tmp;
    spProj5m(tmp,y_v(ss+sp)); 
   tmp = a*x_v(ss+s)+b*tmp;
    coalescedWrite(z_v[ss+s],tmp);
  });
}

template<class vobj,class Coeff> 
void axpby_ssp_pplus(Lattice<vobj> &z,Coeff a,const Lattice<vobj> &x,Coeff b,const Lattice<vobj> &y,int s,int sp)
{
  z.Checkerboard() = x.Checkerboard();
  conformable(x,y);
  conformable(x,z);
  GridBase *grid=x.Grid();
  int Ls = grid->_rdimensions[0];
  autoView( x_v, x, AcceleratorRead);
  autoView( y_v, y, AcceleratorRead);
  autoView( z_v, z, AcceleratorWrite);
  uint64_t nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,vobj::Nsimd(),{
    uint64_t ss = sss*Ls;
    decltype(coalescedRead(y_v[ss+sp])) tmp;
    spProj5p(tmp,y_v(ss+sp));
    tmp = a*x_v(ss+s)+b*tmp;
    coalescedWrite(z_v[ss+s],tmp);
  });
}

template<class vobj> 
void G5R5(Lattice<vobj> &z,const Lattice<vobj> &x)
{
  GridBase *grid=x.Grid();
  z.Checkerboard() = x.Checkerboard();
  conformable(x,z);
  int Ls = grid->_rdimensions[0];
  autoView( x_v, x, AcceleratorRead);
  autoView( z_v, z, AcceleratorWrite);
  uint64_t nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,vobj::Nsimd(),{
    uint64_t ss = sss*Ls;
    for(int s=0;s<Ls;s++){
      int sp = Ls-1-s;
      auto tmp = x_v(ss+s);
      decltype(tmp) tmp_p;
      decltype(tmp) tmp_m;
      spProj5p(tmp_p,tmp);
      spProj5m(tmp_m,tmp);
      // Use of spProj5m, 5p captures the coarse space too
      coalescedWrite(z_v[ss+sp],tmp_p - tmp_m);
    }
  });
}

template<typename vobj>
void G5C(Lattice<vobj> &z, const Lattice<vobj> &x)
{
  GridBase *grid = x.Grid();
  z.Checkerboard() = x.Checkerboard();
  conformable(x, z);

  autoView( x_v, x, AcceleratorRead);
  autoView( z_v, z, AcceleratorWrite);
  uint64_t nloop = grid->oSites();
  accelerator_for(ss,nloop,vobj::Nsimd(),{
    auto tmp = x_v(ss);
    decltype(tmp) tmp_p;
    decltype(tmp) tmp_m;
    spProj5p(tmp_p,tmp);
    spProj5m(tmp_m,tmp);
    coalescedWrite(z_v[ss],tmp_p - tmp_m);
  });
}

/*
template<class CComplex, int nbasis>
void G5C(Lattice<iVector<CComplex, nbasis>> &z, const Lattice<iVector<CComplex, nbasis>> &x)
{
  GridBase *grid = x.Grid();
  z.Checkerboard() = x.Checkerboard();
  conformable(x, z);

  static_assert(nbasis % 2 == 0, "");
  int nb = nbasis / 2;

  autoView( z_v, z, AcceleratorWrite);
  autoView( x_v, x, AcceleratorRead);
  accelerator_for(ss,grid->oSites(),CComplex::Nsimd(),
  {
    for(int n = 0; n < nb; ++n) {
      coalescedWrite(z_v[ss](n), x_v(ss)(n));
    }
    for(int n = nb; n < nbasis; ++n) {
      coalescedWrite(z_v[ss](n), -x_v(ss)(n));
    }
  });
}
*/

NAMESPACE_END(Grid);

