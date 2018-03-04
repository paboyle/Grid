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
  auto x_v = x.View();
  auto z_v = z.View();
  accelerator_loop( ss, x_v,{
    vobj tmp;
    tmp = a*x_v[ss];
    tmp = tmp + G5*(b*timesI(x_v[ss]));
    vstream(z_v[ss],tmp);
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
  auto x_v = x.View();
  auto y_v = y.View();
  auto z_v = z.View();
  // FIXME -- need a new class of accelerator_loop to implement this
  thread_loop( (int ss=0;ss<grid->oSites();ss+=Ls),{ // adds Ls
    vobj tmp = a*x_v[ss+s]+b*y_v[ss+sp];
    vstream(z_v[ss+s],tmp);
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
  auto x_v = x.View();
  auto y_v = y.View();
  auto z_v = z.View();
  thread_loop((int ss=0;ss<grid->oSites();ss+=Ls),{ // adds Ls
    vobj tmp;
    tmp = G5*x_v[ss+s]*a;
    tmp = tmp + b*y_v[ss+sp];
    vstream(z_v[ss+s],tmp);
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
  auto x_v = x.View();
  auto y_v = y.View();
  auto z_v = z.View();
  Gamma G5(Gamma::Algebra::Gamma5);
  thread_loop((int ss=0;ss<grid->oSites();ss+=Ls),{ // adds Ls
    vobj tmp;
    tmp = G5*y_v[ss+sp]*b;
    tmp = tmp + a*x_v[ss+s];
    vstream(z_v[ss+s],tmp);
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

  auto x_v = x.View();
  auto y_v = y.View();
  auto z_v = z.View();
  Gamma G5(Gamma::Algebra::Gamma5);
  thread_loop((int ss=0;ss<grid->oSites();ss+=Ls),{ // adds Ls
    vobj tmp1;
    vobj tmp2;
    tmp1 = a*x_v[ss+s]+b*y_v[ss+sp];
    tmp2 = G5*tmp1;
    vstream(z_v[ss+s],tmp2);
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

  auto x_v = x.View();
  auto y_v = y.View();
  auto z_v = z.View();
  thread_loop((int ss=0;ss<grid->oSites();ss+=Ls),{ // adds Ls
    vobj tmp;
    spProj5m(tmp,y_v[ss+sp]);
    tmp = a*x_v[ss+s]+b*tmp;
    vstream(z_v[ss+s],tmp);
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
  auto x_v = x.View();
  auto y_v = y.View();
  auto z_v = z.View();
  thread_loop((int ss=0;ss<grid->oSites();ss+=Ls),{ // adds Ls
    vobj tmp;
    spProj5p(tmp,y_v[ss+sp]);
    tmp = a*x_v[ss+s]+b*tmp;
    vstream(z_v[ss+s],tmp);
  });
}

template<class vobj> 
void G5R5(Lattice<vobj> &z,const Lattice<vobj> &x)
{
  GridBase *grid=x.Grid();
  z.Checkerboard() = x.Checkerboard();
  conformable(x,z);
  int Ls = grid->_rdimensions[0];
  Gamma G5(Gamma::Algebra::Gamma5);
  auto x_v = x.View();
  auto z_v = z.View();
  thread_loop((int ss=0;ss<grid->oSites();ss+=Ls) {
    vobj tmp;
    for(int s=0;s<Ls;s++){
      int sp = Ls-1-s;
      tmp = G5*x_v[ss+s];
      vstream(z_v[ss+sp],tmp);
    }
  });
}

NAMESPACE_END(Grid);
#endif 
