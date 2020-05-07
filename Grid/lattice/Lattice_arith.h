/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_arith.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Christoph Lehner <christoph@lhnr.de>

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
#ifndef GRID_LATTICE_ARITH_H
#define GRID_LATTICE_ARITH_H

NAMESPACE_BEGIN(Grid);

//////////////////////////////////////////////////////////////////////////////////////////////////////
//  avoid copy back routines for mult, mac, sub, add
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<class obj1,class obj2,class obj3> inline
void mult(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  auto ret_v = ret.AcceleratorView(ViewWrite);
  auto lhs_v = lhs.AcceleratorView(ViewRead);
  auto rhs_v = rhs.AcceleratorView(ViewRead);
  conformable(ret,rhs);
  conformable(lhs,rhs);
  accelerator_for(ss,lhs_v.size(),obj1::Nsimd(),{
    decltype(coalescedRead(obj1())) tmp;
    auto lhs_t = lhs_v(ss);
    auto rhs_t = rhs_v(ss);
    mult(&tmp,&lhs_t,&rhs_t);
    coalescedWrite(ret_v[ss],tmp);
  });
}
  
template<class obj1,class obj2,class obj3> inline
void mac(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(ret,rhs);
  conformable(lhs,rhs);
  auto ret_v = ret.AcceleratorView(ViewWrite);
  auto lhs_v = lhs.AcceleratorView(ViewRead);
  auto rhs_v = rhs.AcceleratorView(ViewRead);
  accelerator_for(ss,lhs_v.size(),obj1::Nsimd(),{
    decltype(coalescedRead(obj1())) tmp;
    auto lhs_t=lhs_v(ss);
    auto rhs_t=rhs_v(ss);
    mac(&tmp,&lhs_t,&rhs_t);
    coalescedWrite(ret_v[ss],tmp);
  });
}
  
template<class obj1,class obj2,class obj3> inline
void sub(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(ret,rhs);
  conformable(lhs,rhs);
  auto ret_v = ret.AcceleratorView(ViewWrite);
  auto lhs_v = lhs.AcceleratorView(ViewRead);
  auto rhs_v = rhs.AcceleratorView(ViewRead);
  accelerator_for(ss,lhs_v.size(),obj1::Nsimd(),{
    decltype(coalescedRead(obj1())) tmp;
    auto lhs_t=lhs_v(ss);
    auto rhs_t=rhs_v(ss);
    sub(&tmp,&lhs_t,&rhs_t);
    coalescedWrite(ret_v[ss],tmp);
  });
}
template<class obj1,class obj2,class obj3> inline
void add(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(ret,rhs);
  conformable(lhs,rhs);
  auto ret_v = ret.AcceleratorView(ViewWrite);
  auto lhs_v = lhs.AcceleratorView(ViewRead);
  auto rhs_v = rhs.AcceleratorView(ViewRead);
  accelerator_for(ss,lhs_v.size(),obj1::Nsimd(),{
    decltype(coalescedRead(obj1())) tmp;
    auto lhs_t=lhs_v(ss);
    auto rhs_t=rhs_v(ss);
    add(&tmp,&lhs_t,&rhs_t);
    coalescedWrite(ret_v[ss],tmp);
  });
}
  
//////////////////////////////////////////////////////////////////////////////////////////////////////
//  avoid copy back routines for mult, mac, sub, add
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<class obj1,class obj2,class obj3> inline
void mult(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(lhs,ret);
  auto ret_v = ret.AcceleratorView(ViewWrite);
  auto lhs_v = lhs.AcceleratorView(ViewRead);
  accelerator_for(ss,lhs_v.size(),obj1::Nsimd(),{
    decltype(coalescedRead(obj1())) tmp;
    mult(&tmp,&lhs_v(ss),&rhs);
    coalescedWrite(ret_v[ss],tmp);
  });
}
  
template<class obj1,class obj2,class obj3> inline
void mac(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(ret,lhs);
  auto ret_v = ret.AcceleratorView(ViewWrite);
  auto lhs_v = lhs.AcceleratorView(ViewRead);
  accelerator_for(ss,lhs_v.size(),obj1::Nsimd(),{
    decltype(coalescedRead(obj1())) tmp;
    auto lhs_t=lhs_v(ss);
    mac(&tmp,&lhs_t,&rhs);
    coalescedWrite(ret_v[ss],tmp);
  });
}
  
template<class obj1,class obj2,class obj3> inline
void sub(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(ret,lhs);
  auto ret_v = ret.AcceleratorView(ViewWrite);
  auto lhs_v = lhs.AcceleratorView(ViewRead);
  accelerator_for(ss,lhs_v.size(),obj1::Nsimd(),{
    decltype(coalescedRead(obj1())) tmp;
    auto lhs_t=lhs_v(ss);
    sub(&tmp,&lhs_t,&rhs);
    coalescedWrite(ret_v[ss],tmp);
  });
}
template<class obj1,class obj2,class obj3> inline
void add(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(lhs,ret);
  auto ret_v = ret.AcceleratorView(ViewWrite);
  auto lhs_v = lhs.AcceleratorView(ViewRead);
  accelerator_for(ss,lhs_v.size(),obj1::Nsimd(),{
    decltype(coalescedRead(obj1())) tmp;
    auto lhs_t=lhs_v(ss);
    add(&tmp,&lhs_t,&rhs);
    coalescedWrite(ret_v[ss],tmp);
  });
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//  avoid copy back routines for mult, mac, sub, add
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<class obj1,class obj2,class obj3> inline
void mult(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = rhs.Checkerboard();
  conformable(ret,rhs);
  auto ret_v = ret.AcceleratorView(ViewWrite);
  auto rhs_v = lhs.AcceleratorView(ViewRead);
  accelerator_for(ss,rhs_v.size(),obj1::Nsimd(),{
    decltype(coalescedRead(obj1())) tmp;
    auto rhs_t=rhs_v(ss);
    mult(&tmp,&lhs,&rhs_t);
    coalescedWrite(ret_v[ss],tmp);
  });
}
  
template<class obj1,class obj2,class obj3> inline
void mac(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = rhs.Checkerboard();
  conformable(ret,rhs);
  auto ret_v = ret.AcceleratorView(ViewWrite);
  auto rhs_v = lhs.AcceleratorView(ViewRead);
  accelerator_for(ss,rhs_v.size(),obj1::Nsimd(),{
    decltype(coalescedRead(obj1())) tmp;
    auto rhs_t=rhs_v(ss);
    mac(&tmp,&lhs,&rhs_t);
    coalescedWrite(ret_v[ss],tmp);
  });
}
  
template<class obj1,class obj2,class obj3> inline
void sub(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = rhs.Checkerboard();
  conformable(ret,rhs);
  auto ret_v = ret.AcceleratorView(ViewWrite);
  auto rhs_v = lhs.AcceleratorView(ViewRead);
  accelerator_for(ss,rhs_v.size(),obj1::Nsimd(),{
    decltype(coalescedRead(obj1())) tmp;
    auto rhs_t=rhs_v(ss);
    sub(&tmp,&lhs,&rhs_t);
    coalescedWrite(ret_v[ss],tmp);
  });
}
template<class obj1,class obj2,class obj3> inline
void add(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = rhs.Checkerboard();
  conformable(ret,rhs);
  auto ret_v = ret.AcceleratorView(ViewWrite);
  auto rhs_v = lhs.AcceleratorView(ViewRead);
  accelerator_for(ss,rhs_v.size(),obj1::Nsimd(),{
    decltype(coalescedRead(obj1())) tmp;
    auto rhs_t=rhs_v(ss);
    add(&tmp,&lhs,&rhs_t);
    coalescedWrite(ret_v[ss],tmp);
  });
}
  
template<class sobj,class vobj> inline
void axpy(Lattice<vobj> &ret,sobj a,const Lattice<vobj> &x,const Lattice<vobj> &y){
  ret.Checkerboard() = x.Checkerboard();
  conformable(ret,x);
  conformable(x,y);
  auto ret_v = ret.AcceleratorView(ViewWrite);
  auto x_v = x.AcceleratorView(ViewRead);
  auto y_v = y.AcceleratorView(ViewRead);
  accelerator_for(ss,x_v.size(),vobj::Nsimd(),{
    auto tmp = a*x_v(ss)+y_v(ss);
    coalescedWrite(ret_v[ss],tmp);
  });
}
template<class sobj,class vobj> inline
void axpby(Lattice<vobj> &ret,sobj a,sobj b,const Lattice<vobj> &x,const Lattice<vobj> &y){
  ret.Checkerboard() = x.Checkerboard();
  conformable(ret,x);
  conformable(x,y);
  auto ret_v = ret.AcceleratorView(ViewWrite);
  auto x_v = x.AcceleratorView(ViewRead);
  auto y_v = y.AcceleratorView(ViewRead);
  accelerator_for(ss,x_v.size(),vobj::Nsimd(),{
    auto tmp = a*x_v(ss)+b*y_v(ss);
    coalescedWrite(ret_v[ss],tmp);
  });
}

template<class sobj,class vobj> inline
RealD axpy_norm(Lattice<vobj> &ret,sobj a,const Lattice<vobj> &x,const Lattice<vobj> &y)
{
    return axpy_norm_fast(ret,a,x,y);
}
template<class sobj,class vobj> inline
RealD axpby_norm(Lattice<vobj> &ret,sobj a,sobj b,const Lattice<vobj> &x,const Lattice<vobj> &y)
{
    return axpby_norm_fast(ret,a,b,x,y);
}

NAMESPACE_END(Grid);
#endif
