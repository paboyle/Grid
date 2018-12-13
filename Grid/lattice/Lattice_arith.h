/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_arith.h

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
#ifndef GRID_LATTICE_ARITH_H
#define GRID_LATTICE_ARITH_H

NAMESPACE_BEGIN(Grid);

//////////////////////////////////////////////////////////////////////////////////////////////////////
//  avoid copy back routines for mult, mac, sub, add
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<class obj1,class obj2,class obj3> inline
void mult(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  auto ret_v = ret.View();
  auto lhs_v = lhs.View();
  auto rhs_v = rhs.View();
  conformable(ret,rhs);
  conformable(lhs,rhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs_v,{
    obj1 tmp;
    mult(&tmp,&lhs_v[ss],&rhs_v[ss]);
    vstream(ret_v[ss],tmp);
  });
#else
  accelerator_loop(ss,lhs_v,{
    mult(&ret_v[ss],&lhs_v[ss],&rhs_v[ss]);
  });
#endif
}
  
template<class obj1,class obj2,class obj3> inline
void mac(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(ret,rhs);
  conformable(lhs,rhs);
  auto ret_v = ret.View();
  auto lhs_v = lhs.View();
  auto rhs_v = rhs.View();
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs_v,{
    obj1 tmp;
    mac(&tmp,&lhs_v[ss],&rhs_v[ss]);
    vstream(ret_v[ss],tmp);
  });
#else
  accelerator_loop(ss,lhs_v,{
    mac(&ret_v[ss],&lhs_v[ss],&rhs_v[ss]);
  });
#endif
}
  
template<class obj1,class obj2,class obj3> inline
void sub(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(ret,rhs);
  conformable(lhs,rhs);
  auto ret_v = ret.View();
  auto lhs_v = lhs.View();
  auto rhs_v = rhs.View();
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs_v,{
    obj1 tmp;
    sub(&tmp,&lhs_v[ss],&rhs_v[ss]);
    vstream(ret_v[ss],tmp);
  });
#else
  accelerator_loop(ss,lhs_v,{
    sub(&ret[ss],&lhs_v[ss],&rhs_v[ss]);
  });
#endif
}
template<class obj1,class obj2,class obj3> inline
void add(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(ret,rhs);
  conformable(lhs,rhs);
  auto ret_v = ret.View();
  auto lhs_v = lhs.View();
  auto rhs_v = rhs.View();
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs_v,{
    obj1 tmp;
    add(&tmp,&lhs_v[ss],&rhs_v[ss]);
    vstream(ret_v[ss],tmp);
  });
#else
  accelerator_loop(ss,lhs_v,{
    add(&ret_v[ss],&lhs_v[ss],&rhs_v[ss]);
  });
#endif
}
  
//////////////////////////////////////////////////////////////////////////////////////////////////////
//  avoid copy back routines for mult, mac, sub, add
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<class obj1,class obj2,class obj3> inline
void mult(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(lhs,ret);
  auto ret_v = ret.View();
  auto lhs_v = lhs.View();
  accelerator_loop(ss,lhs_v,{
    obj1 tmp;
    mult(&tmp,&lhs_v[ss],&rhs);
    vstream(ret_v[ss],tmp);
  });
}
  
template<class obj1,class obj2,class obj3> inline
void mac(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(ret,lhs);
  auto ret_v = ret.View();
  auto lhs_v = lhs.View();
  accelerator_loop(ss,lhs_v,{
    obj1 tmp;
    mac(&tmp,&lhs_v[ss],&rhs);
    vstream(ret_v[ss],tmp);
  });
}
  
template<class obj1,class obj2,class obj3> inline
void sub(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(ret,lhs);
  auto ret_v = ret.View();
  auto lhs_v = lhs.View();
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs_v,{
    obj1 tmp;
    sub(&tmp,&lhs_v[ss],&rhs);
    vstream(ret_v[ss],tmp);
  });
#else 
  accelerator_loop(ss,lhs_v,{
    sub(&ret_v[ss],&lhs_v[ss],&rhs);
  });
#endif
}
template<class obj1,class obj2,class obj3> inline
void add(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(lhs,ret);
  auto ret_v = ret.View();
  auto lhs_v = lhs.View();
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs_v,{
    obj1 tmp;
    add(&tmp,&lhs_v[ss],&rhs);
    vstream(ret_v[ss],tmp);
  });
#else 
  accelerator_loop(ss,lhs_v,{
    add(&ret_v[ss],&lhs_v[ss],&rhs);
  });
#endif
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//  avoid copy back routines for mult, mac, sub, add
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<class obj1,class obj2,class obj3> inline
void mult(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = rhs.Checkerboard();
  conformable(ret,rhs);
  auto ret_v = ret.View();
  auto rhs_v = lhs.View();
#ifdef STREAMING_STORES
  accelerator_loop(ss,rhs_v,{
    obj1 tmp;
    mult(&tmp,&lhs,&rhs_v[ss]);
    vstream(ret_v[ss],tmp);
  });
#else 
  accelerator_loop(ss,rhs_v,{
    mult(&ret_v[ss],&lhs,&rhs_v[ss]);
  });
#endif
}
  
template<class obj1,class obj2,class obj3> inline
void mac(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = rhs.Checkerboard();
  conformable(ret,rhs);
  auto ret_v = ret.View();
  auto rhs_v = lhs.View();
#ifdef STREAMING_STORES
  accelerator_loop(ss,rhs_v,{
    obj1 tmp;
    mac(&tmp,&lhs,&rhs_v[ss]);
    vstream(ret_v[ss],tmp);
  });
#else 
  accelerator_loop(ss,rhs_v,{
    mac(&ret_v[ss],&lhs,&rhs_v[ss]);
  });
#endif
}
  
template<class obj1,class obj2,class obj3> inline
void sub(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = rhs.Checkerboard();
  conformable(ret,rhs);
  auto ret_v = ret.View();
  auto rhs_v = lhs.View();
#ifdef STREAMING_STORES
  accelerator_loop(ss,rhs_v,{
    obj1 tmp;
    sub(&tmp,&lhs,&rhs_v[ss]);
    vstream(ret_v[ss],tmp);
  });
#else 
  accelerator_loop(ss,rhs_v,{
    sub(&ret_v[ss],&lhs,&rhs_v[ss]);
  });
#endif
}
template<class obj1,class obj2,class obj3> inline
void add(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = rhs.Checkerboard();
  conformable(ret,rhs);
  auto ret_v = ret.View();
  auto rhs_v = lhs.View();
#ifdef STREAMING_STORES
  accelerator_loop(ss,rhs_v,{
    obj1 tmp;
    add(&tmp,&lhs,&rhs_v[ss]);
    vstream(ret_v[ss],tmp);
  });
#else 
  accelerator_loop(ss,rhs_v,{
    add(&ret_v[ss],&lhs,&rhs_v[ss]);
  });
#endif
}
  
template<class sobj,class vobj> inline
void axpy(Lattice<vobj> &ret,sobj a,const Lattice<vobj> &x,const Lattice<vobj> &y){
  ret.Checkerboard() = x.Checkerboard();
  conformable(ret,x);
  conformable(x,y);
  auto ret_v = ret.View();
  auto x_v = x.View();
  auto y_v = y.View();
#ifdef STREAMING_STORES
  accelerator_loop(ss,x_v,{
    vobj tmp = a*x_v[ss]+y_v[ss];
    vstream(ret_v[ss],tmp);
  });
#else
  accelerator_loop(ss,x_v,{
    ret_v[ss]=a*x_v[ss]+y_v[ss];
  });
#endif
}
template<class sobj,class vobj> inline
void axpby(Lattice<vobj> &ret,sobj a,sobj b,const Lattice<vobj> &x,const Lattice<vobj> &y){
  ret.Checkerboard() = x.Checkerboard();
  conformable(ret,x);
  conformable(x,y);
  auto ret_v = ret.View();
  auto x_v = x.View();
  auto y_v = y.View();
#ifdef STREAMING_STORES
  accelerator_loop(ss,x_v,{
    vobj tmp = a*x_v[ss]+b*y_v[ss];
    vstream(ret_v[ss],tmp);
  });
#else
  accelerator_loop(ss,x_v,{
    ret_v[ss]=a*x_v[ss]+b*y_v[ss];
  });
#endif
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
