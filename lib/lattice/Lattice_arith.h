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
  conformable(ret,rhs);
  conformable(lhs,rhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs,{
    obj1 tmp;
    mult(&tmp,&lhs[ss],&rhs[ss]);
    vstream(ret[ss],tmp);
  });
#else
  accelerator_loop(ss,lhs,{
    mult(&ret[ss],&lhs[ss],&rhs[ss]);
  });
#endif
}
  
template<class obj1,class obj2,class obj3> inline
void mac(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(ret,rhs);
  conformable(lhs,rhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs,{
    obj1 tmp;
    mac(&tmp,&lhs[ss],&rhs[ss]);
    vstream(ret[ss],tmp);
  });
#else
  accelerator_loop(ss,lhs,{
    mac(&ret[ss],&lhs[ss],&rhs[ss]);
  });
#endif
}
  
template<class obj1,class obj2,class obj3> inline
void sub(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(ret,rhs);
  conformable(lhs,rhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs,{
    obj1 tmp;
    sub(&tmp,&lhs[ss],&rhs[ss]);
    vstream(ret[ss],tmp);
  });
#else
  accelerator_loop(ss,lhs,{
    sub(&ret[ss],&lhs[ss],&rhs[ss]);
  });
#endif
}
template<class obj1,class obj2,class obj3> inline
void add(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(ret,rhs);
  conformable(lhs,rhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs,{
    obj1 tmp;
    add(&tmp,&lhs[ss],&rhs[ss]);
    vstream(ret[ss],tmp);
  });
#else
  accelerator_loop(ss,lhs,{
    add(&ret[ss],&lhs[ss],&rhs[ss]);
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
  accelerator_loop(ss,lhs,{
    obj1 tmp;
    mult(&tmp,&lhs[ss],&rhs);
    vstream(ret[ss],tmp);
  });
}
  
template<class obj1,class obj2,class obj3> inline
void mac(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(ret,lhs);
  accelerator_loop(ss,lhs,{
    obj1 tmp;
    mac(&tmp,&lhs[ss],&rhs);
    vstream(ret[ss],tmp);
  });
}
  
template<class obj1,class obj2,class obj3> inline
void sub(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(ret,lhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs,{
    obj1 tmp;
    sub(&tmp,&lhs[ss],&rhs);
    vstream(ret[ss],tmp);
  });
#else 
  accelerator_loop(ss,lhs,{
    sub(&ret[ss],&lhs[ss],&rhs);
  });
#endif
}
template<class obj1,class obj2,class obj3> inline
void add(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
  ret.Checkerboard() = lhs.Checkerboard();
  conformable(lhs,ret);
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs,{
    obj1 tmp;
    add(&tmp,&lhs[ss],&rhs);
    vstream(ret[ss],tmp);
  });
#else 
  accelerator_loop(ss,lhs,{
    add(&ret[ss],&lhs[ss],&rhs);
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
#ifdef STREAMING_STORES
  accelerator_loop(ss,rhs,{
    obj1 tmp;
    mult(&tmp,&lhs,&rhs[ss]);
    vstream(ret[ss],tmp);
  });
#else 
  accelerator_loop(ss,rhs,{
    mult(&ret[ss],&lhs,&rhs[ss]);
  });
#endif
}
  
template<class obj1,class obj2,class obj3> inline
void mac(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = rhs.Checkerboard();
  conformable(ret,rhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,rhs,{
    obj1 tmp;
    mac(&tmp,&lhs,&rhs[ss]);
    vstream(ret[ss],tmp);
  });
#else 
  accelerator_loop(ss,rhs,{
    mac(&ret[ss],&lhs,&rhs[ss]);
  });
#endif
}
  
template<class obj1,class obj2,class obj3> inline
void sub(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = rhs.Checkerboard();
  conformable(ret,rhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,rhs,{
    obj1 tmp;
    sub(&tmp,&lhs,&rhs[ss]);
    vstream(ret[ss],tmp);
  });
#else 
  accelerator_loop(ss,rhs,{
    sub(&ret[ss],&lhs,&rhs[ss]);
  });
#endif
}
template<class obj1,class obj2,class obj3> inline
void add(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
  ret.Checkerboard() = rhs.Checkerboard();
  conformable(ret,rhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,rhs,{
    obj1 tmp;
    add(&tmp,&lhs,&rhs[ss]);
    vstream(ret[ss],tmp);
  });
#else 
  accelerator_loop(ss,rhs,{
    add(&ret[ss],&lhs,&rhs[ss]);
  });
#endif
}
  
template<class sobj,class vobj> inline
void axpy(Lattice<vobj> &ret,sobj a,const Lattice<vobj> &x,const Lattice<vobj> &y){
  ret.Checkerboard() = x.Checkerboard();
  conformable(ret,x);
  conformable(x,y);
#ifdef STREAMING_STORES
  accelerator_loop(ss,x,{
    vobj tmp = a*x[ss]+y[ss];
    vstream(ret[ss],tmp);
  });
#else
  accelerator_loop(ss,x,{
    ret[ss]=a*x[ss]+y[ss];
  });
#endif
}
template<class sobj,class vobj> inline
void axpby(Lattice<vobj> &ret,sobj a,sobj b,const Lattice<vobj> &x,const Lattice<vobj> &y){
  ret.Checkerboard() = x.Checkerboard();
  conformable(ret,x);
  conformable(x,y);
#ifdef STREAMING_STORES
  accelerator_loop(ss,x,{
    vobj tmp = a*x[ss]+b*y[ss];
    vstream(ret[ss],tmp);
  });
#else
  accelerator_loop(ss,x,{
    ret[ss]=a*x[ss]+b*y[ss];
  });
#endif
}

template<class sobj,class vobj> inline
RealD axpy_norm(Lattice<vobj> &ret,sobj a,const Lattice<vobj> &x,const Lattice<vobj> &y){
  ret.Checkerboard() = x.Checkerboard();
  conformable(ret,x);
  conformable(x,y);
  axpy(ret,a,x,y);
  return norm2(ret);
}
template<class sobj,class vobj> inline
RealD axpby_norm(Lattice<vobj> &ret,sobj a,sobj b,const Lattice<vobj> &x,const Lattice<vobj> &y){
  ret.Checkerboard() = x.Checkerboard();
  conformable(ret,x);
  conformable(x,y);
  axpby(ret,a,b,x,y);
  return norm2(ret); // FIXME implement parallel norm in ss loop
}

NAMESPACE_END(Grid);
#endif
