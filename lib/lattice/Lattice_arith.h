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
  ret.checkerboard = lhs.checkerboard;
  conformable(ret,rhs);
  conformable(lhs,rhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs,{
    obj1 tmp;
    mult(&tmp,&lhs._odata[ss],&rhs._odata[ss]);
    vstream(ret._odata[ss],tmp);
  });
#else
  accelerator_loop(ss,lhs,{
    mult(&ret._odata[ss],&lhs._odata[ss],&rhs._odata[ss]);
  });
#endif
}
  
template<class obj1,class obj2,class obj3> inline
void mac(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
  ret.checkerboard = lhs.checkerboard;
  conformable(ret,rhs);
  conformable(lhs,rhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs,{
    obj1 tmp;
    mac(&tmp,&lhs._odata[ss],&rhs._odata[ss]);
    vstream(ret._odata[ss],tmp);
  });
#else
  accelerator_loop(ss,lhs,{
    mac(&ret._odata[ss],&lhs._odata[ss],&rhs._odata[ss]);
  });
#endif
}
  
template<class obj1,class obj2,class obj3> inline
void sub(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
  ret.checkerboard = lhs.checkerboard;
  conformable(ret,rhs);
  conformable(lhs,rhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs,{
    obj1 tmp;
    sub(&tmp,&lhs._odata[ss],&rhs._odata[ss]);
    vstream(ret._odata[ss],tmp);
  });
#else
  accelerator_loop(ss,lhs,{
    sub(&ret._odata[ss],&lhs._odata[ss],&rhs._odata[ss]);
  });
#endif
}
template<class obj1,class obj2,class obj3> inline
void add(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const Lattice<obj3> &rhs){
  ret.checkerboard = lhs.checkerboard;
  conformable(ret,rhs);
  conformable(lhs,rhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs,{
    obj1 tmp;
    add(&tmp,&lhs._odata[ss],&rhs._odata[ss]);
    vstream(ret._odata[ss],tmp);
  });
#else
  accelerator_loop(ss,lhs,{
    add(&ret._odata[ss],&lhs._odata[ss],&rhs._odata[ss]);
  });
#endif
}
  
//////////////////////////////////////////////////////////////////////////////////////////////////////
//  avoid copy back routines for mult, mac, sub, add
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<class obj1,class obj2,class obj3> inline
void mult(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
  ret.checkerboard = lhs.checkerboard;
  conformable(lhs,ret);
  accelerator_loop(ss,lhs,{
    obj1 tmp;
    mult(&tmp,&lhs._odata[ss],&rhs);
    vstream(ret._odata[ss],tmp);
  });
}
  
template<class obj1,class obj2,class obj3> inline
void mac(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
  ret.checkerboard = lhs.checkerboard;
  conformable(ret,lhs);
  accelerator_loop(ss,lhs,{
    obj1 tmp;
    mac(&tmp,&lhs._odata[ss],&rhs);
    vstream(ret._odata[ss],tmp);
  });
}
  
template<class obj1,class obj2,class obj3> inline
void sub(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
  ret.checkerboard = lhs.checkerboard;
  conformable(ret,lhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs,{
    obj1 tmp;
    sub(&tmp,&lhs._odata[ss],&rhs);
    vstream(ret._odata[ss],tmp);
  });
#else 
  accelerator_loop(ss,lhs,{
    sub(&ret._odata[ss],&lhs._odata[ss],&rhs);
  });
#endif
}
template<class obj1,class obj2,class obj3> inline
void add(Lattice<obj1> &ret,const Lattice<obj2> &lhs,const obj3 &rhs){
  ret.checkerboard = lhs.checkerboard;
  conformable(lhs,ret);
#ifdef STREAMING_STORES
  accelerator_loop(ss,lhs,{
    obj1 tmp;
    add(&tmp,&lhs._odata[ss],&rhs);
    vstream(ret._odata[ss],tmp);
  });
#else 
  accelerator_loop(ss,lhs,{
    add(&ret._odata[ss],&lhs._odata[ss],&rhs);
  });
#endif
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//  avoid copy back routines for mult, mac, sub, add
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<class obj1,class obj2,class obj3> inline
void mult(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
  ret.checkerboard = rhs.checkerboard;
  conformable(ret,rhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,rhs,{
    obj1 tmp;
    mult(&tmp,&lhs,&rhs._odata[ss]);
    vstream(ret._odata[ss],tmp);
  });
#else 
  accelerator_loop(ss,rhs,{
    mult(&ret._odata[ss],&lhs,&rhs._odata[ss]);
  });
#endif
}
  
template<class obj1,class obj2,class obj3> inline
void mac(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
  ret.checkerboard = rhs.checkerboard;
  conformable(ret,rhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,rhs,{
    obj1 tmp;
    mac(&tmp,&lhs,&rhs._odata[ss]);
    vstream(ret._odata[ss],tmp);
  });
#else 
  accelerator_loop(ss,rhs,{
    mac(&ret._odata[ss],&lhs,&rhs._odata[ss]);
  });
#endif
}
  
template<class obj1,class obj2,class obj3> inline
void sub(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
  ret.checkerboard = rhs.checkerboard;
  conformable(ret,rhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,rhs,{
    obj1 tmp;
    sub(&tmp,&lhs,&rhs._odata[ss]);
    vstream(ret._odata[ss],tmp);
  });
#else 
  accelerator_loop(ss,rhs,{
    sub(&ret._odata[ss],&lhs,&rhs._odata[ss]);
  });
#endif
}
template<class obj1,class obj2,class obj3> inline
void add(Lattice<obj1> &ret,const obj2 &lhs,const Lattice<obj3> &rhs){
  ret.checkerboard = rhs.checkerboard;
  conformable(ret,rhs);
#ifdef STREAMING_STORES
  accelerator_loop(ss,rhs,{
    obj1 tmp;
    add(&tmp,&lhs,&rhs._odata[ss]);
    vstream(ret._odata[ss],tmp);
  });
#else 
  accelerator_loop(ss,rhs,{
    add(&ret._odata[ss],&lhs,&rhs._odata[ss]);
  });
#endif
}
  
template<class sobj,class vobj> inline
void axpy(Lattice<vobj> &ret,sobj a,const Lattice<vobj> &x,const Lattice<vobj> &y){
  ret.checkerboard = x.checkerboard;
  conformable(ret,x);
  conformable(x,y);
#ifdef STREAMING_STORES
  accelerator_loop(ss,x,{
    vobj tmp = a*x._odata[ss]+y._odata[ss];
    vstream(ret._odata[ss],tmp);
  });
#else
  accelerator_loop(ss,x,{
    ret._odata[ss]=a*x._odata[ss]+y._odata[ss];
  });
#endif
}
template<class sobj,class vobj> inline
void axpby(Lattice<vobj> &ret,sobj a,sobj b,const Lattice<vobj> &x,const Lattice<vobj> &y){
  ret.checkerboard = x.checkerboard;
  conformable(ret,x);
  conformable(x,y);
#ifdef STREAMING_STORES
  accelerator_loop(ss,x,{
    vobj tmp = a*x._odata[ss]+b*y._odata[ss];
    vstream(ret._odata[ss],tmp);
  });
#else
  accelerator_loop(ss,x,{
    ret._odata[ss]=a*x._odata[ss]+b*y._odata[ss];
  });
#endif
}

template<class sobj,class vobj> inline
RealD axpy_norm(Lattice<vobj> &ret,sobj a,const Lattice<vobj> &x,const Lattice<vobj> &y){
  ret.checkerboard = x.checkerboard;
  conformable(ret,x);
  conformable(x,y);
  axpy(ret,a,x,y);
  return norm2(ret);
}
template<class sobj,class vobj> inline
RealD axpby_norm(Lattice<vobj> &ret,sobj a,sobj b,const Lattice<vobj> &x,const Lattice<vobj> &y){
  ret.checkerboard = x.checkerboard;
  conformable(ret,x);
  conformable(x,y);
  axpby(ret,a,b,x,y);
  return norm2(ret); // FIXME implement parallel norm in ss loop
}

NAMESPACE_END(Grid);
#endif
