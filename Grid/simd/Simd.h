/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/Simd.h

Copyright (C) 2015

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef GRID_SIMD_H
#define GRID_SIMD_H

#if defined(GRID_CUDA) || defined(GRID_HIP)
#include <thrust/complex.h>
#endif

////////////////////////////////////////////////////////////////////////
// Define scalar and vector floating point types
//
// Scalar:   RealF, RealD, ComplexF, ComplexD
//
// Vector:  vRealF, vRealD, vComplexF, vComplexD
//
// Vector types are arch dependent
////////////////////////////////////////////////////////////////////////

#define _MM_SELECT_FOUR_FOUR(A,B,C,D) ((A<<6)|(B<<4)|(C<<2)|(D))
#define _MM_SELECT_FOUR_FOUR_STRING(A,B,C,D) "((" #A "<<6)|(" #B "<<4)|(" #C "<<2)|(" #D "))"
#define _MM_SELECT_EIGHT_TWO(A,B,C,D,E,F,G,H) ((A<<7)|(B<<6)|(C<<5)|(D<<4)|(E<<3)|(F<<2)|(G<<4)|(H))
#define _MM_SELECT_FOUR_TWO (A,B,C,D) _MM_SELECT_EIGHT_TWO(0,0,0,0,A,B,C,D)
#define _MM_SELECT_TWO_TWO  (A,B)     _MM_SELECT_FOUR_TWO(0,0,A,B)

#define RotateBit (0x100)

NAMESPACE_BEGIN(Grid);

typedef uint32_t Integer;

typedef  float  RealF;
typedef  double RealD;
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
typedef RealD   Real;
#else
typedef RealF  Real;
#endif

#if defined(GRID_CUDA) || defined(GRID_HIP)
typedef thrust::complex<RealF> ComplexF;
typedef thrust::complex<RealD> ComplexD;
typedef thrust::complex<Real>  Complex;
typedef thrust::complex<uint16_t>  ComplexH;
template<class T> using complex = thrust::complex<T>;

accelerator_inline ComplexD pow(const ComplexD& r,RealD y){ return(thrust::pow(r,(double)y)); }
accelerator_inline ComplexF pow(const ComplexF& r,RealF y){ return(thrust::pow(r,(float)y)); }
#else 
typedef std::complex<RealF> ComplexF;
typedef std::complex<RealD> ComplexD;
typedef std::complex<Real>  Complex;
typedef std::complex<uint16_t>  ComplexH; // Hack
template<class T> using complex = std::complex<T>;

accelerator_inline ComplexD pow(const ComplexD& r,RealD y){ return(std::pow(r,y)); }
accelerator_inline ComplexF pow(const ComplexF& r,RealF y){ return(std::pow(r,y)); }
#endif

//accelerator_inline RealD pow(const RealD& r,RealD y){ return(std::pow(r,y)); }
//accelerator_inline RealD sqrt(const RealD  & r){ return std::sqrt(r); }

// This comes from ::pow already from math.h and CUDA
// Calls either Grid::pow for complex, or std::pow for real
// Problem is CUDA math_functions is exposing ::pow, and I can't define

using std::abs;
using std::pow;
using std::sqrt;
using std::log;
using std::exp;
using std::sin;
using std::cos;
using std::asin;
using std::acos;


accelerator_inline RealF    conjugate(const RealF  & r){ return r; }
accelerator_inline RealD    conjugate(const RealD  & r){ return r; }
accelerator_inline ComplexD conjugate(const ComplexD& r){ return(conj(r)); }
accelerator_inline ComplexF conjugate(const ComplexF& r ){ return(conj(r)); }

accelerator_inline RealF    adj(const RealF  & r){ return r; }
accelerator_inline RealD    adj(const RealD  & r){ return r; }
accelerator_inline ComplexD adj(const ComplexD& r){ return(conjugate(r)); }
accelerator_inline ComplexF adj(const ComplexF& r ){ return(conjugate(r)); }

accelerator_inline RealF real(const RealF  & r){ return r; }
accelerator_inline RealD real(const RealD  & r){ return r; }
accelerator_inline RealF real(const ComplexF  & r){ return r.real(); }
accelerator_inline RealD real(const ComplexD  & r){ return r.real(); }

accelerator_inline RealF imag(const ComplexF  & r){ return r.imag(); }
accelerator_inline RealD imag(const ComplexD  & r){ return r.imag(); }

accelerator_inline ComplexD innerProduct(const ComplexD & l, const ComplexD & r) { return conjugate(l)*r; }
accelerator_inline ComplexF innerProduct(const ComplexF & l, const ComplexF & r) { return conjugate(l)*r; }
accelerator_inline RealD innerProduct(const RealD & l, const RealD & r) { return l*r; }
accelerator_inline RealF innerProduct(const RealF & l, const RealF & r) { return l*r; }

accelerator_inline ComplexD Reduce(const ComplexD& r){ return r; }
accelerator_inline ComplexF Reduce(const ComplexF& r){ return r; }
accelerator_inline RealD Reduce(const RealD& r){ return r; }
accelerator_inline RealF Reduce(const RealF& r){ return r; }

accelerator_inline RealD toReal(const ComplexD& r){ return r.real(); }
accelerator_inline RealF toReal(const ComplexF& r){ return r.real(); }
accelerator_inline RealD toReal(const RealD& r){ return r; }
accelerator_inline RealF toReal(const RealF& r){ return r; }
  
////////////////////////////////////////////////////////////////////////////////
//Provide support functions for basic real and complex data types required by Grid
//Single and double precision versions. Should be able to template this once only.
////////////////////////////////////////////////////////////////////////////////
accelerator_inline void mac (ComplexD * __restrict__ y,const ComplexD * __restrict__ a,const ComplexD *__restrict__ x){ *y = (*a) * (*x)+(*y); };
accelerator_inline void mult(ComplexD * __restrict__ y,const ComplexD * __restrict__ l,const ComplexD *__restrict__ r){ *y = (*l) * (*r);}
accelerator_inline void sub (ComplexD * __restrict__ y,const ComplexD * __restrict__ l,const ComplexD *__restrict__ r){ *y = (*l) - (*r);}
accelerator_inline void add (ComplexD * __restrict__ y,const ComplexD * __restrict__ l,const ComplexD *__restrict__ r){ *y = (*l) + (*r);}
// conjugate already supported for complex
  
accelerator_inline void mac (ComplexF * __restrict__ y,const ComplexF * __restrict__ a,const ComplexF *__restrict__ x){ *y = (*a) * (*x)+(*y); }
accelerator_inline void mult(ComplexF * __restrict__ y,const ComplexF * __restrict__ l,const ComplexF *__restrict__ r){ *y = (*l) * (*r); }
accelerator_inline void sub (ComplexF * __restrict__ y,const ComplexF * __restrict__ l,const ComplexF *__restrict__ r){ *y = (*l) - (*r); }
accelerator_inline void add (ComplexF * __restrict__ y,const ComplexF * __restrict__ l,const ComplexF *__restrict__ r){ *y = (*l) + (*r); }
  
//conjugate already supported for complex
accelerator_inline ComplexF timesI(const ComplexF &r)     { return(ComplexF(-r.imag(),r.real()));}
accelerator_inline ComplexD timesI(const ComplexD &r)     { return(ComplexD(-r.imag(),r.real()));}
accelerator_inline ComplexF timesMinusI(const ComplexF &r){ return(ComplexF(r.imag(),-r.real()));}
accelerator_inline ComplexD timesMinusI(const ComplexD &r){ return(ComplexD(r.imag(),-r.real()));}
//accelerator_inline ComplexF timesI(const ComplexF &r)     { return(r*ComplexF(0.0,1.0));}
//accelerator_inline ComplexD timesI(const ComplexD &r)     { return(r*ComplexD(0.0,1.0));}
//accelerator_inline ComplexF timesMinusI(const ComplexF &r){ return(r*ComplexF(0.0,-1.0));}
//accelerator_inline ComplexD timesMinusI(const ComplexD &r){ return(r*ComplexD(0.0,-1.0));}

// define projections to real and imaginay parts
accelerator_inline ComplexF projReal(const ComplexF &r){return( ComplexF(r.real(), 0.0));}
accelerator_inline ComplexD projReal(const ComplexD &r){return( ComplexD(r.real(), 0.0));}
accelerator_inline ComplexF projImag(const ComplexF &r){return (ComplexF(r.imag(), 0.0 ));}
accelerator_inline ComplexD projImag(const ComplexD &r){return (ComplexD(r.imag(), 0.0));}

// define auxiliary functions for complex computations
accelerator_inline void timesI(ComplexF &ret,const ComplexF &r)     { ret = timesI(r);}
accelerator_inline void timesI(ComplexD &ret,const ComplexD &r)     { ret = timesI(r);}
accelerator_inline void timesMinusI(ComplexF &ret,const ComplexF &r){ ret = timesMinusI(r);}
accelerator_inline void timesMinusI(ComplexD &ret,const ComplexD &r){ ret = timesMinusI(r);}
  
accelerator_inline void mac (RealD * __restrict__ y,const RealD * __restrict__ a,const RealD *__restrict__ x){ *y = (*a) * (*x)+(*y);}
accelerator_inline void mult(RealD * __restrict__ y,const RealD * __restrict__ l,const RealD *__restrict__ r){ *y = (*l) * (*r);}
accelerator_inline void sub (RealD * __restrict__ y,const RealD * __restrict__ l,const RealD *__restrict__ r){ *y = (*l) - (*r);}
accelerator_inline void add (RealD * __restrict__ y,const RealD * __restrict__ l,const RealD *__restrict__ r){ *y = (*l) + (*r);}
  
accelerator_inline void mac (RealF * __restrict__ y,const RealF * __restrict__ a,const RealF *__restrict__ x){  *y = (*a) * (*x)+(*y); }
accelerator_inline void mult(RealF * __restrict__ y,const RealF * __restrict__ l,const RealF *__restrict__ r){ *y = (*l) * (*r); }
accelerator_inline void sub (RealF * __restrict__ y,const RealF * __restrict__ l,const RealF *__restrict__ r){ *y = (*l) - (*r); }
accelerator_inline void add (RealF * __restrict__ y,const RealF * __restrict__ l,const RealF *__restrict__ r){ *y = (*l) + (*r); }
  
accelerator_inline void vstream(ComplexF &l, const ComplexF &r){ l=r;}
accelerator_inline void vstream(ComplexD &l, const ComplexD &r){ l=r;}
accelerator_inline void vstream(RealF &l, const RealF &r){ l=r;}
accelerator_inline void vstream(RealD &l, const RealD &r){ l=r;}
  
accelerator_inline ComplexD toComplex(const RealD &in) { return ComplexD(in);}
accelerator_inline ComplexF toComplex(const RealF &in) { return ComplexF(in);}
  
class Zero{};
//static Zero Zero();
template<class itype> accelerator_inline void zeroit(itype &arg)   { arg=Zero();};
template<>            accelerator_inline void zeroit(ComplexF &arg){ arg=0; };
template<>            accelerator_inline void zeroit(ComplexD &arg){ arg=0; };
template<>            accelerator_inline void zeroit(RealF &arg)   { arg=0; };
template<>            accelerator_inline void zeroit(RealD &arg)   { arg=0; };

// More limited Integer support  
accelerator_inline Integer Reduce(const Integer& r){ return r; }
accelerator_inline void mac (Integer * __restrict__ y,const Integer * __restrict__ a,const Integer *__restrict__ x){  *y = (*a) * (*x)+(*y); }
accelerator_inline void mult(Integer * __restrict__ y,const Integer * __restrict__ l,const Integer *__restrict__ r){ *y = (*l) * (*r); }
accelerator_inline void sub (Integer * __restrict__ y,const Integer * __restrict__ l,const Integer *__restrict__ r){ *y = (*l) - (*r); }
accelerator_inline void add (Integer * __restrict__ y,const Integer * __restrict__ l,const Integer *__restrict__ r){ *y = (*l) + (*r); }
accelerator_inline void vstream(Integer &l, const RealD &r){ l=r;}
template<>            accelerator_inline void zeroit(Integer &arg)   { arg=0; };

accelerator_inline Integer mod (Integer a,Integer y) { return a%y;}
accelerator_inline Integer div (Integer a,Integer y) { return a/y;}
//accelerator_inline Integer abs (Integer &a) { return a%y;}

//////////////////////////////////////////////////////////
// Permute
// Permute 0 every ABCDEFGH -> BA DC FE HG
// Permute 1 every ABCDEFGH -> CD AB GH EF
// Permute 2 every ABCDEFGH -> EFGH ABCD
// Permute 3 possible on longer iVector lengths (512bit = 8 double = 16 single)
// Permute 4 possible on half precision @512bit vectors.
//
// Defined inside SIMD specialization files
//////////////////////////////////////////////////////////
template<class VectorSIMD>
accelerator_inline void Gpermute(VectorSIMD &y,const VectorSIMD &b,int perm);

NAMESPACE_END(Grid);

#include <Grid/simd/Grid_vector_types.h>
#include <Grid/simd/Grid_doubled_vector.h>
#include <Grid/simd/Grid_vector_unops.h>

NAMESPACE_BEGIN(Grid);

// Default precision is wired to double
typedef vRealD vReal;
typedef vComplexD vComplex;
 
inline std::ostream& operator<< (std::ostream& stream, const vComplexF &o){
  int nn=vComplexF::Nsimd();
  std::vector<ComplexF,alignedAllocator<ComplexF> > buf(nn);
  vstore(o,&buf[0]);
  stream<<"<";
  for(int i=0;i<nn;i++){
    stream<<buf[i];
    if(i<nn-1) stream<<",";
  }
  stream<<">";
  return stream;
}
 
inline std::ostream& operator<< (std::ostream& stream, const vComplexD &o){
  int nn=vComplexD::Nsimd();
  std::vector<ComplexD,alignedAllocator<ComplexD> > buf(nn);
  vstore(o,&buf[0]);
  stream<<"<";
  for(int i=0;i<nn;i++){
    stream<<buf[i];
    if(i<nn-1) stream<<",";
  }
  stream<<">";
  return stream;
}
inline std::ostream& operator<< (std::ostream& stream, const vComplexD2 &o){
  stream<<"<";
  stream<<o.v[0];
  stream<<o.v[1];
  stream<<">";
  return stream;
}

inline std::ostream& operator<< (std::ostream& stream, const vRealF &o){
  int nn=vRealF::Nsimd();
  std::vector<RealF,alignedAllocator<RealF> > buf(nn);
  vstore(o,&buf[0]);
  stream<<"<";
  for(int i=0;i<nn;i++){
    stream<<buf[i];
    if(i<nn-1) stream<<",";
  }
  stream<<">";
  return stream;
}

inline std::ostream& operator<< (std::ostream& stream, const vRealD &o){
  int nn=vRealD::Nsimd();
  std::vector<RealD,alignedAllocator<RealD> > buf(nn);
  vstore(o,&buf[0]);
  stream<<"<";
  for(int i=0;i<nn;i++){
    stream<<buf[i];
    if(i<nn-1) stream<<",";
  }
  stream<<">";
  return stream;
}
inline std::ostream& operator<< (std::ostream& stream, const vInteger &o){
  int nn=vInteger::Nsimd();
  std::vector<Integer,alignedAllocator<Integer> > buf(nn);
  vstore(o,&buf[0]);
  stream<<"<";
  for(int i=0;i<nn;i++){
    stream<<buf[i];
    if(i<nn-1) stream<<",";
  }
  stream<<">";
  return stream;
}

NAMESPACE_END(Grid)
#endif
