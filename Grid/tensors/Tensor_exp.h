/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_exp.h

    Copyright (C) 2015

Author: neo <cossu@post.kek.jp>

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
#ifndef GRID_MATH_EXP_H
#define GRID_MATH_EXP_H

#define DEFAULT_MAT_EXP 20

NAMESPACE_BEGIN(Grid);

/////////////////////////////////////////////// 
// Exponentiate function for scalar, vector, matrix
/////////////////////////////////////////////// 


template<class vtype> accelerator_inline iScalar<vtype> Exponentiate(const iScalar<vtype>&r, RealD alpha ,  Integer Nexp = DEFAULT_MAT_EXP)
{
  iScalar<vtype> ret;
  ret._internal = Exponentiate(r._internal, alpha, Nexp);
  return ret;
}

template<class vtype, int N> accelerator_inline iVector<vtype, N> Exponentiate(const iVector<vtype,N>&r, RealD alpha ,  Integer Nexp = DEFAULT_MAT_EXP)
{
  iVector<vtype, N> ret;
  for (int i = 0; i < N; i++)
    ret._internal[i] = Exponentiate(r._internal[i], alpha, Nexp);
  return ret;
}



// Specialisation: Cayley-Hamilton exponential for SU(3)
#ifndef GRID_CUDA
template<class vtype, typename std::enable_if< GridTypeMapper<vtype>::TensorLevel == 0>::type * =nullptr> 
accelerator_inline iMatrix<vtype,3> Exponentiate(const iMatrix<vtype,3> &arg, RealD alpha  , Integer Nexp = DEFAULT_MAT_EXP )
{
  // for SU(3) 2x faster than the std implementation using Nexp=12
  // notice that it actually computes
  // exp ( input matrix )
  // the i sign is coming from outside
  // input matrix is anti-hermitian NOT hermitian
  typedef iMatrix<vtype,3> mat;
  typedef iScalar<vtype> scalar;
  mat unit(1.0);
  const Complex one_over_three = 1.0 / 3.0;
  const Complex one_over_two = 1.0 / 2.0;

  scalar c0, c1, tmp, c0max, theta, u, w;
  scalar xi0, u2, w2, cosw;
  scalar fden, h0, h1, h2;
  scalar e2iu, emiu, ixi0;
  scalar f0, f1, f2;
  scalar unity(1.0);
      
  mat iQ2 = arg*arg*alpha*alpha;
  mat iQ3 = arg*iQ2*alpha;   
  // sign in c0 from the conventions on the Ta
  scalar imQ3, reQ2;
  imQ3 = imag( trace(iQ3) );
  reQ2 = real( trace(iQ2) );
  c0 = -imQ3 * one_over_three;  
  c1 = -reQ2 * one_over_two;

  // Cayley Hamilton checks to machine precision, tested
  tmp = c1 * one_over_three;
  c0max = 2.0 * pow(tmp, 1.5);

  theta = acos(c0 / c0max) * one_over_three;
  u = sqrt(tmp) * cos(theta);
  w = sqrt(c1) * sin(theta);

  xi0 = sin(w) / w;
  u2 = u * u;
  w2 = w * w;
  cosw = cos(w);

  ixi0 = timesI(xi0);
  emiu = cos(u) - timesI(sin(u));
  e2iu = cos(2.0 * u) + timesI(sin(2.0 * u));

  h0 = e2iu * (u2 - w2) +
    emiu * ((8.0 * u2 * cosw) + (2.0 * u * (3.0 * u2 + w2) * ixi0));
  h1 = e2iu * (2.0 * u) - emiu * ((2.0 * u * cosw) - (3.0 * u2 - w2) * ixi0);
  h2 = e2iu - emiu * (cosw + (3.0 * u) * ixi0);

  fden = unity / (9.0 * u2 - w2);  // reals
  f0 = h0 * fden;
  f1 = h1 * fden;
  f2 = h2 * fden;

  return (f0 * unit + timesMinusI(f1) * arg*alpha - f2 * iQ2);
}
#endif


// General exponential
template<class vtype,int N, typename std::enable_if< GridTypeMapper<vtype>::TensorLevel == 0 >::type * =nullptr> 
accelerator_inline iMatrix<vtype,N> Exponentiate(const iMatrix<vtype,N> &arg, RealD alpha  , Integer Nexp = DEFAULT_MAT_EXP )
{
  // notice that it actually computes
  // exp ( input matrix )
  // the i sign is coming from outside
  // input matrix is anti-hermitian NOT hermitian
  typedef iMatrix<vtype,N> mat;
  mat unit(1.0);
  mat temp(unit);
  for(int i=Nexp; i>=1;--i){
    temp *= alpha/RealD(i);
    temp = unit + temp*arg;
  }
  return temp;

}

NAMESPACE_END(Grid);

#endif
