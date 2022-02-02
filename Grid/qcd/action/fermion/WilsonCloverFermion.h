/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/WilsonCloverFermion.h

    Copyright (C) 2017 - 2022

    Author: Guido Cossu <guido.cossu@ed.ac.uk>
    Author: David Preti <>
    Author: Daniel Richtmann <daniel.richtmann@gmail.com>

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

#include <Grid/qcd/action/fermion/WilsonCloverTypes.h>
#include <Grid/qcd/action/fermion/WilsonCloverHelpers.h>

NAMESPACE_BEGIN(Grid);

///////////////////////////////////////////////////////////////////
// Wilson Clover
//
// Operator ( with anisotropy coefficients):
//
// Q =   1 + (Nd-1)/xi_0 + m
//     + W_t + (nu/xi_0) * W_s
//     - 1/2*[ csw_t * sum_s (sigma_ts F_ts) + (csw_s/xi_0) * sum_ss (sigma_ss F_ss)  ]
//
// s spatial, t temporal directions.
// where W_t and W_s are the temporal and spatial components of the
// Wilson Dirac operator
//
// csw_r = csw_t to recover the isotropic version
//////////////////////////////////////////////////////////////////

template <class Impl>
class WilsonCloverFermion : public WilsonFermion<Impl>,
                            public WilsonCloverHelpers<Impl>
{
public:
  INHERIT_IMPL_TYPES(Impl);
  INHERIT_CLOVER_TYPES(Impl);

  typedef WilsonFermion<Impl>       WilsonBase;
  typedef WilsonCloverHelpers<Impl> Helpers;

  virtual int    ConstEE(void)     { return 0; };
  virtual void Instantiatable(void){};
  // Constructors
  WilsonCloverFermion(GaugeField &_Umu, GridCartesian &Fgrid,
                      GridRedBlackCartesian &Hgrid,
                      const RealD _mass,
                      const RealD _csw_r = 0.0,
                      const RealD _csw_t = 0.0,
                      const WilsonAnisotropyCoefficients &clover_anisotropy = WilsonAnisotropyCoefficients(),
                      const ImplParams &impl_p = ImplParams());

  virtual void M(const FermionField &in, FermionField &out);
  virtual void Mdag(const FermionField &in, FermionField &out);
  virtual void Mooee(const FermionField &in, FermionField &out);
  virtual void MooeeDag(const FermionField &in, FermionField &out);
  virtual void MooeeInv(const FermionField &in, FermionField &out);
  virtual void MooeeInvDag(const FermionField &in, FermionField &out);
  virtual void MooeeInternal(const FermionField &in, FermionField &out, int dag, int inv);

  //virtual void MDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag);
  virtual void MooDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag);
  virtual void MeeDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag);

  void ImportGauge(const GaugeField &_Umu);

  // Derivative parts unpreconditioned pseudofermions
  void MDeriv(GaugeField &force, const FermionField &X, const FermionField &Y, int dag);

public:
  // here fixing the 4 dimensions, make it more general?

  RealD csw_r;                                               // Clover coefficient - spatial
  RealD csw_t;                                               // Clover coefficient - temporal
  RealD diag_mass;                                           // Mass term
  CloverField CloverTerm, CloverTermInv;                     // Clover term
  CloverField CloverTermEven, CloverTermOdd;                 // Clover term EO
  CloverField CloverTermInvEven, CloverTermInvOdd;           // Clover term Inv EO
  CloverField CloverTermDagEven, CloverTermDagOdd;           // Clover term Dag EO
  CloverField CloverTermInvDagEven, CloverTermInvDagOdd;     // Clover term Inv Dag EO
};

NAMESPACE_END(Grid);



