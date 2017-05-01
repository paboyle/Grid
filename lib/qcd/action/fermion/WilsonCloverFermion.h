/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/WilsonCloverFermion.h

    Copyright (C) 2017

    Author: paboyle <paboyle@ph.ed.ac.uk>
    Author: Guido Cossu <guido.cossu@ed.ac.uk>

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
#ifndef GRID_QCD_WILSON_CLOVER_FERMION_H
#define GRID_QCD_WILSON_CLOVER_FERMION_H

#include <Grid/Grid.h>

namespace Grid
{
namespace QCD
{

template <class Impl>
class WilsonCloverFermion : public WilsonFermion<Impl>
{
public:
  // Types definitions
  INHERIT_IMPL_TYPES(Impl);
  template <typename vtype> using iImplClover        = iScalar<iMatrix<iMatrix<vtype, Impl::Dimension>, Ns> >;
  typedef iImplClover<Simd>        SiteCloverType;
  typedef Lattice<SiteCloverType>        CloverFieldType;
public:
  typedef WilsonFermion<Impl> WilsonBase;

  virtual void Instantiatable(void){};
  // Constructors
  WilsonCloverFermion(GaugeField &_Umu, GridCartesian &Fgrid,
                      GridRedBlackCartesian &Hgrid,
                      RealD _mass,
                      RealD _csw,
                      const ImplParams &p = ImplParams()) : WilsonFermion<Impl>(_Umu,
                                                                                Fgrid,
                                                                                Hgrid,
                                                                                _mass, p),
                                                                                CloverTerm(&Fgrid),
                                                                                CloverTermInv(&Fgrid),
                                                                                CloverTermEven(&Hgrid),
                                                                                CloverTermOdd(&Hgrid),
                                                                                CloverTermInvEven(&Hgrid),
                                                                                CloverTermInvOdd(&Hgrid)                                                                                
  {
    csw = _csw;
    assert(Nd == 4); // require 4 dimensions
  }

  virtual RealD M(const FermionField &in, FermionField &out);
  virtual RealD Mdag(const FermionField &in, FermionField &out);

  virtual void Mooee(const FermionField &in, FermionField &out);
  virtual void MooeeDag(const FermionField &in, FermionField &out);
  virtual void MooeeInv(const FermionField &in, FermionField &out);
  virtual void MooeeInvDag(const FermionField &in, FermionField &out);
  virtual void MooeeInternal(const FermionField &in, FermionField &out, int dag, int inv);

  virtual void MDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag);
  virtual void MooDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag);
  virtual void MeeDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag);

  void ImportGauge(const GaugeField &_Umu);

private:
  // here fixing the 4 dimensions, make it more general?

  RealD csw;                                         // Clover coefficient
  CloverFieldType CloverTerm, CloverTermInv; // Clover term
  CloverFieldType CloverTermEven, CloverTermOdd;
  CloverFieldType CloverTermInvEven, CloverTermInvOdd; // Clover term
  // eventually these two can be compressed into 6x6 blocks instead of the 12x12
  // using the DeGrand-Rossi basis for the gamma matrices

  CloverFieldType fillClover(const GaugeLinkField& F){
    CloverFieldType T(F._grid);
    PARALLEL_FOR_LOOP
    for (int i = 0; i < CloverTerm._grid->oSites(); i++){
      for (int s1 = 0; s1 < Nc; s1++)
      for (int s2 = 0; s2 < Nc; s2++)
      T._odata[i]()(s1,s2) = F._odata[i]()();
    }
  return T;
  }
  
};
}
}

#endif  // GRID_QCD_WILSON_CLOVER_FERMION_H
