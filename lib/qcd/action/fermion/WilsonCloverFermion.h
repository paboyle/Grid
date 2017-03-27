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

namespace Grid {
namespace QCD {

template <class Impl>
class WilsonCloverFermion : public WilsonFermion<Impl> {
public:
  INHERIT_IMPL_TYPES(Impl);

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
                                                                                Bx(_Umu._grid),
                                                                                By(_Umu._grid),
                                                                                Bz(_Umu._grid),
                                                                                Ex(_Umu._grid),
                                                                                Ey(_Umu._grid),
                                                                                Ez(_Umu._grid)
  {
    csw = _csw;
    assert(Nd == 4); // require 4 dimensions

  }

  virtual RealD M(const FermionField& in, FermionField& out);
  virtual RealD Mdag(const FermionField& in, FermionField& out);

  virtual void Mooee(const FermionField &in, FermionField &out);
  virtual void MooeeDag(const FermionField &in, FermionField &out);
  virtual void MooeeInv(const FermionField &in, FermionField &out);
  virtual void MooeeInvDag(const FermionField &in, FermionField &out);

  void ImportGauge(const GaugeField &_Umu);
private:
  // here fixing the 4 dimensions, make it more general?

  // Field strengths
  GaugeLinkField Bx, By, Bz, Ex, Ey, Ez;

  RealD csw; // Clover coefficient


  // Methods
  void AddCloverTerm(const FermionField& in, FermionField& out);
};
}
}

#endif  // GRID_QCD_WILSON_CLOVER_FERMION_H
