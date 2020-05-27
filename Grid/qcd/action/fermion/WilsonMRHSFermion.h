/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/WilsonMRHSFermion.h

    Copyright (C) 2015 - 2020

Author: Daniel Richtmann <daniel.richtmann@ur.de>

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

#include <Grid/qcd/action/fermion/WilsonFermion5D.h>

NAMESPACE_BEGIN(Grid);

template<class Impl>
class WilsonMRHSFermion : public WilsonFermion5D<Impl>
{
public:
  INHERIT_IMPL_TYPES(Impl);
  typedef WilsonFermion5D<Impl> WilsonFermion5DBase;

  virtual void Instantiatable(void){};

  WilsonMRHSFermion(GaugeField            &_Umu,
                    GridCartesian         &FiveDimGrid,
                    GridRedBlackCartesian &FiveDimRedBlackGrid,
                    GridCartesian         &FourDimGrid,
                    GridRedBlackCartesian &FourDimRedBlackGrid,
                    RealD                 _mass,
                    const ImplParams      &p = ImplParams())
    : WilsonFermion5D<Impl>(_Umu,
                            FiveDimGrid,
                            FiveDimRedBlackGrid,
                            FourDimGrid,
                            FourDimRedBlackGrid,
                            -_mass,
                            p)
    , mass(_mass)
    , diag_mass(4.0 + mass)
  {}

  using WilsonFermion5DBase::Dhop;
  using WilsonFermion5DBase::DhopOE;
  using WilsonFermion5DBase::DhopEO;
  using WilsonFermion5DBase::DhopDir;
  using WilsonFermion5DBase::DhopDirAll;
  using WilsonFermion5DBase::DhopDirComms;
  using WilsonFermion5DBase::DhopDirCalc;

  virtual void M(const FermionField &in, FermionField &out) {
    out.Checkerboard() = in.Checkerboard();
    Dhop(in, out, DaggerNo);
    axpy(out, diag_mass, in, out);
  }

  virtual void Mdag(const FermionField &in, FermionField &out) {
    out.Checkerboard() = in.Checkerboard();
    Dhop(in, out, DaggerYes);
    axpy(out, diag_mass, in, out);
  }

  virtual void Meooe(const FermionField &in, FermionField &out) {
    if (in.Checkerboard() == Odd) {
      DhopEO(in, out, DaggerNo);
    } else {
      DhopOE(in, out, DaggerNo);
    }
  }

  virtual void Mooee(const FermionField &in, FermionField &out) {
    out.Checkerboard() = in.Checkerboard();
    typename FermionField::scalar_type scal(diag_mass);
    out = scal * in;
  }

  virtual void MooeeInv(const FermionField &in, FermionField &out) {
    out.Checkerboard() = in.Checkerboard();
    out = (1.0/(diag_mass))*in;
  }

  virtual void MeooeDag(const FermionField &in, FermionField &out) {
    if (in.Checkerboard() == Odd) {
      DhopEO(in, out, DaggerYes);
    } else {
      DhopOE(in, out, DaggerYes);
    }
  }

  virtual void MooeeDag(const FermionField &in, FermionField &out) {
    out.Checkerboard() = in.Checkerboard();
    Mooee(in, out);
  }

  virtual void MooeeInvDag(const FermionField &in, FermionField &out) {
    out.Checkerboard() = in.Checkerboard();
    MooeeInv(in,out);
  }

  virtual void Mdir(const FermionField &in, FermionField &out, int dir, int disp) {
    DhopDir(in, out, dir+1, disp); // retain conventions of 4d version
  }

  virtual void MdirAll(const FermionField &in, std::vector<FermionField> &out) {
    DhopDirAll(in, out);
  }

private:
  RealD mass;
  RealD diag_mass;
};

NAMESPACE_END(Grid);
