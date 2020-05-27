/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/WilsonCloverMRHSFermion.h

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

#include <Grid/qcd/action/fermion/WilsonCloverFermion.h>
#include <Grid/qcd/action/fermion/WilsonMRHSFermion.h>

NAMESPACE_BEGIN(Grid);

template<class Impl>
class WilsonCloverMRHSFermion : public WilsonMRHSFermion<Impl>, public WilsonCloverFermion<Impl>
{
public:
  INHERIT_IMPL_TYPES(Impl);
  typedef WilsonMRHSFermion<Impl> WilsonMRHSBase;
  typedef WilsonCloverFermion<Impl> WilsonCloverBase;
  typedef typename WilsonCloverBase::SiteCloverType SiteCloverType;
  typedef typename WilsonCloverBase::CloverFieldType CloverFieldType;

  virtual void Instantiatable(void){};

  WilsonCloverMRHSFermion(GaugeField                         &_Umu,
                          GridCartesian                      &FiveDimGrid,
                          GridRedBlackCartesian              &FiveDimRedBlackGrid,
                          GridCartesian                      &FourDimGrid,
                          GridRedBlackCartesian              &FourDimRedBlackGrid,
                          RealD                               _mass,
                          RealD                               _csw_r = 0.,
                          RealD                               _csw_t = 0.,
                          const WilsonAnisotropyCoefficients &clover_anisotropy = WilsonAnisotropyCoefficients(),
                          const ImplParams                   &p = ImplParams())
    : WilsonMRHSBase(_Umu,
                     FiveDimGrid,
                     FiveDimRedBlackGrid,
                     FourDimGrid,
                     FourDimRedBlackGrid,
                     _mass,
                     p)
    , WilsonCloverBase(_Umu,
                       FourDimGrid,
                       FourDimRedBlackGrid,
                       _mass,
                       _csw_r,
                       _csw_t,
                       clover_anisotropy,
                       p)
  {}

  using WilsonMRHSBase::Dhop;
  using WilsonMRHSBase::DhopOE;
  using WilsonMRHSBase::DhopEO;
  using WilsonMRHSBase::DhopDir;
  using WilsonMRHSBase::DhopDirAll;
  using WilsonMRHSBase::DhopDirComms;
  using WilsonMRHSBase::DhopDirCalc;
  using WilsonMRHSBase::Meooe;
  using WilsonMRHSBase::MeooeDag;
  using WilsonMRHSBase::Mdir;
  using WilsonMRHSBase::MdirAll;

  using WilsonCloverBase::Mooee;
  using WilsonCloverBase::MooeeInv;
  using WilsonCloverBase::MooeeDag;
  using WilsonCloverBase::MooeeInvDag;

  virtual void M(const FermionField &in, FermionField &out) {
    FermionField tmp(out.Grid());

    // Wilson term
    out.Checkerboard() = in.Checkerboard();
    Dhop(in, out, DaggerNo);

    // Clover term
    Mooee(in, tmp);

    out += tmp;
  }

  virtual void Mdag(const FermionField &in, FermionField &out) {
    FermionField tmp(out.Grid());

    // Wilson term
    out.Checkerboard() = in.Checkerboard();
    Dhop(in, out, DaggerYes);

    // Clover term
    MooeeDag(in, tmp);

    out += tmp;
  }

  // emulate behavior in Clover, but do MRHS
  virtual void MooeeInternal(const FermionField &in, FermionField &out, int dag, int inv) {
    out.Checkerboard() = in.Checkerboard();
    assert(in.Checkerboard() == Odd || in.Checkerboard() == Even);

    const CloverFieldType *Clover = WilsonCloverBase::GetCompatibleCloverField(in, inv);
    assert(Clover != nullptr);

    if(dag == DaggerYes)
      WilsonCloverBase::MultClovDagInternal(*Clover,this->Ls,Clover->oSites(),in,out);
    else
      WilsonCloverBase::MultClovInternal(*Clover,this->Ls,Clover->oSites(),in,out);
  }
};

NAMESPACE_END(Grid);
