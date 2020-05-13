/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/DomainWallEOFAFermion.cc

Copyright (C) 2017

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: David Murphy <dmurphy@phys.columbia.edu>

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

#include <Grid/Grid_Eigen_Dense.h>
#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/DomainWallEOFAFermion.h>

NAMESPACE_BEGIN(Grid);

template<class Impl>
DomainWallEOFAFermion<Impl>::DomainWallEOFAFermion(
						   GaugeField            &_Umu,
						   GridCartesian         &FiveDimGrid,
						   GridRedBlackCartesian &FiveDimRedBlackGrid,
						   GridCartesian         &FourDimGrid,
						   GridRedBlackCartesian &FourDimRedBlackGrid,
						   RealD _mq1, RealD _mq2, RealD _mq3,
						   RealD _shift, int _pm, RealD _M5, const ImplParams &p) :
  AbstractEOFAFermion<Impl>(_Umu, FiveDimGrid, FiveDimRedBlackGrid,
			    FourDimGrid, FourDimRedBlackGrid, _mq1, _mq2, _mq3,
			    _shift, _pm, _M5, 1.0, 0.0, p)
{
  RealD eps = 1.0;
  Approx::zolotarev_data *zdata = Approx::higham(eps,this->Ls);
  assert(zdata->n == this->Ls);

  std::cout << GridLogMessage << "DomainWallEOFAFermion with Ls=" << this->Ls << std::endl;
  this->SetCoefficientsTanh(zdata, 1.0, 0.0);

  Approx::zolotarev_free(zdata);
}

/***************************************************************
 * Additional EOFA operators only called outside the inverter.
 * Since speed is not essential, simple axpby-style
 * implementations should be fine.
 ***************************************************************/
template<class Impl>
void DomainWallEOFAFermion<Impl>::Omega(const FermionField& psi, FermionField& Din, int sign, int dag)
{
  int Ls = this->Ls;

  Din = Zero();
  if((sign == 1) && (dag == 0)){ axpby_ssp(Din, 0.0, psi, 1.0, psi, Ls-1, 0); }
  else if((sign == -1) && (dag == 0)){ axpby_ssp(Din, 0.0, psi, 1.0, psi, 0, 0); }
  else if((sign == 1 ) && (dag == 1)){ axpby_ssp(Din, 0.0, psi, 1.0, psi, 0, Ls-1); }
  else if((sign == -1) && (dag == 1)){ axpby_ssp(Din, 0.0, psi, 1.0, psi, 0, 0); }
}

// This is just the identity for DWF
template<class Impl>
void DomainWallEOFAFermion<Impl>::Dtilde(const FermionField& psi, FermionField& chi){ chi = psi; }

// This is just the identity for DWF
template<class Impl>
void DomainWallEOFAFermion<Impl>::DtildeInv(const FermionField& psi, FermionField& chi){ chi = psi; }

/*****************************************************************************************************/

template<class Impl>
void DomainWallEOFAFermion<Impl>::M(const FermionField& psi, FermionField& chi)
{
  FermionField Din(psi.Grid());

  this->Meooe5D(psi, Din);
  this->DW(Din, chi, DaggerNo);
  axpby(chi, 1.0, 1.0, chi, psi);
  this->M5D(psi, chi);
}

template<class Impl>
void DomainWallEOFAFermion<Impl>::Mdag(const FermionField& psi, FermionField& chi)
{
  FermionField Din(psi.Grid());

  this->DW(psi, Din, DaggerYes);
  this->MeooeDag5D(Din, chi);
  this->M5Ddag(psi, chi);
  axpby(chi, 1.0, 1.0, chi, psi);
}

/********************************************************************
 * Performance critical fermion operators called inside the inverter
 ********************************************************************/

template<class Impl>
void DomainWallEOFAFermion<Impl>::M5D(const FermionField& psi, FermionField& chi)
{
  int   Ls    = this->Ls;
  int   pm    = this->pm;
  RealD shift = this->shift;
  RealD mq1   = this->mq1;
  RealD mq2   = this->mq2;
  RealD mq3   = this->mq3;

  // coefficients for shift operator ( = shift*\gamma_{5}*R_{5}*\Delta_{\pm}(mq2,mq3)*P_{\pm} )
  Coeff_t shiftp(0.0), shiftm(0.0);
  if(shift != 0.0){
    if(pm == 1){ shiftp = shift*(mq3-mq2); }
    else{ shiftm = -shift*(mq3-mq2); }
  }

  Vector<Coeff_t> diag(Ls,1.0);
  Vector<Coeff_t> upper(Ls,-1.0); upper[Ls-1] = mq1 + shiftm;
  Vector<Coeff_t> lower(Ls,-1.0); lower[0]    = mq1 + shiftp;

#if(0)
  std::cout << GridLogMessage << "DomainWallEOFAFermion::M5D(FF&,FF&):" << std::endl;
  for(int i=0; i<diag.size(); ++i){
    std::cout << GridLogMessage << "diag[" << i << "] =" << diag[i] << std::endl;
  }
  for(int i=0; i<upper.size(); ++i){
    std::cout << GridLogMessage << "upper[" << i << "] =" << upper[i] << std::endl;
  }
  for(int i=0; i<lower.size(); ++i){
    std::cout << GridLogMessage << "lower[" << i << "] =" << lower[i] << std::endl;
  }
#endif

  this->M5D(psi, chi, chi, lower, diag, upper);
}

template<class Impl>
void DomainWallEOFAFermion<Impl>::M5Ddag(const FermionField& psi, FermionField& chi)
{
  int   Ls    = this->Ls;
  int   pm    = this->pm;
  RealD shift = this->shift;
  RealD mq1   = this->mq1;
  RealD mq2   = this->mq2;
  RealD mq3   = this->mq3;

  // coefficients for shift operator ( = shift*\gamma_{5}*R_{5}*\Delta_{\pm}(mq2,mq3)*P_{\pm} )
  Coeff_t shiftp(0.0), shiftm(0.0);
  if(shift != 0.0){
    if(pm == 1){ shiftp = shift*(mq3-mq2); }
    else{ shiftm = -shift*(mq3-mq2); }
  }

  Vector<Coeff_t> diag(Ls,1.0);
  Vector<Coeff_t> upper(Ls,-1.0); upper[Ls-1] = mq1 + shiftp;
  Vector<Coeff_t> lower(Ls,-1.0); lower[0]    = mq1 + shiftm;

  this->M5Ddag(psi, chi, chi, lower, diag, upper);
}

// half checkerboard operations
template<class Impl>
void DomainWallEOFAFermion<Impl>::Mooee(const FermionField& psi, FermionField& chi)
{
  int Ls = this->Ls;

  Vector<Coeff_t> diag = this->bee;
  Vector<Coeff_t> upper(Ls);
  Vector<Coeff_t> lower(Ls);

  for(int s=0; s<Ls; s++){
    upper[s] = -this->cee[s];
    lower[s] = -this->cee[s];
  }
  upper[Ls-1] = this->dm;
  lower[0]    = this->dp;

  this->M5D(psi, psi, chi, lower, diag, upper);
}

template<class Impl>
void DomainWallEOFAFermion<Impl>::MooeeDag(const FermionField& psi, FermionField& chi)
{
  int Ls = this->Ls;

  Vector<Coeff_t> diag = this->bee;
  Vector<Coeff_t> upper(Ls);
  Vector<Coeff_t> lower(Ls);

  for(int s=0; s<Ls; s++){
    upper[s] = -this->cee[s];
    lower[s] = -this->cee[s];
  }
  upper[Ls-1] = this->dp;
  lower[0]    = this->dm;

  this->M5Ddag(psi, psi, chi, lower, diag, upper);
}

/****************************************************************************************/

//Zolo
template<class Impl>
void DomainWallEOFAFermion<Impl>::SetCoefficientsInternal(RealD zolo_hi, Vector<Coeff_t>& gamma, RealD b, RealD c)
{
  int   Ls    = this->Ls;
  int   pm    = this->pm;
  RealD mq1   = this->mq1;
  RealD mq2   = this->mq2;
  RealD mq3   = this->mq3;
  RealD shift = this->shift;

  ////////////////////////////////////////////////////////
  // Constants for the preconditioned matrix Cayley form
  ////////////////////////////////////////////////////////
  this->bs.resize(Ls);
  this->cs.resize(Ls);
  this->aee.resize(Ls);
  this->aeo.resize(Ls);
  this->bee.resize(Ls);
  this->beo.resize(Ls);
  this->cee.resize(Ls);
  this->ceo.resize(Ls);

  for(int i=0; i<Ls; ++i){
    this->bee[i] = 4.0 - this->M5 + 1.0;
    this->cee[i] = 1.0;
  }

  for(int i=0; i<Ls; ++i){
    this->aee[i] = this->cee[i];
    this->bs[i] = this->beo[i] = 1.0;
    this->cs[i] = this->ceo[i] = 0.0;
  }

  //////////////////////////////////////////
  // EOFA shift terms
  //////////////////////////////////////////
  if(pm == 1){
    this->dp = mq1*this->cee[0] + shift*(mq3-mq2);
    this->dm = mq1*this->cee[Ls-1];
  } else if(this->pm == -1) {
    this->dp = mq1*this->cee[0];
    this->dm = mq1*this->cee[Ls-1] - shift*(mq3-mq2);
  } else {
    this->dp = mq1*this->cee[0];
    this->dm = mq1*this->cee[Ls-1];
  }

  //////////////////////////////////////////
  // LDU decomposition of eeoo
  //////////////////////////////////////////
  this->dee.resize(Ls+1);
  this->lee.resize(Ls);
  this->leem.resize(Ls);
  this->uee.resize(Ls);
  this->ueem.resize(Ls);

  for(int i=0; i<Ls; ++i){

    if(i < Ls-1){

      this->lee[i] = -this->cee[i+1]/this->bee[i]; // sub-diag entry on the ith column

      this->leem[i] = this->dm/this->bee[i];
      for(int j=0; j<i; j++){ this->leem[i] *= this->aee[j]/this->bee[j]; }

      this->dee[i] = this->bee[i];

      this->uee[i] = -this->aee[i]/this->bee[i];   // up-diag entry on the ith row

      this->ueem[i] = this->dp / this->bee[0];
      for(int j=1; j<=i; j++){ this->ueem[i] *= this->cee[j]/this->bee[j]; }

    } else {

      this->lee[i]  = 0.0;
      this->leem[i] = 0.0;
      this->uee[i]  = 0.0;
      this->ueem[i] = 0.0;

    }
  }

  {
    Coeff_t delta_d = 1.0 / this->bee[0];
    for(int j=1; j<Ls-1; j++){ delta_d *= this->cee[j] / this->bee[j]; }
    this->dee[Ls-1] = this->bee[Ls-1] + this->cee[0] * this->dm * delta_d;
    this->dee[Ls] = this->bee[Ls-1] + this->cee[Ls-1] * this->dp * delta_d;
  }
}

// Recompute Cayley-form coefficients for different shift
template<class Impl>
void DomainWallEOFAFermion<Impl>::RefreshShiftCoefficients(RealD new_shift)
{
  this->shift = new_shift;
  Approx::zolotarev_data *zdata = Approx::higham(1.0, this->Ls);
  this->SetCoefficientsTanh(zdata, 1.0, 0.0);
}

NAMESPACE_END(Grid);
