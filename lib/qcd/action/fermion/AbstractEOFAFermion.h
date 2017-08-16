/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/AbstractEOFAFermion.h

Copyright (C) 2017

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#ifndef  GRID_QCD_ABSTRACT_EOFA_FERMION_H
#define  GRID_QCD_ABSTRACT_EOFA_FERMION_H

#include <Grid/qcd/action/fermion/CayleyFermion5D.h>

namespace Grid {
namespace QCD {

  // DJM: Abstract base class for EOFA fermion types.
  // Defines layout of additional EOFA-specific parameters and operators.
  // Use to construct EOFA pseudofermion actions that are agnostic to
  // Shamir / Mobius / etc., and ensure that no one can construct EOFA
  // pseudofermion action with non-EOFA fermion type.
  template<class Impl>
  class AbstractEOFAFermion : public CayleyFermion5D<Impl> {
    public:
      INHERIT_IMPL_TYPES(Impl);

    public:
      // Fermion operator: D(mq1) + shift*\gamma_{5}*R_{5}*\Delta_{\pm}(mq2,mq3)*P_{\pm}
      RealD mq1;
      RealD mq2;
      RealD mq3;
      RealD shift;
      int pm;

      RealD alpha; // Mobius scale
      RealD k;     // EOFA normalization constant

      virtual void Instantiatable(void) = 0;

      // EOFA-specific operations
      // Force user to implement in derived classes
      virtual void  Omega    (const FermionField& in, FermionField& out, int sign, int dag) = 0;
      virtual void  Dtilde   (const FermionField& in, FermionField& out) = 0;
      virtual void  DtildeInv(const FermionField& in, FermionField& out) = 0;

      // Implement derivatives in base class:
      // for EOFA both DWF and Mobius just need d(Dw)/dU
      virtual void MDeriv(GaugeField& mat, const FermionField& U, const FermionField& V, int dag){
        this->DhopDeriv(mat, U, V, dag);
      };
      virtual void MoeDeriv(GaugeField& mat, const FermionField& U, const FermionField& V, int dag){
        this->DhopDerivOE(mat, U, V, dag);
      };
      virtual void MeoDeriv(GaugeField& mat, const FermionField& U, const FermionField& V, int dag){
        this->DhopDerivEO(mat, U, V, dag);
      };

      // Recompute 5D coefficients for different value of shift constant
      // (needed for heatbath loop over poles)
      virtual void RefreshShiftCoefficients(RealD new_shift) = 0;

      // Constructors
      AbstractEOFAFermion(GaugeField& _Umu, GridCartesian& FiveDimGrid, GridRedBlackCartesian& FiveDimRedBlackGrid,
        GridCartesian& FourDimGrid, GridRedBlackCartesian& FourDimRedBlackGrid,
        RealD _mq1, RealD _mq2, RealD _mq3, RealD _shift, int _pm,
        RealD _M5, RealD _b, RealD _c, const ImplParams& p=ImplParams())
        : CayleyFermion5D<Impl>(_Umu, FiveDimGrid, FiveDimRedBlackGrid, FourDimGrid, FourDimRedBlackGrid,
          _mq1, _M5, p), mq1(_mq1), mq2(_mq2), mq3(_mq3), shift(_shift), pm(_pm)
      {
        int Ls = this->Ls;
        this->alpha = _b + _c;
        this->k = this->alpha * (_mq3-_mq2) * std::pow(this->alpha+1.0,2*Ls) /
                    ( std::pow(this->alpha+1.0,Ls) + _mq2*std::pow(this->alpha-1.0,Ls) ) /
                    ( std::pow(this->alpha+1.0,Ls) + _mq3*std::pow(this->alpha-1.0,Ls) );
      };
  };
}}

#endif
