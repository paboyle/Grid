/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/MobiusEOFAFermion.h

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
#ifndef  GRID_QCD_MOBIUS_EOFA_FERMION_H
#define  GRID_QCD_MOBIUS_EOFA_FERMION_H

#include <Grid/qcd/action/fermion/AbstractEOFAFermion.h>

namespace Grid {
namespace QCD {

  template<class Impl>
  class MobiusEOFAFermion : public AbstractEOFAFermion<Impl>
  {
    public:
      INHERIT_IMPL_TYPES(Impl);

    public:
      // Shift operator coefficients for red-black preconditioned Mobius EOFA
      std::vector<Coeff_t> Mooee_shift;
      std::vector<Coeff_t> MooeeInv_shift_lc;
      std::vector<Coeff_t> MooeeInv_shift_norm;
      std::vector<Coeff_t> MooeeInvDag_shift_lc;
      std::vector<Coeff_t> MooeeInvDag_shift_norm;

      virtual void Instantiatable(void) {};

      // EOFA-specific operations
      virtual void  Omega            (const FermionField& in, FermionField& out, int sign, int dag);
      virtual void  Dtilde           (const FermionField& in, FermionField& out);
      virtual void  DtildeInv        (const FermionField& in, FermionField& out);

      // override multiply
      virtual RealD M                (const FermionField& in, FermionField& out);
      virtual RealD Mdag             (const FermionField& in, FermionField& out);

      // half checkerboard operations
      virtual void  Mooee            (const FermionField& in, FermionField& out);
      virtual void  MooeeDag         (const FermionField& in, FermionField& out);
      virtual void  MooeeInv         (const FermionField& in, FermionField& out);
      virtual void  MooeeInv_shift   (const FermionField& in, FermionField& out);
      virtual void  MooeeInvDag      (const FermionField& in, FermionField& out);
      virtual void  MooeeInvDag_shift(const FermionField& in, FermionField& out);

      virtual void   M5D             (const FermionField& psi, FermionField& chi);
      virtual void   M5Ddag          (const FermionField& psi, FermionField& chi);

      /////////////////////////////////////////////////////
      // Instantiate different versions depending on Impl
      /////////////////////////////////////////////////////
      void M5D(const FermionField& psi, const FermionField& phi, FermionField& chi,
        std::vector<Coeff_t>& lower, std::vector<Coeff_t>& diag, std::vector<Coeff_t>& upper);

      void M5D_shift(const FermionField& psi, const FermionField& phi, FermionField& chi,
        std::vector<Coeff_t>& lower, std::vector<Coeff_t>& diag, std::vector<Coeff_t>& upper,
        std::vector<Coeff_t>& shift_coeffs);

      void M5Ddag(const FermionField& psi, const FermionField& phi, FermionField& chi,
        std::vector<Coeff_t>& lower, std::vector<Coeff_t>& diag, std::vector<Coeff_t>& upper);

      void M5Ddag_shift(const FermionField& psi, const FermionField& phi, FermionField& chi,
        std::vector<Coeff_t>& lower, std::vector<Coeff_t>& diag, std::vector<Coeff_t>& upper,
        std::vector<Coeff_t>& shift_coeffs);

      void MooeeInternal(const FermionField& in, FermionField& out, int dag, int inv);

      void MooeeInternalCompute(int dag, int inv, Vector<iSinglet<Simd>>& Matp, Vector<iSinglet<Simd>>& Matm);

      void MooeeInternalAsm(const FermionField& in, FermionField& out, int LLs, int site,
        Vector<iSinglet<Simd>>& Matp, Vector<iSinglet<Simd>>& Matm);

      void MooeeInternalZAsm(const FermionField& in, FermionField& out, int LLs, int site,
        Vector<iSinglet<Simd>>& Matp, Vector<iSinglet<Simd>>& Matm);

      virtual void RefreshShiftCoefficients(RealD new_shift);

      // Constructors
      MobiusEOFAFermion(GaugeField& _Umu, GridCartesian& FiveDimGrid, GridRedBlackCartesian& FiveDimRedBlackGrid,
        GridCartesian& FourDimGrid, GridRedBlackCartesian& FourDimRedBlackGrid,
        RealD _mq1, RealD _mq2, RealD _mq3, RealD _shift, int pm,
        RealD _M5, RealD _b, RealD _c, const ImplParams& p=ImplParams());

    protected:
      void SetCoefficientsPrecondShiftOps(void);
  };
}}

#define INSTANTIATE_DPERP_MOBIUS_EOFA(A)\
template void MobiusEOFAFermion<A>::M5D(const FermionField& psi, const FermionField& phi, FermionField& chi, \
  std::vector<Coeff_t>& lower, std::vector<Coeff_t>& diag, std::vector<Coeff_t>& upper); \
template void MobiusEOFAFermion<A>::M5D_shift(const FermionField& psi, const FermionField& phi, FermionField& chi, \
  std::vector<Coeff_t>& lower, std::vector<Coeff_t>& diag, std::vector<Coeff_t>& upper, std::vector<Coeff_t>& shift_coeffs); \
template void MobiusEOFAFermion<A>::M5Ddag(const FermionField& psi, const FermionField& phi, FermionField& chi, \
  std::vector<Coeff_t>& lower, std::vector<Coeff_t>& diag, std::vector<Coeff_t>& upper); \
template void MobiusEOFAFermion<A>::M5Ddag_shift(const FermionField& psi, const FermionField& phi, FermionField& chi, \
  std::vector<Coeff_t>& lower, std::vector<Coeff_t>& diag, std::vector<Coeff_t>& upper, std::vector<Coeff_t>& shift_coeffs); \
template void MobiusEOFAFermion<A>::MooeeInv(const FermionField& psi, FermionField& chi); \
template void MobiusEOFAFermion<A>::MooeeInv_shift(const FermionField& psi, FermionField& chi); \
template void MobiusEOFAFermion<A>::MooeeInvDag(const FermionField& psi, FermionField& chi); \
template void MobiusEOFAFermion<A>::MooeeInvDag_shift(const FermionField& psi, FermionField& chi);

#undef  MOBIUS_EOFA_DPERP_DENSE
#define MOBIUS_EOFA_DPERP_CACHE
#undef  MOBIUS_EOFA_DPERP_LINALG
#define MOBIUS_EOFA_DPERP_VEC

#endif
