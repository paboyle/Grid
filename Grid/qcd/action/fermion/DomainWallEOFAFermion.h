/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/DomainWallEOFAFermion.h

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
#ifndef  GRID_QCD_DOMAIN_WALL_EOFA_FERMION_H
#define  GRID_QCD_DOMAIN_WALL_EOFA_FERMION_H

#include <Grid/qcd/action/fermion/AbstractEOFAFermion.h>

namespace Grid {
namespace QCD {

  template<class Impl>
  class DomainWallEOFAFermion : public AbstractEOFAFermion<Impl>
  {
    public:
      INHERIT_IMPL_TYPES(Impl);

    public:
      // Modified (0,Ls-1) and (Ls-1,0) elements of Mooee
      // for red-black preconditioned Shamir EOFA
      Coeff_t dm;
      Coeff_t dp;

      virtual void Instantiatable(void) {};

      // EOFA-specific operations
      virtual void  Omega      (const FermionField& in, FermionField& out, int sign, int dag);
      virtual void  Dtilde     (const FermionField& in, FermionField& out);
      virtual void  DtildeInv  (const FermionField& in, FermionField& out);

      // override multiply
      virtual RealD M          (const FermionField& in, FermionField& out);
      virtual RealD Mdag       (const FermionField& in, FermionField& out);

      // half checkerboard operations
      virtual void  Mooee      (const FermionField& in, FermionField& out);
      virtual void  MooeeDag   (const FermionField& in, FermionField& out);
      virtual void  MooeeInv   (const FermionField& in, FermionField& out);
      virtual void  MooeeInvDag(const FermionField& in, FermionField& out);

      virtual void   M5D       (const FermionField& psi, FermionField& chi);
      virtual void   M5Ddag    (const FermionField& psi, FermionField& chi);

      /////////////////////////////////////////////////////
      // Instantiate different versions depending on Impl
      /////////////////////////////////////////////////////
      void M5D(const FermionField& psi, const FermionField& phi, FermionField& chi,
        std::vector<Coeff_t>& lower, std::vector<Coeff_t>& diag, std::vector<Coeff_t>& upper);

      void M5Ddag(const FermionField& psi, const FermionField& phi, FermionField& chi,
        std::vector<Coeff_t>& lower, std::vector<Coeff_t>& diag, std::vector<Coeff_t>& upper);

      void MooeeInternal(const FermionField& in, FermionField& out, int dag, int inv);

      void MooeeInternalCompute(int dag, int inv, Vector<iSinglet<Simd>>& Matp, Vector<iSinglet<Simd>>& Matm);

      void MooeeInternalAsm(const FermionField& in, FermionField& out, int LLs, int site,
        Vector<iSinglet<Simd>>& Matp, Vector<iSinglet<Simd>>& Matm);

      void MooeeInternalZAsm(const FermionField& in, FermionField& out, int LLs, int site,
        Vector<iSinglet<Simd>>& Matp, Vector<iSinglet<Simd>>& Matm);

      virtual void RefreshShiftCoefficients(RealD new_shift);

      // Constructors
      DomainWallEOFAFermion(GaugeField& _Umu, GridCartesian& FiveDimGrid, GridRedBlackCartesian& FiveDimRedBlackGrid,
        GridCartesian& FourDimGrid, GridRedBlackCartesian& FourDimRedBlackGrid,
        RealD _mq1, RealD _mq2, RealD _mq3, RealD _shift, int pm,
        RealD _M5, const ImplParams& p=ImplParams());

    protected:
      void SetCoefficientsInternal(RealD zolo_hi, std::vector<Coeff_t>& gamma, RealD b, RealD c);
  };
}}

#define INSTANTIATE_DPERP_DWF_EOFA(A)\
template void DomainWallEOFAFermion<A>::M5D(const FermionField& psi, const FermionField& phi, FermionField& chi, \
  std::vector<Coeff_t>& lower, std::vector<Coeff_t>& diag, std::vector<Coeff_t>& upper); \
template void DomainWallEOFAFermion<A>::M5Ddag(const FermionField& psi, const FermionField& phi, FermionField& chi, \
  std::vector<Coeff_t>& lower, std::vector<Coeff_t>& diag, std::vector<Coeff_t>& upper); \
template void DomainWallEOFAFermion<A>::MooeeInv(const FermionField& psi, FermionField& chi); \
template void DomainWallEOFAFermion<A>::MooeeInvDag(const FermionField& psi, FermionField& chi);

#undef  DOMAIN_WALL_EOFA_DPERP_DENSE
#define DOMAIN_WALL_EOFA_DPERP_CACHE
#undef  DOMAIN_WALL_EOFA_DPERP_LINALG
#define DOMAIN_WALL_EOFA_DPERP_VEC

#endif
