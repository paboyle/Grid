/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/DomainWallEOFAFermionssp.cc

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

#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/DomainWallEOFAFermion.h>

namespace Grid {
namespace QCD {

    // FIXME -- make a version of these routines with site loop outermost for cache reuse.
    // Pminus fowards
    // Pplus  backwards
    template<class Impl>
    void DomainWallEOFAFermion<Impl>::M5D(const FermionField& psi, const FermionField& phi,
        FermionField& chi, std::vector<Coeff_t>& lower, std::vector<Coeff_t>& diag, std::vector<Coeff_t>& upper)
    {
        Coeff_t one(1.0);
        int Ls = this->Ls;
        for(int s=0; s<Ls; s++){
            if(s==0) {
              axpby_ssp_pminus(chi, diag[s], phi, upper[s], psi, s, s+1);
              axpby_ssp_pplus (chi, one, chi, lower[s], psi, s, Ls-1);
            } else if (s==(Ls-1)) {
              axpby_ssp_pminus(chi, diag[s], phi, upper[s], psi, s, 0);
              axpby_ssp_pplus (chi, one, chi, lower[s], psi, s, s-1);
            } else {
              axpby_ssp_pminus(chi, diag[s], phi, upper[s], psi, s, s+1);
              axpby_ssp_pplus(chi, one, chi, lower[s], psi, s, s-1);
            }
        }
    }

    template<class Impl>
    void DomainWallEOFAFermion<Impl>::M5Ddag(const FermionField& psi, const FermionField& phi,
        FermionField& chi, std::vector<Coeff_t>& lower, std::vector<Coeff_t>& diag, std::vector<Coeff_t>& upper)
    {
        Coeff_t one(1.0);
        int Ls = this->Ls;
        for(int s=0; s<Ls; s++){
            if(s==0) {
              axpby_ssp_pplus (chi, diag[s], phi, upper[s], psi, s, s+1);
              axpby_ssp_pminus(chi, one, chi, lower[s], psi, s, Ls-1);
            } else if (s==(Ls-1)) {
              axpby_ssp_pplus (chi, diag[s], phi, upper[s], psi, s, 0);
              axpby_ssp_pminus(chi, one, chi, lower[s], psi, s, s-1);
            } else {
              axpby_ssp_pplus (chi, diag[s], phi, upper[s], psi, s, s+1);
              axpby_ssp_pminus(chi, one, chi, lower[s], psi, s, s-1);
            }
        }
    }

    template<class Impl>
    void DomainWallEOFAFermion<Impl>::MooeeInv(const FermionField& psi, FermionField& chi)
    {
        Coeff_t one(1.0);
        Coeff_t czero(0.0);
        chi.checkerboard = psi.checkerboard;
        int Ls = this->Ls;

        FermionField tmp(psi._grid);

        // Apply (L^{\prime})^{-1}
        axpby_ssp(chi, one, psi, czero, psi, 0, 0);      // chi[0]=psi[0]
        for(int s=1; s<Ls; s++){
            axpby_ssp_pplus(chi, one, psi, -this->lee[s-1], chi, s, s-1);// recursion Psi[s] -lee P_+ chi[s-1]
        }

        // L_m^{-1}
        for(int s=0; s<Ls-1; s++){ // Chi[ee] = 1 - sum[s<Ls-1] -leem[s]P_- chi
            axpby_ssp_pminus(chi, one, chi, -this->leem[s], chi, Ls-1, s);
        }

        // U_m^{-1} D^{-1}
        for(int s=0; s<Ls-1; s++){
            axpby_ssp_pplus(chi, one/this->dee[s], chi, -this->ueem[s]/this->dee[Ls], chi, s, Ls-1);
        }
        axpby_ssp_pminus(tmp, czero, chi, one/this->dee[Ls-1], chi, Ls-1, Ls-1);
        axpby_ssp_pplus(chi, one, tmp, one/this->dee[Ls], chi, Ls-1, Ls-1);

        // Apply U^{-1}
        for(int s=Ls-2; s>=0; s--){
            axpby_ssp_pminus(chi, one, chi, -this->uee[s], chi, s, s+1);  // chi[Ls]
        }
    }

    template<class Impl>
    void DomainWallEOFAFermion<Impl>::MooeeInvDag(const FermionField& psi, FermionField& chi)
    {
        Coeff_t one(1.0);
        Coeff_t czero(0.0);
        chi.checkerboard = psi.checkerboard;
        int Ls = this->Ls;

        FermionField tmp(psi._grid);

        // Apply (U^{\prime})^{-dagger}
        axpby_ssp(chi, one, psi, czero, psi, 0, 0);      // chi[0]=psi[0]
        for(int s=1; s<Ls; s++){
            axpby_ssp_pminus(chi, one, psi, -conjugate(this->uee[s-1]), chi, s, s-1);
        }

        // U_m^{-\dagger}
        for(int s=0; s<Ls-1; s++){
            axpby_ssp_pplus(chi, one, chi, -conjugate(this->ueem[s]), chi, Ls-1, s);
        }

        // L_m^{-\dagger} D^{-dagger}
        for(int s=0; s<Ls-1; s++){
            axpby_ssp_pminus(chi, one/conjugate(this->dee[s]), chi, -conjugate(this->leem[s]/this->dee[Ls-1]), chi, s, Ls-1);
        }
        axpby_ssp_pminus(tmp, czero, chi, one/conjugate(this->dee[Ls-1]), chi, Ls-1, Ls-1);
        axpby_ssp_pplus(chi, one, tmp, one/conjugate(this->dee[Ls]), chi, Ls-1, Ls-1);

        // Apply L^{-dagger}
        for(int s=Ls-2; s>=0; s--){
            axpby_ssp_pplus(chi, one, chi, -conjugate(this->lee[s]), chi, s, s+1);  // chi[Ls]
        }
    }

    #ifdef DOMAIN_WALL_EOFA_DPERP_LINALG

        INSTANTIATE_DPERP_DWF_EOFA(WilsonImplF);
        INSTANTIATE_DPERP_DWF_EOFA(WilsonImplD);
        INSTANTIATE_DPERP_DWF_EOFA(GparityWilsonImplF);
        INSTANTIATE_DPERP_DWF_EOFA(GparityWilsonImplD);
        INSTANTIATE_DPERP_DWF_EOFA(ZWilsonImplF);
        INSTANTIATE_DPERP_DWF_EOFA(ZWilsonImplD);

        INSTANTIATE_DPERP_DWF_EOFA(WilsonImplFH);
        INSTANTIATE_DPERP_DWF_EOFA(WilsonImplDF);
        INSTANTIATE_DPERP_DWF_EOFA(GparityWilsonImplFH);
        INSTANTIATE_DPERP_DWF_EOFA(GparityWilsonImplDF);
        INSTANTIATE_DPERP_DWF_EOFA(ZWilsonImplFH);
        INSTANTIATE_DPERP_DWF_EOFA(ZWilsonImplDF);

    #endif

}}
