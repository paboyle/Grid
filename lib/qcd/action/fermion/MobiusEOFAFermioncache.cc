/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/MobiusEOFAFermioncache.cc

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
#include <Grid/qcd/action/fermion/MobiusEOFAFermion.h>

namespace Grid {
namespace QCD {

  // FIXME -- make a version of these routines with site loop outermost for cache reuse.

  template<class Impl>
  void MobiusEOFAFermion<Impl>::M5D(const FermionField &psi, const FermionField &phi, FermionField &chi,
    std::vector<Coeff_t> &lower, std::vector<Coeff_t> &diag, std::vector<Coeff_t> &upper)
  {
    int Ls = this->Ls;
    GridBase *grid = psi._grid;

    assert(phi.checkerboard == psi.checkerboard);
    chi.checkerboard = psi.checkerboard;

    // Flops = 6.0*(Nc*Ns) *Ls*vol
    this->M5Dcalls++;
    this->M5Dtime -= usecond();

    parallel_for(int ss=0; ss<grid->oSites(); ss+=Ls){
      for(int s=0; s<Ls; s++){
        auto tmp = psi._odata[0];
        if(s==0){
          spProj5m(tmp, psi._odata[ss+s+1]);
          chi[ss+s] = diag[s]*phi[ss+s] + upper[s]*tmp;
          spProj5p(tmp, psi._odata[ss+Ls-1]);
          chi[ss+s] = chi[ss+s] + lower[s]*tmp;
        } else if(s==(Ls-1)) {
          spProj5m(tmp, psi._odata[ss+0]);
          chi[ss+s] = diag[s]*phi[ss+s] + upper[s]*tmp;
          spProj5p(tmp, psi._odata[ss+s-1]);
          chi[ss+s] = chi[ss+s] + lower[s]*tmp;
        } else {
          spProj5m(tmp, psi._odata[ss+s+1]);
          chi[ss+s] = diag[s]*phi[ss+s] + upper[s]*tmp;
          spProj5p(tmp, psi._odata[ss+s-1]);
          chi[ss+s] = chi[ss+s] + lower[s]*tmp;
        }
      }
    }

    this->M5Dtime += usecond();
  }

  template<class Impl>
  void MobiusEOFAFermion<Impl>::M5D_shift(const FermionField &psi, const FermionField &phi, FermionField &chi,
    std::vector<Coeff_t> &lower, std::vector<Coeff_t> &diag, std::vector<Coeff_t> &upper,
    std::vector<Coeff_t> &shift_coeffs)
  {
    int Ls = this->Ls;
    int shift_s = (this->pm == 1) ? (Ls-1) : 0; // s-component modified by shift operator
    GridBase *grid = psi._grid;

    assert(phi.checkerboard == psi.checkerboard);
    chi.checkerboard = psi.checkerboard;

    // Flops = 6.0*(Nc*Ns) *Ls*vol
    this->M5Dcalls++;
    this->M5Dtime -= usecond();

    parallel_for(int ss=0; ss<grid->oSites(); ss+=Ls){
      for(int s=0; s<Ls; s++){
        auto tmp = psi._odata[0];
        if(s==0){
          spProj5m(tmp, psi._odata[ss+s+1]);
          chi[ss+s] = diag[s]*phi[ss+s] + upper[s]*tmp;
          spProj5p(tmp, psi._odata[ss+Ls-1]);
          chi[ss+s] = chi[ss+s] + lower[s]*tmp;
        } else if(s==(Ls-1)) {
          spProj5m(tmp, psi._odata[ss+0]);
          chi[ss+s] = diag[s]*phi[ss+s] + upper[s]*tmp;
          spProj5p(tmp, psi._odata[ss+s-1]);
          chi[ss+s] = chi[ss+s] + lower[s]*tmp;
        } else {
          spProj5m(tmp, psi._odata[ss+s+1]);
          chi[ss+s] = diag[s]*phi[ss+s] + upper[s]*tmp;
          spProj5p(tmp, psi._odata[ss+s-1]);
          chi[ss+s] = chi[ss+s] + lower[s]*tmp;
        }
        if(this->pm == 1){ spProj5p(tmp, psi._odata[ss+shift_s]); }
        else{ spProj5m(tmp, psi._odata[ss+shift_s]); }
        chi[ss+s] = chi[ss+s] + shift_coeffs[s]*tmp;
      }
    }

    this->M5Dtime += usecond();
  }

  template<class Impl>
  void MobiusEOFAFermion<Impl>::M5Ddag(const FermionField &psi, const FermionField &phi, FermionField &chi,
    std::vector<Coeff_t> &lower, std::vector<Coeff_t> &diag, std::vector<Coeff_t> &upper)
  {
    int Ls = this->Ls;
    GridBase *grid = psi._grid;

    assert(phi.checkerboard == psi.checkerboard);
    chi.checkerboard = psi.checkerboard;

    // Flops = 6.0*(Nc*Ns) *Ls*vol
    this->M5Dcalls++;
    this->M5Dtime -= usecond();

    parallel_for(int ss=0; ss<grid->oSites(); ss+=Ls){
      auto tmp = psi._odata[0];
      for(int s=0; s<Ls; s++){
        if(s==0) {
          spProj5p(tmp, psi._odata[ss+s+1]);
          chi[ss+s] = diag[s]*phi[ss+s] + upper[s]*tmp;
          spProj5m(tmp, psi._odata[ss+Ls-1]);
          chi[ss+s] = chi[ss+s] + lower[s]*tmp;
        } else if(s==(Ls-1)) {
          spProj5p(tmp, psi._odata[ss+0]);
          chi[ss+s] = diag[s]*phi[ss+s] + upper[s]*tmp;
          spProj5m(tmp, psi._odata[ss+s-1]);
          chi[ss+s] = chi[ss+s] + lower[s]*tmp;
        } else {
          spProj5p(tmp, psi._odata[ss+s+1]);
          chi[ss+s] = diag[s]*phi[ss+s] + upper[s]*tmp;
          spProj5m(tmp, psi._odata[ss+s-1]);
          chi[ss+s] = chi[ss+s] + lower[s]*tmp;
        }
      }
    }

    this->M5Dtime += usecond();
  }

  template<class Impl>
  void MobiusEOFAFermion<Impl>::M5Ddag_shift(const FermionField &psi, const FermionField &phi, FermionField &chi,
    std::vector<Coeff_t> &lower, std::vector<Coeff_t> &diag, std::vector<Coeff_t> &upper,
    std::vector<Coeff_t> &shift_coeffs)
  {
    int Ls = this->Ls;
    int shift_s = (this->pm == 1) ? (Ls-1) : 0; // s-component modified by shift operator
    GridBase *grid = psi._grid;

    assert(phi.checkerboard == psi.checkerboard);
    chi.checkerboard = psi.checkerboard;

    // Flops = 6.0*(Nc*Ns) *Ls*vol
    this->M5Dcalls++;
    this->M5Dtime -= usecond();

    parallel_for(int ss=0; ss<grid->oSites(); ss+=Ls){
      chi[ss+Ls-1] = zero;
      auto tmp = psi._odata[0];
      for(int s=0; s<Ls; s++){
        if(s==0) {
          spProj5p(tmp, psi._odata[ss+s+1]);
          chi[ss+s] = diag[s]*phi[ss+s] + upper[s]*tmp;
          spProj5m(tmp, psi._odata[ss+Ls-1]);
          chi[ss+s] = chi[ss+s] + lower[s]*tmp;
        } else if(s==(Ls-1)) {
          spProj5p(tmp, psi._odata[ss+0]);
          chi[ss+s] = chi[ss+s] + diag[s]*phi[ss+s] + upper[s]*tmp;
          spProj5m(tmp, psi._odata[ss+s-1]);
          chi[ss+s] = chi[ss+s] + lower[s]*tmp;
        } else {
          spProj5p(tmp, psi._odata[ss+s+1]);
          chi[ss+s] = diag[s]*phi[ss+s] + upper[s]*tmp;
          spProj5m(tmp, psi._odata[ss+s-1]);
          chi[ss+s] = chi[ss+s] + lower[s]*tmp;
        }
        if(this->pm == 1){ spProj5p(tmp, psi._odata[ss+s]); }
        else{ spProj5m(tmp, psi._odata[ss+s]); }
        chi[ss+shift_s] = chi[ss+shift_s] + shift_coeffs[s]*tmp;
      }
    }

    this->M5Dtime += usecond();
  }

  template<class Impl>
  void MobiusEOFAFermion<Impl>::MooeeInv(const FermionField &psi, FermionField &chi)
  {
    if(this->shift != 0.0){ MooeeInv_shift(psi,chi); return; }

    GridBase *grid = psi._grid;
    int Ls = this->Ls;

    chi.checkerboard = psi.checkerboard;

    this->MooeeInvCalls++;
    this->MooeeInvTime -= usecond();

    parallel_for(int ss=0; ss<grid->oSites(); ss+=Ls){

      auto tmp = psi._odata[0];

      // Apply (L^{\prime})^{-1}
      chi[ss] = psi[ss]; // chi[0]=psi[0]
      for(int s=1; s<Ls; s++){
        spProj5p(tmp, chi[ss+s-1]);
        chi[ss+s] = psi[ss+s] - this->lee[s-1]*tmp;
      }

      // L_m^{-1}
      for(int s=0; s<Ls-1; s++){ // Chi[ee] = 1 - sum[s<Ls-1] -leem[s]P_- chi
        spProj5m(tmp, chi[ss+s]);
        chi[ss+Ls-1] = chi[ss+Ls-1] - this->leem[s]*tmp;
      }

      // U_m^{-1} D^{-1}
      for(int s=0; s<Ls-1; s++){ // Chi[s] + 1/d chi[s]
        spProj5p(tmp, chi[ss+Ls-1]);
        chi[ss+s] = (1.0/this->dee[s])*chi[ss+s] - (this->ueem[s]/this->dee[Ls-1])*tmp;
      }
      chi[ss+Ls-1] = (1.0/this->dee[Ls-1])*chi[ss+Ls-1];

      // Apply U^{-1}
      for(int s=Ls-2; s>=0; s--){
        spProj5m(tmp, chi[ss+s+1]);
        chi[ss+s] = chi[ss+s] - this->uee[s]*tmp;
      }
    }

    this->MooeeInvTime += usecond();
  }

  template<class Impl>
  void MobiusEOFAFermion<Impl>::MooeeInv_shift(const FermionField &psi, FermionField &chi)
  {
    GridBase *grid = psi._grid;
    int Ls = this->Ls;

    chi.checkerboard = psi.checkerboard;

    this->MooeeInvCalls++;
    this->MooeeInvTime -= usecond();

    parallel_for(int ss=0; ss<grid->oSites(); ss+=Ls){

      auto tmp1        = psi._odata[0];
      auto tmp2        = psi._odata[0];
      auto tmp2_spProj = psi._odata[0];

      // Apply (L^{\prime})^{-1} and accumulate MooeeInv_shift_lc[j]*psi[j] in tmp2
      chi[ss] = psi[ss]; // chi[0]=psi[0]
      tmp2 = MooeeInv_shift_lc[0]*psi[ss];
      for(int s=1; s<Ls; s++){
        spProj5p(tmp1, chi[ss+s-1]);
        chi[ss+s] = psi[ss+s] - this->lee[s-1]*tmp1;
        tmp2 = tmp2 + MooeeInv_shift_lc[s]*psi[ss+s];
      }
      if(this->pm == 1){ spProj5p(tmp2_spProj, tmp2);}
      else{ spProj5m(tmp2_spProj, tmp2); }

      // L_m^{-1}
      for(int s=0; s<Ls-1; s++){ // Chi[ee] = 1 - sum[s<Ls-1] -leem[s]P_- chi
        spProj5m(tmp1, chi[ss+s]);
        chi[ss+Ls-1] = chi[ss+Ls-1] - this->leem[s]*tmp1;
      }

      // U_m^{-1} D^{-1}
      for(int s=0; s<Ls-1; s++){ // Chi[s] + 1/d chi[s]
        spProj5p(tmp1, chi[ss+Ls-1]);
        chi[ss+s] = (1.0/this->dee[s])*chi[ss+s] - (this->ueem[s]/this->dee[Ls-1])*tmp1;
      }
      // chi[ss+Ls-1] = (1.0/this->dee[Ls-1])*chi[ss+Ls-1] + MooeeInv_shift_norm[Ls-1]*tmp2_spProj;
      chi[ss+Ls-1] = (1.0/this->dee[Ls-1])*chi[ss+Ls-1];
      spProj5m(tmp1, chi[ss+Ls-1]);
      chi[ss+Ls-1] = chi[ss+Ls-1] + MooeeInv_shift_norm[Ls-1]*tmp2_spProj;

      // Apply U^{-1} and add shift term
      for(int s=Ls-2; s>=0; s--){
        chi[ss+s] = chi[ss+s] - this->uee[s]*tmp1;
        spProj5m(tmp1, chi[ss+s]);
        chi[ss+s] = chi[ss+s] + MooeeInv_shift_norm[s]*tmp2_spProj;
      }
    }

    this->MooeeInvTime += usecond();
  }

  template<class Impl>
  void MobiusEOFAFermion<Impl>::MooeeInvDag(const FermionField &psi, FermionField &chi)
  {
    if(this->shift != 0.0){ MooeeInvDag_shift(psi,chi); return; }

    GridBase *grid = psi._grid;
    int Ls = this->Ls;

    chi.checkerboard = psi.checkerboard;

    this->MooeeInvCalls++;
    this->MooeeInvTime -= usecond();

    parallel_for(int ss=0; ss<grid->oSites(); ss+=Ls){

      auto tmp = psi._odata[0];

      // Apply (U^{\prime})^{-dag}
      chi[ss] = psi[ss];
      for(int s=1; s<Ls; s++){
        spProj5m(tmp, chi[ss+s-1]);
        chi[ss+s] = psi[ss+s] - this->uee[s-1]*tmp;
      }

      // U_m^{-\dag}
      for(int s=0; s<Ls-1; s++){
        spProj5p(tmp, chi[ss+s]);
        chi[ss+Ls-1] = chi[ss+Ls-1] - this->ueem[s]*tmp;
      }

      // L_m^{-\dag} D^{-dag}
      for(int s=0; s<Ls-1; s++){
        spProj5m(tmp, chi[ss+Ls-1]);
        chi[ss+s] = (1.0/this->dee[s])*chi[ss+s] - (this->leem[s]/this->dee[Ls-1])*tmp;
      }
      chi[ss+Ls-1] = (1.0/this->dee[Ls-1])*chi[ss+Ls-1];

      // Apply L^{-dag}
      for(int s=Ls-2; s>=0; s--){
        spProj5p(tmp, chi[ss+s+1]);
        chi[ss+s] = chi[ss+s] - this->lee[s]*tmp;
      }
    }

    this->MooeeInvTime += usecond();
  }

  template<class Impl>
  void MobiusEOFAFermion<Impl>::MooeeInvDag_shift(const FermionField &psi, FermionField &chi)
  {
    GridBase *grid = psi._grid;
    int Ls = this->Ls;

    chi.checkerboard = psi.checkerboard;

    this->MooeeInvCalls++;
    this->MooeeInvTime -= usecond();

    parallel_for(int ss=0; ss<grid->oSites(); ss+=Ls){

      auto tmp1        = psi._odata[0];
      auto tmp2        = psi._odata[0];
      auto tmp2_spProj = psi._odata[0];

      // Apply (U^{\prime})^{-dag} and accumulate MooeeInvDag_shift_lc[j]*psi[j] in tmp2
      chi[ss] = psi[ss];
      tmp2 = MooeeInvDag_shift_lc[0]*psi[ss];
      for(int s=1; s<Ls; s++){
        spProj5m(tmp1, chi[ss+s-1]);
        chi[ss+s] = psi[ss+s] - this->uee[s-1]*tmp1;
        tmp2 = tmp2 + MooeeInvDag_shift_lc[s]*psi[ss+s];
      }
      if(this->pm == 1){ spProj5p(tmp2_spProj, tmp2);}
      else{ spProj5m(tmp2_spProj, tmp2); }

      // U_m^{-\dag}
      for(int s=0; s<Ls-1; s++){
        spProj5p(tmp1, chi[ss+s]);
        chi[ss+Ls-1] = chi[ss+Ls-1] - this->ueem[s]*tmp1;
      }

      // L_m^{-\dag} D^{-dag}
      for(int s=0; s<Ls-1; s++){
        spProj5m(tmp1, chi[ss+Ls-1]);
        chi[ss+s] = (1.0/this->dee[s])*chi[ss+s] - (this->leem[s]/this->dee[Ls-1])*tmp1;
      }
      chi[ss+Ls-1] = (1.0/this->dee[Ls-1])*chi[ss+Ls-1];
      spProj5p(tmp1, chi[ss+Ls-1]);
      chi[ss+Ls-1] = chi[ss+Ls-1] + MooeeInvDag_shift_norm[Ls-1]*tmp2_spProj;

      // Apply L^{-dag}
      for(int s=Ls-2; s>=0; s--){
        chi[ss+s] = chi[ss+s] - this->lee[s]*tmp1;
        spProj5p(tmp1, chi[ss+s]);
        chi[ss+s] = chi[ss+s] + MooeeInvDag_shift_norm[s]*tmp2_spProj;
      }
    }

    this->MooeeInvTime += usecond();
  }

  #ifdef MOBIUS_EOFA_DPERP_CACHE

    INSTANTIATE_DPERP_MOBIUS_EOFA(WilsonImplF);
    INSTANTIATE_DPERP_MOBIUS_EOFA(WilsonImplD);
    INSTANTIATE_DPERP_MOBIUS_EOFA(GparityWilsonImplF);
    INSTANTIATE_DPERP_MOBIUS_EOFA(GparityWilsonImplD);
    INSTANTIATE_DPERP_MOBIUS_EOFA(ZWilsonImplF);
    INSTANTIATE_DPERP_MOBIUS_EOFA(ZWilsonImplD);

    INSTANTIATE_DPERP_MOBIUS_EOFA(WilsonImplFH);
    INSTANTIATE_DPERP_MOBIUS_EOFA(WilsonImplDF);
    INSTANTIATE_DPERP_MOBIUS_EOFA(GparityWilsonImplFH);
    INSTANTIATE_DPERP_MOBIUS_EOFA(GparityWilsonImplDF);
    INSTANTIATE_DPERP_MOBIUS_EOFA(ZWilsonImplFH);
    INSTANTIATE_DPERP_MOBIUS_EOFA(ZWilsonImplDF);

  #endif

}}
