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

NAMESPACE_BEGIN(Grid);

 
template<class Impl>
void MobiusEOFAFermion<Impl>::M5D(const FermionField &psi_i, const FermionField &phi_i, FermionField &chi_i,
				  Vector<Coeff_t> &lower, Vector<Coeff_t> &diag, Vector<Coeff_t> &upper)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase *grid = psi_i.Grid();
  int Ls = this->Ls;
  auto psi = psi_i.View();
  auto phi = phi_i.View();
  auto chi = chi_i.View();

  assert(phi.Checkerboard() == psi.Checkerboard());

  // Flops = 6.0*(Nc*Ns) *Ls*vol
  this->M5Dcalls++;
  this->M5Dtime -= usecond();

  int nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    uint64_t ss = sss*Ls;
    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp1;
    spinor tmp2;
    for(int s=0; s<Ls; s++){
      uint64_t idx_u = ss+((s+1)%Ls);
      uint64_t idx_l = ss+((s+Ls-1)%Ls);
      spProj5m(tmp1, psi(idx_u));
      spProj5p(tmp2, psi(idx_l));
      coalescedWrite(chi[ss+s], diag[s]*phi(ss+s) + upper[s]*tmp1 + lower[s]*tmp2);
    }
  });

  this->M5Dtime += usecond();
}

template<class Impl>
void MobiusEOFAFermion<Impl>::M5D_shift(const FermionField &psi_i, const FermionField &phi_i, FermionField &chi_i,
					Vector<Coeff_t> &lower, Vector<Coeff_t> &diag, Vector<Coeff_t> &upper,
					Vector<Coeff_t> &shift_coeffs)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase *grid = psi_i.Grid();
  int Ls = this->Ls;
  auto psi = psi_i.View();
  auto phi = phi_i.View();
  auto chi = chi_i.View();

  auto pm  = this->pm;
  int shift_s = (pm == 1) ? (Ls-1) : 0; // s-component modified by shift operator

  assert(phi.Checkerboard() == psi.Checkerboard());

  // Flops = 6.0*(Nc*Ns) *Ls*vol
  this->M5Dcalls++;
  this->M5Dtime -= usecond();

  int nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    uint64_t ss = sss*Ls;
    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp1;
    spinor tmp2;
    spinor tmp;
    for(int s=0; s<Ls; s++){
      uint64_t idx_u = ss+((s+1)%Ls);
      uint64_t idx_l = ss+((s+Ls-1)%Ls);
      spProj5m(tmp1, psi(idx_u));
      spProj5p(tmp2, psi(idx_l));

      if(pm == 1){ spProj5p(tmp, psi(ss+shift_s)); }
      else       { spProj5m(tmp, psi(ss+shift_s)); }

      coalescedWrite(chi[ss+s], diag[s]*phi(ss+s) + upper[s]*tmp1 +lower[s]*tmp2 + shift_coeffs[s]*tmp);
    }
  });

  this->M5Dtime += usecond();
}

template<class Impl>
void MobiusEOFAFermion<Impl>::M5Ddag(const FermionField &psi_i, const FermionField &phi_i, FermionField &chi_i,
				     Vector<Coeff_t> &lower, Vector<Coeff_t> &diag, Vector<Coeff_t> &upper)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase *grid = psi_i.Grid();
  int Ls = this->Ls;
  auto psi = psi_i.View();
  auto phi = phi_i.View();
  auto chi = chi_i.View();

  assert(phi.Checkerboard() == psi.Checkerboard());

  // Flops = 6.0*(Nc*Ns) *Ls*vol
  this->M5Dcalls++;
  this->M5Dtime -= usecond();

  int nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(), {
    uint64_t ss = sss*Ls;

    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp1, tmp2;

    for(int s=0; s<Ls; s++){
      uint64_t idx_u = ss+((s+1)%Ls);
      uint64_t idx_l = ss+((s+Ls-1)%Ls);
      spProj5p(tmp1, psi(idx_u));
      spProj5m(tmp2, psi(idx_l));
      coalescedWrite(chi[ss+s], diag[s]*phi(ss+s) + upper[s]*tmp1 + lower[s]*tmp2);
    }
  });

  this->M5Dtime += usecond();
}

template<class Impl>
void MobiusEOFAFermion<Impl>::M5Ddag_shift(const FermionField &psi_i, const FermionField &phi_i, FermionField &chi_i,
					   Vector<Coeff_t> &lower, Vector<Coeff_t> &diag, Vector<Coeff_t> &upper,
					   Vector<Coeff_t> &shift_coeffs)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase *grid = psi_i.Grid();
  int Ls = this->Ls;
  int shift_s = (this->pm == 1) ? (Ls-1) : 0; // s-component modified by shift operator
  auto psi = psi_i.View();
  auto phi = phi_i.View();
  auto chi = chi_i.View();

  assert(phi.Checkerboard() == psi.Checkerboard());

  // Flops = 6.0*(Nc*Ns) *Ls*vol
  this->M5Dcalls++;
  this->M5Dtime -= usecond();

  auto pm = this->pm;

  int nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    uint64_t ss = sss*Ls;

    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp1, tmp2, tmp;
    tmp1=Zero();
    coalescedWrite(chi[ss+Ls-1],tmp1);

    for(int s=0; s<Ls; s++){

      uint64_t idx_u = ss+((s+1)%Ls);
      uint64_t idx_l = ss+((s+Ls-1)%Ls);

      spProj5p(tmp1, psi(idx_u));
      spProj5m(tmp2, psi(idx_l));

      if(s==(Ls-1)) coalescedWrite(chi[ss+s], chi(ss+s)+ diag[s]*phi(ss+s) + upper[s]*tmp1 + lower[s]*tmp2);
      else          coalescedWrite(chi[ss+s], diag[s]*phi(ss+s) + upper[s]*tmp1 + lower[s]*tmp2);
      if(pm == 1){ spProj5p(tmp, psi(ss+s)); }
      else       { spProj5m(tmp, psi(ss+s)); }

      coalescedWrite(chi[ss+shift_s],chi(ss+shift_s)+shift_coeffs[s]*tmp);
    }
  });

  this->M5Dtime += usecond();
}

template<class Impl>
void MobiusEOFAFermion<Impl>::MooeeInv(const FermionField &psi_i, FermionField &chi_i)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase *grid = psi_i.Grid();
  int Ls = this->Ls;
  auto psi = psi_i.View();
  auto chi = chi_i.View();

  auto plee = & this->lee [0];
  auto pdee = & this->dee [0];
  auto puee = & this->uee [0];
  auto pleem= & this->leem[0];
  auto pueem= & this->ueem[0];

  if(this->shift != 0.0){ MooeeInv_shift(psi_i,chi_i); return; }

  this->MooeeInvCalls++;
  this->MooeeInvTime -= usecond();

  int nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{

    uint64_t ss = sss*Ls;

    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp;

    // Apply (L^{\prime})^{-1}
    coalescedWrite(chi[ss], psi(ss)); // chi[0]=psi[0]
    for(int s=1; s<Ls; s++){
      spProj5p(tmp, chi(ss+s-1));
      coalescedWrite(chi[ss+s], psi(ss+s) - plee[s-1]*tmp);
    }

    // L_m^{-1}
    for(int s=0; s<Ls-1; s++){ // Chi[ee] = 1 - sum[s<Ls-1] -leem[s]P_- chi
      spProj5m(tmp, chi(ss+s));
      coalescedWrite(chi[ss+Ls-1], chi(ss+Ls-1) - pleem[s]*tmp);
    }

    // U_m^{-1} D^{-1}
    for(int s=0; s<Ls-1; s++){ // Chi[s] + 1/d chi[s]
      spProj5p(tmp, chi(ss+Ls-1));
      coalescedWrite(chi[ss+s], (1.0/pdee[s])*chi(ss+s) - (pueem[s]/pdee[Ls-1])*tmp);
    }
    coalescedWrite(chi[ss+Ls-1], (1.0/pdee[Ls-1])*chi(ss+Ls-1));

    // Apply U^{-1}
    for(int s=Ls-2; s>=0; s--){
      spProj5m(tmp, chi(ss+s+1));
      coalescedWrite(chi[ss+s], chi(ss+s) - puee[s]*tmp);
    }
  });
   
  this->MooeeInvTime += usecond();
}

template<class Impl>
void MobiusEOFAFermion<Impl>::MooeeInv_shift(const FermionField &psi_i, FermionField &chi_i)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase *grid = psi_i.Grid();
  int Ls = this->Ls;
  auto psi = psi_i.View();
  auto chi = chi_i.View();

  auto pm = this->pm;
  auto plee = & this->lee [0];
  auto pdee = & this->dee [0];
  auto puee = & this->uee [0];
  auto pleem= & this->leem[0];
  auto pueem= & this->ueem[0];
  auto pMooeeInv_shift_lc   = &MooeeInv_shift_lc[0];
  auto pMooeeInv_shift_norm = &MooeeInv_shift_norm[0];
  this->MooeeInvCalls++;
  this->MooeeInvTime -= usecond();

  int nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{

    uint64_t ss = sss*Ls;

    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp1,tmp2,tmp2_spProj;

    // Apply (L^{\prime})^{-1} and accumulate MooeeInv_shift_lc[j]*psi[j] in tmp2
    coalescedWrite(chi[ss], psi(ss)); // chi[0]=psi[0]
    tmp2 = pMooeeInv_shift_lc[0]*psi(ss);
    for(int s=1; s<Ls; s++){
      spProj5p(tmp1, chi(ss+s-1));
      coalescedWrite(chi[ss+s], psi(ss+s) - plee[s-1]*tmp1);
      tmp2 = tmp2 + pMooeeInv_shift_lc[s]*psi(ss+s);
    }
    if(pm == 1){ spProj5p(tmp2_spProj, tmp2);}
    else       { spProj5m(tmp2_spProj, tmp2); }

    // L_m^{-1}
    for(int s=0; s<Ls-1; s++){ // Chi[ee] = 1 - sum[s<Ls-1] -leem[s]P_- chi
      spProj5m(tmp1, chi(ss+s));
      coalescedWrite(chi[ss+Ls-1], chi(ss+Ls-1) - pleem[s]*tmp1);
    }

    // U_m^{-1} D^{-1}
    for(int s=0; s<Ls-1; s++){ // Chi[s] + 1/d chi[s]
      spProj5p(tmp1, chi(ss+Ls-1));
      coalescedWrite(chi[ss+s], (1.0/pdee[s])*chi(ss+s) - (pueem[s]/pdee[Ls-1])*tmp1);
    }
    // chi[ss+Ls-1] = (1.0/pdee[Ls-1])*chi[ss+Ls-1] + MooeeInv_shift_norm[Ls-1]*tmp2_spProj;
    coalescedWrite(chi[ss+Ls-1], (1.0/pdee[Ls-1])*chi(ss+Ls-1));
    spProj5m(tmp1, chi(ss+Ls-1));
    coalescedWrite(chi[ss+Ls-1], chi(ss+Ls-1) + pMooeeInv_shift_norm[Ls-1]*tmp2_spProj);

    // Apply U^{-1} and add shift term
    for(int s=Ls-2; s>=0; s--){
      coalescedWrite(chi[ss+s] , chi(ss+s) - puee[s]*tmp1);
      spProj5m(tmp1, chi(ss+s));
      coalescedWrite(chi[ss+s], chi(ss+s) + pMooeeInv_shift_norm[s]*tmp2_spProj);
    }
  });

  this->MooeeInvTime += usecond();
}

template<class Impl>
void MobiusEOFAFermion<Impl>::MooeeInvDag(const FermionField &psi_i, FermionField &chi_i)
{
  if(this->shift != 0.0){ MooeeInvDag_shift(psi_i,chi_i); return; }

  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase *grid = psi_i.Grid();
  int Ls = this->Ls;
  auto psi = psi_i.View();
  auto chi = chi_i.View();

  auto plee = & this->lee [0];
  auto pdee = & this->dee [0];
  auto puee = & this->uee [0];
  auto pleem= & this->leem[0];
  auto pueem= & this->ueem[0];

  this->MooeeInvCalls++;
  this->MooeeInvTime -= usecond();

  int nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{

    uint64_t ss = sss*Ls;

    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp;

    // Apply (U^{\prime})^{-dag}
    coalescedWrite(chi[ss], psi(ss));
    for(int s=1; s<Ls; s++){
      spProj5m(tmp, chi(ss+s-1));
      coalescedWrite(chi[ss+s], psi(ss+s) - puee[s-1]*tmp);
    }
    
    // U_m^{-\dag}
    for(int s=0; s<Ls-1; s++){
      spProj5p(tmp, chi(ss+s));
      coalescedWrite(chi[ss+Ls-1], chi(ss+Ls-1) - pueem[s]*tmp);
    }

    // L_m^{-\dag} D^{-dag}
    for(int s=0; s<Ls-1; s++){
      spProj5m(tmp, chi(ss+Ls-1));
      coalescedWrite(chi[ss+s], (1.0/pdee[s])*chi(ss+s) - (pleem[s]/pdee[Ls-1])*tmp);
    }
    coalescedWrite(chi[ss+Ls-1], (1.0/pdee[Ls-1])*chi(ss+Ls-1));

    // Apply L^{-dag}
    for(int s=Ls-2; s>=0; s--){
      spProj5p(tmp, chi(ss+s+1));
      coalescedWrite(chi[ss+s], chi(ss+s) - plee[s]*tmp);
    }
  });

  this->MooeeInvTime += usecond();
}

template<class Impl>
void MobiusEOFAFermion<Impl>::MooeeInvDag_shift(const FermionField &psi_i, FermionField &chi_i)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase *grid = psi_i.Grid();
  auto psi = psi_i.View();
  auto chi = chi_i.View();
  int Ls = this->Ls;

  auto pm = this->pm;
  auto plee = & this->lee [0];
  auto pdee = & this->dee [0];
  auto puee = & this->uee [0];
  auto pleem= & this->leem[0];
  auto pueem= & this->ueem[0];
  auto pMooeeInvDag_shift_lc   = &MooeeInvDag_shift_lc[0];
  auto pMooeeInvDag_shift_norm = &MooeeInvDag_shift_norm[0];

  this->MooeeInvCalls++;
  this->MooeeInvTime -= usecond();

  int nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{

    uint64_t ss = sss*Ls;

    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp1,tmp2,tmp2_spProj;

    // Apply (U^{\prime})^{-dag} and accumulate MooeeInvDag_shift_lc[j]*psi[j] in tmp2
    coalescedWrite(chi[ss], psi(ss));
    tmp2 = pMooeeInvDag_shift_lc[0]*psi(ss);
    for(int s=1; s<Ls; s++){
      spProj5m(tmp1, chi(ss+s-1));
      coalescedWrite(chi[ss+s],psi(ss+s) - puee[s-1]*tmp1);
      tmp2 = tmp2 + pMooeeInvDag_shift_lc[s]*psi(ss+s);
    }

    if(pm == 1){ spProj5p(tmp2_spProj, tmp2);}
    else       { spProj5m(tmp2_spProj, tmp2);}

    // U_m^{-\dag}
    for(int s=0; s<Ls-1; s++){
      spProj5p(tmp1, chi(ss+s));
      coalescedWrite(chi[ss+Ls-1], chi(ss+Ls-1) - pueem[s]*tmp1);
    }

    // L_m^{-\dag} D^{-dag}
    for(int s=0; s<Ls-1; s++){
      spProj5m(tmp1, chi(ss+Ls-1));
      coalescedWrite(chi[ss+s], (1.0/pdee[s])*chi(ss+s) - (pleem[s]/pdee[Ls-1])*tmp1);
    }
    coalescedWrite(chi[ss+Ls-1], (1.0/pdee[Ls-1])*chi(ss+Ls-1));
    spProj5p(tmp1, chi(ss+Ls-1));
    coalescedWrite(chi[ss+Ls-1], chi(ss+Ls-1) + pMooeeInvDag_shift_norm[Ls-1]*tmp2_spProj);

    // Apply L^{-dag}
    for(int s=Ls-2; s>=0; s--){
      coalescedWrite(chi[ss+s], chi(ss+s) - plee[s]*tmp1);
      spProj5p(tmp1, chi(ss+s));
      coalescedWrite(chi[ss+s], chi(ss+s) + pMooeeInvDag_shift_norm[s]*tmp2_spProj);
    }
  });

  this->MooeeInvTime += usecond();
}

NAMESPACE_END(Grid);
