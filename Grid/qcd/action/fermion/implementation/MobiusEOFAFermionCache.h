/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/MobiusEOFAFermioncache.cc

Copyright (C) 2017

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: David Murphy <dmurphy@phys.columbia.edu>
Author: Gianluca Filaci <g.filaci@ed.ac.uk>

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
				  std::vector<Coeff_t> &lower, std::vector<Coeff_t> &diag, std::vector<Coeff_t> &upper)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase *grid = psi_i.Grid();
  int Ls = this->Ls;
  autoView(psi , psi_i, AcceleratorRead);
  autoView(phi , phi_i, AcceleratorRead);
  autoView(chi , chi_i, AcceleratorWrite);

  assert(phi.Checkerboard() == psi.Checkerboard());

  auto pdiag  = &this->d_diag[0];
  auto pupper = &this->d_upper[0];
  auto plower = &this->d_lower[0];

  acceleratorCopyToDevice(&diag[0],&pdiag[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&upper[0],&pupper[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&lower[0],&plower[0],Ls*sizeof(Coeff_t));
  
  // Flops = 6.0*(Nc*Ns) *Ls*vol
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
      coalescedWrite(chi[ss+s], pdiag[s]*phi(ss+s) + pupper[s]*tmp1 + plower[s]*tmp2);
    }
  });

}

template<class Impl>
void MobiusEOFAFermion<Impl>::M5D_shift(const FermionField &psi_i, const FermionField &phi_i, FermionField &chi_i,
					std::vector<Coeff_t> &lower, std::vector<Coeff_t> &diag, std::vector<Coeff_t> &upper,
					std::vector<Coeff_t> &shift_coeffs)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase *grid = psi_i.Grid();
  int Ls = this->Ls;
  autoView(psi , psi_i, AcceleratorRead);
  autoView(phi , phi_i, AcceleratorRead);
  autoView(chi , chi_i, AcceleratorWrite);

  auto pm  = this->pm;
  int shift_s = (pm == 1) ? (Ls-1) : 0; // s-component modified by shift operator
  
  assert(phi.Checkerboard() == psi.Checkerboard());

  auto pdiag  = &this->d_diag[0];
  auto pupper = &this->d_upper[0];
  auto plower = &this->d_lower[0];
  auto pshift_coeffs = &this->d_shift_coefficients[0];

  acceleratorCopyToDevice(&diag[0],&pdiag[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&upper[0],&pupper[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&lower[0],&plower[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&shift_coeffs[0],&pshift_coeffs[0],Ls*sizeof(Coeff_t));

  // Flops = 6.0*(Nc*Ns) *Ls*vol
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

      coalescedWrite(chi[ss+s], pdiag[s]*phi(ss+s) + pupper[s]*tmp1 +plower[s]*tmp2 + pshift_coeffs[s]*tmp);
    }
  });

}

template<class Impl>
void MobiusEOFAFermion<Impl>::M5Ddag(const FermionField &psi_i, const FermionField &phi_i, FermionField &chi_i,
				     std::vector<Coeff_t> &lower, std::vector<Coeff_t> &diag, std::vector<Coeff_t> &upper)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase *grid = psi_i.Grid();
  int Ls = this->Ls;
  autoView(psi , psi_i, AcceleratorRead);
  autoView(phi , phi_i, AcceleratorRead);
  autoView(chi , chi_i, AcceleratorWrite);

  assert(phi.Checkerboard() == psi.Checkerboard());
  
  auto pdiag  = &this->d_diag[0];
  auto pupper = &this->d_upper[0];
  auto plower = &this->d_lower[0];

  acceleratorCopyToDevice(&diag[0],&pdiag[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&upper[0],&pupper[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&lower[0],&plower[0],Ls*sizeof(Coeff_t));

  // Flops = 6.0*(Nc*Ns) *Ls*vol
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
      coalescedWrite(chi[ss+s], pdiag[s]*phi(ss+s) + pupper[s]*tmp1 + plower[s]*tmp2);
    }
  });
}

template<class Impl>
void MobiusEOFAFermion<Impl>::M5Ddag_shift(const FermionField &psi_i, const FermionField &phi_i, FermionField &chi_i,
					   std::vector<Coeff_t> &lower, std::vector<Coeff_t> &diag, std::vector<Coeff_t> &upper,
					   std::vector<Coeff_t> &shift_coeffs)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase *grid = psi_i.Grid();
  int Ls = this->Ls;
  int shift_s = (this->pm == 1) ? (Ls-1) : 0; // s-component modified by shift operator
  autoView(psi , psi_i, AcceleratorRead);
  autoView(phi , phi_i, AcceleratorRead);
  autoView(chi , chi_i, AcceleratorWrite);

  assert(phi.Checkerboard() == psi.Checkerboard());

  auto pdiag  = &this->d_diag[0];
  auto pupper = &this->d_upper[0];
  auto plower = &this->d_lower[0];
  auto pshift_coeffs = &this->d_shift_coefficients[0];

  acceleratorCopyToDevice(&diag[0],&pdiag[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&upper[0],&pupper[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&lower[0],&plower[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&shift_coeffs[0],&pshift_coeffs[0],Ls*sizeof(Coeff_t));
  
  // Flops = 6.0*(Nc*Ns) *Ls*vol
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

      if(s==(Ls-1)) coalescedWrite(chi[ss+s], chi(ss+s)+ pdiag[s]*phi(ss+s) + pupper[s]*tmp1 + plower[s]*tmp2);
      else          coalescedWrite(chi[ss+s], pdiag[s]*phi(ss+s) + pupper[s]*tmp1 + plower[s]*tmp2);
      if(pm == 1){ spProj5p(tmp, psi(ss+s)); }
      else       { spProj5m(tmp, psi(ss+s)); }

      coalescedWrite(chi[ss+shift_s],chi(ss+shift_s)+pshift_coeffs[s]*tmp);
    }
  });

}

template<class Impl>
void MobiusEOFAFermion<Impl>::MooeeInv(const FermionField &psi_i, FermionField &chi_i)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase *grid = psi_i.Grid();
  int Ls = this->Ls;
  autoView(psi , psi_i, AcceleratorRead);
  autoView(chi , chi_i, AcceleratorWrite);

  auto plee  = & this->d_lee [0];
  auto pdee  = & this->d_dee [0];
  auto puee  = & this->d_uee [0];
  auto pleem = & this->d_leem[0];
  auto pueem = & this->d_ueem[0];

  acceleratorCopyToDevice(&this->lee[0],&plee[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&this->dee[0],&pdee[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&this->uee[0],&puee[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&this->leem[0],&pleem[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&this->ueem[0],&pueem[0],Ls*sizeof(Coeff_t));

  if(this->shift != 0.0){ MooeeInv_shift(psi_i,chi_i); return; }

  int nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    uint64_t ss=sss*Ls;
    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp, acc, res;

    // X = Nc*Ns
    // flops = 2X + (Ls-2)(4X + 4X) + 6X + 1 + 2X + (Ls-1)(10X + 1) = -16X + Ls(1+18X) = -192 + 217*Ls flops
    // Apply (L^{\prime})^{-1} L_m^{-1}
    res = psi(ss);
    spProj5m(tmp,res);
    acc = pleem[0]*tmp;
    spProj5p(tmp,res);
    coalescedWrite(chi[ss],res);
    
    for(int s=1;s<Ls-1;s++){
      res = psi(ss+s);
      res -= plee[s-1]*tmp;
      spProj5m(tmp,res);
      acc += pleem[s]*tmp;
      spProj5p(tmp,res);
      coalescedWrite(chi[ss+s],res);
    }
    res = psi(ss+Ls-1) - plee[Ls-2]*tmp - acc;
    
    // Apply U_m^{-1} D^{-1} U^{-1}
    res = (1.0/pdee[Ls-1])*res;
    coalescedWrite(chi[ss+Ls-1],res);
    spProj5p(acc,res);
    spProj5m(tmp,res);
    for (int s=Ls-2;s>=0;s--){
      res = (1.0/pdee[s])*chi(ss+s) - puee[s]*tmp - pueem[s]*acc;
      spProj5m(tmp,res);
      coalescedWrite(chi[ss+s],res);
    }
  });
   
}

template<class Impl>
void MobiusEOFAFermion<Impl>::MooeeInv_shift(const FermionField &psi_i, FermionField &chi_i)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase *grid = psi_i.Grid();
  int Ls = this->Ls;
  autoView(psi , psi_i, AcceleratorRead);
  autoView(chi , chi_i, AcceleratorWrite);

  // Move into object and constructor
  auto pm = this->pm;
  auto plee  = & this->d_lee [0];
  auto pdee  = & this->d_dee [0];
  auto puee  = & this->d_uee [0];
  auto pleem = & this->d_leem[0];
  auto pueem = & this->d_ueem[0];
  auto pMooeeInv_shift_lc   = &this->d_MooeeInv_shift_lc[0];
  auto pMooeeInv_shift_norm = &this->d_MooeeInv_shift_norm[0];

  acceleratorCopyToDevice(&this->lee[0],&plee[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&this->dee[0],&pdee[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&this->uee[0],&puee[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&this->leem[0],&pleem[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&this->ueem[0],&pueem[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&MooeeInv_shift_lc[0],&pMooeeInv_shift_lc[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&MooeeInv_shift_norm[0],&pMooeeInv_shift_norm[0],Ls*sizeof(Coeff_t));

  int nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
      uint64_t ss=sss*Ls;
      typedef decltype(coalescedRead(psi[0])) spinor;
      spinor tmp, acc, res, tmp_spProj;

      // Apply (L^{\prime})^{-1} L_m^{-1}
      res = psi(ss);
      spProj5m(tmp,res);
      acc = pleem[0]*tmp;
      spProj5p(tmp,res);
      coalescedWrite(chi[ss],res);
      tmp_spProj = pMooeeInv_shift_lc[0]*res;

      for(int s=1;s<Ls-1;s++){
	res = psi(ss+s);
	tmp_spProj += pMooeeInv_shift_lc[s]*res;
	res -= plee[s-1]*tmp;
	spProj5m(tmp,res);
	acc += pleem[s]*tmp;
	spProj5p(tmp,res);
	coalescedWrite(chi[ss+s],res);
      }
      res = psi(ss+Ls-1);

      tmp_spProj += pMooeeInv_shift_lc[Ls-1]*res;
      if(pm == 1){ spProj5p(tmp_spProj, tmp_spProj);}
      else       { spProj5m(tmp_spProj, tmp_spProj); }

      res = res - plee[Ls-2]*tmp - acc;

      // Apply U_m^{-1} D^{-1} U^{-1}
      res = (1.0/pdee[Ls-1])*res;
      spProj5p(acc,res);
      spProj5m(tmp,res);
      coalescedWrite(chi[ss+Ls-1], res + pMooeeInv_shift_norm[Ls-1]*tmp_spProj);
      for (int s=Ls-2;s>=0;s--){
	res = (1.0/pdee[s])*chi(ss+s) - puee[s]*tmp - pueem[s]*acc;
	spProj5m(tmp,res);
	coalescedWrite(chi[ss+s], res + pMooeeInv_shift_norm[s]*tmp_spProj);
      }
  });

}

template<class Impl>
void MobiusEOFAFermion<Impl>::MooeeInvDag(const FermionField &psi_i, FermionField &chi_i)
{
  if(this->shift != 0.0){ MooeeInvDag_shift(psi_i,chi_i); return; }

  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase *grid = psi_i.Grid();
  int Ls = this->Ls;
  autoView(psi , psi_i, AcceleratorRead);
  autoView(chi , chi_i, AcceleratorWrite);

  auto plee  = &this->d_lee [0];
  auto pdee  = &this->d_dee [0];
  auto puee  = &this->d_uee [0];
  auto pleem = &this->d_leem[0];
  auto pueem = &this->d_ueem[0];

  acceleratorCopyToDevice(&this->lee[0],&plee[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&this->dee[0],&pdee[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&this->uee[0],&puee[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&this->leem[0],&pleem[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&this->ueem[0],&pueem[0],Ls*sizeof(Coeff_t));

  int nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    uint64_t ss=sss*Ls;
    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp, acc, res;

    // X = Nc*Ns
    // flops = 2X + (Ls-2)(4X + 4X) + 6X + 1 + 2X + (Ls-1)(10X + 1) = -16X + Ls(1+18X) = -192 + 217*Ls flops
    // Apply (U^{\prime})^{-dagger} U_m^{-\dagger}
    res = psi(ss);
    spProj5p(tmp,res);
    acc = pueem[0]*tmp;
    spProj5m(tmp,res);
    coalescedWrite(chi[ss],res);
    
    for(int s=1;s<Ls-1;s++){
      res = psi(ss+s);
      res -= puee[s-1]*tmp;
      spProj5p(tmp,res);
      acc += pueem[s]*tmp;
      spProj5m(tmp,res);
      coalescedWrite(chi[ss+s],res);
    }
    res = psi(ss+Ls-1) - puee[Ls-2]*tmp - acc;
    
    // Apply L_m^{-\dagger} D^{-dagger} L^{-dagger}
    res = (1.0/pdee[Ls-1])*res;
    coalescedWrite(chi[ss+Ls-1],res);
    spProj5m(acc,res);
    spProj5p(tmp,res);
    for (int s=Ls-2;s>=0;s--){
      res = (1.0/pdee[s])*chi(ss+s) - plee[s]*tmp - pleem[s]*acc;
      spProj5p(tmp,res);
      coalescedWrite(chi[ss+s],res);
    }
  });
}

template<class Impl>
void MobiusEOFAFermion<Impl>::MooeeInvDag_shift(const FermionField &psi_i, FermionField &chi_i)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase *grid = psi_i.Grid();
  autoView(psi , psi_i, AcceleratorRead);
  autoView(chi , chi_i, AcceleratorWrite);
  int Ls = this->Ls;

  auto pm = this->pm;
  auto plee  = & this->d_lee [0];
  auto pdee  = & this->d_dee [0];
  auto puee  = & this->d_uee [0];
  auto pleem = & this->d_leem[0];
  auto pueem = & this->d_ueem[0];

  auto pMooeeInvDag_shift_lc   = &this->d_MooeeInv_shift_lc[0];
  auto pMooeeInvDag_shift_norm = &this->d_MooeeInv_shift_norm[0];

  acceleratorCopyToDevice(&this->lee[0],&plee[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&this->dee[0],&pdee[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&this->uee[0],&puee[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&this->leem[0],&pleem[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&this->ueem[0],&pueem[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&MooeeInvDag_shift_lc[0],&pMooeeInvDag_shift_lc[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&MooeeInvDag_shift_norm[0],&pMooeeInvDag_shift_norm[0],Ls*sizeof(Coeff_t));

  //  auto pMooeeInvDag_shift_lc   = &MooeeInvDag_shift_lc[0];
  //  auto pMooeeInvDag_shift_norm = &MooeeInvDag_shift_norm[0];

  int nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
      uint64_t ss=sss*Ls;
      typedef decltype(coalescedRead(psi[0])) spinor;
      spinor tmp, acc, res, tmp_spProj;

      // Apply (U^{\prime})^{-dagger} U_m^{-\dagger}
      res = psi(ss);
      spProj5p(tmp,res);
      acc = pueem[0]*tmp;
      spProj5m(tmp,res);
      coalescedWrite(chi[ss],res);
      tmp_spProj = pMooeeInvDag_shift_lc[0]*res;

      for(int s=1;s<Ls-1;s++){
	res = psi(ss+s);
	tmp_spProj += pMooeeInvDag_shift_lc[s]*res;
	res -= puee[s-1]*tmp;
	spProj5p(tmp,res);
	acc += pueem[s]*tmp;
	spProj5m(tmp,res);
	coalescedWrite(chi[ss+s],res);
      }
      res = psi(ss+Ls-1);

      tmp_spProj += pMooeeInvDag_shift_lc[Ls-1]*res;
      if(pm == 1){ spProj5p(tmp_spProj, tmp_spProj); }
      else       { spProj5m(tmp_spProj, tmp_spProj); }

      res = res - puee[Ls-2]*tmp - acc;

      // Apply L_m^{-\dagger} D^{-dagger} L^{-dagger}
      res = (1.0/pdee[Ls-1])*res;
      spProj5m(acc,res);
      spProj5p(tmp,res);
      coalescedWrite(chi[ss+Ls-1], res + pMooeeInvDag_shift_norm[Ls-1]*tmp_spProj);
      for (int s=Ls-2;s>=0;s--){
	res = (1.0/pdee[s])*chi(ss+s) - plee[s]*tmp - pleem[s]*acc;
	spProj5p(tmp,res);
	coalescedWrite(chi[ss+s], res + pMooeeInvDag_shift_norm[s]*tmp_spProj);
      }
  });

}

NAMESPACE_END(Grid);
