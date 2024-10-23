/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/DomainWallEOFAFermioncache.cc

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
#include <Grid/qcd/action/fermion/DomainWallEOFAFermion.h>

NAMESPACE_BEGIN(Grid);

// FIXME -- make a version of these routines with site loop outermost for cache reuse.
// Pminus fowards
// Pplus  backwards..
template<class Impl>
void DomainWallEOFAFermion<Impl>::M5D(const FermionField& psi_i, const FermionField& phi_i,FermionField& chi_i, 
				      std::vector<Coeff_t>& lower, std::vector<Coeff_t>& diag, std::vector<Coeff_t>& upper)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  int Ls = this->Ls;
  GridBase* grid = psi_i.Grid();
  autoView( phi , phi_i, AcceleratorRead);
  autoView( psi , psi_i, AcceleratorRead);
  autoView( chi , chi_i, AcceleratorWrite);
  assert(phi.Checkerboard() == psi.Checkerboard());

  auto pdiag  = &this->d_diag[0];
  auto pupper = &this->d_upper[0];
  auto plower = &this->d_lower[0];

  acceleratorCopyToDevice(&diag[0],&pdiag[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&upper[0],&pupper[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&lower[0],&plower[0],Ls*sizeof(Coeff_t));

  // Flops = 6.0*(Nc*Ns) *Ls*vol
  
  auto nloop=grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    auto ss=sss*Ls;
    typedef decltype(coalescedRead(psi[0])) spinor;
    for(int s=0; s<Ls; s++){
      spinor tmp1, tmp2;
      uint64_t idx_u = ss+((s+1)%Ls);
      uint64_t idx_l = ss+((s+Ls-1)%Ls);
      spProj5m(tmp1, psi(idx_u));
      spProj5p(tmp2, psi(idx_l));
      coalescedWrite(chi[ss+s], pdiag[s]*phi(ss+s) + pupper[s]*tmp1 + plower[s]*tmp2);
    }
  });

}

template<class Impl>
void DomainWallEOFAFermion<Impl>::M5Ddag(const FermionField& psi_i, const FermionField& phi_i, FermionField& chi_i, 
					 std::vector<Coeff_t>& lower, std::vector<Coeff_t>& diag, std::vector<Coeff_t>& upper)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase* grid = psi_i.Grid();
  int Ls = this->Ls;

  autoView( psi , psi_i, AcceleratorRead);
  autoView( phi , phi_i, AcceleratorRead);
  autoView( chi , chi_i, AcceleratorWrite);
  assert(phi.Checkerboard() == psi.Checkerboard());
  
  auto pdiag  = &this->d_diag[0];
  auto pupper = &this->d_upper[0];
  auto plower = &this->d_lower[0];

  acceleratorCopyToDevice(&diag[0] ,&pdiag[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&upper[0],&pupper[0],Ls*sizeof(Coeff_t));
  acceleratorCopyToDevice(&lower[0],&plower[0],Ls*sizeof(Coeff_t));

  // Flops = 6.0*(Nc*Ns) *Ls*vol

  auto nloop=grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    typedef decltype(coalescedRead(psi[0])) spinor;
    auto ss=sss*Ls;
    for(int s=0; s<Ls; s++){
      spinor tmp1, tmp2;
      uint64_t idx_u = ss+((s+1)%Ls);
      uint64_t idx_l = ss+((s+Ls-1)%Ls);
      spProj5p(tmp1, psi(idx_u));
      spProj5m(tmp2, psi(idx_l));
      coalescedWrite(chi[ss+s], pdiag[s]*phi(ss+s) + pupper[s]*tmp1 + plower[s]*tmp2);
    }
  });

}

template<class Impl>
void DomainWallEOFAFermion<Impl>::MooeeInv(const FermionField& psi_i, FermionField& chi_i)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase* grid = psi_i.Grid();
  autoView( psi, psi_i, AcceleratorRead);
  autoView( chi, chi_i, AcceleratorWrite);
  int Ls = this->Ls;

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

  uint64_t nloop=grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    uint64_t ss=sss*Ls;
    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp, acc, res;

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
    acc = (1.0/pdee[Ls  ])*res;
    tmp = (1.0/pdee[Ls-1])*res;
    spProj5p(acc,acc);
    spProj5m(tmp,tmp);
    coalescedWrite(chi[ss+Ls-1], acc + tmp);
    for (int s=Ls-2;s>=0;s--){
      res = (1.0/pdee[s])*chi(ss+s) - puee[s]*tmp - pueem[s]*acc;
      spProj5m(tmp,res);
      coalescedWrite(chi[ss+s],res);
    }
  });
}

template<class Impl>
void DomainWallEOFAFermion<Impl>::MooeeInvDag(const FermionField& psi_i, FermionField& chi_i)
{
  chi_i.Checkerboard() = psi_i.Checkerboard();
  GridBase* grid = psi_i.Grid();
  autoView( psi, psi_i, AcceleratorRead);
  autoView( chi, chi_i, AcceleratorWrite);
  int Ls = this->Ls;

  auto plee  = & this->lee[0];
  auto pdee  = & this->dee[0];
  auto puee  = & this->uee[0];

  auto pleem = & this->leem[0];
  auto pueem = & this->ueem[0];

  assert(psi.Checkerboard() == psi.Checkerboard());

  auto nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    uint64_t ss=sss*Ls;
    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp, acc, res;

    // Apply (U^{\prime})^{-dagger} U_m^{-\dagger} 
    res = psi(ss);
    spProj5p(tmp,res);
    acc = conjugate(pueem[0])*tmp;
    spProj5m(tmp,res);
    coalescedWrite(chi[ss],res);
    
    for(int s=1;s<Ls-1;s++){
      res = psi(ss+s);
      res -= conjugate(puee[s-1])*tmp;
      spProj5p(tmp,res);
      acc += conjugate(pueem[s])*tmp;
      spProj5m(tmp,res);
      coalescedWrite(chi[ss+s],res);
    }
    res = psi(ss+Ls-1) - conjugate(puee[Ls-2])*tmp - acc;
    
    // Apply L_m^{-\dagger} D^{-dagger} L^{-dagger}
    acc = conjugate(1.0/pdee[Ls-1])*res;
    tmp = conjugate(1.0/pdee[Ls  ])*res;
    spProj5m(acc,acc);
    spProj5p(tmp,tmp);
    coalescedWrite(chi[ss+Ls-1], acc + tmp);
    for (int s=Ls-2;s>=0;s--){
      res = conjugate(1.0/pdee[s])*chi(ss+s) - conjugate(plee[s])*tmp - conjugate(pleem[s])*acc;
      spProj5p(tmp,res);
      coalescedWrite(chi[ss+s],res);
    }
  });

}

NAMESPACE_END(Grid);
