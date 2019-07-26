/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/CayleyFermion5D.cc

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/qcd/action/fermion/CayleyFermion5D.h>


NAMESPACE_BEGIN(Grid);

// Pminus fowards
// Pplus  backwards..
template<class Impl>  
void
CayleyFermion5D<Impl>::M5D(const FermionField &psi_i,
			       const FermionField &phi_i, 
			       FermionField &chi_i,
			       Vector<Coeff_t> &lower,
			       Vector<Coeff_t> &diag,
			       Vector<Coeff_t> &upper)
{
  
  chi_i.Checkerboard()=psi_i.Checkerboard();
  GridBase *grid=psi_i.Grid();
  auto psi = psi_i.View();
  auto phi = phi_i.View();
  auto chi = chi_i.View();
  assert(phi.Checkerboard() == psi.Checkerboard());

  int Ls =this->Ls;

  // 10 = 3 complex mult + 2 complex add
  // Flops = 10.0*(Nc*Ns) *Ls*vol (/2 for red black counting)
  M5Dcalls++;
  M5Dtime-=usecond();

  uint64_t nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    uint64_t ss= sss*Ls;
    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp1, tmp2;
    for(int s=0;s<Ls;s++){
      uint64_t idx_u = ss+((s+1)%Ls);
      uint64_t idx_l = ss+((s+Ls-1)%Ls);
      spProj5m(tmp1,psi(idx_u));
      spProj5p(tmp2,psi(idx_l));
      coalescedWrite(chi[ss+s],diag[s]*phi(ss+s)+upper[s]*tmp1+lower[s]*tmp2);
    }
  });
  M5Dtime+=usecond();
}

template<class Impl>  
void
CayleyFermion5D<Impl>::M5Ddag(const FermionField &psi_i,
			      const FermionField &phi_i, 
			      FermionField &chi_i,
			      Vector<Coeff_t> &lower,
			      Vector<Coeff_t> &diag,
			      Vector<Coeff_t> &upper)
{
  chi_i.Checkerboard()=psi_i.Checkerboard();
  GridBase *grid=psi_i.Grid();
  auto psi = psi_i.View();
  auto phi = phi_i.View();
  auto chi = chi_i.View();
  assert(phi.Checkerboard() == psi.Checkerboard());

  int Ls=this->Ls;

  // Flops = 6.0*(Nc*Ns) *Ls*vol
  M5Dcalls++;
  M5Dtime-=usecond();

  uint64_t nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    uint64_t ss=sss*Ls;
    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp1,tmp2;
    for(int s=0;s<Ls;s++){
      uint64_t idx_u = ss+((s+1)%Ls);
      uint64_t idx_l = ss+((s+Ls-1)%Ls);
      spProj5p(tmp1,psi(idx_u));
      spProj5m(tmp2,psi(idx_l));
      coalescedWrite(chi[ss+s],diag[s]*phi(ss+s)+upper[s]*tmp1+lower[s]*tmp2);
    }
  });
  M5Dtime+=usecond();
}

template<class Impl>
void
CayleyFermion5D<Impl>::MooeeInv    (const FermionField &psi_i, FermionField &chi_i)
{
  chi_i.Checkerboard()=psi_i.Checkerboard();
  GridBase *grid=psi_i.Grid();

  auto psi = psi_i.View();
  auto chi = chi_i.View();

  int Ls=this->Ls;

  auto plee  = & lee [0];
  auto pdee  = & dee [0];
  auto puee  = & uee [0];
  auto pleem = & leem[0];
  auto pueem = & ueem[0];

  MooeeInvCalls++;
  MooeeInvTime-=usecond();
  uint64_t nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    uint64_t ss=sss*Ls;
    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp;

    // flops = 12*2*Ls + 12*2*Ls + 3*12*Ls + 12*2*Ls  = 12*Ls * (9) = 108*Ls flops
    // Apply (L^{\prime})^{-1}
    coalescedWrite(chi[ss],psi(ss)); // chi[0]=psi[0]
    for(int s=1;s<Ls;s++){
      spProj5p(tmp,chi(ss+s-1));  
      coalescedWrite(chi[ss+s] , psi(ss+s)-plee[s-1]*tmp);
    }

    // L_m^{-1} 
    for (int s=0;s<Ls-1;s++){ // Chi[ee] = 1 - sum[s<Ls-1] -pleem[s]P_- chi
      spProj5m(tmp,chi(ss+s));    
      coalescedWrite(chi[ss+Ls-1], chi(ss+Ls-1) - pleem[s]*tmp);
    }

    // U_m^{-1} D^{-1}
    for (int s=0;s<Ls-1;s++){
      // Chi[s] + 1/d chi[s] 
      spProj5p(tmp,chi(ss+Ls-1)); 
      coalescedWrite(chi[ss+s], (1.0/pdee[s])*chi(ss+s)-(pueem[s]/pdee[Ls-1])*tmp);
    }	
    coalescedWrite(chi[ss+Ls-1], (1.0/pdee[Ls-1])*chi(ss+Ls-1));
      
    // Apply U^{-1}
    for (int s=Ls-2;s>=0;s--){
      spProj5m(tmp,chi(ss+s+1));  
      coalescedWrite(chi[ss+s], chi(ss+s) - puee[s]*tmp);
    }
  });

  MooeeInvTime+=usecond();

}

template<class Impl>
void
CayleyFermion5D<Impl>::MooeeInvDag (const FermionField &psi_i, FermionField &chi_i)
{
  chi_i.Checkerboard()=psi_i.Checkerboard();
  GridBase *grid=psi_i.Grid();
  int Ls=this->Ls;

  auto psi = psi_i.View();
  auto chi = chi_i.View();

  auto plee  = & lee [0];
  auto pdee  = & dee [0];
  auto puee  = & uee [0];
  auto pleem = & leem[0];
  auto pueem = & ueem[0];

  assert(psi.Checkerboard() == psi.Checkerboard());

  MooeeInvCalls++;
  MooeeInvTime-=usecond();


  uint64_t nloop = grid->oSites()/Ls;
  accelerator_for(sss,nloop,Simd::Nsimd(),{
    uint64_t ss=sss*Ls;
    typedef decltype(coalescedRead(psi[0])) spinor;
    spinor tmp;

    // Apply (U^{\prime})^{-dagger}
    coalescedWrite(chi[ss],psi(ss));
    for (int s=1;s<Ls;s++){
      spProj5m(tmp,chi(ss+s-1));
      coalescedWrite(chi[ss+s], psi(ss+s)-conjugate(puee[s-1])*tmp);
    }
    // U_m^{-\dagger} 
    for (int s=0;s<Ls-1;s++){
      spProj5p(tmp,chi(ss+s));
      coalescedWrite(chi[ss+Ls-1], chi(ss+Ls-1) - conjugate(pueem[s])*tmp);
    }

    // L_m^{-\dagger} D^{-dagger}
    for (int s=0;s<Ls-1;s++){
      spProj5m(tmp,chi(ss+Ls-1));
      coalescedWrite(chi[ss+s], conjugate(1.0/pdee[s])*chi(ss+s)-conjugate(pleem[s]/pdee[Ls-1])*tmp);
    }	
    coalescedWrite(chi[ss+Ls-1], conjugate(1.0/pdee[Ls-1])*chi(ss+Ls-1));
  
    // Apply L^{-dagger}
    for (int s=Ls-2;s>=0;s--){
      spProj5p(tmp,chi(ss+s+1));
      coalescedWrite(chi[ss+s], chi(ss+s) - conjugate(plee[s])*tmp);
    }
  });
  MooeeInvTime+=usecond();

}

NAMESPACE_END(Grid);
