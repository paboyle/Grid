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
void CayleyFermion5D<Impl>::M5D(const FermionField &psi_i,
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
  Coeff_t *lower_v = &lower[0];
  Coeff_t *diag_v  = &diag[0];
  Coeff_t *upper_v = &upper[0];
  int Ls =this->Ls;
  assert(phi.Checkerboard() == psi.Checkerboard());
  
  const uint64_t nsimd = grid->Nsimd();
  const uint64_t sites4d = nsimd * grid->oSites() / Ls;
  
  // 10 = 3 complex mult + 2 complex add
  // Flops = 10.0*(Nc*Ns) *Ls*vol (/2 for red black counting)
  M5Dcalls++;
  M5Dtime-=usecond();
  
  accelerator_loopN( sss, sites4d ,{
    uint64_t lane = sss % nsimd;
    uint64_t ss   = Ls * (sss / nsimd);
    
    for(int s=0;s<Ls;s++){
      auto res = extractLane(lane,phi[ss+s]);
      res = diag_v[s]*res;
      
      auto tmp = extractLane(lane,psi[ss+(s+1)%Ls]);
      spProj5m(tmp,tmp);
      res += upper_v[s]*tmp;
      
      tmp = extractLane(lane,psi[ss+(s+Ls-1)%Ls]);
      spProj5p(tmp,tmp);
      res += lower_v[s]*tmp;
      
      insertLane(lane,chi[ss+s],res);
    }
  });
  M5Dtime+=usecond();
}

template<class Impl>  
void CayleyFermion5D<Impl>::M5Ddag(const FermionField &psi_i,
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
  Coeff_t *lower_v = &lower[0];
  Coeff_t *diag_v  = &diag[0];
  Coeff_t *upper_v = &upper[0];
  int Ls =this->Ls;
  assert(phi.Checkerboard() == psi.Checkerboard());

  const uint64_t nsimd = grid->Nsimd();
  const uint64_t sites4d = nsimd * grid->oSites() / Ls;
  
  // 10 = 3 complex mult + 2 complex add
  // Flops = 10.0*(Nc*Ns) *Ls*vol (/2 for red black counting)
  M5Dcalls++;
  M5Dtime-=usecond();

  accelerator_loopN( sss, sites4d ,{
    uint64_t lane = sss % nsimd;
    uint64_t ss   = Ls * (sss / nsimd);
    
    for(int s=0;s<Ls;s++){
      auto res = extractLane(lane,phi[ss+s]);
      res = diag_v[s]*res;
      
      auto tmp = extractLane(lane,psi[ss+(s+1)%Ls]);
      spProj5p(tmp,tmp);
      res += upper_v[s]*tmp;
      
      tmp = extractLane(lane,psi[ss+(s+Ls-1)%Ls]);
      spProj5m(tmp,tmp);
      res += lower_v[s]*tmp;
      
      insertLane(lane,chi[ss+s],res);
    }
  });
  M5Dtime+=usecond();
}

template<class Impl>
void CayleyFermion5D<Impl>::MooeeInv    (const FermionField &psi_i, FermionField &chi_i)
{
  chi_i.Checkerboard()=psi_i.Checkerboard();
  GridBase *grid=psi_i.Grid();

  auto psi = psi_i.View();
  auto chi = chi_i.View();

  int Ls=this->Ls;

  MooeeInvCalls++;
  MooeeInvTime-=usecond();

  thread_loop((int ss=0;ss<grid->oSites();ss+=Ls),{ // adds Ls
    auto tmp = psi[0];

    // flops = 12*2*Ls + 12*2*Ls + 3*12*Ls + 12*2*Ls  = 12*Ls * (9) = 108*Ls flops
    // Apply (L^{\prime})^{-1}
    chi[ss]=psi[ss]; // chi[0]=psi[0]
    for(int s=1;s<Ls;s++){
      spProj5p(tmp,chi[ss+s-1]);  
      chi[ss+s] = psi[ss+s]-lee[s-1]*tmp;
    }
    // L_m^{-1} 
    for (int s=0;s<Ls-1;s++){ // Chi[ee] = 1 - sum[s<Ls-1] -leem[s]P_- chi
      spProj5m(tmp,chi[ss+s]);    
      chi[ss+Ls-1] = chi[ss+Ls-1] - leem[s]*tmp;
    }
    // U_m^{-1} D^{-1}
    for (int s=0;s<Ls-1;s++){
      // Chi[s] + 1/d chi[s] 
      spProj5p(tmp,chi[ss+Ls-1]); 
      chi[ss+s] = (1.0/dee[s])*chi[ss+s]-(ueem[s]/dee[Ls-1])*tmp;
    }	
    chi[ss+Ls-1]= (1.0/dee[Ls-1])*chi[ss+Ls-1];
      
    // Apply U^{-1}
    for (int s=Ls-2;s>=0;s--){
      spProj5m(tmp,chi[ss+s+1]);  
      chi[ss+s] = chi[ss+s] - uee[s]*tmp;
    }
  });

  MooeeInvTime+=usecond();

}

template<class Impl>
void CayleyFermion5D<Impl>::MooeeInvDag (const FermionField &psi_i, FermionField &chi_i)
{
  chi_i.Checkerboard()=psi_i.Checkerboard();
  GridBase *grid=psi_i.Grid();
  int Ls=this->Ls;

  auto psi = psi_i.View();
  auto chi = chi_i.View();

  assert(psi.Checkerboard() == psi.Checkerboard());

  MooeeInvCalls++;
  MooeeInvTime-=usecond();

  thread_loop((int ss=0;ss<grid->oSites();ss+=Ls),{ // adds Ls

    auto tmp = psi[0];

    // Apply (U^{\prime})^{-dagger}
    chi[ss]=psi[ss];
    for (int s=1;s<Ls;s++){
      spProj5m(tmp,chi[ss+s-1]);
      chi[ss+s] = psi[ss+s]-conjugate(uee[s-1])*tmp;
    }
    // U_m^{-\dagger} 
    for (int s=0;s<Ls-1;s++){
      spProj5p(tmp,chi[ss+s]);
      chi[ss+Ls-1] = chi[ss+Ls-1] - conjugate(ueem[s])*tmp;
    }

    // L_m^{-\dagger} D^{-dagger}
    for (int s=0;s<Ls-1;s++){
      spProj5m(tmp,chi[ss+Ls-1]);
      chi[ss+s] = conjugate(1.0/dee[s])*chi[ss+s]-conjugate(leem[s]/dee[Ls-1])*tmp;
    }	
    chi[ss+Ls-1]= conjugate(1.0/dee[Ls-1])*chi[ss+Ls-1];
  
    // Apply L^{-dagger}
    for (int s=Ls-2;s>=0;s--){
      spProj5p(tmp,chi[ss+s+1]);
      chi[ss+s] = chi[ss+s] - conjugate(lee[s])*tmp;
    }
  });

  MooeeInvTime+=usecond();

}

#ifdef CAYLEY_DPERP_GPU
INSTANTIATE_DPERP(WilsonImplF);
INSTANTIATE_DPERP(WilsonImplD);
INSTANTIATE_DPERP(GparityWilsonImplF);
INSTANTIATE_DPERP(GparityWilsonImplD);
INSTANTIATE_DPERP(ZWilsonImplF);
INSTANTIATE_DPERP(ZWilsonImplD);

INSTANTIATE_DPERP(WilsonImplFH);
INSTANTIATE_DPERP(WilsonImplDF);
INSTANTIATE_DPERP(GparityWilsonImplFH);
INSTANTIATE_DPERP(GparityWilsonImplDF);
INSTANTIATE_DPERP(ZWilsonImplFH);
INSTANTIATE_DPERP(ZWilsonImplDF);
#endif

NAMESPACE_END(Grid);
