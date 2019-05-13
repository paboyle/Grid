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
  Coeff_t *lee_v  = &lee[0];
  Coeff_t *leem_v = &leem[0];
  Coeff_t *uee_v  = &uee[0];
  Coeff_t *ueem_v = &ueem[0];
  Coeff_t *dee_v  = &dee[0];
  
  int Ls=this->Ls;
  const uint64_t nsimd = grid->Nsimd();
  const uint64_t sites4d = nsimd * grid->oSites() / Ls;
  
  typedef typename SiteSpinor::scalar_object ScalarSiteSpinor;
  
  MooeeInvCalls++;
  MooeeInvTime-=usecond();
  
  accelerator_loopN( sss, sites4d ,{
    uint64_t lane = sss % nsimd;
    uint64_t ss   = Ls * (sss / nsimd);
    ScalarSiteSpinor res, tmp, acc;
    
    // X = Nc*Ns
    // flops = 2X + (Ls-2)(4X + 4X) + 6X + 1 + 2X + (Ls-1)(10X + 1) = -16X + Ls(1+18X) = -192 + 217*Ls flops
    // Apply (L^{\prime})^{-1} L_m^{-1}
    res = extractLane(lane,psi[ss]);
    spProj5m(tmp,res);
    acc = leem_v[0]*tmp;
    spProj5p(tmp,res);
    insertLane(lane,chi[ss],res);
    
    for(int s=1;s<Ls-1;s++){
      res = extractLane(lane,psi[ss+s]);
      res -= lee_v[s-1]*tmp;
      spProj5m(tmp,res);
      acc += leem_v[s]*tmp;
      spProj5p(tmp,res);
      insertLane(lane,chi[ss+s],res);
    }
    res = extractLane(lane,psi[ss+Ls-1]);
    res = res - lee_v[Ls-2]*tmp - acc;
    
    // Apply U_m^{-1} D^{-1} U^{-1}
    res = (1.0/dee_v[Ls-1])*res;
    insertLane(lane,chi[ss+Ls-1],res);
    spProj5p(acc,res);
    spProj5m(tmp,res);
    for (int s=Ls-2;s>=0;s--){
      res = extractLane(lane,chi[ss+s]);
      res = (1.0/dee_v[s])*res - uee_v[s]*tmp - ueem_v[s]*acc;
      spProj5m(tmp,res);
      insertLane(lane,chi[ss+s],res);
    }
  });
  
  MooeeInvTime+=usecond();

}

template<class Impl>
void CayleyFermion5D<Impl>::MooeeInvDag (const FermionField &psi_i, FermionField &chi_i)
{
  chi_i.Checkerboard()=psi_i.Checkerboard();
  GridBase *grid=psi_i.Grid();
  
  auto psi = psi_i.View();
  auto chi = chi_i.View();
  Coeff_t *lee_v  = &lee[0];
  Coeff_t *leem_v = &leem[0];
  Coeff_t *uee_v  = &uee[0];
  Coeff_t *ueem_v = &ueem[0];
  Coeff_t *dee_v  = &dee[0];
  
  int Ls=this->Ls;
  const uint64_t nsimd = grid->Nsimd();
  const uint64_t sites4d = nsimd * grid->oSites() / Ls;
  
  typedef typename SiteSpinor::scalar_object ScalarSiteSpinor;
  
  MooeeInvCalls++;
  MooeeInvTime-=usecond();
  
  accelerator_loopN( sss, sites4d ,{
    uint64_t lane = sss % nsimd;
    uint64_t ss   = Ls * (sss / nsimd);
    ScalarSiteSpinor res, tmp, acc;
    
    // X = Nc*Ns
    // flops = 2X + (Ls-2)(4X + 4X) + 6X + 1 + 2X + (Ls-1)(10X + 1) = -16X + Ls(1+18X) = -192 + 217*Ls flops
    // Apply (U^{\prime})^{-dagger} U_m^{-\dagger}
    res = extractLane(lane,psi[ss]);
    spProj5p(tmp,res);
    acc = conjugate(ueem_v[0])*tmp;
    spProj5m(tmp,res);
    insertLane(lane,chi[ss],res);
    
    for(int s=1;s<Ls-1;s++){
      res = extractLane(lane,psi[ss+s]);
      res -= conjugate(uee_v[s-1])*tmp;
      spProj5p(tmp,res);
      acc += conjugate(ueem_v[s])*tmp;
      spProj5m(tmp,res);
      insertLane(lane,chi[ss+s],res);
    }
    res = extractLane(lane,psi[ss+Ls-1]);
    res = res - conjugate(uee_v[Ls-2])*tmp - acc;
    
    // Apply L_m^{-\dagger} D^{-dagger} L^{-dagger}
    res = conjugate(1.0/dee_v[Ls-1])*res;
    insertLane(lane,chi[ss+Ls-1],res);
    spProj5m(acc,res);
    spProj5p(tmp,res);
    for (int s=Ls-2;s>=0;s--){
      res = extractLane(lane,chi[ss+s]);
      res = conjugate(1.0/dee_v[s])*res - conjugate(lee_v[s])*tmp - conjugate(leem_v[s])*acc;
      spProj5p(tmp,res);
      insertLane(lane,chi[ss+s],res);
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
