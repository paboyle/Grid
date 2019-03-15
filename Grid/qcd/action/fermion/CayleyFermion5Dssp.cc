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


namespace Grid {
namespace QCD {

  // FIXME -- make a version of these routines with site loop outermost for cache reuse.
  // Pminus fowards
  // Pplus  backwards
template<class Impl>  
void CayleyFermion5D<Impl>::M5D(const FermionField &psi,
				const FermionField &phi, 
				FermionField &chi,
				std::vector<Coeff_t> &lower,
				std::vector<Coeff_t> &diag,
				std::vector<Coeff_t> &upper)
{
  Coeff_t one(1.0);
  int Ls=this->Ls;
  for(int s=0;s<Ls;s++){
    if ( s==0 ) {
      axpby_ssp_pminus(chi,diag[s],phi,upper[s],psi,s,s+1);
      axpby_ssp_pplus (chi,one,chi,lower[s],psi,s,Ls-1);
    } else if ( s==(Ls-1)) { 
      axpby_ssp_pminus(chi,diag[s],phi,upper[s],psi,s,0);
      axpby_ssp_pplus (chi,one,chi,lower[s],psi,s,s-1);
    } else {
      axpby_ssp_pminus(chi,diag[s],phi,upper[s],psi,s,s+1);
      axpby_ssp_pplus(chi,one,chi,lower[s],psi,s,s-1);
    }
  }
}
template<class Impl>  
void CayleyFermion5D<Impl>::M5Ddag(const FermionField &psi,
				   const FermionField &phi, 
				   FermionField &chi,
				   std::vector<Coeff_t> &lower,
				   std::vector<Coeff_t> &diag,
				   std::vector<Coeff_t> &upper)
{
  Coeff_t one(1.0);
  int Ls=this->Ls;
  for(int s=0;s<Ls;s++){
    if ( s==0 ) {
      axpby_ssp_pplus (chi,diag[s],phi,upper[s],psi,s,s+1);
      axpby_ssp_pminus(chi,one,chi,lower[s],psi,s,Ls-1);
    } else if ( s==(Ls-1)) { 
      axpby_ssp_pplus (chi,diag[s],phi,upper[s],psi,s,0);
      axpby_ssp_pminus(chi,one,chi,lower[s],psi,s,s-1);
    } else {
      axpby_ssp_pplus (chi,diag[s],phi,upper[s],psi,s,s+1);
      axpby_ssp_pminus(chi,one,chi,lower[s],psi,s,s-1);
    }
  }
}

template<class Impl>
void CayleyFermion5D<Impl>::MooeeInv    (const FermionField &psi, FermionField &chi)
{
  Coeff_t one(1.0);
  Coeff_t czero(0.0);
  chi.checkerboard=psi.checkerboard;
  int Ls=this->Ls;
  // Apply (L^{\prime})^{-1}
  axpby_ssp (chi,one,psi,     czero,psi,0,0);      // chi[0]=psi[0]
  for (int s=1;s<Ls;s++){
    axpby_ssp_pplus(chi,one,psi,-lee[s-1],chi,s,s-1);// recursion Psi[s] -lee P_+ chi[s-1]
  }
  // L_m^{-1} 
  for (int s=0;s<Ls-1;s++){ // Chi[ee] = 1 - sum[s<Ls-1] -leem[s]P_- chi
    axpby_ssp_pminus(chi,one,chi,-leem[s],chi,Ls-1,s);
  }
  // U_m^{-1} D^{-1}
  for (int s=0;s<Ls-1;s++){
    // Chi[s] + 1/d chi[s] 
    axpby_ssp_pplus(chi,one/dee[s],chi,-ueem[s]/dee[Ls-1],chi,s,Ls-1);
  }	
  axpby_ssp(chi,one/dee[Ls-1],chi,czero,chi,Ls-1,Ls-1); // Modest avoidable 
  
  // Apply U^{-1}
  for (int s=Ls-2;s>=0;s--){
    axpby_ssp_pminus (chi,one,chi,-uee[s],chi,s,s+1);  // chi[Ls]
  }
}

template<class Impl>
void CayleyFermion5D<Impl>::MooeeInvDag (const FermionField &psi, FermionField &chi)
{
  Coeff_t one(1.0);
  Coeff_t czero(0.0);
  chi.checkerboard=psi.checkerboard;
  int Ls=this->Ls;
  // Apply (U^{\prime})^{-dagger}
  axpby_ssp (chi,one,psi,     czero,psi,0,0);      // chi[0]=psi[0]
  for (int s=1;s<Ls;s++){
    axpby_ssp_pminus(chi,one,psi,-conjugate(uee[s-1]),chi,s,s-1);
  }
  // U_m^{-\dagger} 
  for (int s=0;s<Ls-1;s++){
    axpby_ssp_pplus(chi,one,chi,-conjugate(ueem[s]),chi,Ls-1,s);
  }
  // L_m^{-\dagger} D^{-dagger}
  for (int s=0;s<Ls-1;s++){
    axpby_ssp_pminus(chi,one/conjugate(dee[s]),chi,-conjugate(leem[s]/dee[Ls-1]),chi,s,Ls-1);
  }	
  axpby_ssp(chi,one/conjugate(dee[Ls-1]),chi,czero,chi,Ls-1,Ls-1); // Modest avoidable 
  
  // Apply L^{-dagger}
  for (int s=Ls-2;s>=0;s--){
    axpby_ssp_pplus (chi,one,chi,-conjugate(lee[s]),chi,s,s+1);  // chi[Ls]
  }
}


#ifdef CAYLEY_DPERP_LINALG
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

}
}
