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

#include <Grid/Grid.h>


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
  int Ls=this->Ls;
  for(int s=0;s<Ls;s++){
    if ( s==0 ) {
      axpby_ssp_pminus(chi,diag[s],phi,upper[s],psi,s,s+1);
      axpby_ssp_pplus (chi,1.0,chi,lower[s],psi,s,Ls-1);
    } else if ( s==(Ls-1)) { 
      axpby_ssp_pminus(chi,diag[s],phi,upper[s],psi,s,0);
      axpby_ssp_pplus (chi,1.0,chi,lower[s],psi,s,s-1);
    } else {
      axpby_ssp_pminus(chi,diag[s],phi,upper[s],psi,s,s+1);
      axpby_ssp_pplus(chi,1.0,chi,lower[s],psi,s,s-1);
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
  int Ls=this->Ls;
  for(int s=0;s<Ls;s++){
    if ( s==0 ) {
      axpby_ssp_pplus (chi,diag[s],phi,upper[s],psi,s,s+1);
      axpby_ssp_pminus(chi,1.0,chi,lower[s],psi,s,Ls-1);
    } else if ( s==(Ls-1)) { 
      axpby_ssp_pplus (chi,diag[s],phi,upper[s],psi,s,0);
      axpby_ssp_pminus(chi,1.0,chi,lower[s],psi,s,s-1);
    } else {
      axpby_ssp_pplus (chi,diag[s],phi,upper[s],psi,s,s+1);
      axpby_ssp_pminus(chi,1.0,chi,lower[s],psi,s,s-1);
    }
  }
}

template<class Impl>
void CayleyFermion5D<Impl>::MooeeInv    (const FermionField &psi, FermionField &chi)
{
  chi.checkerboard=psi.checkerboard;
  int Ls=this->Ls;
  // Apply (L^{\prime})^{-1}
  axpby_ssp (chi,1.0,psi,     0.0,psi,0,0);      // chi[0]=psi[0]
  for (int s=1;s<Ls;s++){
    axpby_ssp_pplus(chi,1.0,psi,-lee[s-1],chi,s,s-1);// recursion Psi[s] -lee P_+ chi[s-1]
  }
  // L_m^{-1} 
  for (int s=0;s<Ls-1;s++){ // Chi[ee] = 1 - sum[s<Ls-1] -leem[s]P_- chi
    axpby_ssp_pminus(chi,1.0,chi,-leem[s],chi,Ls-1,s);
  }
  // U_m^{-1} D^{-1}
  for (int s=0;s<Ls-1;s++){
    // Chi[s] + 1/d chi[s] 
    axpby_ssp_pplus(chi,1.0/dee[s],chi,-ueem[s]/dee[Ls-1],chi,s,Ls-1);
  }	
  axpby_ssp(chi,1.0/dee[Ls-1],chi,0.0,chi,Ls-1,Ls-1); // Modest avoidable 
  
  // Apply U^{-1}
  for (int s=Ls-2;s>=0;s--){
    axpby_ssp_pminus (chi,1.0,chi,-uee[s],chi,s,s+1);  // chi[Ls]
  }
}

template<class Impl>
void CayleyFermion5D<Impl>::MooeeInvDag (const FermionField &psi, FermionField &chi)
{
  chi.checkerboard=psi.checkerboard;
  int Ls=this->Ls;
  // Apply (U^{\prime})^{-dagger}
  axpby_ssp (chi,1.0,psi,     0.0,psi,0,0);      // chi[0]=psi[0]
  for (int s=1;s<Ls;s++){
    axpby_ssp_pminus(chi,1.0,psi,-uee[s-1],chi,s,s-1);
  }
  // U_m^{-\dagger} 
  for (int s=0;s<Ls-1;s++){
    axpby_ssp_pplus(chi,1.0,chi,-ueem[s],chi,Ls-1,s);
  }
  // L_m^{-\dagger} D^{-dagger}
  for (int s=0;s<Ls-1;s++){
    axpby_ssp_pminus(chi,1.0/dee[s],chi,-leem[s]/dee[Ls-1],chi,s,Ls-1);
  }	
  axpby_ssp(chi,1.0/dee[Ls-1],chi,0.0,chi,Ls-1,Ls-1); // Modest avoidable 
  
  // Apply L^{-dagger}
  for (int s=Ls-2;s>=0;s--){
    axpby_ssp_pplus (chi,1.0,chi,-lee[s],chi,s,s+1);  // chi[Ls]
  }
}


#ifdef CAYLEY_DPERP_LINALG
  INSTANTIATE(WilsonImplF);
  INSTANTIATE(WilsonImplD);
  INSTANTIATE(GparityWilsonImplF);
  INSTANTIATE(GparityWilsonImplD);
#endif

}
}
