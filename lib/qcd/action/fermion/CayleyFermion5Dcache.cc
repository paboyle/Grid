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
  // Pplus  backwards..
template<class Impl>  
void CayleyFermion5D<Impl>::M5D(const FermionField &psi,
				const FermionField &phi, 
				FermionField &chi,
				std::vector<Coeff_t> &lower,
				std::vector<Coeff_t> &diag,
				std::vector<Coeff_t> &upper)
{
  int Ls =this->Ls;
  GridBase *grid=psi._grid;
  assert(phi.checkerboard == psi.checkerboard);
  chi.checkerboard=psi.checkerboard;
  // Flops = 6.0*(Nc*Ns) *Ls*vol
  M5Dcalls++;
  M5Dtime-=usecond();
PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss+=Ls){ // adds Ls
    for(int s=0;s<Ls;s++){
      auto tmp = psi._odata[0];
      if ( s==0 ) {
 	                            spProj5m(tmp,psi._odata[ss+s+1]);
	chi[ss+s]=diag[s]*phi[ss+s]+upper[s]*tmp;

	                    spProj5p(tmp,psi._odata[ss+Ls-1]);
	chi[ss+s]=chi[ss+s]+lower[s]*tmp;
      } else if ( s==(Ls-1)) {
	                            spProj5m(tmp,psi._odata[ss+0]);
	chi[ss+s]=diag[s]*phi[ss+s]+upper[s]*tmp;

 	                    spProj5p(tmp,psi._odata[ss+s-1]);
	chi[ss+s]=chi[ss+s]+lower[s]*tmp;
      } else { 
	                            spProj5m(tmp,psi._odata[ss+s+1]);
	chi[ss+s]=diag[s]*phi[ss+s]+upper[s]*tmp;

	                    spProj5p(tmp,psi._odata[ss+s-1]);
	chi[ss+s]=chi[ss+s]+lower[s]*tmp;
      }
    }
  }
  M5Dtime+=usecond();
}

template<class Impl>  
void CayleyFermion5D<Impl>::M5Ddag(const FermionField &psi,
				   const FermionField &phi, 
				   FermionField &chi,
				   std::vector<Coeff_t> &lower,
				   std::vector<Coeff_t> &diag,
				   std::vector<Coeff_t> &upper)
{
  int Ls =this->Ls;
  GridBase *grid=psi._grid;
  assert(phi.checkerboard == psi.checkerboard);
  chi.checkerboard=psi.checkerboard;

  // Flops = 6.0*(Nc*Ns) *Ls*vol
  M5Dcalls++;
  M5Dtime-=usecond();
PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss+=Ls){ // adds Ls
    auto tmp = psi._odata[0];
    for(int s=0;s<Ls;s++){
      if ( s==0 ) {
	spProj5p(tmp,psi._odata[ss+s+1]);
	chi[ss+s]=diag[s]*phi[ss+s]+upper[s]*tmp;

	spProj5m(tmp,psi._odata[ss+Ls-1]);
	chi[ss+s]=chi[ss+s]+lower[s]*tmp;
      } else if ( s==(Ls-1)) {
	spProj5p(tmp,psi._odata[ss+0]);
	chi[ss+s]=diag[s]*phi[ss+s]+upper[s]*tmp;

	spProj5m(tmp,psi._odata[ss+s-1]);
	chi[ss+s]=chi[ss+s]+lower[s]*tmp;
      } else { 
	spProj5p(tmp,psi._odata[ss+s+1]);
	chi[ss+s]=diag[s]*phi[ss+s]+upper[s]*tmp;

	spProj5m(tmp,psi._odata[ss+s-1]);
	chi[ss+s]=chi[ss+s]+lower[s]*tmp;
      }
    }
  }
  M5Dtime+=usecond();
}

template<class Impl>
void CayleyFermion5D<Impl>::MooeeInv    (const FermionField &psi, FermionField &chi)
{
  GridBase *grid=psi._grid;
  int Ls=this->Ls;

  chi.checkerboard=psi.checkerboard;

  MooeeInvCalls++;
  MooeeInvTime-=usecond();

PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss+=Ls){ // adds Ls
    auto tmp = psi._odata[0];

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
  }

  MooeeInvTime+=usecond();

}

template<class Impl>
void CayleyFermion5D<Impl>::MooeeInvDag (const FermionField &psi, FermionField &chi)
{
  GridBase *grid=psi._grid;
  int Ls=this->Ls;

  assert(psi.checkerboard == psi.checkerboard);
  chi.checkerboard=psi.checkerboard;

  MooeeInvCalls++;
  MooeeInvTime-=usecond();

PARALLEL_FOR_LOOP
  for(int ss=0;ss<grid->oSites();ss+=Ls){ // adds Ls

    auto tmp = psi._odata[0];

    // Apply (U^{\prime})^{-dagger}
    chi[ss]=psi[ss];
    for (int s=1;s<Ls;s++){
                            spProj5m(tmp,chi[ss+s-1]);
      chi[ss+s] = psi[ss+s]-uee[s-1]*tmp;
    }
    // U_m^{-\dagger} 
    for (int s=0;s<Ls-1;s++){
                                   spProj5p(tmp,chi[ss+s]);
      chi[ss+Ls-1] = chi[ss+Ls-1] - ueem[s]*tmp;
    }

    // L_m^{-\dagger} D^{-dagger}
    for (int s=0;s<Ls-1;s++){
      spProj5m(tmp,chi[ss+Ls-1]);
      chi[ss+s] = (1.0/dee[s])*chi[ss+s]-(leem[s]/dee[Ls-1])*tmp;
    }	
    chi[ss+Ls-1]= (1.0/dee[Ls-1])*chi[ss+Ls-1];
  
    // Apply L^{-dagger}
    for (int s=Ls-2;s>=0;s--){
      spProj5p(tmp,chi[ss+s+1]);
      chi[ss+s] = chi[ss+s] - lee[s]*tmp;
    }
  }

  MooeeInvTime+=usecond();

}

#ifdef CAYLEY_DPERP_CACHE
  INSTANTIATE_DPERP(WilsonImplF);
  INSTANTIATE_DPERP(WilsonImplD);
  INSTANTIATE_DPERP(GparityWilsonImplF);
  INSTANTIATE_DPERP(GparityWilsonImplD);
  INSTANTIATE_DPERP(ZWilsonImplF);
  INSTANTIATE_DPERP(ZWilsonImplD);
#endif

}}
