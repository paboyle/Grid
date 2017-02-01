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

#include <Grid/Eigen/Dense>
#include <Grid/Grid.h>


namespace Grid {
namespace QCD {
  /*
   * Dense matrix versions of routines
   */

  /*
template<class Impl>
void CayleyFermion5D<Impl>::MooeeInvDag (const FermionField &psi, FermionField &chi)
{
  this->MooeeInternal(psi,chi,DaggerYes,InverseYes);
}
  
template<class Impl>
void CayleyFermion5D<Impl>::MooeeInv(const FermionField &psi, FermionField &chi)
{
  this->MooeeInternal(psi,chi,DaggerNo,InverseYes);
}
  */
template<class Impl>
void CayleyFermion5D<Impl>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv)
{
  int Ls=this->Ls;
  int LLs = psi._grid->_rdimensions[0];
  int vol = psi._grid->oSites()/LLs;
  
  chi.checkerboard=psi.checkerboard;
  
  assert(Ls==LLs);
  
  Eigen::MatrixXd Pplus  = Eigen::MatrixXd::Zero(Ls,Ls);
  Eigen::MatrixXd Pminus = Eigen::MatrixXd::Zero(Ls,Ls);
  
  for(int s=0;s<Ls;s++){
    Pplus(s,s) = bee[s];
    Pminus(s,s)= bee[s];
  }
  
  for(int s=0;s<Ls-1;s++){
    Pminus(s,s+1) = -cee[s];
  }
  
  for(int s=0;s<Ls-1;s++){
    Pplus(s+1,s) = -cee[s+1];
  }
  Pplus (0,Ls-1) = mass*cee[0];
  Pminus(Ls-1,0) = mass*cee[Ls-1];
  
  Eigen::MatrixXd PplusMat ;
  Eigen::MatrixXd PminusMat;
  
  if ( inv ) {
    PplusMat =Pplus.inverse();
    PminusMat=Pminus.inverse();
  } else { 
    PplusMat =Pplus;
    PminusMat=Pminus;
  }
  
  if(dag){
    PplusMat.adjointInPlace();
    PminusMat.adjointInPlace();
  }

  // For the non-vectorised s-direction this is simple
  
  for(auto site=0;site<vol;site++){
    
    SiteSpinor     SiteChi;
    SiteHalfSpinor SitePplus;
    SiteHalfSpinor SitePminus;
    
    for(int s1=0;s1<Ls;s1++){
      SiteChi =zero;
      for(int s2=0;s2<Ls;s2++){
	int lex2 = s2+Ls*site;
	
	if ( PplusMat(s1,s2) != 0.0 ) {
	  spProj5p(SitePplus,psi[lex2]);
	  accumRecon5p(SiteChi,PplusMat (s1,s2)*SitePplus);
	}
	
	if ( PminusMat(s1,s2) != 0.0 ) {
	  spProj5m(SitePminus,psi[lex2]);
	  accumRecon5m(SiteChi,PminusMat(s1,s2)*SitePminus);
	}
      }
      chi[s1+Ls*site] = SiteChi*0.5;
    }
  }
}

template void CayleyFermion5D<GparityWilsonImplF>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<GparityWilsonImplD>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<WilsonImplF>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);
template void CayleyFermion5D<WilsonImplD>::MooeeInternal(const FermionField &psi, FermionField &chi,int dag, int inv);

}}
