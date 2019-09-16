/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/DomainWallEOFAFermiondense.cc

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

#include <Grid/Grid_Eigen_Dense.h>
#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/DomainWallEOFAFermion.h>

NAMESPACE_BEGIN(Grid);

/*
 * Dense matrix versions of routines
 */
template<class Impl>
void DomainWallEOFAFermion<Impl>::MooeeInvDag(const FermionField& psi, FermionField& chi)
{
  this->MooeeInternal(psi, chi, DaggerYes, InverseYes);
}

template<class Impl>
void DomainWallEOFAFermion<Impl>::MooeeInv(const FermionField& psi, FermionField& chi)
{
  this->MooeeInternal(psi, chi, DaggerNo, InverseYes);
}

template<class Impl>
void DomainWallEOFAFermion<Impl>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv)
{
  int Ls = this->Ls;
  int LLs = psi.Grid()->_rdimensions[0];
  int vol = psi.Grid()->oSites()/LLs;

  chi.Checkerboard() = psi.Checkerboard();

  assert(Ls==LLs);

  Eigen::MatrixXd Pplus  = Eigen::MatrixXd::Zero(Ls,Ls);
  Eigen::MatrixXd Pminus = Eigen::MatrixXd::Zero(Ls,Ls);

  for(int s=0;s<Ls;s++){
    Pplus(s,s)  = this->bee[s];
    Pminus(s,s) = this->bee[s];
  }

  for(int s=0; s<Ls-1; s++){
    Pminus(s,s+1) = -this->cee[s];
  }

  for(int s=0; s<Ls-1; s++){
    Pplus(s+1,s) = -this->cee[s+1];
  }

  Pplus (0,Ls-1) = this->dp;
  Pminus(Ls-1,0) = this->dm;

  Eigen::MatrixXd PplusMat ;
  Eigen::MatrixXd PminusMat;

  if(inv) {
    PplusMat  = Pplus.inverse();
    PminusMat = Pminus.inverse();
  } else {
    PplusMat  = Pplus;
    PminusMat = Pminus;
  }

  if(dag){
    PplusMat.adjointInPlace();
    PminusMat.adjointInPlace();
  }

  // For the non-vectorised s-direction this is simple

  for(auto site=0; site<vol; site++){

    SiteSpinor     SiteChi;
    SiteHalfSpinor SitePplus;
    SiteHalfSpinor SitePminus;

    for(int s1=0; s1<Ls; s1++){
      SiteChi = Zero();
      for(int s2=0; s2<Ls; s2++){
	int lex2 = s2 + Ls*site;
	if(PplusMat(s1,s2) != 0.0){
	  spProj5p(SitePplus,psi[lex2]);
	  accumRecon5p(SiteChi, PplusMat(s1,s2)*SitePplus);
	}
	if(PminusMat(s1,s2) != 0.0){
	  spProj5m(SitePminus, psi[lex2]);
	  accumRecon5m(SiteChi, PminusMat(s1,s2)*SitePminus);
	}
      }
      chi[s1+Ls*site] = SiteChi*0.5;
    }
  }
}

#ifdef DOMAIN_WALL_EOFA_DPERP_DENSE

INSTANTIATE_DPERP_DWF_EOFA(GparityWilsonImplF);
INSTANTIATE_DPERP_DWF_EOFA(GparityWilsonImplD);
INSTANTIATE_DPERP_DWF_EOFA(WilsonImplF);
INSTANTIATE_DPERP_DWF_EOFA(WilsonImplD);
INSTANTIATE_DPERP_DWF_EOFA(ZWilsonImplF);
INSTANTIATE_DPERP_DWF_EOFA(ZWilsonImplD);

template void DomainWallEOFAFermion<GparityWilsonImplF>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
template void DomainWallEOFAFermion<GparityWilsonImplD>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
template void DomainWallEOFAFermion<WilsonImplF>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
template void DomainWallEOFAFermion<WilsonImplD>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
template void DomainWallEOFAFermion<ZWilsonImplF>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
template void DomainWallEOFAFermion<ZWilsonImplD>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);

INSTANTIATE_DPERP_DWF_EOFA(GparityWilsonImplFH);
INSTANTIATE_DPERP_DWF_EOFA(GparityWilsonImplDF);
INSTANTIATE_DPERP_DWF_EOFA(WilsonImplFH);
INSTANTIATE_DPERP_DWF_EOFA(WilsonImplDF);
INSTANTIATE_DPERP_DWF_EOFA(ZWilsonImplFH);
INSTANTIATE_DPERP_DWF_EOFA(ZWilsonImplDF);

template void DomainWallEOFAFermion<GparityWilsonImplFH>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
template void DomainWallEOFAFermion<GparityWilsonImplDF>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
template void DomainWallEOFAFermion<WilsonImplFH>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
template void DomainWallEOFAFermion<WilsonImplDF>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
template void DomainWallEOFAFermion<ZWilsonImplFH>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);
template void DomainWallEOFAFermion<ZWilsonImplDF>::MooeeInternal(const FermionField& psi, FermionField& chi, int dag, int inv);

#endif

NAMESPACE_END(Grid);
