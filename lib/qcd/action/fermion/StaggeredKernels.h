/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/StaggeredKernels.h

Copyright (C) 2015

Author: Azusa Yamaguchi, Peter Boyle

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef GRID_QCD_STAGGERED_KERNELS_H
#define GRID_QCD_STAGGERED_KERNELS_H

namespace Grid {
namespace QCD {

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Helper routines that implement Staggered stencil for a single site.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class StaggeredKernelsStatic { 
 public:
  enum { OptGeneric, OptHandUnroll, OptInlineAsm };
  // S-direction is INNERMOST and takes no part in the parity.
  static int Opt;  // these are a temporary hack
};
 
template<class Impl> class StaggeredKernels : public FermionOperator<Impl> , public StaggeredKernelsStatic { 
 public:
   
  INHERIT_IMPL_TYPES(Impl);
  typedef FermionOperator<Impl> Base;
   
public:
    
   void DhopDir(StencilImpl &st, DoubledGaugeField &U, DoubledGaugeField &UUU, SiteSpinor * buf,
		      int sF, int sU, const FermionField &in, FermionField &out, int dir,int disp);

   void DhopSiteDepth(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteSpinor * buf,
		     int sF, int sU, const FermionField &in, SiteSpinor &out,int threeLink);


   void DhopSiteDepthHand(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, SiteSpinor * buf,
		     int sF, int sU, const FermionField &in, SiteSpinor&out,int threeLink);

   void DhopSiteHand(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, DoubledGaugeField &UUU,SiteSpinor * buf,
		     int LLs, int sU, const FermionField &in, FermionField &out, int dag);

   void DhopSiteAsm(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U,DoubledGaugeField &UUU, SiteSpinor * buf,
			 int LLs, int sU, const FermionField &in, FermionField &out);
      
   void DhopSite(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, DoubledGaugeField &UUU, SiteSpinor * buf,
		int sF, int sU, const FermionField &in, FermionField &out);

   void DhopSiteDag(StencilImpl &st, LebesgueOrder &lo, DoubledGaugeField &U, DoubledGaugeField &UUU, SiteSpinor *buf, 
                   int LLs, int sU, const FermionField &in, FermionField &out);
  
public:

  StaggeredKernels(const ImplParams &p = ImplParams());

};
    
}}

#endif
