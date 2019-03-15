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
  enum { CommsAndCompute, CommsThenCompute };
  static int Opt;
  static int Comms;
};
 
template<class Impl> class StaggeredKernels : public FermionOperator<Impl> , public StaggeredKernelsStatic { 
 public:
   
  INHERIT_IMPL_TYPES(Impl);
  typedef FermionOperator<Impl> Base;
   
public:
    
   void DhopDir(StencilImpl &st, DoubledGaugeField &U, DoubledGaugeField &UUU, SiteSpinor * buf,
		      int sF, int sU, const FermionField &in, FermionField &out, int dir,int disp);

   ///////////////////////////////////////////////////////////////////////////////////////
   // Generic Nc kernels
   ///////////////////////////////////////////////////////////////////////////////////////
   void DhopSiteGeneric(StencilImpl &st, LebesgueOrder &lo, 
			DoubledGaugeField &U, DoubledGaugeField &UUU, 
			SiteSpinor * buf, int LLs, int sU, 
			const FermionField &in, FermionField &out,int dag);
   void DhopSiteGenericInt(StencilImpl &st, LebesgueOrder &lo, 
			   DoubledGaugeField &U, DoubledGaugeField &UUU, 
			   SiteSpinor * buf, int LLs, int sU, 
			   const FermionField &in, FermionField &out,int dag);
   void DhopSiteGenericExt(StencilImpl &st, LebesgueOrder &lo, 
			   DoubledGaugeField &U, DoubledGaugeField &UUU,
			   SiteSpinor * buf, int LLs, int sU, 
			   const FermionField &in, FermionField &out,int dag);

   ///////////////////////////////////////////////////////////////////////////////////////
   // Nc=3 specific kernels
   ///////////////////////////////////////////////////////////////////////////////////////
   void DhopSiteHand(StencilImpl &st, LebesgueOrder &lo, 
		     DoubledGaugeField &U,DoubledGaugeField &UUU, 
		     SiteSpinor * buf, int LLs, int sU, 
		     const FermionField &in, FermionField &out,int dag);
   void DhopSiteHandInt(StencilImpl &st, LebesgueOrder &lo, 
			DoubledGaugeField &U,DoubledGaugeField &UUU, 
			SiteSpinor * buf, int LLs, int sU, 
			const FermionField &in, FermionField &out,int dag);
   void DhopSiteHandExt(StencilImpl &st, LebesgueOrder &lo, 
			DoubledGaugeField &U,DoubledGaugeField &UUU, 
			SiteSpinor * buf, int LLs, int sU, 
			const FermionField &in, FermionField &out,int dag);

   ///////////////////////////////////////////////////////////////////////////////////////
   // Asm Nc=3 specific kernels
   ///////////////////////////////////////////////////////////////////////////////////////
   void DhopSiteAsm(StencilImpl &st, LebesgueOrder &lo, 
		    DoubledGaugeField &U,DoubledGaugeField &UUU, 
		    SiteSpinor * buf, int LLs, int sU, 
		    const FermionField &in, FermionField &out,int dag);
   ///////////////////////////////////////////////////////////////////////////////////////////////////
   // Generic interface; fan out to right routine
   ///////////////////////////////////////////////////////////////////////////////////////////////////
   void DhopSite(StencilImpl &st, LebesgueOrder &lo, 
		 DoubledGaugeField &U, DoubledGaugeField &UUU, 
		 SiteSpinor * buf, int LLs, int sU,
		 const FermionField &in, FermionField &out, int interior=1,int exterior=1);

   void DhopSiteDag(StencilImpl &st, LebesgueOrder &lo, 
		    DoubledGaugeField &U, DoubledGaugeField &UUU, 
		    SiteSpinor * buf, int LLs, int sU,
		    const FermionField &in, FermionField &out, int interior=1,int exterior=1);

   void DhopSite(StencilImpl &st, LebesgueOrder &lo, 
		 DoubledGaugeField &U, DoubledGaugeField &UUU, 
		 SiteSpinor * buf, int LLs, int sU,
		 const FermionField &in, FermionField &out, int dag, int interior,int exterior);
  
public:

  StaggeredKernels(const ImplParams &p = ImplParams());

};
    
}}

#endif
