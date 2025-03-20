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
#pragma once

NAMESPACE_BEGIN(Grid)

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

  void DhopImproved(StencilImpl &st,
		    DoubledGaugeField &U, DoubledGaugeField &UUU, 
		    const FermionField &in, FermionField &out, int dag, int interior,int exterior);
  void DhopNaive(StencilImpl &st,
		 DoubledGaugeField &U,
		 const FermionField &in, FermionField &out, int dag, int interior,int exterior);
  
  void DhopDirKernel(StencilImpl &st, DoubledGaugeFieldView &U, DoubledGaugeFieldView &UUU, SiteSpinor * buf,
		     int sF, int sU, const FermionFieldView &in, FermionFieldView &out, int dir,int disp);
 protected:    

   ///////////////////////////////////////////////////////////////////////////////////////
   // Generic Nc kernels
   ///////////////////////////////////////////////////////////////////////////////////////
   template<int Naik> 
   static accelerator_inline
   void DhopSiteGeneric(StencilView &st, 
			DoubledGaugeFieldView &U, DoubledGaugeFieldView &UUU, 
			SiteSpinor * buf, int LLs, int sU, 
			const FermionFieldView &in, FermionFieldView &out,int dag);
   
   template<int Naik> static accelerator_inline
   void DhopSiteGenericInt(StencilView &st, 
			   DoubledGaugeFieldView &U, DoubledGaugeFieldView &UUU, 
			   SiteSpinor * buf, int LLs, int sU, 
			   const FermionFieldView &in, FermionFieldView &out,int dag);
   
   template<int Naik> static accelerator_inline
   void DhopSiteGenericExt(StencilView &st, 
			   DoubledGaugeFieldView &U, DoubledGaugeFieldView &UUU,
			   SiteSpinor * buf, int LLs, int sU, 
			   const FermionFieldView &in, FermionFieldView &out,int dag);

   ///////////////////////////////////////////////////////////////////////////////////////
   // Nc=3 specific kernels
   ///////////////////////////////////////////////////////////////////////////////////////
   
   template<int Naik> static accelerator_inline
   void DhopSiteHand(StencilView &st, 
		     DoubledGaugeFieldView &U,DoubledGaugeFieldView &UUU, 
		     SiteSpinor * buf, int LLs, int sU, 
		     const FermionFieldView &in, FermionFieldView &out,int dag);
   
   template<int Naik> static accelerator_inline
   void DhopSiteHandInt(StencilView &st, 
			DoubledGaugeFieldView &U,DoubledGaugeFieldView &UUU, 
			SiteSpinor * buf, int LLs, int sU, 
			const FermionFieldView &in, FermionFieldView &out,int dag);
   
   template<int Naik> static accelerator_inline
   void DhopSiteHandExt(StencilView &st, 
			DoubledGaugeFieldView &U,DoubledGaugeFieldView &UUU, 
			SiteSpinor * buf, int LLs, int sU, 
			const FermionFieldView &in, FermionFieldView &out,int dag);

   ///////////////////////////////////////////////////////////////////////////////////////
   // Asm Nc=3 specific kernels
   ///////////////////////////////////////////////////////////////////////////////////////
   
   void DhopSiteAsm(StencilView &st, 
		    DoubledGaugeFieldView &U,DoubledGaugeFieldView &UUU, 
		    SiteSpinor * buf, int LLs, int sU, 
		    const FermionFieldView &in, FermionFieldView &out,int dag);
  
public:

  StaggeredKernels(const ImplParams &p = ImplParams());

};
NAMESPACE_END(Grid);    
